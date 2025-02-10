#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <type_traits>
#include "bigint.h"

static constexpr int KARATSUBA_THRESHOLD = 30;

// https://graphics.stanford.edu/~seander/bithacks.html
static ::uint32_t de_bruijn(::uint32_t n)
{
    static ::uint32_t multiply_de_bruijn_bit_position[32] =
    {
        0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
        8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
    };
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    return multiply_de_bruijn_bit_position[(n * 0x07C4ACDDU) >> 27];
}

namespace bigint
{
    static limb_t g_zero = 0;
    static void * (*g_allocate)(size_t) = ::malloc;
    static void (*g_deallocate)(void *) = ::free;

    static int mul_with_scratch_memory(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size, limb_t * memory);

    static int long_multiplication(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size)
    {
        constexpr limb_double_t base = static_cast<limb_double_t>(1) << (sizeof(limb_t) * 8);

        for (int i = 0; i < (lhs_size + rhs_size); ++i)
            result[i] = 0;

        for (int i = 0; i < lhs_size; ++i)
        {
            limb_t carry = 0;
            for (int j = 0; j < rhs_size; ++j)
            {
                limb_double_t product = static_cast<limb_double_t>(result[i + j]);
                product += carry;
                product += static_cast<limb_double_t>(lhs[i]) * rhs[j];
                carry = static_cast<limb_t>(product / base);
                result[i + j] = product % base;
            }
            result[i + rhs_size] = carry;
        }
        return result[lhs_size + rhs_size - 1] == 0 ? lhs_size + rhs_size - 1 : lhs_size + rhs_size;
    }

    static int karatsuba(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size, limb_t * memory)
    {
        int const m = std::max(lhs_size, rhs_size);
        int const m2 = m / 2;

        limb_t const * const low1 = lhs;
        limb_t const * const low2 = rhs;
        limb_t const * high1 = &g_zero;
        limb_t const * high2 = &g_zero;
        int low1_size = std::min(m2, lhs_size);
        int low2_size = std::min(m2, rhs_size);
        int high1_size = std::min(m - m2, lhs_size - m2);
        int high2_size = std::min(m - m2, rhs_size - m2);

        while (low1_size > 1 && low1[low1_size - 1] == 0)
            low1_size--;
        while (low2_size > 1 && low2[low2_size - 1] == 0)
            low2_size--;

        if (high1_size <= 0)
            high1_size = 1;
        else
            high1 = lhs + m2;

        if (high2_size <= 0)
            high2_size = 1;
        else
            high2 = rhs + m2;

        limb_t * const z1_lhs = result;
        limb_t * const z1_rhs = result + m2 + 2;
        int const z1_lhs_size = add(z1_lhs, low1, high1, low1_size, high1_size);
        assert(z1_lhs_size <= (m2 + 2));
        int const z1_rhs_size = add(z1_rhs, low2, high2, low2_size, high2_size);
        assert((m2 + 2 + z1_rhs_size) <= (lhs_size + rhs_size));

        limb_t * const z1 = memory;
        int z1_size = mul_with_scratch_memory(z1, z1_lhs, z1_rhs, z1_lhs_size, z1_rhs_size, memory + z1_lhs_size + z1_rhs_size);
        memory += z1_size;

        limb_t * const z0 = result;
        limb_t * const z2 = result + m2 * 2;

        int const z0_size = mul_with_scratch_memory(z0, low1, low2, low1_size, low2_size, memory);
        assert(z0_size <= (m2 * 2));
        int const z2_size = mul_with_scratch_memory(z2, high1, high2, high1_size, high2_size, memory);

        z1_size = sub(z1, z1, z2, z1_size, z2_size);
        z1_size = sub(z1, z1, z0, z1_size, z0_size);

        int result_size = z0_size;
        if (!(z2_size == 1 && z2[0] == 0))
        {
            result_size = m2 * 2 + z2_size;
            for (limb_t * z = z0 + z0_size; z != z2; ++z)
                *z = 0;
        }
        assert(result_size <= (lhs_size + rhs_size));
        if (!(z1_size == 1 && z1[0] == 0))
        {
            result_size = add(result + m2, result + m2, z1, result_size - m2, z1_size) + m2;
        }

        return result_size;
    }

    static int mul_with_scratch_memory(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size, limb_t * memory)
    {
        if ((lhs_size == 1 && *lhs == 0) || (rhs_size == 1 && *rhs == 0))
        {
            *result = 0;
            return 1;
        }
        if (lhs_size == 1 && rhs_size == 1)
        {
            limb_double_t const product = static_cast<limb_double_t>(*lhs) * static_cast<limb_double_t>(*rhs);
            result[0] = static_cast<limb_t>(product);
            result[1] = static_cast<limb_t>(product >> (sizeof(limb_t) * 8));
            return result[1] == 0 ? 1 : 2;
        }

        if (lhs_size > KARATSUBA_THRESHOLD && rhs_size > KARATSUBA_THRESHOLD)
        {
            return karatsuba(result, lhs, rhs, lhs_size, rhs_size, memory);
        }
        return long_multiplication(result, lhs, rhs, lhs_size, rhs_size);
    }

    // https://skanthak.hier-im-netz.de/division.html
    template <bool ignore_quotient, bool ignore_remainder>
    static int algorithm_d(limb_t * quotient, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size, limb_t * remainder, int * remainder_size)
    {
        constexpr limb_t limb_bits = sizeof(limb_t) * 8;
        constexpr limb_double_t base = static_cast<limb_double_t>(1) << (sizeof(limb_t) * 8);

        /* Normalize by shifting rhs left just enough so that its high-order
        bit is on, and shift lhs left the same amount. We may have to append a
        high-order digit on the dividend; we do that unconditionally. */

        limb_t const shift = limb_bits - de_bruijn(rhs[rhs_size - 1]) - 1;
        limb_t * const rhs_normalized = static_cast<limb_t *>(g_allocate(sizeof(limb_t) * rhs_size));
        for (int i = rhs_size - 1; i > 0; i--)
            rhs_normalized[i] = (rhs[i] << shift) | (static_cast<limb_double_t>(rhs[i - 1]) >> (limb_bits - shift));
        rhs_normalized[0] = rhs[0] << shift;

        limb_t * const lhs_normalized = static_cast<limb_t *>(g_allocate(sizeof(limb_t) * (lhs_size + 1)));
        lhs_normalized[lhs_size] = static_cast<limb_double_t>(lhs[lhs_size - 1]) >> (limb_bits - shift);
        for (int i = lhs_size - 1; i > 0; i--)
            lhs_normalized[i] = (lhs[i] << shift) | (static_cast<limb_double_t>(lhs[i - 1]) >> (limb_bits - shift));
        lhs_normalized[0] = lhs[0] << shift;

        for (int j = lhs_size - rhs_size; j >= 0; j--) {
            // Compute estimate quotient_digit of quotient[j].
            limb_double_t quotient_digit = (lhs_normalized[j + rhs_size] * base + lhs_normalized[j + rhs_size - 1]) / rhs_normalized[rhs_size - 1];
            limb_double_t remainder_digit = (lhs_normalized[j + rhs_size] * base + lhs_normalized[j + rhs_size - 1]) % rhs_normalized[rhs_size - 1];

            while (quotient_digit >= base || static_cast<limb_t>(quotient_digit) * static_cast<limb_double_t>(rhs_normalized[rhs_size - 2]) > base * remainder_digit + lhs_normalized[j + rhs_size - 2])
            {
                quotient_digit = quotient_digit - 1;
                remainder_digit = remainder_digit + rhs_normalized[rhs_size - 1];
                if (remainder_digit >= base)
                {
                    break;
                }
            }

            // Multiply and subtract.
            typename std::make_signed<limb_double_t>::type diff = 0;
            typename std::make_signed<limb_double_t>::type sum;
            for (int i = 0; i < rhs_size; i++) {
                limb_double_t const product = static_cast<limb_t>(quotient_digit) * static_cast<limb_double_t>(rhs_normalized[i]);
                sum = lhs_normalized[i + j] - diff - (product & (base - 1));
                lhs_normalized[i + j] = static_cast<limb_t>(sum);
                diff = (product >> limb_bits) - (sum >> limb_bits);
            }
            sum = lhs_normalized[j + rhs_size] - diff;
            lhs_normalized[j + rhs_size] = static_cast<limb_t>(sum);

            if (!ignore_quotient)
                quotient[j] = static_cast<limb_t>(quotient_digit);

            // If we subtracted too much, add back
            if (sum < 0) {
                if (!ignore_quotient)
                    quotient[j] = quotient[j] - 1;
                diff = 0;
                for (int i = 0; i < rhs_size; i++) {
                    sum = static_cast<limb_double_t>(lhs_normalized[i + j]) + rhs_normalized[i] + diff;
                    lhs_normalized[i + j] = static_cast<limb_t>(sum);
                    diff = sum >> limb_bits;
                }
                lhs_normalized[j + rhs_size] += static_cast<limb_t>(diff);
            }
        }

        // If the caller wants the remainder, unnormalize it and pass it back.
        if (!ignore_remainder) {
            for (int i = 0; i < rhs_size - 1; i++)
                remainder[i] = (lhs_normalized[i] >> shift) | (static_cast<limb_double_t>(lhs_normalized[i + 1]) << (limb_bits - shift));
            remainder[rhs_size - 1] = lhs_normalized[rhs_size - 1] >> shift;
            *remainder_size = rhs_size;
            while (*remainder_size > 1 && remainder[*remainder_size - 1] == 0)
                (*remainder_size)--;
        }
        g_deallocate(lhs_normalized);
        g_deallocate(rhs_normalized);
        if (!ignore_quotient)
        {
            int quotient_size = lhs_size - rhs_size + 1;
            if (quotient[quotient_size - 1] == 0 && quotient_size > 1)
                quotient_size--;
            return quotient_size;
        }
        return 0;
    }

    static int divide_by_one_digit(limb_t * quotient, limb_t const * lhs, limb_t rhs, int lhs_size, limb_t * remainder)
    {
        limb_t lhs_digit = lhs[lhs_size - 1];
        limb_double_t num;

        quotient[lhs_size - 1] = lhs_digit / rhs;
        lhs_digit = lhs_digit % rhs;
        for (int i = lhs_size - 2; i >= 0; --i)
        {
            num = static_cast<limb_double_t>(lhs_digit) << (sizeof(limb_t) * 8);
            num |= lhs[i];
            quotient[i] = static_cast<limb_t>(num / rhs);
            lhs_digit = num % rhs;
            if (lhs_digit >= rhs)
            {
                quotient[i] = lhs_digit / rhs;
                lhs_digit = lhs_digit % rhs;
            }
        }
        *remainder = lhs_digit;
        if (lhs_size > 1 && quotient[lhs_size - 1] == 0)
            return lhs_size - 1;
        return lhs_size;
    }

    static limb_t mod_by_one_digit(limb_t const * lhs, limb_t rhs, int lhs_size)
    {
        limb_t lhs_digit = lhs[lhs_size - 1];
        limb_double_t num;

        lhs_digit = lhs_digit % rhs;
        for (int i = lhs_size - 2; i >= 0; --i)
        {
            num = static_cast<limb_double_t>(lhs_digit) << (sizeof(limb_t) * 8);
            num |= lhs[i];
            lhs_digit = num % rhs;
            if (lhs_digit >= rhs)
            {
                lhs_digit = lhs_digit % rhs;
            }
        }
        return lhs_digit;
    }

    static constexpr limb_t largest_base10_numerator_fitting_in_limb()
    {
        constexpr limb_double_t limit = static_cast<limb_double_t>(1) << 32;
        limb_double_t n = 1;
        while (n < limit)
        {
            n *= 10;
        }
        return static_cast<limb_t>(n / 10);
    }

    static constexpr int largest_base10_numerator_fitting_in_limb_size()
    {
        constexpr limb_double_t limit = static_cast<limb_double_t>(1) << 32;
        limb_double_t n = 1;
        int result = 0;
        while (n < limit)
        {
            n *= 10;
            result++;
        }
        return result - 1;
    }

    void set_allocator(void * (*alloc)(size_t), void (*dealloc)(void *))
    {
        g_allocate = alloc == nullptr ? ::malloc : alloc;
        g_deallocate = dealloc == nullptr ? ::free : dealloc;
    }

    int from_base10(limb_t * digits, char const * str, int length)
    {
        int result_size = 1;

        *digits = 0;
        for (int i = 0; i < length; ++i)
        {
            constexpr limb_double_t base = static_cast<limb_double_t>(1) << (sizeof(limb_t) * 8);

            limb_t carry = str[i] - '0';
            if (carry >= 10)
                return -1;
            for (int j = 0; j < result_size; ++j)
            {
                limb_double_t const product = carry + static_cast<limb_double_t>(digits[j]) * 10;
                carry = static_cast<limb_t>(product / base);
                digits[j] = product % base;
            }
            if (carry > 0)
            {
                digits[result_size] = carry;
                result_size++;
            }
        }
        return result_size;
    }

    int from_base16(limb_t * digits, char const * str, int length)
    {
        int index = -1;
        limb_t bitshift = 0;

        while (length > 1 && *str == '0')
        {
            length--;
            str++;
        }

        while (length > 0)
        {
            if ((bitshift % (sizeof(limb_t) * 8)) == 0)
            {
                bitshift = 0;
                index++;
                digits[index] = 0;
            }
            --length;
            if (str[length] >= 'a' && str[length] <= 'f')
            {
                digits[index] |= ((str[length] - 'a' + 10) << bitshift);
            }
            else if (str[length] >= 'A' && str[length] <= 'F')
            {
                digits[index] |= ((str[length] - 'A' + 10) << bitshift);
            }
            else if (str[length] >= '0' && str[length] <= '9')
            {
                digits[index] |= ((str[length] - '0') << bitshift);
            }
            else
                return -1;
            bitshift += 4;
        }
        return index + 1;
    }

    int to_base10(char * str, int capacity, limb_t const * digits, int size)
    {
        constexpr limb_t numerator = largest_base10_numerator_fitting_in_limb();
        constexpr int numerator_size = largest_base10_numerator_fitting_in_limb_size();
        limb_t * const copy = static_cast<limb_t *>(g_allocate(sizeof(limb_t) * size));
        int len = 0;

        for (int i = 0; i < size; ++i)
            copy[i] = digits[i];
        while (size > 1 || *copy >= numerator)
        {
            limb_t remainder;
            size = divide_by_one_digit(copy, copy, numerator, size, &remainder);
            for (int i = 0; i < numerator_size; ++i)
            {
                if (len < capacity)
                    str[len] = (remainder % 10) + '0';
                ++len;
                remainder /= 10;
            }
        }
        while (*copy > 0)
        {
            if (len < capacity)
                str[len] = (*copy % 10) + '0';
            ++len;
            *copy /= 10;
        }
        g_deallocate(copy);
        
        int i = 0;
        int j = std::min(capacity, len - 1);
        while (i < j)
        {
            char const tmp = str[i];
            str[i++] = str[j];
            str[j--] = tmp;
        }
        return len;
    }

    int to_base16(char * str, int capacity, limb_t const * digits, int size)
    {
        int len = 0;
        limb_t shift = sizeof(limb_t) * 8;
        
        while (((digits[size - 1] >> (shift - 4)) & 0xf) == 0 && shift > 4)
        {
            shift -= 4;
        }
        for (int i = size - 1; i >= 0; --i)
        {
            while (shift > 0)
            {
                shift -= 4;
                if (len < capacity)
                {
                    limb_t const hexdigit = (digits[i] >> shift) & 0xf;
                    if (hexdigit >= 0 && hexdigit < 10)
                        str[len] = '0' + static_cast<char>(hexdigit);
                    else
                        str[len] = 'a' + static_cast<char>(hexdigit - 10);
                }
                ++len;
            }
            shift = sizeof(limb_t) * 8;
        }
        return len;
    }

    int compare(limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size)
    {
        if (lhs_size > rhs_size)
            return 1;
        if (lhs_size < rhs_size)
            return -1;
        for (int i = lhs_size - 1; i >= 0; --i)
        {
            if (lhs[i] > rhs[i])
                return 1;
            if (lhs[i] < rhs[i])
                return -1;
        }
        return 0;
    }

    int add(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size)
    {
        if (lhs_size < rhs_size)
        {
            limb_t const * const tmp = lhs;
            lhs = rhs;
            rhs = tmp;
            int const tmp_size = lhs_size;
            lhs_size = rhs_size;
            rhs_size = tmp_size;
        }

        int size = 0;
        limb_t carry = 0;

        while (size < rhs_size)
        {
            limb_t digit = lhs[size];
            digit += carry;
            carry = (digit < carry) ? 1 : 0;
            digit += rhs[size];
            carry = (digit < rhs[size]) ? 1 : carry;
            result[size] = digit;
            ++size;
        }
        while (size < lhs_size)
        {
            limb_t digit = lhs[size];
            digit += carry;
            carry = (digit < carry) ? 1 : 0;
            result[size] = digit;
            ++size;
        }
        if (carry > 0)
            result[size] = carry;
        size += carry;
        return size;
    }

    int sub(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size)
    {
        assert(lhs_size >= rhs_size);
        int size = 0;
        limb_t carry = 0;

        while (size < rhs_size)
        {
            limb_t digit = lhs[size];
            digit -= carry;
            carry = (digit > lhs[size]) ? 1 : 0;
            limb_t const prev = digit;
            digit -= rhs[size];
            carry = (digit > prev) ? 1 : carry;
            result[size] = digit;
            ++size;
        }
        while (size < lhs_size)
        {
            limb_t digit = lhs[size];
            digit -= carry;
            carry = (digit > lhs[size]) ? 1 : 0;
            result[size] = digit;
            ++size;
        }
        assert(carry == 0);
        while (size > 1 && result[size - 1] == 0)
        {
            --size;
        }
        return size;
    }

    int mul(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size)
    {
        if (lhs_size > KARATSUBA_THRESHOLD && rhs_size > KARATSUBA_THRESHOLD)
        {
            int const capacity = (lhs_size + rhs_size) * 2;
            limb_t * const memory = static_cast<limb_t *>(g_allocate(sizeof(limb_t) * capacity));
            int const result_size = mul_with_scratch_memory(result, lhs, rhs, lhs_size, rhs_size, memory);
            g_deallocate(memory);
            return result_size;
        }
        else
        {
            return mul_with_scratch_memory(result, lhs, rhs, lhs_size, rhs_size, nullptr);
        }
    }

    int div(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size)
    {
        assert(!(rhs_size == 1 && *rhs == 0));

        if (lhs_size < rhs_size)
        {
            *result = 0;
            return 1;
        }

        if (rhs_size > 1)
        {
            return algorithm_d<false, true>(result, lhs, rhs, lhs_size, rhs_size, nullptr, nullptr);
        }
        else
        {
            limb_t remainder_digit;
            return divide_by_one_digit(result, lhs, rhs[0], lhs_size, &remainder_digit);
        }
    }

    int mod(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size)
    {
        assert(!(rhs_size == 1 && *rhs == 0));

        if (lhs_size < rhs_size)
        {
            for (int i = 0; i < lhs_size; ++i)
                result[i] = lhs[i];
            return lhs_size;
        }

        if (rhs_size > 1)
        {
            int result_size;
            algorithm_d<true, false>(nullptr, lhs, rhs, lhs_size, rhs_size, result, &result_size);
            return result_size;
        }
        else
        {
            *result = mod_by_one_digit(lhs, rhs[0], lhs_size);
            return 1;
        }
    }

    number::number()
    {
        m_size = 1;
        m_capacity = 0;
        m_negative = false;
        m_digits = &g_zero;
    }

    number::number(number const & other)
    {
        m_size = other.m_size;
        m_capacity = other.m_size;
        m_negative = other.m_negative;
        m_digits = static_cast<limb_t *>(g_allocate(m_capacity * sizeof(limb_t)));
        for (int i = 0; i < m_size; ++i)
            m_digits[i] = other.m_digits[i];
    }

    number::number(number && other)
    {
        m_size = other.m_size;
        m_capacity = other.m_capacity;
        m_negative = other.m_negative;
        m_digits = other.m_digits;
        other.m_size = 1;
        other.m_capacity = 0;
        other.m_negative = false;
        other.m_digits = &g_zero;
    }

    number & number::operator=(number const & other)
    {
        if (this != &other)
        {
            if (m_digits != &g_zero)
                g_deallocate(m_digits);
            m_size = other.m_size;
            m_capacity = other.m_size;
            m_negative = other.m_negative;
            m_digits = static_cast<limb_t *>(g_allocate(m_capacity * sizeof(limb_t)));
            for (int i = 0; i < m_size; ++i)
                m_digits[i] = other.m_digits[i];
        }
        return *this;
    }

    number & number::operator=(number && other)
    {
        if (this != &other)
        {
            if (m_digits != &g_zero)
                g_deallocate(m_digits);
            m_size = other.m_size;
            m_capacity = other.m_capacity;
            m_negative = other.m_negative;
            m_digits = other.m_digits;
            other.m_size = 1;
            other.m_capacity = 0;
            other.m_negative = false;
            other.m_digits = &g_zero;
        }
        return *this;
    }

    number::~number()
    {
        if (m_digits != &g_zero)
            g_deallocate(m_digits);
    }

    void number::from_int(limb_t const * digits, int size)
    {
        number result;
        
        result.m_size = size;
        result.m_capacity = size;
        result.m_negative = false;
        result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
        for (int i = 0; i < size; ++i)
            result.m_digits[i] = digits[i];
        *this = std::move(result);
    }

    bool number::from_base10(char const * str, int length)
    {
        constexpr int numerator = largest_base10_numerator_fitting_in_limb_size();
        number result;

        if (length > 0 && str[0] == '-')
        {
            result.m_negative = true;
            str++;
            length--;
        }
        if (length <= 0)
            return false;
        result.m_capacity = length / numerator;
        if ((length % numerator) > 0)
            result.m_capacity++;
        result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
        result.m_size = bigint::from_base10(result.m_digits, str, length);
        if (result.m_size > 0)
        {
            *this = std::move(result);
            return true;
        }
        return false;
    }

    bool number::from_base16(char const * str, int length)
    {
        number result;
        
        if (length > 0 && str[0] == '-')
        {
            result.m_negative = true;
            str++;
            length--;
        }
        if (length <= 0)
            return false;
        result.m_capacity = length / (sizeof(limb_t) * 2);
        if ((length % (sizeof(limb_t) * 2)) > 0)
            result.m_capacity++;
        result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
        result.m_size = bigint::from_base16(result.m_digits, str, length);
        if (result.m_size > 0)
        {
            *this = std::move(result);
            return true;
        }
        return false;
    }

    int number::to_base10(char * str, int capacity)
    {
        if (m_negative)
        {
            if (capacity > 0)
                str[0] = '-';
            return bigint::to_base10(str + 1, capacity - 1, m_digits, m_size) + 1;
        }
        return bigint::to_base10(str, capacity, m_digits, m_size);
    }

    int number::to_base16(char * str, int capacity)
    {
        if (m_negative)
        {
            if (capacity > 0)
                str[0] = '-';
            return bigint::to_base16(str + 1, capacity - 1, m_digits, m_size) + 1;
        }
        return bigint::to_base16(str, capacity, m_digits, m_size);
    }

    bool operator<(number const & lhs, number const & rhs)
    {
        return compare(lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size) < 0;
    }

    bool operator>(number const & lhs, number const & rhs)
    {
        return compare(lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size) > 0;
    }

    bool operator<=(number const & lhs, number const & rhs)
    {
        return compare(lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size) <= 0;
    }

    bool operator>=(number const & lhs, number const & rhs)
    {
        return compare(lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size) >= 0;
    }

    bool operator==(number const & lhs, number const & rhs)
    {
        return compare(lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size) == 0;
    }

    bool operator!=(number const & lhs, number const & rhs)
    {
        return compare(lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size) != 0;
    }

    number operator+(number const & lhs, number const & rhs)
    {
        number result;

        if (lhs.m_negative == rhs.m_negative)
        {
            result.m_capacity = std::max(lhs.m_size, rhs.m_size) + 1;
            result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
            result.m_negative = lhs.m_negative;
            result.m_size = add(result.m_digits, lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size);
        }
        else
        {
            switch (compare(lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size))
            {
            case -1: // rhs > lhs
                result.m_capacity = rhs.m_size;
                result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
                result.m_negative = rhs.m_negative;
                result.m_size = sub(result.m_digits, rhs.m_digits, lhs.m_digits, rhs.m_size, lhs.m_size);
                break;
            case 1: // lhs > rhs
                result.m_capacity = lhs.m_size;
                result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
                result.m_negative = lhs.m_negative;
                result.m_size = sub(result.m_digits, lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size);
                break;
            default:
                // result is already default initialized to zero
                break;
            }
        }
        return result;
    }

    number operator-(number const & lhs, number const & rhs)
    {
        number result;

        if (lhs.m_negative != rhs.m_negative)
        {
            result.m_capacity = std::max(lhs.m_size, rhs.m_size) + 1;
            result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
            result.m_negative = lhs.m_negative;
            result.m_size = add(result.m_digits, lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size);
        }
        else
        {
            switch (compare(lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size))
            {
            case -1: // rhs > lhs
                result.m_capacity = rhs.m_size;
                result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
                result.m_negative = !rhs.m_negative;
                result.m_size = sub(result.m_digits, rhs.m_digits, lhs.m_digits, rhs.m_size, lhs.m_size);
                break;
            case 1: // lhs > rhs
                result.m_capacity = lhs.m_size;
                result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
                result.m_negative = lhs.m_negative;
                result.m_size = sub(result.m_digits, lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size);
                break;
            default:
                // result is already default initialized to zero
                break;
            }
        }
        return result;
    }

    number operator*(number const & lhs, number const & rhs)
    {
        number result;

        result.m_capacity = lhs.m_size + rhs.m_size;
        result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
        result.m_negative = lhs.m_negative != rhs.m_negative;
        result.m_size = mul(result.m_digits, lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size);
        return result;
    }

    number operator/(number const & lhs, number const & rhs)
    {
        number result;

        result.m_capacity = lhs.m_size;
        result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
        result.m_negative = lhs.m_negative != rhs.m_negative;
        result.m_size = div(result.m_digits, lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size);
        return result;
    }

    number operator%(number const & lhs, number const & rhs)
    {
        number result;

        result.m_capacity = rhs.m_size;
        result.m_digits = static_cast<limb_t *>(g_allocate(result.m_capacity * sizeof(limb_t)));
        result.m_negative = lhs.m_negative;
        result.m_size = mod(result.m_digits, lhs.m_digits, rhs.m_digits, lhs.m_size, rhs.m_size);
        return result;
    }
}
