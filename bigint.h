#pragma once

#include <cstdint>

namespace bigint
{
    using limb_t = ::uint32_t;
    using limb_double_t = ::uint64_t;

    // passing nullptr will set the allocator to malloc/free (the default allocator)
    void set_allocator(void * (*alloc)(size_t), void (*dealloc)(void *));

    // returns the number of digits written, or -1 if an error occured
    // digits capcity must be big enough to store the result
    // slower than from_base16
    int from_base10(limb_t * digits, char const * str, int length);

    // returns the number of digits written, or -1 if an error occured
    // digits capcity must be big enough to store the result
    int from_base16(limb_t * digits, char const * str, int length);

    // create a string representation of digits in base 10, returns the numbers of characters written in str
    // slower than to_base16
    int to_base10(char * str, int capacity, limb_t const * digits, int size);

    // create a string representation of digits in base 16, returns the numbers of characters written in str
    int to_base16(char * str, int capacity, limb_t const * digits, int size);

    // returns -1 if lhs < rhs, 1 if rhs > lhs, or 0 if lhs == rhs
    int compare(limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size);

    // result capacity must be at least max(lhs_size, rhs_size) + 1
    // returns the size of the result
    int add(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size);

    // result capacity must be at least max(lhs_size, rhs_size)
    // lhs must be greater than rhs
    // returns the size of the result
    int sub(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size);

    // result capacity must be at least lhs_size + rhs_size
    // returns the size of the result
    int mul(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size);

    // result capacity must be at least lhs_size
    // rhs must be greater than zero
    // returns the size of the result
    int div(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size);

    // result capacity must be at least rhs_size
    // rhs must be greater than zero
    // returns the size of the result
    int mod(limb_t * result, limb_t const * lhs, limb_t const * rhs, int lhs_size, int rhs_size);

    class number
    {
    public:
        number();
        number(number const & other);
        number(number && other);
        number & operator=(number const & other);
        number & operator=(number && other);
        ~number();

        void from_int(limb_t const * digits, int size);
        bool from_base10(char const * str, int length);
        bool from_base16(char const * str, int length);
        int to_base10(char * str, int capacity);
        int to_base16(char * str, int capacity);

        friend bool operator<(number const & lhs, number const & rhs);
        friend bool operator>(number const & lhs, number const & rhs);
        friend bool operator<=(number const & lhs, number const & rhs);
        friend bool operator>=(number const & lhs, number const & rhs);
        friend bool operator==(number const & lhs, number const & rhs);
        friend bool operator!=(number const & lhs, number const & rhs);

        friend number operator+(number const & lhs, number const & rhs);
        friend number operator-(number const & lhs, number const & rhs);
        friend number operator*(number const & lhs, number const & rhs);
        friend number operator/(number const & lhs, number const & rhs);
        friend number operator%(number const & lhs, number const & rhs);

    private:
        int m_size;
        int m_capacity;
        bool m_negative;
        limb_t * m_digits;
    };
}
