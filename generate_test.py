import random
from math import ceil


FILE_TEMPLATE_HEADER = """#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "bigint.h"

struct bigint_data
{
    bigint::limb_t * digits;
    int size;
};

static bigint_data const * g_inputs;
static bigint_data const * g_expected_results;

void test_op(int (*op)(bigint::limb_t *, bigint::limb_t const *, bigint::limb_t const *, int, int), bigint_data const & lhs, bigint_data const & rhs, bigint_data const & expected_result)
{
    static bigint::limb_t result[0x1000];
    int const size = op(result, lhs.digits, rhs.digits, lhs.size, rhs.size);
    assert(bigint::compare(result, expected_result.digits, size, expected_result.size) == 0);
}

"""

FILE_TEMPLATE_FOOTER = """
int main()
{
    generate_test_numbers();
    test_add();
    test_sub();
    test_mul();
    test_div();
    test_mod();
    free_test_numbers();
    return 0;
}
"""

def gen_numbers():
    return ['0', '1', 'f', 'ff', 'ffffffffffffffff', '10000000000000000'] + [''.join(random.choices('0123456789abcdef', k=random.randrange(1000, 8000))) for _ in range(50)]

numbers_outputs_map = {}
numbers_outputs_array = []

def gen_test_op(numbers, name, cond, op):
    n = 0
    lines = []
    lines.append(f'void test_{name}()')
    lines.append('{')
    lines.append('    auto const chrono_start = std::chrono::high_resolution_clock::now();')
    for i, lhs in enumerate(numbers):
        for j, rhs in enumerate(numbers):
            if not cond(int(lhs, 16), int(rhs, 16)):
                continue
            expected_result = hex(op(int(lhs, 16), int(rhs, 16)))[2:]
            if expected_result in numbers_outputs_map:
                k = numbers_outputs_map[expected_result]
            else:
                k = len(numbers_outputs_array)
                numbers_outputs_map[expected_result] = k
                numbers_outputs_array.append(expected_result)
            lines.append(f'    test_op(bigint::{name}, g_inputs[{i}], g_inputs[{j}], g_expected_results[{k}]);')
            n += 1
    lines.append('    auto const chrono_end = std::chrono::high_resolution_clock::now();')
    lines.append('    auto const time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end - chrono_start).count();')
    lines.append(f'    std::cout << "{name}: " << (time_ns / {n}) << "ns/op" << std::endl;')
    lines.append('}')
    lines.append('')
    return '\n'.join(lines)

numbers = gen_numbers()

test_add = gen_test_op(numbers, 'add', lambda lhs, rhs: True, lambda lhs, rhs: lhs + rhs)
test_sub = gen_test_op(numbers, 'sub', lambda lhs, rhs: lhs >= rhs, lambda lhs, rhs: lhs - rhs)
test_mul = gen_test_op(numbers, 'mul', lambda lhs, rhs: True, lambda lhs, rhs: lhs * rhs)
test_div = gen_test_op(numbers, 'div', lambda lhs, rhs: lhs >= rhs and rhs != 0, lambda lhs, rhs: lhs // rhs)
test_mod = gen_test_op(numbers, 'mod', lambda lhs, rhs: lhs >= rhs and rhs != 0, lambda lhs, rhs: lhs % rhs)

test_numbers_arrays_gen = []
test_numbers_arrays_gen.append('void generate_test_numbers()')
test_numbers_arrays_gen.append('{')
test_numbers_arrays_gen.append(f'    bigint_data * const inputs = new bigint_data[{len(numbers)}];')
for i, n in enumerate(numbers):
    test_numbers_arrays_gen.append(f'    inputs[{i}].digits = new bigint::limb_t[{max(1, int(ceil(len(n) / 8)))}];')
    test_numbers_arrays_gen.append(f'    inputs[{i}].size = bigint::from_base16(inputs[{i}].digits, "{n}", {len(n)});')
test_numbers_arrays_gen.append('    g_inputs = inputs;')
test_numbers_arrays_gen.append('')
test_numbers_arrays_gen.append(f'    bigint_data * const expected_results = new bigint_data[{len(numbers_outputs_array)}];')
for i, n in enumerate(numbers_outputs_array):
    test_numbers_arrays_gen.append(f'    expected_results[{i}].digits = new bigint::limb_t[{max(1, int(ceil(len(n) / 8)))}];')
    test_numbers_arrays_gen.append(f'    expected_results[{i}].size = bigint::from_base16(expected_results[{i}].digits, "{n}", {len(n)});')
test_numbers_arrays_gen.append('    g_expected_results = expected_results;')
test_numbers_arrays_gen.append('}')
test_numbers_arrays_gen.append('')
test_numbers_arrays_gen.append('void free_test_numbers()')
test_numbers_arrays_gen.append('{')
test_numbers_arrays_gen.append(f'    for (int i = 0; i < {len(numbers)}; ++i)')
test_numbers_arrays_gen.append('        delete[] g_inputs[i].digits;')
test_numbers_arrays_gen.append(f'    for (int i = 0; i < {len(numbers_outputs_array)}; ++i)')
test_numbers_arrays_gen.append('        delete[] g_expected_results[i].digits;')
test_numbers_arrays_gen.append('}')
test_numbers_arrays_gen.append('')
test_numbers_arrays_gen = '\n'.join(test_numbers_arrays_gen)

with open('test.cpp', 'w+') as f:
    f.write(FILE_TEMPLATE_HEADER)
    f.write(test_numbers_arrays_gen)
    f.write(test_add)
    f.write(test_sub)
    f.write(test_mul)
    f.write(test_div)
    f.write(test_mod)
    f.write(FILE_TEMPLATE_FOOTER)
