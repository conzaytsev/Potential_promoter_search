"""
Copyright (C) 2022 Konstantin Zaytsev

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from math import sqrt
from scipy.optimize import fsolve
import numpy as np
import json


def p(i, j, probability_of_bases):
    return probability_of_bases[i] * probability_of_bases[j] / 600


def R(matrix_f, bases):
    R_f = 0
    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                R_f += matrix_f[i][base1][base2] ** 2
    return R_f


def D(matrix_f, bases, probability_of_bases):
    D_f = 0
    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                D_f += matrix_f[i][base1][base2] * p(base1, base2, probability_of_bases)
    return D_f


def matrix_transformation(old_matrix, k1_f, k2_f, bases, probability_of_bases):
    new_matrix = [{'A': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'T': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                   'C': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'G': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}}]
    for i in range(599):
        new_matrix.append({'A': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'T': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                           'C': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                           'G': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}})
    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                new_matrix[i][base1][base2] = old_matrix[i][base1][base2] * k1_f + p(base1, base2,
                                                                                     probability_of_bases) * k2_f
    return new_matrix


def kd_transformation(old_matrix, k, bases, probability_of_bases):
    new_matrix = [{'A': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'T': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                   'C': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'G': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}}]
    for i in range(599):
        new_matrix.append({'A': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'T': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                           'C': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                           'G': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}})
    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                new_matrix[i][base1][base2] = old_matrix[i][base1][base2] + probability_of_bases[base1] * \
                                              probability_of_bases[base2] * k
    return new_matrix


def r2_transformation(old_matrix, r2, bases):
    old_r2 = R(old_matrix, bases)
    k = sqrt(r2 / old_r2)

    new_matrix = [{'A': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'T': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                   'C': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'G': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}}]
    for i in range(599):
        new_matrix.append({'A': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'T': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                           'C': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                           'G': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}})
    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                new_matrix[i][base1][base2] = old_matrix[i][base1][base2] * k
    return new_matrix


def create(matrix_number, R2, D2, number, start_coordinate, end_coordinate):
    filename = str(1000 + number) + '.txt'
    bases = ['A', 'T', 'C', 'G']
    probability_of_bases = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    base_pair_counter = [0 for i in range(600)]

    matrix = []
    for i in range(600):
        matrix.append({'A': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'T': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                       'C': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'G': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}})

    file = open('bin/' + filename, 'r')
    for line in file:
        if line[0:4] == 'CORM':
            if int(line[26]) == 1:
                base = 'A'
            elif int(line[26]) == 2:
                base = 'T'
            elif int(line[26]) == 3:
                base = 'C'
            else:
                base = 'G'

            matrix[int(line[11:14]) - 1][base]['A'] = int(line[31:34])
            matrix[int(line[11:14]) - 1][base]['T'] = int(line[36:39])
            matrix[int(line[11:14]) - 1][base]['C'] = int(line[41:44])
            matrix[int(line[11:14]) - 1][base]['G'] = int(line[46:49])

            probability_of_bases['A'] += int(line[31:34])
            probability_of_bases['T'] += int(line[36:39])
            probability_of_bases['C'] += int(line[41:44])
            probability_of_bases['G'] += int(line[46:49])

    file.close()

    sum_of_bases = 0
    for base in bases:
        sum_of_bases += probability_of_bases[base]
    for base in bases:
        probability_of_bases[base] /= sum_of_bases

    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                base_pair_counter[i] += matrix[i][base1][base2]

    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                matrix[i][base1][base2] = (matrix[i][base1][base2] - base_pair_counter[i] * probability_of_bases[
                    base1] * probability_of_bases[base2]) / sqrt(
                    base_pair_counter[i] * probability_of_bases[base1] * probability_of_bases[base2] * (
                            1 - probability_of_bases[base1] * probability_of_bases[base2]))

    R1 = R(matrix, bases)
    D1 = D(matrix, bases, probability_of_bases)
    print('pair matrix:')
    print(R1, D1)


    def coeff_funk(k):
        return [D(kd_transformation(matrix, k[0], bases, probability_of_bases), bases, probability_of_bases) - D2]

    root = fsolve(coeff_funk, [1])
    print(root[0])

    matrix = kd_transformation(matrix, root[0], bases, probability_of_bases)
    R1 = R(matrix, bases)
    D1 = D(matrix, bases, probability_of_bases)
    print(R1, D1)

    matrix = r2_transformation(matrix, R2, bases)
    R1 = R(matrix, bases)
    D1 = D(matrix, bases, probability_of_bases)
    print(R1, D1)

    for i in range(600):
        if i < start_coordinate or i > end_coordinate:
            for base1 in bases:
                for base2 in bases:
                    matrix[i][base2][base1] = -1

    np.save('pair_matrix/' + str(matrix_number) + '.pair_matrix.npy', matrix)

    single_matrix = [{'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0} for i in range(600)]
    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                single_matrix[i][base1] += matrix[i][base2][base1] / 4

    for i in range(600):
        if i < start_coordinate or i > end_coordinate:
            for base1 in bases:
                single_matrix[i][base1] = -1

    np.save('single_matrix/' + str(matrix_number) + '.single_matrix.npy', single_matrix)

    with open('single_matrix/' + str(matrix_number) + '.single_matrix.json', 'w') as f:
        # indent=2 is not needed but makes the file human-readable
        json.dump(single_matrix, f)

    with open('pair_matrix/' + str(matrix_number) + '.pair_matrix.json', 'w') as f:
        # indent=2 is not needed but makes the file human-readable
        json.dump(matrix, f)

    return 0


def old_create(matrix_number, R2, D2):
    filename = '600.100'
    bases = ['A', 'T', 'C', 'G']
    probability_of_bases = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    base_pair_counter = [0 for i in range(600)]

    matrix = []
    for i in range(600):
        matrix.append({'A': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'T': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                       'C': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'G': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}})

    file = open('bin/' + filename, 'r')
    for line in file:
        if line[0:4] == 'CORM':
            if int(line[26]) == 1:
                base = 'A'
            elif int(line[26]) == 2:
                base = 'T'
            elif int(line[26]) == 3:
                base = 'C'
            else:
                base = 'G'

            matrix[int(line[11:14]) - 1][base]['A'] = int(line[31:34])
            matrix[int(line[11:14]) - 1][base]['T'] = int(line[36:39])
            matrix[int(line[11:14]) - 1][base]['C'] = int(line[41:44])
            matrix[int(line[11:14]) - 1][base]['G'] = int(line[46:49])

            probability_of_bases['A'] += int(line[31:34])
            probability_of_bases['T'] += int(line[36:39])
            probability_of_bases['C'] += int(line[41:44])
            probability_of_bases['G'] += int(line[46:49])

    file.close()

    sum_of_bases = 0
    for base in bases:
        sum_of_bases += probability_of_bases[base]
    for base in bases:
        probability_of_bases[base] /= sum_of_bases

    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                base_pair_counter[i] += matrix[i][base1][base2]

    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                matrix[i][base1][base2] = (matrix[i][base1][base2] - base_pair_counter[i] * probability_of_bases[
                    base1] * probability_of_bases[base2]) / sqrt(
                    base_pair_counter[i] * probability_of_bases[base1] * probability_of_bases[base2] * (
                            1 - probability_of_bases[base1] * probability_of_bases[base2]))

    R1 = R(matrix, bases)
    D1 = D(matrix, bases, probability_of_bases)
    print('pair matrix:')
    print(R1, D1)

    def coeff_funk(k):
        return [D(kd_transformation(matrix, k[0], bases, probability_of_bases), bases, probability_of_bases) - D2]

    root = fsolve(coeff_funk, [1])
    print(root[0])

    matrix = kd_transformation(matrix, root[0], bases, probability_of_bases)
    R1 = R(matrix, bases)
    D1 = D(matrix, bases, probability_of_bases)
    print(R1, D1)

    matrix = r2_transformation(matrix, R2, bases)
    R1 = R(matrix, bases)
    D1 = D(matrix, bases, probability_of_bases)
    print(R1, D1)

    np.save('pair_matrix/' + str(matrix_number) + '.pair_matrix.npy', matrix)

    single_matrix = [{'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0} for i in range(600)]
    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                single_matrix[i][base1] += matrix[i][base2][base1] / 4

    np.save('single_matrix/' + str(matrix_number) + '.single_matrix.npy', single_matrix)

    with open('single_matrix/' + str(matrix_number) + '.single_matrix.json', 'w') as f:
        # indent=2 is not needed but makes the file human-readable
        json.dump(single_matrix, f)

    with open('pair_matrix/' + str(matrix_number) + '.pair_matrix.json', 'w') as f:
        # indent=2 is not needed but makes the file human-readable
        json.dump(matrix, f)

    return 0


def create_by_name(matrix_number, R2, D2, filename):
    bases = ['A', 'T', 'C', 'G']
    probability_of_bases = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    base_pair_counter = [0 for i in range(600)]

    matrix = []
    for i in range(600):
        matrix.append({'A': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'T': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0},
                       'C': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}, 'G': {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}})

    file = open('bin/' + filename, 'r')
    for line in file:
        if line[0:4] == 'CORM':
            if int(line[26]) == 1:
                base = 'A'
            elif int(line[26]) == 2:
                base = 'T'
            elif int(line[26]) == 3:
                base = 'C'
            else:
                base = 'G'

            matrix[int(line[11:14]) - 1][base]['A'] = int(line[31:34])
            matrix[int(line[11:14]) - 1][base]['T'] = int(line[36:39])
            matrix[int(line[11:14]) - 1][base]['C'] = int(line[41:44])
            matrix[int(line[11:14]) - 1][base]['G'] = int(line[46:49])

            probability_of_bases['A'] += int(line[31:34])
            probability_of_bases['T'] += int(line[36:39])
            probability_of_bases['C'] += int(line[41:44])
            probability_of_bases['G'] += int(line[46:49])

    file.close()

    sum_of_bases = 0
    for base in bases:
        sum_of_bases += probability_of_bases[base]
    for base in bases:
        probability_of_bases[base] /= sum_of_bases

    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                base_pair_counter[i] += matrix[i][base1][base2]

    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                matrix[i][base1][base2] = (matrix[i][base1][base2] - base_pair_counter[i] * probability_of_bases[
                    base1] * probability_of_bases[base2]) / sqrt(
                    base_pair_counter[i] * probability_of_bases[base1] * probability_of_bases[base2] * (
                            1 - probability_of_bases[base1] * probability_of_bases[base2]))

    R1 = R(matrix, bases)
    D1 = D(matrix, bases, probability_of_bases)
    print('pair matrix:')
    print(R1, D1)

    def coeff_funk(k):
        return [D(kd_transformation(matrix, k[0], bases, probability_of_bases), bases, probability_of_bases) - D2]

    root = fsolve(coeff_funk, [1])
    print(root[0])

    matrix = kd_transformation(matrix, root[0], bases, probability_of_bases)
    R1 = R(matrix, bases)
    D1 = D(matrix, bases, probability_of_bases)
    print(R1, D1)

    matrix = r2_transformation(matrix, R2, bases)
    R1 = R(matrix, bases)
    D1 = D(matrix, bases, probability_of_bases)
    print(R1, D1)

    for i in range(600):
        if i > 520:
            for base1 in bases:
                for base2 in bases:
                    matrix[i][base2][base1] = -1

    np.save('pair_matrix/' + str(matrix_number) + '.pair_matrix.npy', matrix)

    single_matrix = [{'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0} for i in range(600)]
    for i in range(600):
        for base1 in bases:
            for base2 in bases:
                single_matrix[i][base1] += matrix[i][base2][base1] / 4

    for i in range(600):
        if i > 520:
            for base1 in bases:
                single_matrix[i][base1] = -1

    np.save('single_matrix/' + str(matrix_number) + '.single_matrix.npy', single_matrix)

    with open('single_matrix/' + str(matrix_number) + '.single_matrix.json', 'w') as f:
        # indent=2 is not needed but makes the file human-readable
        json.dump(single_matrix, f)

    with open('pair_matrix/' + str(matrix_number) + '.pair_matrix.json', 'w') as f:
        # indent=2 is not needed but makes the file human-readable
        json.dump(matrix, f)

    return 0
