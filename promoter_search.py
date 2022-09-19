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

from random import randint
from local_alignment import align
from multiprocessing import Pool
import json
import sys
import statistics
from math import sqrt
import os
import subprocess
import re


def initializer(single_import_matrix, pair_import_matrix, g_start, g_continue, corridor_w, thrshld, query_l, flg,
                mn, var, mtrx):
    global gap_start
    gap_start = g_start
    global gap_continue
    gap_continue = g_continue
    global corridor_width
    corridor_width = corridor_w
    global threshold
    threshold = thrshld
    global query_length
    query_length = query_l

    global single_matrix
    single_matrix = single_import_matrix
    global pair_matrix
    pair_matrix = pair_import_matrix

    global flag
    flag = flg

    global z_mean
    z_mean = mn
    global z_variance
    z_variance = var

    global matrix_number
    matrix_number = mtrx


def replace(query, old, new):
    for i in range(query.count(old)):
        index = query.find(old)
        query = query[:index] + new[randint(0, len(new) - 1)] + query[index + 1:]
    return query


def check(query, position):
    best_result = ()
    if query_length - query.count('A') - query.count('T') - query.count('G') - query.count('C') > 100:
        return best_result

    if query_length - query.count('A') - query.count('T') - query.count('G') - query.count('C') > 0:
        query = replace(query, 'N', ['A', 'T', 'G', 'C'])

    if query_length - query.count('A') - query.count('T') - query.count('G') - query.count('C') > 0:
        query = replace(query, 'U', ['T'])
        query = replace(query, 'R', ['A', 'G'])
        query = replace(query, 'Y', ['T', 'C'])
        query = replace(query, 'K', ['T', 'G'])
        query = replace(query, 'M', ['A', 'C'])
        query = replace(query, 'S', ['G', 'C'])
        query = replace(query, 'W', ['A', 'T'])
        query = replace(query, 'B', ['T', 'G', 'C'])
        query = replace(query, 'D', ['A', 'T', 'G'])
        query = replace(query, 'H', ['A', 'T', 'C'])
        query = replace(query, 'V', ['A', 'G', 'C'])

    if query_length - query.count('A') - query.count('T') - query.count('G') - query.count('C') > 0:
        return best_result
    if flag == 1:
        # standard order, base spiral
        result = align(query, single_matrix, pair_matrix, gap_start, gap_continue, corridor_width)
        z = (result[0] - z_mean) / z_variance
        if z > threshold:
            if best_result and z > best_result[2] or not best_result:
                best_result = (
                    1, result[0], z, matrix_number, result[3] + position, result[4] + position, result[1], result[2])
        # spiral aligned, score, z-value, matrix number, promoter start position, promoter end position, aligned sequence, aligned matrix

    elif flag == 2:
        reversed_query = query[::-1]
        print(query == reversed_query[::-1])
        # reversed order, base spiral
        result = align(reversed_query, single_matrix, pair_matrix, gap_start, gap_continue,
                       corridor_width)
        z = (result[0] - z_mean) / z_variance
        if z > threshold:
            if best_result and z > best_result[2] or not best_result:
                best_result = (2, result[0], z, matrix_number, position + query_length - 1 - result[4],
                               position + query_length - 1 - result[3], result[1], result[2])

    elif flag == 3:
        complementary_query = query.replace('A', '1').replace('T', '2').replace('G', '3').replace('C', '4')
        complementary_query = complementary_query.replace('1', 'T').replace('2', 'A').replace('3', 'C').replace('4',
                                                                                                                'G')
        # standard order, complementary spiral
        result = align(complementary_query, single_matrix, pair_matrix, gap_start, gap_continue,
                       corridor_width)
        z = (result[0] - z_mean) / z_variance
        if z > threshold:
            if best_result and z > best_result[2] or not best_result:
                best_result = (
                    3, result[0], z, matrix_number, result[3] + position, result[4] + position, result[1], result[2])

    else:
        complementary_query = query.replace('A', '1').replace('T', '2').replace('G', '3').replace('C', '4')
        complementary_query = complementary_query.replace('1', 'T').replace('2', 'A').replace('3', 'C').replace('4',
                                                                                                                'G')
        reversed_complementary_query = complementary_query[::-1]
        # reversed order, complementary spiral
        result = align(reversed_complementary_query, single_matrix, pair_matrix, gap_start,
                       gap_continue, corridor_width)
        z = (result[0] - z_mean) / z_variance
        if z > threshold:
            if best_result and z > best_result[2] or not best_result:
                best_result = (4, result[0], z, matrix_number, position + query_length - 1 - result[4],
                               position + query_length - 1 - result[3], result[1], result[2])
    return best_result


def main():
    from config import gap_start_fine, gap_continue_fine, chromosome_corridor_width, query_length, z_threshold, cores
    file_name = 'chromosome.' + str(sys.argv[1])
    node = int(sys.argv[2])
    matrix = int(sys.argv[3])

    step = int(chromosome_corridor_width * query_length / 600) - 20
    position = 0

    if not os.path.isdir('results/' + str(re.findall(r'\d+', file_name)[0]) + '/' + str(matrix)):
        subprocess.call('mkdir -p -- results/' + str(re.findall(r'\d+', file_name)[0]) + '/' + str(matrix),
                        shell=True, stdout=subprocess.DEVNULL)

    file = open('results/' + str(re.findall(r'\d+', file_name)[0]) + '/' + str(matrix) + '/results' + str(node) + '.' + file_name + '.' + str(matrix), 'w')
    file.close()

    scores = []
    file = open('z/' + str(re.findall(r'\d+', file_name)[0]) + '/z.' + str(matrix), 'r')
    for line in file:
        scores.append(float(line[:-1]))
    file.close()

    mean = statistics.mean(scores)
    vari = sqrt(statistics.variance(scores))

    file = open('edited_chromosomes/' + file_name, 'r')
    chromosome = ''
    for line in file:
        chromosome += line
    file.close()

    with open('single_matrix/' + str(matrix) + '.single_matrix.json', 'r') as f:
        single_import_matrix = json.load(f)
    with open('pair_matrix/' + str(matrix) + '.pair_matrix.json', 'r') as f:
        pair_import_matrix = json.load(f)

    pool = Pool(cores, initializer, (single_import_matrix, pair_import_matrix, gap_start_fine, gap_continue_fine,
                                     chromosome_corridor_width, z_threshold, query_length, node, mean, vari, matrix))

    while position < (len(chromosome) - query_length):
        request = []
        for i in range(cores):
            if position < (len(chromosome) - query_length):
                query = chromosome[position: position + query_length]
                request.append((query, position))  # position + begin if chromosome splits
                position += step

        try:
            request_result = pool.starmap(check, request)
        except Exception:
            continue

        for core_result in request_result:
            if core_result:
                file = open('results/' + str(re.findall(r'\d+', file_name)[0]) + '/' + str(matrix) + '/results' + str(
                    node) + '.' + file_name + '.' + str(matrix), 'a')
                file.write(str(core_result) + '\n')
                file.close()


main()
