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


def initializer(single_import_matrix, pair_import_matrix, g_start, g_continue, corridor_w, thrshld, mn, vr, sc, ec, ml, il):
    global gap_start
    gap_start = g_start
    global gap_continue
    gap_continue = g_continue
    global corridor_width
    corridor_width = corridor_w
    global threshold
    threshold = thrshld

    global single_matrix
    single_matrix = single_import_matrix
    global pair_matrix
    pair_matrix = pair_import_matrix

    global z_mean
    z_mean = mn
    global z_variance
    z_variance = vr

    global start
    start = sc
    global end
    end = ec
    global minimum_length
    minimum_length = ml
    global intersection
    intersection = il


def replace(query, old, new):
    for i in range(query.count(old)):
        index = query.find(old)
        query = query[:index] + new[randint(0, len(new) - 1)] + query[index + 1:]
    return query


def check(query, position):
    best_result = ()
    query_length = len(query)
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

    # standard order, base spiral
    result = align(query, single_matrix, pair_matrix, gap_start, gap_continue, corridor_width)
    z = (result[0] - z_mean) / z_variance
    if z > threshold:
        length = result[4] - result[3]
        if length > minimum_length:
            if length * (1 - intersection) + result[3] > start and end + length * (1 - intersection) > result[4]:
                best_result = (position, result[0], z, result[3], result[4], result[1], result[2])
    # position, score, z-value, promoter start position, promoter end position, aligned sequence, aligned matrix
    return best_result


def scan():
    matrix_number = int(sys.argv[1])  # 0
    threshold = int(sys.argv[2])  # 5  # 0 z-value
    gap_start = int(sys.argv[3])  # -30
    gap_continue = int(sys.argv[4])  # -6
    corridor_width = int(sys.argv[5])  # 80
    start_coordinate = int(sys.argv[6])
    end_coordinate = int(sys.argv[7])
    min_length = int(sys.argv[8])
    intersection_level = float(sys.argv[9])
    cores = int(sys.argv[10])
    
    position = 0

    scores = []
    file = open('bin/z/' + str(matrix_number) + '.z.txt', 'r')
    for line in file:
        scores.append(float(line[:-1]))
    file.close()

    mean = statistics.mean(scores)
    vari = sqrt(statistics.variance(scores))


    database = []
    file = open('bin/database_work_copy.fasta', 'r')
    for line in file:
        if line[0] == '>':
            database.append(line)
        else:
            database[-1] = database[-1] + line
    file.close()


    with open('single_matrix/' + str(matrix_number) + '.single_matrix.json', 'r') as f:
        single_import_matrix = json.load(f)
    with open('pair_matrix/' + str(matrix_number) + '.pair_matrix.json', 'r') as f:
        pair_import_matrix = json.load(f)


    pool = Pool(cores, initializer, (single_import_matrix, pair_import_matrix, gap_start, gap_continue, corridor_width,
                                     threshold, mean, vari, start_coordinate, end_coordinate, min_length,
                                     intersection_level))

    file = open('classes/' + str(matrix_number) + '.scan_result.txt', 'w')
    file.close()
    file = open('classes/' + str(matrix_number) + '.scan_result_positions.txt', 'w')
    file.close()

    while position < len(database):
        request = []
        for i in range(cores):
            if position < len(database):
                query = database[position][database[position].find('\n'):].replace('\n', '')
                request.append((query, position))
                position += 1

        #start = timeit.default_timer()

        try:
            request_result = pool.starmap(check, request)
        except Exception:
            continue

        #stop = timeit.default_timer()
        #print('Run time:', stop - start)

        for core_result in request_result:
            if core_result:
                file = open('classes/' + str(matrix_number) + '.scan_result.txt', 'a')
                file.write(database[core_result[0]][:9] + '\n' + str(core_result[1:]) + '\n')
                file.close()

                file = open('classes/' + str(matrix_number) + '.scan_result_positions.txt', 'a')
                file.write(str(core_result[0]) + '\n')
                file.close()
                

scan()
