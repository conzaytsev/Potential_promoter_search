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

from local_alignment import align
from multiprocessing import Pool, cpu_count
import json
import sys
from random import randint


def initializer(single_import_matrix, pair_import_matrix, g_start, g_continue, corridor_w):
    global gap_start
    gap_start = g_start
    global gap_continue
    gap_continue = g_continue
    global corridor_width
    corridor_width = corridor_w

    global single_matrix
    single_matrix = single_import_matrix
    global pair_matrix
    pair_matrix = pair_import_matrix


def mixer(line):
    mixed_line = ''
    while line != '':
        position = randint(0, len(line) - 1)
        mixed_line += line[position]
        line = line[:position] + line[position + 1:]
    return mixed_line


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
    best_result = (result[0])
    return best_result


def scan():
    matrix_number = int(sys.argv[1])  # 0
    sample_size = int(sys.argv[2])  # 100
    gap_start = int(sys.argv[3])  # -30
    gap_continue = int(sys.argv[4])  # -6
    corridor_width = int(sys.argv[5])  # 80

    cores = int(cpu_count() / 2)
    position = 0

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

    pool = Pool(cores, initializer, (single_import_matrix, pair_import_matrix, gap_start, gap_continue, corridor_width))

    file = open('bin/z/' + str(matrix_number) + '.z.txt', 'w')
    file.close()

    while position < sample_size:
        request = []
        for i in range(cores):
            if position < sample_size:
                random_number = randint(0, len(database) - 1)
                query = database[random_number][database[random_number].find('\n'):].replace('\n', '')
                request.append((mixer(query), position))
                position += 1

        # start = timeit.default_timer()

        try:
            request_result = pool.starmap(check, request)
        except Exception:
            continue

        # stop = timeit.default_timer()
        # print('Run time:', stop - start)

        for core_result in request_result:
            if core_result:
                file = open('bin/z/' + str(matrix_number) + '.z.txt', 'a')
                file.write(str(core_result) + '\n')
                file.close()


def chromosome_scan():
    matrix_number = int(sys.argv[1])  # 0
    sample_size = int(sys.argv[2])  # 100
    gap_start = int(sys.argv[3])  # -30
    gap_continue = int(sys.argv[4])  # -6
    corridor_width = int(sys.argv[5])  # 80
    z_filename = str(sys.argv[6])
    cores = int(sys.argv[7])
    step = 600

    position = 0

    file = open(z_filename, 'r')
    chromosome = ''
    for line in file:
        chromosome += line
    file.close()

    with open('single_matrix/' + str(matrix_number) + '.single_matrix.json', 'r') as f:
        single_import_matrix = json.load(f)
    with open('pair_matrix/' + str(matrix_number) + '.pair_matrix.json', 'r') as f:
        pair_import_matrix = json.load(f)

    pool = Pool(cores, initializer, (single_import_matrix, pair_import_matrix, gap_start, gap_continue, corridor_width))

    file = open('bin/z/' + str(matrix_number) + '.z.txt', 'w')
    file.close()

    while position < sample_size:
        request = []
        for i in range(cores):
            if position < sample_size:
                request.append((chromosome[position * step:(position + 1) * step + 1], position))
                position += 1
                # start = timeit.default_timer()

        try:
            request_result = pool.starmap(check, request)
        except Exception:
            continue

        # stop = timeit.default_timer()
        # print('Run time:', stop - start)

        for core_result in request_result:
            if core_result:
                file = open('bin/z/' + str(matrix_number) + '.z.txt', 'a')
                file.write(str(core_result) + '\n')
                file.close()


chromosome_scan()
