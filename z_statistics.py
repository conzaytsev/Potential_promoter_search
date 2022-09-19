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
import re
import os
import subprocess


def initializer(single_import_matrixes, pair_import_matrixes, g_start, g_continue, corridor_w, query_l, mtrxs):
    global gap_start
    gap_start = g_start
    global gap_continue
    gap_continue = g_continue
    global corridor_width
    corridor_width = corridor_w
    global query_length
    query_length = query_l
    global matrix_n
    matrix_n = mtrxs
    global single_matrix
    single_matrix = single_import_matrixes
    global pair_matrix
    pair_matrix = pair_import_matrixes


def replace(query, old, new):
    for i in range(query.count(old)):
        index = query.find(old)
        query = query[:index] + new[randint(0, len(new) - 1)] + query[index + 1:]
    return query


def check(query, position):
    results = [0 for i in range(matrix_n)]
    # standard order, base spiral
    for number in range(matrix_n):
        result = align(query, single_matrix[number], pair_matrix[number], gap_start, gap_continue, corridor_width)
        results[number] = result[0]
        # spiral aligned, score, matrix number, promoter start position, promoter end position, aligned sequence, aligned matrix

    return results


def calculate(filename, gap_start, gap_continue, corridor_width, query_length, step, max_length, cores, matrices):
    '''filename = 'chromosome.mix.' + str(sys.argv[1])
    gap_start = -25  # -59
    gap_continue = -6  # -22
    corridor_width = 80
    query_length = 700  # 700
    step = int(corridor_width * query_length / 600) - 20  # 60
    max_length = 1000000
    cores = int(cpu_count() / 2)
    matrices = 1  # 25'''

    min_matrix_number = 0
    position = 0

    if not os.path.isdir('z/' + str(re.findall(r'\d+', filename)[0])):
        subprocess.call('mkdir -p -- z/' + str(re.findall(r'\d+', filename)[0]), shell=True, stdout=subprocess.DEVNULL)

    for number in range(matrices):
        file = open('z/' + str(re.findall(r'\d+', filename)[0]) + '/z.' + str(number + min_matrix_number), 'w')
        file.close()

    chromosome = ''
    file = open('edited_chromosomes/' + filename, 'r')
    for line in file:
        chromosome += line
    file.close()

    single_import_matrixes = []
    pair_import_matrixes = []

    for number in range(matrices):
        with open('single_matrix/' + str(number + min_matrix_number) + '.single_matrix.json', 'r') as f:
            single_import_matrixes.append(json.load(f))
        with open('pair_matrix/' + str(number + min_matrix_number) + '.pair_matrix.json', 'r') as f:
            pair_import_matrixes.append(json.load(f))

    pool = Pool(cores, initializer, (single_import_matrixes, pair_import_matrixes,
                                     gap_start, gap_continue, corridor_width, query_length, matrices))

    while position < min((len(chromosome) - query_length), max_length):
        request = []
        for i in range(cores):
            if position < (len(chromosome) - query_length):
                query = chromosome[position: position + query_length]
                request.append((query, position))  # position + begin if chromosome splits
                position += step

        results = pool.starmap(check, request)
        print(results)

        for number in range(matrices):
            file = open('z/' + str(re.findall(r'\d+', filename)[0]) + '/z.' + str(number + min_matrix_number), 'a')
            for line in results:
                file.write(str(line[number]) + '\n')
            file.close()


def main():
    from config import min_chromosome, max_chromosome, gap_start_fine, gap_continue_fine, chromosome_corridor_width, query_length, z_statistics_sample_size, cores, classes
    step = int(chromosome_corridor_width * query_length / 600) - 20
    max_length = z_statistics_sample_size * step + query_length
    for number in range(min_chromosome, max_chromosome + 1):
        filename = 'chromosome.mix.' + str(number)
        calculate(filename, gap_start_fine, gap_continue_fine, chromosome_corridor_width, query_length, step, max_length, cores, classes)


main()
