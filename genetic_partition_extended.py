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

import json
import database_sampler
import matrix_creator
from random import randint
import subprocess
import matrix_exporter
import os
import glob


def work_copy_rewrite(filename):
    with open('bin/' + filename, 'r') as file1:
        with open('bin/database_work_copy.fasta', 'w') as file2:
            for line in file1:
                file2.write(line)


def repeats_remove(child):
    database = []
    database_index = []
    file = open('bin/database_work_copy.fasta', 'r')
    for line in file:
        if line[0] == '>':
            database.append('')
            database_index.append(line[:9])
        else:
            database[-1] = database[-1] + line
    file.close()

    positions = []
    repeated_positions = []
    for position in child:
        if position not in positions:
            positions.append(position)
        else:
            repeated_positions.append(position)

    for position in repeated_positions:
        number = child.index(position)
        new_position = randint(0, len(database) - 1)
        while new_position in positions:
            new_position = randint(0, len(database) - 1)
        child[number] = new_position

    return child


def crossingover(set_scores, position_sets):
    set_size = len(position_sets[0])

    min_score = min(set_scores)
    max_score = max(set_scores)

    adjusted_scores = []
    for score in set_scores:
        adjusted_scores.append(score + max_score - 2 * min_score)

    sum_score = sum(adjusted_scores)

    rnd_number = randint(0, sum_score - 1)
    parent_1 = 0
    for i in range(len(adjusted_scores)):
        rnd_number -= adjusted_scores[i]
        if rnd_number <= 0:
            parent_1 = i
            break

    parent_2 = parent_1
    while parent_2 == parent_1:
        rnd_number = randint(0, sum_score - 1)
        for i in range(len(adjusted_scores)):
            rnd_number -= adjusted_scores[i]
            if rnd_number <= 0:
                parent_2 = i
                break

    point_1 = randint(0, set_size - 1)
    point_2 = randint(0, set_size - 1)
    if point_1 > point_2:
        point = point_1
        point_1 = point_2
        point_2 = point

    child_1 = position_sets[parent_1][0: point_1]
    child_1.extend(position_sets[parent_2][point_1: point_2])
    child_1.extend(position_sets[parent_1][point_2:])

    child_2 = position_sets[parent_2][0: point_1]
    child_2.extend(position_sets[parent_1][point_1: point_2])
    child_2.extend(position_sets[parent_2][point_2:])

    return [child_1, child_2]


def mutation(child, mut_per, mut_bal):
    database = []
    database_index = []
    file = open('bin/database_work_copy.fasta', 'r')
    for line in file:
        if line[0] == '>':
            database.append('')
            database_index.append(line[:9])
        else:
            database[-1] = database[-1] + line
    file.close()

    mut_n = int(len(child) * mut_per / 100)
    for i in range(mut_n):
        mutation_type = randint(0, 100)
        if mutation_type < mut_bal:
            number_1 = randint(0, len(child) - 1)
            number_2 = randint(0, len(database) - 1)
            while number_2 in child:
                number_2 = randint(0, len(database) - 1)
            child[number_1] = number_2
        else:
            number_1 = randint(0, len(child) - 1)
            number_2 = randint(0, len(child) - 1)
            position = child[number_1]
            child[number_1] = child[number_2]
            child[number_2] = position
    return child


def scan(position_set, number, class_number, R2, Kd, z_statistics_sample_size, gap_start_fine,
         gap_continue_fine, corridor_width, z_threshold, sample, start_coordinate, end_coordinate, extension,
         min_length, intersection_level, z_filename, cores, login, password):
    matrix_creator.web_create(position_set, number, login, password)
    matrix_exporter.create(class_number, R2, Kd, number, start_coordinate, end_coordinate)

    if sample != 'learn_sample.fasta':
        work_copy_rewrite(sample)

    if sample == 'extended_full_sample.fasta':
        start_coordinate += extension
        end_coordinate += extension

    return_code = subprocess.call(
        ['pypy3', 'class_z_statistics.py', str(class_number), str(z_statistics_sample_size),
         str(gap_start_fine), str(gap_continue_fine), str(corridor_width), z_filename, str(cores)], shell=False,
        stdout=subprocess.DEVNULL)
    return_code += subprocess.call(
        ['pypy3', 'database_scanner.py', str(class_number), str(z_threshold),
         str(gap_start_fine),
         str(gap_continue_fine), str(corridor_width), str(start_coordinate),
         str(end_coordinate), str(min_length), str(intersection_level), str(cores)], shell=False,
        stdout=subprocess.DEVNULL)

    print('z + scan exit code ', return_code)

    if sample != 'learn_sample.fasta':
        work_copy_rewrite('learn_sample.fasta')

    if sample == 'extended_full_sample.fasta':
        start_coordinate -= extension
        end_coordinate -= extension

    result_indexes = []
    file = open('classes/' + str(class_number) + '.scan_result_positions.txt', 'r')
    for line in file:
        if line != '':
            result_indexes.append(int(line))
    file.close()
    return len(result_indexes)


def main():
    from config import promoter_database_filename as database_filename, classifier_mixed_chromosome, n_sets, set_size, \
        n_children, mut_per, mut_bal, \
        z_threshold, z_statistics_sample_size, gap_start_fine, gap_continue_fine, \
        database_corridor_width, extended_corridor_width, Kd, R2, learn_partition, min_set, iteration_limit, \
        start_coordinate, \
        end_coordinate, extension, min_length, intersection_level, cores, mahds_login as login, \
        mahds_password as password

    extended_database_filename = 'extended_sorted_' + database_filename
    z_filename = 'edited_chromosomes/chromosome.mix.' + str(classifier_mixed_chromosome)
    number = 0
    restart = 1
    class_number = 0
    best_score = 0
    best_check_score = 0

    if restart == 1:
        file = open('journal.txt', 'w')
        file.close()
        file = open('overfit.txt', 'w')
        file.close()

        for directory in ['classes', 'pair_matrix', 'single_matrix']:
            files = glob.glob(directory + '/*')
            for f in files:
                os.remove(f)

        with open(database_filename, 'r') as file1:
            with open('bin/database_work_copy.fasta', 'w') as file2:
                for line in file1:
                    file2.write(line)

        with open(extended_database_filename, 'r') as file1:
            with open('bin/extended_full_sample.fasta', 'w') as file2:
                for line in file1:
                    file2.write(line)

    else:
        names = [f[-6:-5] for f in os.listdir('classes') if f[-4:] == 'json']
        class_number = int(names[-1])
        work_copy_rewrite('full_sample.fasta')

    database_len = 0
    file = open('bin/database_work_copy.fasta', 'r')
    for line in file:
        if line[0] == '>':
            database_len += 1

    while database_len >= min_set:
        if restart == 1:
            if database_len * learn_partition < min_set:
                learn_partition = min_set / database_len

            database_sampler.database_divider(learn_partition)
            work_copy_rewrite('learn_sample.fasta')

            with open('journal.txt', 'a') as file:
                file.write('\n\n\n' + str(class_number) + '\n')
            with open('overfit.txt', 'a') as file:
                file.write('\n\n\n' + str(class_number) + '\n')

            position_sets = []
            for i in range(n_sets):
                position_sets.append(database_sampler.create_index(set_size))

            set_scores = []
            for position_set in position_sets:
                set_scores.append(
                    scan(position_set, number, class_number, R2, Kd, z_statistics_sample_size, gap_start_fine,
                         gap_continue_fine, database_corridor_width, z_threshold, 'learn_sample.fasta',
                         start_coordinate,
                         end_coordinate, extension, min_length, intersection_level, z_filename, cores, login, password))

                with open('journal.txt', 'a') as file:
                    file.write(str(set_scores) + '\n')

            set_scores, position_sets = (list(t) for t in zip(*sorted(zip(set_scores, position_sets), reverse=True)))

        else:
            restart = 1
            work_copy_rewrite('learn_sample.fasta')

            with open('classes/class' + str(class_number) + '.json', 'r') as file:
                position_sets = json.load(file)

            set_scores = []
            for position_set in position_sets:
                set_scores.append(
                    scan(position_set, number, class_number, R2, Kd, z_statistics_sample_size, gap_start_fine,
                         gap_continue_fine, database_corridor_width, z_threshold, 'learn_sample.fasta',
                         start_coordinate,
                         end_coordinate, extension, min_length, intersection_level, z_filename, cores, login, password))

                with open('journal.txt', 'a') as file:
                    file.write(str(set_scores) + '\n')

            set_scores, position_sets = (list(t) for t in zip(*sorted(zip(set_scores, position_sets), reverse=True)))

        iteration_counter = 0
        while True:
            children = []
            for i in range(int(n_children / 2)):
                children.extend(crossingover(set_scores, position_sets))

            for i in range(len(children)):
                children[i] = repeats_remove(children[i])
                children[i] = mutation(children[i], mut_per, mut_bal)

            position_sets.extend(children)

            for position_set in children:
                set_scores.append(
                    scan(position_set, number, class_number, R2, Kd, z_statistics_sample_size, gap_start_fine,
                         gap_continue_fine, database_corridor_width, z_threshold, 'learn_sample.fasta',
                         start_coordinate,
                         end_coordinate, extension, min_length, intersection_level, z_filename, cores, login, password))

            set_scores, position_sets = (list(t) for t in zip(*sorted(zip(set_scores, position_sets), reverse=True)))

            for i in range(n_children):
                set_scores.pop(-1)
                position_sets.pop(-1)

            if best_score < set_scores[0]:
                best_score = set_scores[0]

                check_score = scan(position_sets[0], number, class_number, R2, Kd, z_statistics_sample_size,
                                   gap_start_fine, gap_continue_fine, extended_corridor_width, z_threshold,
                                   'extended_full_sample.fasta', start_coordinate, end_coordinate, extension,
                                   min_length, intersection_level, z_filename, cores, login, password)

                with open('overfit.txt', 'a') as file:
                    file.write(str(best_score) + '   ' + str(check_score) + '   ' + '{0:.3f}'.format(
                        check_score / best_score) + '\n')
                    file.close()

                if check_score > best_check_score:
                    best_check_score = check_score
                    with open('classes/class' + str(class_number) + '.json', 'w') as file:
                        # indent=2 is not needed but makes the file human-readable
                        json.dump(position_sets, file)
                        file.close()

                    matrix_creator.web_create(position_sets[0], number, login, password)
                    with open('bin/' + str(1000 + number) + '.txt', 'r') as file1:
                        with open('classes/matrix' + str(class_number) + '.txt', 'w') as file2:
                            for line in file1:
                                file2.write(line)

                iteration_counter = 0

            if iteration_counter > iteration_limit:
                scan(position_sets[0], number, class_number, R2, Kd, z_statistics_sample_size,
                     gap_start_fine, gap_continue_fine, extended_corridor_width, z_threshold,
                     'extended_full_sample.fasta', start_coordinate, end_coordinate, extension,
                     min_length, intersection_level, z_filename, cores, login, password)

                result_indexes = []
                with open('classes/' + str(class_number) + '.scan_result_positions.txt', 'r') as file:
                    for line in file:
                        if line != '':
                            result_indexes.append(int(line))
                result_indexes.reverse()

                database = []
                with open('bin/full_sample.fasta', 'r') as file:
                    for line in file:
                        if line[0] == '>':
                            database.append(line)
                        else:
                            database[-1] = database[-1] + line

                extended_database = []
                with open('bin/extended_full_sample.fasta', 'r') as file:
                    for line in file:
                        if line[0] == '>':
                            extended_database.append(line)
                        else:
                            extended_database[-1] = extended_database[-1] + line

                for index in result_indexes:
                    database.pop(index)
                    extended_database.pop(index)

                with open('bin/database_work_copy.fasta', 'w') as file:
                    for entity in database:
                        file.write(entity)

                with open('bin/extended_full_sample.fasta', 'w') as file:
                    for entity in extended_database:
                        file.write(entity)

                database_len = 0
                file = open('bin/database_work_copy.fasta', 'r')
                for line in file:
                    if line[0] == '>':
                        database_len += 1

                best_score = 0
                best_check_score = 0
                class_number += 1
                break

            with open('journal.txt', 'a') as file:
                file.write(str(set_scores) + '\n')
                file.close()

            iteration_counter += 1


main()
