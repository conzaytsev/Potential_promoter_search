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


def create_index(sample_length):
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

    sample_index = []
    sample_position = []
    for i in range(sample_length):
        random_number = randint(0, len(database) - 1)
        while random_number in sample_position:
            random_number = randint(0, len(database) - 1)
        sample_index.append(database_index[random_number])
        sample_position.append(random_number)

    return sample_position


def create_file(sample_position, number):
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

    sample = []
    for position in sample_position:
        sample.append(database[position])

    file = open('bin/' + str(1000 + number) + '.cds', 'w')
    file.write('>\n')
    for entity in sample:
        file.write(entity)
    file.write('>')
    file.close()


def old_create_file(sample_position):
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

    sample = []
    for position in sample_position:
        sample.append(database[position])

    file = open('bin/epdat_600.fasta', 'w')
    file.write('>\n')
    for entity in sample:
        file.write(entity)
    file.write('>')
    file.close()


def create(sample_length, class_number, number):
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

    sample = []
    sample_index = []
    sample_position = []
    for i in range(sample_length):
        random_number = randint(0, len(database) - 1)
        while random_number in sample_position:
            random_number = randint(0, len(database) - 1)
        sample.append(database[random_number])
        sample_index.append(database_index[random_number])
        sample_position.append(random_number)

    file = open('bin/' + str(1000 + number) + '.cds', 'w')
    file.write('>\n')
    for entity in sample:
        file.write(entity)
    file.write('>')
    file.close()

    file = open('classes/' + str(class_number) + '.class.txt', 'w')
    for index in sample_index:
        file.write(index + '\n')
    file.close()

    file = open('classes/' + str(class_number) + '.class_positions.txt', 'w')
    for position in sample_position:
        file.write(str(position) + '\n')
    file.close()


def recreate(class_number, number):
    class_positions = []
    file = open('classes/' + str(class_number) + '.scan_result_positions.txt', 'r')
    for line in file:
        class_positions.append(int(line[:-1]))
    file.close()

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

    sample = []
    sample_index = []
    for position in class_positions:
        sample.append(database[position])
        sample_index.append(database_index[position])

    file = open('bin/' + str(1000 + number) + '.cds', 'w')
    file.write('>\n')
    for entity in sample:
        file.write(entity)
    file.write('>')
    file.close()

    file = open('classes/' + str(class_number) + '.class.txt', 'a')
    for index in sample_index:
        file.write(index + '\n')
    file.close()


def add(class_number):
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

    sample_position = []
    file = open('classes/' + str(class_number) + '.class_positions.txt', 'r')
    for line in file:
        sample_position.append(int(line))
    file.close()

    random_number = randint(0, len(database) - 1)
    while random_number in sample_position:
        random_number = randint(0, len(database) - 1)
    sample_position.append(random_number)

    sample = []
    sample_index = []
    for position in sample_position:
        sample.append(database[position])
        sample_index.append(database_index[position])

    file = open('bin/epdat_600.fasta', 'w')
    file.write('>\n')
    for entity in sample:
        file.write(entity)
    file.write('>')
    file.close()

    file = open('classes/' + str(class_number) + '.class.txt', 'w')
    for index in sample_index:
        file.write(index + '\n')
    file.close()

    file = open('classes/' + str(class_number) + '.class_positions.txt', 'w')
    for position in sample_position:
        file.write(str(position) + '\n')
    file.close()


def replace(class_number):
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

    sample_position = []
    file = open('classes/' + str(class_number) + '.class_positions.txt', 'r')
    for line in file:
        sample_position.append(int(line))
    file.close()

    sample_position.pop(-1)

    random_number = randint(0, len(database) - 1)
    while random_number in sample_position:
        random_number = randint(0, len(database) - 1)
    sample_position.append(random_number)

    sample = []
    sample_index = []
    for position in sample_position:
        sample.append(database[position])
        sample_index.append(database_index[position])

    file = open('bin/epdat_600.fasta', 'w')
    file.write('>\n')
    for entity in sample:
        file.write(entity)
    file.write('>')
    file.close()

    file = open('classes/' + str(class_number) + '.class.txt', 'w')
    for index in sample_index:
        file.write(index + '\n')
    file.close()

    file = open('classes/' + str(class_number) + '.class_positions.txt', 'w')
    for position in sample_position:
        file.write(str(position) + '\n')
    file.close()


def remove(class_number, position):
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

    sample_position = []
    file = open('classes/' + str(class_number) + '.class_positions.txt', 'r')
    for line in file:
        sample_position.append(int(line))
    file.close()

    sample_position.pop(position)

    sample = []
    sample_index = []
    for position in sample_position:
        sample.append(database[position])
        sample_index.append(database_index[position])

    file = open('bin/epdat_600.fasta', 'w')
    file.write('>\n')
    for entity in sample:
        file.write(entity)
    file.write('>')
    file.close()

    file = open('classes/' + str(class_number) + '.class.txt', 'w')
    for index in sample_index:
        file.write(index + '\n')
    file.close()

    file = open('classes/' + str(class_number) + '.class_positions.txt', 'w')
    for position in sample_position:
        file.write(str(position) + '\n')
    file.close()


def database_sample(database_name, divider):
    database = []
    file = open(database_name, 'r')
    for line in file:
        if line[0] == '>':
            database.append(line)
        else:
            database[-1] = database[-1] + line
    file.close()

    size = int(len(database) / divider)

    sample = []
    for i in range(size):
        random_number = randint(0, len(database) - 1)
        sample.append(database[random_number])
        database.pop(random_number)

    file = open('bin/database_work_copy.fasta', 'w')
    for entity in sample:
        file.write(entity)
    file.close()


def database_divider(learn_partition):
    database = []
    file = open('bin/database_work_copy.fasta', 'r')
    for line in file:
        if line[0] == '>':
            database.append(line)
        else:
            database[-1] = database[-1] + line
    file.close()

    size = int(len(database) * learn_partition)

    sample = []
    sample_index = []
    for i in range(size):
        random_number = randint(0, len(database) - 1)
        if random_number in sample_index:
            random_number = randint(0, len(database) - 1)
        sample.append(database[random_number])
        sample_index.append(random_number)

    with open('bin/learn_sample.fasta', 'w') as file:
        for entity in sample:
            file.write(entity)

    with open('bin/full_sample.fasta', 'w') as file:
        for entity in database:
            file.write(entity)
