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

from config import min_chromosome, max_chromosome, chromosome_filename, chromosome_extension


for number in range(min_chromosome, max_chromosome + 1):
    file_name = 'chromosome.' + str(number)
    file = open('chromosomes/' + chromosome_filename + str(number) + chromosome_extension, 'r')
    trigger = 0
    chromosome = ''
    for line in file:
        if trigger == 0:
            trigger = 1
        else:
            chromosome += line[:-1]

    file.close()

    print(number)

    file = open('edited_chromosomes/chromosome' + str(number), 'w')
    file.write(chromosome)
    file.close()
