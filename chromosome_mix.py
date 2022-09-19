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

import re
from random import shuffle
from config import min_chromosome, max_chromosome


for number in range(min_chromosome, max_chromosome + 1):
    file_name = 'chromosome.' + str(number)
    file = open('edited_chromosomes/' + file_name, 'r')
    chromosome = ''

    for line in file:
        chromosome += line
    file.close()

    print(len(chromosome))
    chromosome = re.sub('[NURYKMSWBDHV]', '', chromosome)
    print(len(chromosome))

    ch_list = list(chromosome)
    shuffle(ch_list)
    mixed_chromosome = ''.join(ch_list)
    print(len(mixed_chromosome))

    file = open('edited_chromosomes/chromosome.mix.' + str(number), 'w')
    file.write(mixed_chromosome)
    file.close()
