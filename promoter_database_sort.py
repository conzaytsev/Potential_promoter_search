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

from config import promoter_database_filename, extended_promoter_database_filename


full_database_1 = []
database_index_1 = []
file = open(promoter_database_filename, 'r')
for line in file:
    if line[0] == '>':
        database_index_1.append(line[:9])
        full_database_1.append(line)
    else:
        full_database_1[-1] = full_database_1[-1] + line
file.close()
print(database_index_1)
print(len(database_index_1))
print(len(full_database_1))

full_database_2 = []
database_index_2 = []
file = open(extended_promoter_database_filename, 'r')
for line in file:
    if line[0] == '>':
        database_index_2.append(line[:9])
        full_database_2.append(line)
    else:
        full_database_2[-1] = full_database_2[-1] + line
full_database_2[-1] = full_database_2[-1] + '\n'
file.close()
print(database_index_2)
print(len(database_index_2))

sorted_full_database_2 = []
for index in database_index_1:
    for i in range(len(database_index_2)):
        if database_index_2[i] == index:
            sorted_full_database_2.append(full_database_2[i])
            database_index_2.pop(i)
            full_database_2.pop(i)
            break

print(len(sorted_full_database_2))

with open('extended_sorted_' + promoter_database_filename, 'w') as file:
    for entity in sorted_full_database_2:
        file.write(entity)
