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
import re
import sys


def sort(data):
    sorted_data = []
    while data:
        index = nested_min(data, 4)
        sorted_data.append(data[index])
        data.pop(index)
    return sorted_data


def join(data1, data2):
    joined = []

    for i in range(len(data2)):
        trigger = 0
        for j in range(len(data1)):
            if min(data1[j][5], data2[i][5]) - max(data1[j][4], data2[i][4]) > 0:
                trigger = 1
                break
        if trigger == 0:
            joined.append(data2[i])

    for j in range(len(data1)):
        joined.append(data1[j])

    return joined


def nested_max(data, key):
    max_value = data[0][key]
    index = 0
    for i in range(len(data)):
        line = data[i]
        if line[key] > max_value:
            max_value = line[key]
            index = i
    return index


def nested_min(data, key):
    min_value = data[0][key]
    index = 0
    for i in range(len(data)):
        line = data[i]
        if line[key] < min_value:
            min_value = line[key]
            index = i
    return index


def main():
    from config import classes, z_threshold as threshold, min_length
    chromosome = str(sys.argv[1])
    node = int(sys.argv[2])
    min_class_number = 0

    for matrix in range(min_class_number, min_class_number + classes):
        results = []
        file = open('results/' + str(re.findall(r'\d+', chromosome)[0]) + '/' + str(matrix) + '/results' + str(node) + '.chromosome.' + chromosome + '.' + str(matrix), 'r')
        for line in file:
            string = line[:-1].replace('(', '').replace(')', '').replace('\' ', '').replace('\'', '').split(', ')
            if float(string[2]) >= threshold and int(string[5]) - int(string[4]) > min_length:
                results.append((int(string[0]), float(string[1]), float(string[2]), int(string[3]), int(string[4]), int(string[5]), string[6].split(' '), string[7].split(' ')))
        file.close()

        sort_results = []
        while results:
            index = nested_max(results, 2)
            sort_results.append(results[index])

            if index < len(results) - 1:
                next_start = results[index + 1][4]
                while next_start <= results[index][5]:
                    results.pop(index + 1)
                    if index < len(results) - 1:
                        next_start = results[index + 1][4]
                    else:
                        break
            if index > 0:
                prev_end = results[index - 1][5]
                while prev_end >= results[index][4]:
                    results.pop(index - 1)
                    index -= 1
                    if index > 0:
                        prev_end = results[index - 1][5]
                    else:
                        break

            results.pop(index)

        with open('results/' + str(re.findall(r'\d+', chromosome)[0]) + '/' + str(matrix) + '/sorted_results' + str(node) + '.chromosome.' + chromosome + '.' + str(matrix), 'w') as f:
            json.dump(sort_results, f)

    results = []
    for i in range(min_class_number, min_class_number + classes):
        with open('results/' + str(re.findall(r'\d+', chromosome)[0]) + '/' + str(i) + '/sorted_results' + str(node) + '.chromosome.' + chromosome + '.' + str(i), 'r') as f:
            load = (json.load(f))

        results = join(results, sort(load))

    with open('results/length/' + str(re.findall(r'\d+', chromosome)[0]) + '.txt', 'a') as file:
        file.write(chromosome + ' ' + str(node) + ' ' + str(len(results)) + '\n')

    with open('results/' + str(re.findall(r'\d+', chromosome)[0]) + '/results' + str(node) + '.chromosome.' + chromosome + '.joined', 'w') as f:
        json.dump(results, f)


main()
