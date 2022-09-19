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

from math import ceil


def align(sequence, single_matrix, pair_matrix, gap_start_fine, gap_continue_fine, corridor_width):
    K1 = len(sequence)
    K = 600

    F = {}
    for i in range(K):
        F[i] = {}

    k1 = K / K1
    k2 = K1 / K

    for i in range(K):
        for j in range(K1):
            if k1 * j - corridor_width < i < k1 * j + corridor_width or k2 * i - corridor_width < j < k2 * i + corridor_width:
                F[i][j] = {'score': 0, 'pr_i': 0, 'pr_j': 0, 'pr_pj': 0, 'gap': gap_start_fine}

    F[0][0]['score'] = single_matrix[0][sequence[0]]
    F[0][0]['pr_i'] = -1
    F[0][0]['pr_j'] = -1
    F[0][0]['pr_pj'] = 0

    for i in range(1, min(max(corridor_width, ceil(corridor_width * k1)), K)):
        F[i][0]['score'] = max(single_matrix[i][sequence[0]], F[i - 1][0]['score'] + F[i - 1][0]['gap'])

        if F[i][0]['score'] == single_matrix[i][sequence[0]]:
            F[i][0]['pr_i'] = -1
            F[i][0]['pr_j'] = -1
            F[i][0]['pr_pj'] = -1
        else:
            F[i][0]['pr_i'] = i - 1
            F[i][0]['pr_j'] = 0
            F[i][0]['pr_pj'] = -1
            F[i][0]['gap'] = gap_continue_fine

    for j in range(1, min(max(corridor_width, ceil(corridor_width * k2)), K1)):
        F[0][j]['score'] = max(single_matrix[0][sequence[j]], F[0][j - 1]['score'] + F[0][j - 1]['gap'])

        if F[0][j]['score'] == single_matrix[0][sequence[j]]:
            F[0][j]['pr_i'] = -1
            F[0][j]['pr_j'] = -1
            F[0][j]['pr_pj'] = j
        else:
            F[0][j]['pr_i'] = 0
            F[0][j]['pr_j'] = j - 1
            F[0][j]['pr_pj'] = F[0][j - 1]['pr_pj']
            F[0][j]['gap'] = gap_continue_fine

    for i in range(1, K):
        for j in F[i]:
            if j != 0:
                trigger = 0

                vertical = 0
                horizontal = 0
                diagonal = 0

                if F[i].get(j - 1) is not None:
                    vertical = F[i][j - 1]['score'] + F[i][j - 1]['gap']
                    trigger = 1
                    F[i][j]['score'] = vertical

                if F[i - 1].get(j) is not None:
                    horizontal = F[i - 1][j]['score'] + F[i - 1][j]['gap']
                    if trigger == 1:
                        trigger = 5
                    else:
                        trigger = 2
                        F[i][j]['score'] = horizontal

                if F[i - 1].get(j - 1) is not None:
                    if F[i - 1][j - 1]['pr_pj'] >= 0:
                        diagonal = F[i - 1][j - 1]['score'] + pair_matrix[i][sequence[F[i - 1][j - 1]['pr_pj']]][
                            sequence[j]]
                    else:
                        diagonal = F[i - 1][j - 1]['score'] + single_matrix[i][sequence[j]]
                    if trigger == 1:
                        F[i][j]['score'] = max(vertical, diagonal)
                    elif trigger == 2:
                        F[i][j]['score'] = max(horizontal, diagonal)
                    else:
                        F[i][j]['score'] = max(vertical, horizontal, diagonal)

                if F[i][j]['score'] == diagonal:
                    F[i][j]['pr_i'] = i - 1
                    F[i][j]['pr_j'] = j - 1
                    F[i][j]['pr_pj'] = j

                elif F[i][j]['score'] == horizontal:
                    F[i][j]['pr_i'] = i - 1
                    F[i][j]['pr_j'] = j
                    F[i][j]['pr_pj'] = -1
                    F[i][j]['gap'] = gap_continue_fine

                else:
                    F[i][j]['pr_i'] = i
                    F[i][j]['pr_j'] = j - 1
                    F[i][j]['pr_pj'] = F[i][j - 1]['pr_pj']
                    F[i][j]['gap'] = gap_continue_fine

                if F[i][j]['score'] < single_matrix[i][sequence[j]]:
                    F[i][j]['score'] = single_matrix[i][sequence[j]]
                    F[i][j]['pr_i'] = -1
                    F[i][j]['pr_j'] = -1
                    F[i][j]['pr_pj'] = j

    max_score = F[0][0]['score']
    end_i = 0
    end_j = 0

    for i in range(1, K):
        for j in F[i]:
            if F[i][j]['score'] > max_score:
                max_score = F[i][j]['score']
                end_i = i
                end_j = j

    al_sequence = ''
    al_numbers = ''
    start_j = 0
    finish_j = end_j

    while end_i >= 0 and end_j >= 0:
        pr_i = int(F[end_i][end_j]['pr_i'])
        pr_j = int(F[end_i][end_j]['pr_j'])

        if pr_i != end_i:
            al_numbers = ' ' + str(end_i) + al_numbers
        else:
            al_numbers = ' -' + al_numbers

        if pr_j != end_j:
            al_sequence = ' ' + sequence[end_j] + al_sequence
        else:
            al_sequence = ' -' + al_sequence

        start_j = end_j
        end_i = pr_i
        end_j = pr_j
    return max_score, al_sequence, al_numbers, start_j, finish_j
