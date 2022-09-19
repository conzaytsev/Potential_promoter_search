# Potential promoter search

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

## Introduction

Here are several programs, which can be used for promoter sequence classification and potential promoters search in genome.
Promoter classification is done using genetic algorithm and MAHDS method.

MAHDS calculation is performed on a remote server, which is accessed via web connection.
Login and password for server calculations are required. They can be obtained at http://victoria.biengi.ac.ru/mahds/auth.

Potential promoter search is done using pairwise alignments between created matrices for promoter classes and genome sequences.

Parameters are set in the config.py file.

## Requirements

`Python 3.7` with `numpy`, `scipy`, `requests`, `beautifulsoup4`
`PyPy 3.7`

## Input files

Chromosomes in FASTA format should be placed in chromosomes/ directory.
Promoter database in FASTA format should be placed in / directory.
Promoter database with extended promoter region in FASTA format should be placed in / directory.

## Run order
python3 chromosomes_to_single_line.py

	Converts chromosome files to supported format.


python3 chromosome_mix.py

	Creates files with mixed chromosomes for false positive level estimation.


python3 promoter_database_sort.py

	Sorts Promoter database and Extended promoter database in the same order.


python3 genetic_partition_extended.py

	Creates promoter sequences classes and their matrices.

	Calculated position weight matrices for all the classes are saved to pair_matrix and single_matrix directories. 
	Matrices in pair_matrix directory include correlations between neighbouring nucleotides. 
	Matrices in single_matrix directory are simple position weight matrices. 
	All matrices are saved in JSON and as a numpy dump. 
	When new matrix appears in the folder, previous one is fully calculated. 
	For example, when last matrix in pair_matrix folder is 1.pair_matrix.json, then 0.pair_matrix.json is ready for use.

	journal.txt contains class volumes for each of n_sets on each iteration.
	overfit.txt contains only best class volumes for each of classes. 

	As an example:
	0
	880   7623   8.662
	986   7732   7.842
	1170   9487   8.109
	1173   8833   7.530



	1
	224   1272   5.679
	225   1176   5.227

	Here, 0 is number of class. 
	First column contains best volumes for the 0 class when scanning learn sample from database_filename. 
	Second column contains volumes for the 0 class when scanning full extended_database_filename with the same matrix, which achieved the result from first column. 
	Third column value is fraction between second column value and first column value. 
	When class calculation is finished, calculation of the next class starts.


pypy3 z_statistics.py

	Calculates Monte Carlo statistics for potential promoters Z value estimation.


pypy3 promoter_search.py chromosome spiral matrix

	Scans chromosome for potential promoters.

	chromosome: chromosome number or mixed chromosome number (N or mix.N , min_chromosome ≤ N ≤ max_chromosome).
	spiral: 1 -- + strand, standard order
		2 -- + strand, reversed order
		3 -- - strand, standard order
		4 -- - strand, reversed order
	matrix: class matrix number (0, 1, ...)

	Results are saved to results/chromosome_number/spiral


pypy3 results_rank_intersect.py chromosome spiral

	Removes intersecting results.

	Intersected results are saved to results/chromosome_number


## WARNING
	genetic_partition_extended.py and parallel_promoter_search.py require a lot of resources and lots of time.
	All of these were run on a system with 4 x Intel Xeon Platinum 8270 (104 cores in total). 
	It took about 2 weeks to create promoter classes and 4 weeks to scan all 24 chromosomes.



## PARAMETERS
	min_chromosome -- smallest chromosome number.

	max_chromosome -- largest chromosome number. 
	Only numerical chromosome names can be used. 
	For Human genome X chromosome was designated as 23 and Y as 24.

	chromosome_filename -- name of chromosome files before the number.

	chromosome_extension -- name of chromosome files after the number. 
	(Full name should be chromosome_filename + number + chromosome_extension)

	promoter_database_filename -- file with promoter database in FASTA format. 
	All promoters should have coordinates from -499 to 100.

	extended_promoter_database_filename -- file with FASTA promoter database. 
	Promoters should have coordinates from -499 + start_coordinate - extension to end_coordinate - 500 + extension. 
	Promoters in this database should be sorted in the same order, as promoters from database_filename.

	mahds_login -- login for MAHDS server.

	mahds_password -- password for MAHDS server.

	classifier_mixed_chromosome -- mixed chromosome number. 
	Mixed chromosome is used to calculate statistical significance value. 

	n_sets -- number of promoter sets for genetic algorithm.

	set_size -- number of promoters in each set.

	n_children -- number of new sequence sets generated on each iteration of the algorithm.

	mut_per -- percentage of mutations in each new set.

	mut_bal -- balance between type 1 and 2 mutations in new sets. 
	0 -> all mutations are type 1, 100 -> all mutations are type 2.

	Kd -- normalization parameter for position weight matrices.

	R2 -- normalization parameter for position weight matrices.

	learn_partition -- maximum volume of learn database sample as a fraction of full database volume.

	min_set -- minimum volume of learn database sample as number of sequences.

	iteration_limit -- number of algorithm iterations without improvement for stopping class selection.

	start_coordinate -- coordinate of the first nucleotide for position weight matrix calculation as its position in sequence (0..599).

	end_coordinate -- coordinate of the last nucleotide for position weight matrix calculation as its position in sequence (0..599).

	extension -- half difference between extended_database_filename sequences length and end_coordinate - start_coordinate.

	intersection_level -- minimum intersection between extended aligned sequence and its -499:20 coordinates as a fraction of its length.

	classes -- number of classes, created by genetic_partition_extended.py and used for potential promoters search.

	z_statistics_sample_size -- volume of generated statistics for statistical significance calculation.

	gap_start_fine -- gap for starting a gap in a sequence to position weight matrix alignment.

	gap_continue_fine -- gap for continuing an open gap in a sequence to position weight matrix alignment.

	query_length -- length of the window for chromosome scanning.

	min_length -- minimum alignment length for class selection.

	z_threshold -- minimum statistical significance value for class definition.

	database_corridor_width -- width of the corridor for alignments with sequences from promoter_database_filename. 
	Maximum distance from main diagonal. 
	Used for scanning database_filename.

	extended_corridor_width -- width of the corridor for alignments with sequences from extended_promoter_database_filename.

	chromosome_corridor_width -- width of the corridor for potential promoters search in chromosomes.

	cores -- number of threads used for calculation.
	Used in genetic_partition_extended.py, z_statistics.py, parallel_promoter_search.py.
