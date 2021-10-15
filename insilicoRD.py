import csv
from glob import glob
import pandas as pd
import os
import ray
import itertools
import re
import subprocess
import pybedtools
from multiprocessing import Pool, Manager


ray.init(memory=10**9, object_store_memory=10**9, ignore_reinit_error=True)


def insilicoRD(enzyme, enzyme2=None, x=60, fragment_range_start=None, fragment_range_end=2000):
    # get_coordinates -- output: digested_COORD.bed
    def get_coordinates(enzyme, enzyme2=None):
        enzyme_string = ""
        enzyme_string2 = ""
        with open('RE_table.csv', 'r') as read_obj:
            reader = csv.reader(read_obj)
            for row in reader:
                if enzyme.lower() == row[0].lower():
                    enzyme_string = row[1]
                if enzyme2 is not None:
                    if enzyme2.lower() == row[0].lower():
                        enzyme_string2 = row[1]
        if enzyme_string == "":
            raise Exception("Enzyme 1 is invalid!")
        if enzyme_string2 == "" and enzyme2 is not None:
            raise Exception("Enzyme 2 is invalid!")
        if enzyme2 is None:
            enzyme_to_compile = f'({enzyme_string})'
        else:
            enzyme_to_compile = f'({enzyme_string}|{enzyme_string2})'

        @ray.remote
        def create_chr(chromosome_number):
            chromosome_file = "mm10_chr" + str(chromosome_number) + ".fa"
            f = open(chromosome_file, 'r')
            genome_list = f.read().strip().split("\n")
            genome_list.pop(0)  # delete the first line
            genome_string = ""
            for line in genome_list:
                genome_string += line
            regex = re.compile(enzyme_to_compile)  # finding instances of enzyme in genome
            # output file
            output_file = 'COORD_chr' + str(chromosome_number) + '.csv'
            # write a csv file with chr & start
            with open(output_file, 'w+', newline='') as write_obj:
                headers = ['chr', 'start', 'end', '.', '.', 'strand', 'length']
                writer = csv.DictWriter(write_obj, fieldnames=headers, delimiter='\t')
                writer.writeheader()
                # regular expression
                for match in regex.finditer(genome_string):
                    row_dict = {'chr': chromosome_number, 'start': match.start(), '.': ".", 'strand': "+"}
                    writer.writerow(row_dict)

            # count number of rows
            read_obj = open(output_file, 'r')
            reader = csv.DictReader(read_obj, delimiter='\t')
            number_of_rows = len(list(reader))

            # make a list of the end
            with open(output_file, 'r') as read_obj:
                reader = csv.DictReader(read_obj, delimiter='\t')
                counter = 0
                end_list = []
                for row in reader:
                    counter += 1
                    if not counter == number_of_rows:
                        for next_row in reader:
                            row_end = next_row['start']
                            end_list.append(row_end)
                            row['end'] = row_end

            # now can make a csv with chr, start, end, & length
            with open(output_file, 'w+', newline='') as write_obj:
                headers = ['chr', 'start', 'end', '.', '.', 'strand', 'length']
                writer = csv.DictWriter(write_obj, fieldnames=headers, delimiter='\t')
                writer.writeheader()
                i = 0
                for match in regex.finditer(genome_string):
                    if i < number_of_rows - 1:
                        length = int(end_list[i]) - int(match.start())
                        start = int(match.start())
                        end = int(end_list[i])
                        row_dict = {'chr': "chr" + str(chromosome_number), 'start': start, 'end': end, '.': ".",
                                    'strand': "+", 'length': length}
                        writer.writerow(row_dict)
                    else:
                        length = int(len(genome_string)) - int(match.start())
                        start = int(match.start())
                        end = len(genome_string)
                        row_dict = {'chr': "chr" + str(chromosome_number), 'start': start, 'end': end, '.': ".",
                                    'strand': "+", 'length': length}
                        writer.writerow(row_dict)
                    i += 1
            print("Coordinate file outputted: Chromosome " + str(chromosome_number))

        def create_chr2():
            tasks = [create_chr.remote(i) for i in range(1, 20)]
            taskM = create_chr.remote('M')
            taskX = create_chr.remote('X')
            taskY = create_chr.remote('Y')
            ray.get(tasks)
            ray.get(taskM)
            ray.get(taskX)
            ray.get(taskY)

        create_chr2()

        # genome sequence to a single string
        output_paradigm = "COORD_chr*"
        merged_csv = "digested_COORD1.bed"
        list_of_files = []
        for file in glob(output_paradigm):
            list_of_files.append(file)
        merged_obj = pd.concat([pd.read_csv(file) for file in list_of_files], sort=False, axis='rows')
        merged_obj = merged_obj.iloc[1:] # CHANGED 8/6
        merged_obj.to_csv(merged_csv, index=False, encoding="utf-8")
        with open("digested_COORD1.bed", 'r') as f:
            with open("digested_COORD.bed", 'w') as f1:
                f.readline()  # skip header line
                for line in f:
                    f1.write(line)
        chromosome_numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 'M', 'X', 'Y']
        for chromosome_number in chromosome_numbers:
            output_file = 'COORD_chr' + str(chromosome_number) + '.csv'
            os.remove(output_file)

        print("---------------------------------------finished get_coordinates---------------------------------------")

    def bedtools_mod():

        a = pybedtools.BedTool('features2.bed')
        a1 = a.merge(d=-1).saveas('merged_features.bed')
        b = pybedtools.BedTool('digested_COORD.bed')
        b1 = b.sort().saveas('sorted_digested_COORD.bed')
        c = b1.intersect(a1, wa=True).saveas('intersected.bed')
        d = c.merge(c=[4, 5, 6, 7], o=['distinct'], d=-1).saveas('merged_intersected.bed')

    ray.init(memory=9.28 * 10 ** 9, object_store_memory=111.6 * 10 ** 9, ignore_reinit_error=True)

    # input = merged_intersected.bed --> output intersected_COV

    def cytosines_covered(x=60, fragment_range_start=None, fragment_range_end=2000):
        headers = ['chr', 'start', 'end', '.', '.', 'strand', 'length', 'CG', 'CA', 'CT', 'CC', 'TG', 'AG', 'GG']
        d = {}
        chromosome_numbers = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
                              '14', '15', '16', '17', '18', '19', 'M', 'X', 'Y']
        if fragment_range_start is None:
            fragment_range_start = 2 * x
        # df = pd.read_csv("merged_intersected.bed", sep='\t',
        #                  names=['chr', 'start', 'end', '.', '_', 'strand', 'length'])
        # max_length = df['length'].max()
        # fragment_range_end = max_length

        with open('intersected_COV.bed', 'w') as write_obj:
            writer = csv.DictWriter(write_obj, delimiter='\t', fieldnames=headers)
            writer.writeheader()

        @ray.remote
        def create(n):
            chr_re = f'([0-9]+)|(M)|(X)|(Y)'
            regex = re.compile(chr_re)
            with open('merged_intersected.bed', 'r', newline='') as read_obj:
                reader = csv.DictReader(read_obj, delimiter='\t', fieldnames=headers)
                num_lines = sum(1 for each in open('merged_intersected.bed'))
                factor = num_lines // 7 + 1
                line_start = factor * (n - 1)
                if n == 7:
                    line_end = num_lines
                else:
                    line_end = factor * n
                for row in itertools.islice(reader, line_start, line_end):
                    if fragment_range_start < int(row['length']) < fragment_range_end:
                        # chr1_start = [end, ".", ".", "+", ]
                        for match in regex.finditer(row['chr']):
                            chr_num = match.group()  # finding 1-19 or M or X or Y
                        d[str(row['chr']) + "_" + str(int(row['start']))] = [chr_num, str(row['length'])]  # start 1
                        d[str(row['chr']) + "_" + str(int(row['end']) - x)] = [chr_num, str(row['length'])]  # start 2
                print("Successfully created chr list and start list")
                d_genome = {}  # dictionary of chr: genome_string
                for chromosome in chromosome_numbers:
                    f = open("mm10_chr" + str(chromosome) + ".fa", 'r')
                    genome_list = f.read().strip().split("\n")
                    genome_list.pop(0)  # delete the first line
                    genome_string = ""
                    for line in genome_list:
                        genome_string += line
                    # total_C += genome_string.count('C')
                    d_genome[str(chromosome)] = genome_string
                print("Successfully created genome string dict")
                start_re = f'_([0-9]+)'
                regex = re.compile(start_re)
                for each in d:
                    chr_num = d[each][0]
                    genome_string = d_genome[str(chr_num)]  # genome of chr
                    for match in regex.finditer(each):
                        start = int(match.group(1))
                    end = int(start) + x
                    fragment_string = genome_string[start:end + 1]
                    length = d[each][1]
                    d[each] = [str(chr_num), start, end, length, fragment_string.count('CG'),
                               fragment_string.count('CA'),
                               fragment_string.count('CT'), fragment_string.count('CC'), fragment_string.count('TG'),
                               fragment_string.count('AG'), fragment_string.count('GG')]
                    # now dict is chr1_start : {chr num}, {CG count}, {CA count}, {CT count}, {CC count},
                    # {TG count}, {AG count}, {GG count}
                print("Successfully created fragment list")
                print("Successfully counted Cs")
            # intersected_COV.bed has to be deleted before running
            with open('intersected_COV.bed', 'a+') as write_obj:
                writer = csv.DictWriter(write_obj, delimiter='\t', fieldnames=headers)
                for each in d:
                    row_dict = {'chr': "chr" + str(d[each][0]), 'start': d[each][1], 'end': d[each][2], '.': ".",
                                'strand': "+", 'length': d[each][3], 'CG': d[each][4], 'CA': d[each][5],
                                'CT': d[each][6], 'CC': d[each][7], 'TG': d[each][8], 'AG': d[each][9],
                                'GG': d[each][10]}
                    writer.writerow(row_dict)
            print("Successfully wrote file: " + str(n))

        def create2():
            tasks = [create.remote(i) for i in range(1, 8)]
            ray.get(tasks)

        create2()
        with open('intersected_COV.bed', 'r') as read_obj:
            read_obj.readline()
            C_covered = sum(int(each['CC']) + int(each['CA']) + int(each['CT']) + int(each['CG']) + int(each['TG']) +
                            int(each['AG']) + int(each['GG']) for each in
                            csv.DictReader(read_obj, delimiter='\t', fieldnames=headers))
        total_C = 0
        for chromosome in chromosome_numbers:
            f = open("mm10_chr" + str(chromosome) + ".fa", 'r')
            genome_list = f.read().strip().split("\n")
            genome_list.pop(0)  # delete the first line
            genome_string = ""
            for line in genome_list:
                genome_string += line
            total_C += genome_string.count('C')
        percent_covered = C_covered / total_C
        c_txt = open("C_coverage.txt", "w+")
        c_txt.write(str(C_covered) + '\n')
        c_txt.write(str(total_C) + '\n')
        c_txt.write(str(percent_covered))
        c_txt.close()
        print(
            "------------------------------------------finished cytosines_covered-------------------------------------")
    get_coordinates(enzyme, enzyme2)
    bedtools_mod()
    cytosines_covered(x, fragment_range_start, fragment_range_end)


insilicoRD('mspI', 'BfaI', 60)
