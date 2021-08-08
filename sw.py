import itertools
import csv
from multiprocessing import Process, Manager, Pool
import re
import pandas as pd

fragment_range_start = None
x = 60
fragment_range_end = 2000
enzyme = 'MspI'
enzyme2 = 'BfaI'


##### new sw.py #####
def grab_chunk(length_size, step):
    #  Determine all rows with length between n and n + step # now n + 200
    end_length = length_size + 200  # changed from + step
    rows = df[(df['length'] >= length_size) & (df['length'] < end_length)]

    #  Set up patterns for column names
    patterns = ['CG', 'CA', 'CT', 'CC', 'TG', 'AG', 'GG']

    #  For each of the columns, get the sum and then sum those
    counts = rows[patterns].sum().sum()

    #  Assign result to the global dictionary
    key = "c_" + str(length_size) + "_" + str(end_length)
    d[key] = [counts, len(rows)]


if __name__ == '__main__':
    #  Determine the number of cores to be used
    numCores = 20

    #  Read data
    filename = 'intersected_COV.bed'
    # filename = 'short_COV.bed'
    df = pd.read_csv(filename, sep='\t',
                     names=['chr', 'start', 'end', '.', '_', 'strand', 'length', 'CG', 'CA', 'CT', 'CC', 'TG', 'AG',
                            'GG'], skiprows=[0])
    #  Define the global dictionary
    manager = Manager()
    d = manager.dict()

    if fragment_range_start is None:
        fragment_range_start = 2 * x

    step = 10

    # Iniialize dictionary
    fragment_range = range(fragment_range_start, fragment_range_end - 190, step)
    for x in fragment_range:
        d["c_{0}_{1}".format(x, x + 200)] = [0, 0]  # changed from x + step

    # Parallelize the code that determines the counts for the dictionary
    pool = Pool(numCores)
    results = pool.starmap(grab_chunk, zip(fragment_range, itertools.repeat(step)))  # may have to change to 200?
    
    if enzyme2 is not None:
        output = enzyme + "_" + enzyme2 + "_sw.csv"
    else:
        output = enzyme + "_sw.csv"
    with open(output, 'w+') as write_obj:
        headers = ["Start_Length", "End_Length", "C_per_read", "Total_C", "Total_fragments", "Features_covered"]
        writer = csv.DictWriter(write_obj, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        start_re = f'_([0-9]+)_'
        regex2 = re.compile(start_re)
        maximum = 0
        max_start = 0
        max_end = 0
        for each in d.keys():
            for match in regex2.finditer(each):
                start = int(match.group(1))  # finding numerical start
            if d[each][1] == 0:
                C_per_read = 0
            else:
                C_per_read = d[each][0] / d[each][1]
            if C_per_read > maximum:
                maximum = C_per_read
                max_start = start
                max_end = start + 200
            row_dict = {'Start_Length': start, 'End_Length': start + 200, 'C_per_read': C_per_read,
                        'Total_C': d[each][0], 'Total_fragments': d[each][1]}  # "Features_covered": 0
            writer.writerow(row_dict)
        print("max:", maximum)
        c_txt = open("maximum_c.txt", "w+")
        c_txt.write(
            "Maximum cytosine count per fragment is " + str(maximum) + " in the fragment range " + str(max_start)
            + " to " + str(max_end) + ".")
        c_txt.close()
