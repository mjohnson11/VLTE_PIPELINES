# Removing "_pilon" from fasta titles and writing the new reference ot the right position
old_fasta = '../../Output/WGS/reference_refinement/pilon_output/w303_vlte.fasta'
new_fasta = '../../Output/WGS/reference/w303_vlte.fasta'
with open(new_fasta, 'w') as outfile:
    with open(old_fasta, 'r') as infile:
        for line in infile:
            outfile.write(line.replace("_pilon", ""))

# Updates a gff based on indels in a pilon .changes file

import pandas as pd
import csv
from collections import defaultdict

change_file = '../../Output/WGS/reference_refinement/pilon_output/w303_vlte.changes'
old_gff = '../accessory_files/orig_w303_ref/w303_ref.gff'
new_gff = '../../Output/WGS/reference/w303_vlte.gff'

# Generating map between old reference index and new reference index
base_map = defaultdict(dict)
with open(change_file, 'r') as infile:
    reader = csv.reader(infile, delimiter=' ')
    for row in reader:
        chromo = row[0].split(':')[0].split('_')[0]
        old_loc = int(row[0].split(':')[1].split('-')[0])
        new_loc = int(row[1].split(':')[1].split('-')[0])
        if row[2] == '.':
            base_map[chromo][old_loc] = new_loc+len(row[3])
        elif row[3] == '.':
            base_map[chromo][old_loc+len(row[2])] = new_loc
        elif len(row[2]) != len(row[3]):
            base_map[chromo][old_loc+len(row[2])] = new_loc + len(row[3])


def old_to_new_base(chromo, old_base):  # makes the conversions
    td = base_map[chromo]
    anchors = sorted(td.keys())
    shift = 0
    for anch in anchors:
        if anch <= old_base:
            shift = td[anch] - anch
        else:
            break
    return old_base + shift

## Double check: all simple substitutions in the change file should 
with open(change_file, 'r') as infile:
    reader = csv.reader(infile, delimiter=' ')
    for row in reader:
        chromo = row[0].split(':')[0].split('_')[0]
        old_loc = int(row[0].split(':')[1].split('-')[0])
        new_loc = int(row[1].split(':')[1].split('-')[0])
        if len(row[2]) == len(row[3]) and '.' not in row:
            try:
                assert new_loc == old_to_new_base(chromo, old_loc)
            except AssertionError:
                print('Index change error, for this row:', row, '. Thinks this is the new base:', old_to_new_base(chromo, old_loc))

with open(old_gff, 'r') as infile:
    reader = csv.reader(infile, delimiter='\t')
    next(reader) #skip header
    with open(new_gff, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['##gff-version 3'])
        for row in reader:
            writer.writerow(row[:3] + [old_to_new_base(row[0], int(r)) for r in row[3:5]] + row[5:])

