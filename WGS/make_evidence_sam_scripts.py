import pandas as pd
from glob import glob
from collections import defaultdict
import csv
import time
import subprocess
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('start_well_index', help='key to identify which wells to do')
args = parser.parse_args()

start_well_index = int(args.start_well_index)

sam_cols = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL'] + ['EXTRA'+str(i) for i in range(7)]

# this is lenient, not all of these reads will actually be used, but that's ok, just showing everything close is fine
def is_evidence(row, cmp):
    left = int(row['POS']) 
    right = left+150
    mut_matches = []
    for pos in cmp[row['RNAME']]:
        if pos > left and pos < right:
            mut_matches.append(row['RNAME']+str(pos))
    return ';'.join(mut_matches)

def alignment(r):
    # uses the cigar to format the sequence to show in an alignment
    # adds - for deletions, puts insertions in parentheses (will be shown only by a red line that shows the insertion on hover)
    seq, cigar = r['SEQ'], r['CIGAR']
    new_seq = ''
    current_nums = ''
    current_spot = 0
    for character in cigar:
        if character in ['M', 'D', 'I', 'S', 'H', '=', 'X']:
            if character in ['S', 'H', '=', 'X']:
                return new_seq
            elif character == 'M':
                new_seq += seq[current_spot:current_spot+int(current_nums)]
                current_spot += int(current_nums)
            elif character == 'D':
                new_seq += '-'*int(current_nums)
            else:
                new_seq += '(' + seq[current_spot:current_spot+int(current_nums)] + ')'
                current_spot += int(current_nums)
            current_nums = ''
        else:
            current_nums += character
    return new_seq


wells = ['P1B02', 'P1B03', 'P1B04', 'P1B07', 'P1B11', 'P1C02', 'P1C04', 'P1C05', 'P1C06', 'P1C07', 'P1C08', 'P1C09', 'P1C11', 'P1D03', 'P1D09', 'P1E04', 'P1E09', 
'P1E11', 'P1F05', 'P1F07', 'P1F08', 'P1F10', 'P1F11', 'P1G04', 'P1G05', 'P1G08', 'P1G09', 'P1G10', 'P1G11', 'P1H11', 'P2B04', 'P2B05', 'P2B07', 'P2B08', 
'P2B09', 'P2B10', 'P2B11', 'P2C02', 'P2C04', 'P2C05', 'P2C06', 'P2C10', 'P2C11', 'P2D03', 'P2D06', 'P2D08', 'P2D11', 'P2E06', 'P2E08', 'P2E11', 'P2F02', 
'P2F07', 'P2F09', 'P2F11', 'P2G04', 'P2G05', 'P2G08', 'P2G09', 'P2G10', 'P2G11', 'P3B07', 'P3B08', 'P3B10', 'P3B11', 'P3C03', 'P3C04', 'P3C05', 'P3C07', 
'P3C10', 'P3C11', 'P3D02', 'P3D03', 'P3D05', 'P3D09', 'P3D10', 'P3D11', 'P3E02', 'P3E08', 'P3E11', 'P3F03', 'P3F05', 'P3F07', 'P3F09', 'P3F11', 'P3G02', 
'P3G05', 'P3G06', 'P3G09', 'P3G10', 'P3G11']
gens = [70, 1410, 2640, 5150, 7530, 10150]
chromo_mut_positions = defaultdict(list)
#done_wells = [i.split('/')[-1].split('.')[0] for i in glob('../../Output/Browser/evidence_sams/*.tsv')]
for well in wells[start_well_index*10:start_well_index*10+10]:
    f = '../../Output/WGS/well_output/' + well + '_filtered.tsv'
    otime = time.time()
    #if well not in done_wells:
    print(well)
    td = pd.read_csv(f, delimiter='\t')
    td['pos2'] = td['POS']
    bedfile = '../../Output/Browser/evidence_sams/'+well+'_mut_regions.bed'
    td[['CHROM', 'POS', 'pos2']].to_csv(bedfile, sep='\t', index=False, header=False)
    for entry in td.as_matrix(['CHROM', 'POS']):
        chromo_mut_positions[entry[0]].append(entry[1])
    dats = []
    for gen in gens:
        try:
            outsam = '../../Output/Browser/evidence_sams/G'+str(gen)+well+'.tsv'
            subprocess.call(['samtools view ../../Output/WGS/work/G' + str(gen) + '_' + well + '.sam' + ' -L ' + bedfile + ' -o ' + outsam], shell=True)
            td = pd.read_csv(outsam, delimiter='\t', header=None, names=sam_cols, index_col=False)
            td['Gen'] = [gen]*len(td)
            td['Evidence_for'] = td.apply(lambda r: is_evidence(r, chromo_mut_positions), axis=1)
            td['Aligned'] = td.apply(lambda r: alignment(r), axis=1)
            dats.append(td)
        except FileNotFoundError:
            print('File missing for', str(gen), well)
    pd.concat(dats).to_csv('../../Output/Browser/evidence_sams/'+well+'_expanded_sam.tsv', sep='\t', index=False)
    print('done in', time.time()-otime)
