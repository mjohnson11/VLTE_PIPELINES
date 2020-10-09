from glob import glob
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('outdir', help='directory where chromo dbs etc and future output will go')
args = parser.parse_args()
outdir = args.outdir

alt_allele_thresh = 5

def get_af_traj(r, gu):
    # row, gens_use
    afs = []
    for g in gu:
        ref, alt = r['G' + str(g) + '_ref_counts'], r['G' + str(g) + '_alt_counts']
        tot = ref + alt
        if tot >= 5:
            afs.append(str(g) + '_' + str(alt/tot))
    return ';'.join(afs)

def num_tps_thresh(r, thresh, use_gens):
    return len([g for g in use_gens if r['G'+str(g)+'_alt_counts'] >= thresh])


def annotation_splitter(a):
    # input is full annotation column
    # output is dictionary from ALT alleles (strings) to lists of annotations:
    anns = defaultdict(list)
    for annotation in a.split(','):
        alt_text = annotation.split('|')[0]
        anns[alt_text].append(annotation)
    return anns


def kill_repetitive_bnds(sv_id):
    # takes a row and returns true or false on whether it is a repetitive BND event
    # which will be the same for two paired events, so we can exclude the second event in each case
    # this is indicated by a _2 in the ID
    if '_2' in str(sv_id):
        return False
    else:
        return True

def format_sv_alt(row):
    alt = row['ALT']
    if row['SVTYPE'] == 'BND':
        return alt
    else:
        return alt + str(int(row['END']))

# STEP 0: PARSING SV DATA FROM LUMPY
sv_dat = pd.read_csv('../../Output/WGS/smoove_output/all_samples.smoove.square.ann.tsv', delimiter='\t').fillna(0)
sv_data_cols = [i.split('.AO')[0] for i in sv_dat if '.AO' in i]
sv_renamer = {i+'.AO':i+'_alt_counts' for i in sv_data_cols}
sv_renamer.update({i+'.RO':i+'_ref_counts' for i in sv_data_cols})
sv_dat.rename(columns=sv_renamer, inplace=True)
sv_dat['ALT'] = sv_dat.apply(lambda r: format_sv_alt(r), axis=1)
sv_dat['not_bnd_repetitive'] = sv_dat['ID'].apply(lambda i: kill_repetitive_bnds(i))
sv_dat = sv_dat[sv_dat['not_bnd_repetitive']]

# STEP 1: Combine chromosome VCFs and format annotations:
flist = glob(outdir+'chromo_vcfs/*ann.tsv')
base_cols = ['CHROM', 'POS', 'REF', 'QUAL', 'SVTYPE']
dfs = []

# Getting column name info
nd = pd.read_csv(flist[0], delimiter='\t')
other_cols = [i for i in nd if i not in base_cols+['ALT', 'ANN']]
sample_names = [i.split('.AD')[0] for i in other_cols]
wells = sorted(set([i.split('_')[1] for i in sample_names]))
all_cols = base_cols+['ALT', 'ANN']+other_cols
new_col_names = base_cols + ['ALT', 'ANN'] + [i+'_ref_counts' for i in sample_names] + [i+'_alt_counts' for i in sample_names]

# Ok now looping through files, splitting multi-allelic records
for fin in flist:
    nd = pd.read_csv(fin, delimiter='\t')
    nd['SVTYPE'] = ['NA']*len(nd)
    new_rows = []
    for row in nd.as_matrix(all_cols):
        alt = row[all_cols.index('ALT')]
        if ',' not in alt: # only one alternate allele here
            new_rows.append(list(row[:len(base_cols)+2]) + [int(row[c].split(',')[0]) for c in range(len(base_cols)+2,len(all_cols))] + 
                [int(row[c].split(',')[1]) for c in range(len(base_cols)+2,len(all_cols))])
        else:
            ad = annotation_splitter(row[all_cols.index('ANN')]) # this yeilds a dictionary from alternate alleles to their annotations
            alt_split = alt.split(',')
            base_row = list(row[:len(base_cols)])
            for i in range(len(alt_split)): # for each alternate allele
                current_alt = alt_split[i]
                current_ann = ','.join(ad[current_alt])
                new_rows.append(base_row + [current_alt, current_ann] + [int(row[c].split(',')[0]) for c in range(len(base_cols)+2,len(all_cols))] + 
                [int(row[c].split(',')[i+1]) for c in range(len(base_cols)+2,len(all_cols))])
    dfs.append(pd.DataFrame(new_rows, columns=new_col_names))

nd = pd.concat(dfs+[sv_dat[new_col_names]])

# STEP 2: Tallying up total alternate allele counts
# Tallying up all alt allele counts
nd['total_alt_counts'] = np.nansum(nd[[sn + '_alt_counts' for sn in sample_names]], axis=1)
# Tallying up all Gen 70 alt allele counts
nd['all_G70_alt_counts'] = np.nansum(nd[['G70_' + w + '_alt_counts' for w in wells]], axis=1)
# This is a first filter - every population will need to have at least this many reads, so if they all don't we're totally out of luck
nd = nd[nd['total_alt_counts']>=alt_allele_thresh]

for w in wells:
    acols = [i for i in nd if w in i and '_alt_counts' in i]
    nd[w + '_total_alt_counts'] = np.sum(nd[acols], axis=1)

# Outputting large variant table
nd.to_csv(outdir+'allele_table_expanded.tsv', index=False, sep='\t')
print('done w step 2')

# STEP 3: Filtering and outputting per well
# Jumpstart if the first part only ran
#nd = pd.read_csv(outdir+'allele_table_expanded.tsv', delimiter='\t')

gens = [70, 1410, 2640, 5150, 7530, 10150]
extra_cols = ['af_trajectory', 'well_total_alt_counts', 'total_alt_counts', 'all_G70_alt_counts', 'perc_of_alt', 'filter_1', 'filter_2', 'tps_w_at_least_5']
for w in wells:
    well_cols = [i for i in nd if w in i]
    gens_use = [g for g in gens if 'G' + str(g) + '_' + w + '_alt_counts' in well_cols]
    td = nd[nd[w + '_total_alt_counts']>=alt_allele_thresh].rename(columns={i:i.replace('_'+w, '').replace(w, 'well') for i in well_cols})
    td['af_trajectory'] = td.apply(lambda row: get_af_traj(row, gens_use), axis=1)
    td['perc_of_alt'] = td['well_total_alt_counts'] / td['total_alt_counts']
    td['tps_w_at_least_5'] = td.apply(lambda row: num_tps_thresh(row, 5, gens_use), axis=1)
    td['filter_1'] = (td['perc_of_alt']>0.9)
    td['filter_2'] = ((td['well_total_alt_counts']/(td['all_G70_alt_counts']+td['well_total_alt_counts']-td['G70_alt_counts']) > 0.9) & (td['tps_w_at_least_5']>=2))
    tdfilt = td[td['filter_1'] | td['filter_2']]
    use_cols = base_cols + ['ALT', 'ANN'] + extra_cols + ['G'+str(g)+'_ref_counts' for g in gens_use] + ['G'+str(g)+'_alt_counts' for g in gens_use]
    td[use_cols].to_csv(outdir+'/well_output/' + w + '_full.tsv', index=False, sep='\t')
    tdfilt[use_cols].to_csv(outdir+'/well_output/' + w + '_filtered.tsv', index=False, sep='\t')