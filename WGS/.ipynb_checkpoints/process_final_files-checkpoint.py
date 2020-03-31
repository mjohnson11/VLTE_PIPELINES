import pandas as pd
import numpy as np

def unique_allele_test(row):
    if row['total_total_alt_allele_counts'] > 0:
        top_well = sorted(wells, key=lambda w: row[w + '_total_alt_allele_counts'])[-1]
        if row[top_well + '_total_alt_allele_counts'] / row['total_total_alt_allele_counts'] > 0.9:
                return top_well
    return 'No'

def get_allele_counts(x):
    # just picking out the AD portion of the column
    s = x.split(':')
    if len(s) > 1:
        return s[1]
    else:
        return '0'
    
def get_ref_allele_counts(x):
    return int(x.split(',')[0])

def get_alt_allele_counts(x):
    return sum([int(i) for i in x.split(',')[1:]])

def get_af_traj(r, gu):
    # row, gens_use
    afs = []
    for g in gu:
        if r['G' + str(g) + '_tot_counts'] > 0:
            afs.append(str(g) + '_' + str(r['G'+str(g)+'_alt_allele_counts'] / r['G'+str(g)+'_tot_counts']))
    return ';'.join(afs)


nd = pd.read_csv('../../Output/WGS/FINAL_FULL_VARIANTS.tsv', delimiter='\t')
new_cols = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'ANN', 'FORMAT']
wells = sorted(list(set([i.split('_')[-1] for i in nd if i not in new_cols])))
print('file read')

for i in [j for j in nd.columns if j not in new_cols]:
    nd[i + '_allele_counts'] = nd[i].apply(lambda x: get_allele_counts(x))
    nd[i + '_ref_allele_counts'] = nd[i + '_allele_counts'].apply(lambda x: get_ref_allele_counts(x))
    nd[i + '_alt_allele_counts'] = nd[i + '_allele_counts'].apply(lambda x: get_alt_allele_counts(x))
print('1')
for w in wells:
    acols = [i for i in nd if w in i and '_alt_allele_counts' in i]
    nd[w + '_total_alt_allele_counts'] = np.sum(nd[acols], axis=1)
print('2')
nd['total_total_alt_allele_counts'] = np.nansum(nd[[w + '_total_alt_allele_counts' for w in wells]], axis=1)
nd['unique_allele'] = nd.apply(lambda r: unique_allele_test(r), axis=1)
nd.to_csv('../../Output/WGS/allele_table_expanded.tsv', index=False, sep='\t')
print('3')

nd = pd.read_csv('../Final_output/allele_table_expanded.tsv', delimiter='\t')
gens = [70, 1410, 2640, 5150, 7530, 10150]
for w in wells:
    well_cols = [i for i in nd if w in i]
    gens_use = [g for g in gens if 'G' + str(g) + '_' + w in well_cols]
    td = nd[nd[w + '_total_alt_allele_counts']>2][new_cols + well_cols + ['total_total_alt_allele_counts']].rename(columns={i:i.replace('_'+w, '') for i in well_cols})
    td['percentage_of_total_alt_counts'] = td[w + '_total_alt_allele_counts'] / td['total_total_alt_allele_counts']
    for g in gens_use:
        td['G' + str(g) + '_tot_counts'] = td['G'+str(g)+'_alt_allele_counts'] + td['G'+str(g)+'_ref_allele_counts']
        td['G' + str(g) + '_af'] = td['G'+str(g)+'_alt_allele_counts'] / td['G' + str(g) + '_tot_counts']
    td['af_trajectory'] = td.apply(lambda row: get_af_traj(row, gens_use), axis=1)
    use_cols = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'ANN', 'af_trajectory'] + ['G'+str(g)+'_allele_counts' for g in gens_use]
    td[use_cols + ['percentage_of_total_alt_counts']].to_csv('../../Output/WGS/well_output/' + w + '_full.tsv', index=False, sep='\t')
    td[td['percentage_of_total_alt_counts'] > 0.9][use_cols].to_csv('../../Output/WGS/well_output/' + w + '_unique90.tsv', index=False, sep='\t')
        
