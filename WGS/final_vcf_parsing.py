from glob import glob
import pandas as pd
import numpy as np
from collections import Counter

alt_allele_thresh = 5

def get_allele_counts(x):
    # just picking out the AD portion of the column
    s = x.split(':')
    if len(s) > 1:
        return s[1]
    else:
        return '0'

def get_alt_allele_counts(x):
    return sum([int(i) for i in x.split(',')[1:]])

def get_af_traj(r, gu):
    # row, gens_use
    afs = []
    for g in gu:
        col = 'G' + str(g) + '_allele_counts'
        if col in r:
            ref, alt = int(r[col].split(',')[0]), int(r[col].split(',')[1])
            tot = ref + alt
            if tot >= 5:
                afs.append(str(g) + '_' + str(alt/tot))
    return ';'.join(afs)

def num_tps_thresh(r, thresh):
    c = 0
    for col in ['G'+str(g)+'_allele_counts' for g in gens]:
        if col in r:
            if int(r[col].split(',')[1]) >= thresh:
                c += 1
    return c

def annotation_parser(a):
    # input is a single annotation
    sa = a.split('|')
    allele, effs, gene_hit = sa[0], sa[1], sa[3]
    extra_ann = list(orf2ann.setdefault(gene_hit, ['']*2))
    # excluding intergenic or downstream annotations
    if effs not in ['intergenic_region', 'downstream_gene_variant']:
        if 'upstream' in effs: # if upstream, must be within 100 bp of the feature to keep it
            distance_to_feature = sa[14]
            if distance_to_feature != '':
                if int(distance_to_feature) < 100:
                    return [allele, effs, gene_hit, distance_to_feature] + extra_ann
        else:
            return [allele, effs, gene_hit, ''] + extra_ann
    return []
                    
def format_annotations(info):
    all_anns = []
    for entry in info.split('ANN=')[-1].split(';')[0].split(','):
        tmp = annotation_parser(entry)
        if len(tmp) > 0:
            all_anns.append('|'.join(tmp))
    return ';'.join(all_anns)

def get_sv_allele_counts(v):
    allele_counts = []
    tmp = str(v).split(':')
    for spot in [5, 6]:
        if len(tmp) > spot:
            try:
                allele_counts.append(int(tmp[spot]))
            except ValueError:
                allele_counts.append(0)   
        else:
            allele_counts.append(0)
    return ','.join([str(i) for i in allele_counts])


orig_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
new_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'ANN', 'FILTER', 'INFO', 'FORMAT']
# PARSING SV DATA FROM LUMPY
sv_dat = pd.read_csv('../../Output/WGS/smoove_output/all_samples.smoove.square.vcf', delimiter='\t', skiprows=1159).rename(columns={'#CHROM': 'CHROM'})
assert list(sv_dat['FORMAT'])[0].split(':')[5] == 'RO'
assert list(sv_dat['FORMAT'])[0].split(':')[6] == 'AO'
sv_data_cols = [i for i in sv_dat if i not in new_cols]
for c in sv_data_cols:
    sv_dat[c+'_allele_counts'] = sv_dat[c].apply(lambda val: get_sv_allele_counts(val))
    sv_dat[c+'_alt_allele_counts'] = sv_dat[c + '_allele_counts'].apply(lambda x: get_alt_allele_counts(x))
# Don't have SNPEFF annotation for SVs
sv_dat['ANN'] = ['SV']*len(sv_dat)


# STEP 1: Combine chromosome VCFs and format annotations:
orf2ann = {i[0]: [str(j) for j in i[1:]] for i in pd.read_csv('../accessory_files/yeast_gene_annotations.csv').as_matrix(['ORF', 'Gene', 'briefDescription'])}
flist = glob('../../Output/WGS/chromo_vcfs/*ann.vcf')
nd = pd.concat([pd.read_csv(fin, skiprows=62, delimiter='\t') for fin in flist])
nd['ANN'] = nd['INFO'].apply(lambda i: format_annotations(i))
nd = nd.rename(columns={'#CHROM': 'CHROM'})
nd = nd[new_cols + [i for i in nd if i not in orig_cols and i not in new_cols]]
nd.to_csv('../../Output/WGS/FINAL_FULL_VARIANTS.tsv', index=False, sep='\t')
print("done w step 1")

# STEP 2: Adding filtering columns and splitting multi-allelic records, also adding in SV stuff
# Jumpstart if step 1 is done
#nd = pd.read_csv('../../Output/WGS/FINAL_FULL_VARIANTS.tsv', delimiter='\t')

# check that format is right, then lose the column
assert list(nd['FORMAT'])[0].split(':')[1] == 'AD'
data_cols = [j for j in nd if j not in new_cols]
wells = sorted(list(set([i.split('_')[-1] for i in data_cols])))
print('file read')
for i in data_cols:
    nd[i + '_allele_counts'] = nd[i].apply(lambda x: get_allele_counts(x))
    nd[i + '_alt_allele_counts'] = nd[i + '_allele_counts'].apply(lambda x: get_alt_allele_counts(x))

# Checking that adding SV data will work
assert Counter(data_cols) == Counter(sv_data_cols)
these_columns = new_cols+[i + '_allele_counts' for i in data_cols]+[i + '_alt_allele_counts' for i in data_cols]
nd = pd.concat([nd[these_columns], sv_dat[these_columns]])

# Tallying up all alt allele counts
nd['total_total_alt_allele_counts'] = np.nansum(nd[[i + '_alt_allele_counts' for i in data_cols]], axis=1)
# Tallying up all Gen 70 alt allele counts
nd['all_G70_alt_allele_counts'] = np.nansum(nd[['G70_' + w + '_alt_allele_counts' for w in wells]], axis=1)
# This is a first filter - every population will need to have at least this many reads, so if they all don't we're totally out of luck
nd = nd[nd['total_total_alt_allele_counts']>=alt_allele_thresh]

# Splitting multi-allelic records into separate rows
big_mat = []
new_cols_plus = new_cols + ['total_total_alt_allele_counts', 'all_G70_alt_allele_counts']
for row in nd.as_matrix(new_cols_plus + [c + '_allele_counts' for c in data_cols]):
    alts = str(row[4])
    for i in range(len(alts.split(','))):
        row = list(row)
        ac, refc, altc = [], [], []
        for r in row[12:]:
            rs = r.split(',')
            ac.append(rs[0]+','+rs[i+1])
            refc.append(int(rs[0]))
            altc.append(int(rs[i+1]))
        big_mat.append(row[:4] + [alts.split(',')[i]] + row[5:12] + ac + refc + altc)
nd = pd.DataFrame(big_mat, columns=new_cols_plus + [c + '_allele_counts' for c in data_cols] + [c + '_ref_allele_counts' for c in data_cols] + [c + '_alt_allele_counts' for c in data_cols])

for w in wells:
    acols = [i for i in nd if w in i and '_alt_allele_counts' in i]
    nd[w + '_total_alt_allele_counts'] = np.sum(nd[acols], axis=1)

# Outputting large variant table
nd.to_csv('../../Output/WGS/allele_table_expanded.tsv', index=False, sep='\t')
print('done w step 2')

# STEP 3: Filtering and outputting per well
# Jumpstart if the first part only ran
#nd = pd.read_csv('../../Output/WGS/allele_table_expanded.tsv', delimiter='\t')

gens = [70, 1410, 2640, 5150, 7530, 10150]
extra_cols = ['af_trajectory', 'well_total_alt_allele_counts', 'total_total_alt_allele_counts', 'all_G70_alt_allele_counts', 'perc_of_alt', 'filter_1', 'filter_2', 'tps_w_at_least_5']
for w in wells:
    well_cols = [i for i in nd if w in i]
    gens_use = [g for g in gens if 'G' + str(g) + '_' + w + '_allele_counts' in well_cols]
    td = nd[nd[w + '_total_alt_allele_counts']>=alt_allele_thresh].rename(columns={i:i.replace('_'+w, '').replace(w, 'well') for i in well_cols})
    td['af_trajectory'] = td.apply(lambda row: get_af_traj(row, gens_use), axis=1)
    td['perc_of_alt'] = td['well_total_alt_allele_counts'] / td['total_total_alt_allele_counts']
    td['tps_w_at_least_5'] = td.apply(lambda row: num_tps_thresh(row, 5), axis=1)
    td['filter_1'] = (td['perc_of_alt']>0.9)
    td['filter_2'] = ((td['well_total_alt_allele_counts']/(td['all_G70_alt_allele_counts']+td['well_total_alt_allele_counts']-td['G70_alt_allele_counts']) > 0.9) & (td['tps_w_at_least_5']>=2))
    tdfilt = td[td['filter_1'] | td['filter_2']]
    use_cols = new_cols + extra_cols + ['G'+str(g)+'_allele_counts' for g in gens_use]
    td[use_cols].to_csv('../../Output/WGS/well_output/' + w + '_full.tsv', index=False, sep='\t')
    tdfilt[use_cols].to_csv('../../Output/WGS/well_output/' + w + '_filtered.tsv', index=False, sep='\t')