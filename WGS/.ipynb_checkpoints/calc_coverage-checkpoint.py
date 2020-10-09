import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from hmmlearn import hmm
import itertools

outdir = '../../Output/WGS/coverage/'
window_size=500

# Step 1 STANDARDIZE depths
# a) calculate relative (to whole genome coverage) coverage in 500 bp windows for all gen 70 samples
# b) for each strain, find the median coverage in each window, store this information
# c) calculate relative (to whole genome coverage) coverage in 500 bp windows for all samples, and add a column corrected by median gen 70 coverage
# d) output this file (well_depth_summary.tsv)
# e) output windows of inferred non-1 states (well_cnv_predictions.tsv)
# f) plot in another script


def get_window_depth(fin):
    d = pd.read_csv(fin, delimiter='\t', header=None, names=['CHROM', 'POS', 'DEPTH'])
    d['window'] = d['POS']//window_size
    td = d[['CHROM', 'window', 'DEPTH']].groupby(['CHROM', 'window'], as_index=False).sum()
    td['relative_depth'] = td['DEPTH']/np.median(td['DEPTH'])
    return td

def get_stand_depth(row, sd, sc):
    s_depth = sd[row['CHROM']+'_'+str(row['window'])]
    if s_depth/sc < 0.25:
        return np.nan
    else:
        return s_depth

with open('../accessory_files/Samples.txt', 'r') as infile:
    samps = [l.strip() for l in infile.readlines()]
wells = sorted(set([i.split('_')[1] for i in samps]))
seq_gens = [70, 1410, 2640, 5150, 7530, 10150]
chromos = ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII', 'chrVIII','chrIX', 'chrX', 
           'chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI', 'chrMito', '2-micron']

wellinfo = pd.read_csv('../accessory_files/VLTE_by_well_info.csv')[['plate.well', 'contam', 'strain']]
wellinfo['plate.well'] = wellinfo['plate.well'].apply(lambda p: p[:2]+p[3:]) #reformatting to match for merge
platewell_to_strain = {i[0]: i[1] for i in np.array(wellinfo[['plate.well', 'strain']])}

strains = ['diploid', 'a', 'alpha']

strain_cov_info = {s: Counter() for s in strains}
strain_counts = Counter()

for well in wells:
    strain = platewell_to_strain[well]
    print(well, strain)
    td = get_window_depth(outdir+'G70_'+well+'_depth.tsv')
    strain_counts[strain] += 1
    for entry in np.array(td[['CHROM', 'window', 'relative_depth']]):
        if entry[2] != 0:
            strain_cov_info[strain][entry[0]+'_'+str(entry[1])] += entry[2]

for well in wells:
    strain = platewell_to_strain[well]
    sd = strain_cov_info[strain]
    sc = strain_counts[strain]
    print(well, strain)
    dats = []
    used_gens = []
    for g in range(len(seq_gens)):
        gen = seq_gens[g]
        if 'G'+str(gen)+'_'+well in samps:
            samp = 'G'+str(gen)+'_'+well
            td = get_window_depth(outdir+samp+'_depth.tsv').fillna(0)
            td['standardized_depth'] = td.apply(lambda row: sc*row['relative_depth']/get_stand_depth(row, sd, sc), axis=1)
            td = td[pd.notnull(td['standardized_depth'])]
            # VERY ROUGH HMM FOR FINDING CANDIDATE CNVs
            X = np.array(td['standardized_depth'])
            dvar = np.nanvar(X)
            if dvar > 0:
                model = hmm.GaussianHMM(n_components=7, n_iter=100, params='')
                model.startprob_ = np.array([0.01, 0.01, 0.94, 0.01, 0.01, 0.01, 0.01])
                ## HMM transition matrix
                tp = 0.0001
                tm = np.ones((7,7))*tp
                for i in range(len(tm)):
                    tm[i][i] = 1-6*tp
                model.transmat_ = tm
                model.means_ = np.array([[0], [0.5], [1], [1.5], [2], [3], [4]])
                state_to_mean = {0: 0, 1: 0.5, 2: 1, 3: 1.5, 4: 2, 5: 3, 6: 4}
                # variances scale with state except 0 has half the variance of 1 so it is possible
                model.covars_ = np.array([[0.05*dvar], [0.05*dvar], [1*dvar], [1.5*dvar], [2*dvar], [3*dvar], [4*dvar]])

                inferred_states = dict()
                for chromo in chromos:
                    ttd = td[td['CHROM']==chromo].sort_values(by='window')
                    ttd = ttd[pd.notnull(ttd['standardized_depth'])]
                    window_list = list(ttd['window'])
                    Xtmp = np.reshape(np.array(ttd['standardized_depth']), (len(ttd), 1))
                    hmm_result = model.predict(Xtmp)
                    inferred_states[chromo] = {window_list[i]: hmm_result[i] for i in range(len(window_list))}
                td['inferred_state'] = td.apply(lambda row: state_to_mean.get(inferred_states[row['CHROM']].get(row['window'], np.nan), np.nan), axis=1)
            else:
                td['inferred_state'] = [np.nan]*len(td)
            # AND FINALLY APPEND DATA TO LIST FOR ONE WELL OUTPUT
            dats.append(td)
            used_gens.append(str(gen))
            
    running_dat = dats[0].merge(dats[1], on=['CHROM', 'window'], how='outer', suffixes=('_'+used_gens[0], ''))
    for g in range(2,len(used_gens)-1):
        running_dat = running_dat.merge(dats[g], on=['CHROM', 'window'], how='outer', suffixes=('_'+used_gens[g-1],''))
    running_dat = running_dat.merge(dats[-1], on=['CHROM', 'window'], how='outer', suffixes=('_'+used_gens[-2], '_'+used_gens[-1]))
    running_dat.sort_values(by=['CHROM', 'window']).to_csv(outdir+well+'_depth_summary.tsv', sep='\t', index=False)