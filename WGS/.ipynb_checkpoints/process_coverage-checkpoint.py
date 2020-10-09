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

## FOR CALLING MUTATATIONS IN TELOMERES

chromo_lens = {
    'chrI': 230218,
    'chrII': 813184,
    'chrIII': 316620,
    'chrIV': 1531933,
    'chrV': 576874,
    'chrVI': 270161,
    'chrVII': 1090940,
    'chrVIII': 562643,
    'chrIX': 439888,
    'chrX': 745751,
    'chrXI': 666816,
    'chrXII': 1078177,
    'chrXIII': 924431,
    'chrXIV': 784333,
    'chrXV': 1091291,
    'chrXVI': 948066,
}

romans = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15, 'XVI': 16}

telo_lens = {i[0]: i[1] for i in np.array(pd.read_csv('../accessory_files/yeast_telomere_lengths.tsv', delimiter='\t', header=None))}

exclude_regions = { # excluding regions with very frequent CNVs, which will be looked at in every population separately
    'chrVIII': (210000, 220000), # CUP1-2 
    'chrXII': (450000, 475000), # rDNA
    'chrIV': (1150000, 1170000), # HXTs
}

def same_window(s1, e1, s2, e2):
    threshold = np.min([np.mean([(e1-s1)+1,(e2-s2)+1])/4, 5])
    if np.abs(e1-e2) < threshold and np.abs(s1-s2) < threshold:
        return True
    else:
        return False
    
def likely_good(row):
    # Includes CNVs if they are at least 4 windows (2kb) big AND appear in multiple generations
    if (';' in row['gens']) and (row['length'] > 3):
        return True
    return False

def make_cnv_file(well, td, outfile, use_gens):
    calls = [] # like [chromo, gen, start, end, length, state, description]
    for chromo in chromo_lens:
        # getting telomere regions, with a 1kb buffer
        left_edge = telo_lens['TEL'+str(romans[chromo[3:]]).zfill(2)+'L']
        right_edge = chromo_lens[chromo] - telo_lens['TEL'+str(romans[chromo[3:]]).zfill(2)+'R']
        for gen in use_gens:
            w_to_s = {i[0]: i[1] for i in np.array(td[td['CHROM']==chromo][['window', 'inferred_state_'+str(gen)]])}
            winds = sorted([w for w in w_to_s if pd.notnull(w_to_s[w])])
            # jump through blocks of states
            for k, g in itertools.groupby(winds, key=w_to_s.get):
                if k != 1:
                    tmp = list(g)
                    if not (tmp[-1]*window_size < (left_edge+2000) or tmp[0]*window_size > (right_edge-2000)): # excluding those IN telomeres
                        if not (tmp[0]*window_size < (left_edge+2000) or tmp[-1]*window_size > (right_edge-2000)): # marking those NEAR telomeres
                            near_telomere = False
                        else:
                            near_telomere = True
                        # excluding those in common regions (analyzed separately)
                        common = False
                        if chromo in exclude_regions:
                            if (tmp[-1]*window_size < (exclude_regions[chromo][1]+2000) or tmp[0]*window_size > (exclude_regions[chromo][0]-2000)):
                                common = True
                        if not common:
                            repeat_call = False
                            for call in calls:
                                if call[0] == chromo:
                                    if same_window(call[2], call[3], tmp[0], tmp[-1]):
                                        call[1] += ';'+str(gen)
                                        repeat_call = True
                                        break
                            if not repeat_call:
                                calls.append([chromo, str(gen), tmp[0], tmp[-1], len(tmp), k, near_telomere])
    if len(calls) > 0:
        tmp = pd.DataFrame([[well] + c for c in calls], columns = ['Well', 'CHROM', 'gens', 'start', 'end', 'length', 'state', 'near_telomere'])
        tmp['likely_good'] = tmp.apply(lambda r: likely_good(r), axis=1)
        tmp['CNV_start'] = tmp['start']*window_size
        tmp['CNV_end'] = tmp['end']*window_size
        tmp[tmp['likely_good']][['Well', 'CHROM', 'CNV_start', 'CNV_end', 'gens', 'state', 'near_telomere']].sort_values(by=['CHROM', 'CNV_start']).to_csv(outfile, sep='\t', index=False)
        return tmp[tmp['likely_good']][['Well', 'CHROM', 'CNV_start', 'CNV_end', 'gens', 'state', 'near_telomere']], True
    else:
        return None, False

with open('../accessory_files/Samples.txt', 'r') as infile:
    samps = [l.strip() for l in infile.readlines()]
    
wells = sorted(set([i.split('_')[1] for i in samps]))
seq_gens = [70, 1410, 2640, 5150, 7530, 10150]
all_cnv_dats = []
for well in wells:
    used_gens = [gen for gen in seq_gens if 'G'+str(gen)+'_'+well in samps]
    running_dat = pd.read_csv(outdir+well+'_depth_summary.tsv', delimiter='\t')
    #for g in used_gens:
    #    running_dat['state_is_one_'+str(g)] = running_dat['inferred_state_'+str(g)]==1
    td, worked = make_cnv_file(well, running_dat, outdir+well+'_cnv_predictions.tsv', used_gens)
    if worked:
        print(well, len(td))
        all_cnv_dats.append(td)
    
pd.concat(all_cnv_dats).to_csv(outdir+'all_cnvs_to_check.csv', index=False)