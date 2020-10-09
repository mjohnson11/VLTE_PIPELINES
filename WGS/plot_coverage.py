import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pl
import pandas as pd
import seaborn as sns
import numpy as np
import argparse
cbs = sns.color_palette('colorblind')

parser = argparse.ArgumentParser()
parser.add_argument('whichplate', help='which plate to work on')
args = parser.parse_args()
whichplate = args.whichplate

outdir = '../../Output/WGS/coverage_output/'
figout = '../../Output/Browser/coverage/'

sns.set_style("white")
window_size = 500

with open('../accessory_files/Samples.txt', 'r') as infile:
    samps = [l.strip() for l in infile.readlines()]
wells = sorted(set([i.split('_')[1] for i in samps]))
seq_gens = [70, 1410, 2640, 5150, 7530, 10150]
chromos = ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII', 'chrVIII','chrIX', 
           'chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI', 'chrMito', '2-micron']
colorss = ['k', '#648FFF']

use_wells = [w for w in wells if w[1]==whichplate]
for well in use_wells:
    print(well)
    fig, these_subs = pl.subplots(6, 1, figsize=(13, 5), dpi=100)
    fig_2mic, these_subs_2mic = pl.subplots(6, 1, figsize=(13, 5), dpi=100)
    td = pd.read_csv(outdir+well+'_depth_summary.tsv', delimiter='\t')
    
    for g in range(len(seq_gens)):
        gen = seq_gens[g]
        offsets = [0]
        if 'G'+str(gen)+'_'+well in samps:
            samp = 'G'+str(gen)+'_'+well
            # plotting per generation plots
            f, subps = pl.subplots(4, 4, figsize=(13, 5), dpi=100)
            pl.subplots_adjust(hspace=0.7)
            subs = [subps[i][j] for i in range(4) for j in range(4)]
            for c in range(len(chromos)-2):
                chromo = chromos[c]
                ttd = td[td['CHROM']==chromo]
                subs[c].plot(ttd['window']*(window_size/1000), ttd['standardized_depth_'+str(gen)], lw=0.5, c='#648FFF', alpha=0.6)
                subs[c].plot(ttd['window']*(window_size/1000), ttd['inferred_state_'+str(gen)], lw=1, c='k')
                subs[c].set_ylim([0, 3])
                subs[c].set_yticks([0,1,2,3])
                subs[c].set_title(chromo, y=0.65)
                subs[c].grid(axis='y')
                # plotting combined generation plot
                these_subs[g].plot(ttd['window']+offsets[-1], ttd['standardized_depth_'+str(gen)], lw=0.2, c=colorss[len(offsets)%2], alpha=0.5)
                these_subs[g].plot(ttd['window']+offsets[-1], ttd['inferred_state_'+str(gen)], lw=1, c='k')
                offsets.append(offsets[-1] + max(ttd['window'])+10)
            sns.despine()
            f.savefig(figout+samp+'_depth.png', background='transparent', bbox_inches='tight', pad_inches=0.1)
            ttd = td[td['CHROM']=='2-micron']
            these_subs_2mic[g].plot(ttd['window']+offsets[-1], ttd['standardized_depth_'+str(gen)], lw=2, c=colorss[len(offsets)%2], alpha=1)
            #these_subs_2mic[g].plot(ttd['window']+offsets[-1], ttd['inferred_state_'+str(gen)], lw=1, c='k')
            these_subs[g].set_title(samp.replace('_', ' '), y=0.7)
            these_subs_2mic[g].set_title(samp.replace('_', ' '), y=0.7)
            other_chromo_medians = [np.nanmedian(td[td['CHROM']=='2-micron']['standardized_depth_'+str(gen)]), 
                                     np.nanmedian(td[td['CHROM']=='chrMito']['standardized_depth_'+str(gen)])]
            these_subs[g].scatter([-2200, -800], other_chromo_medians, s=10, c="#444444")
        these_subs[g].set_ylim([0,2.5])
        these_subs[g].set_yticks([0,1,2])
        sns.despine(ax=these_subs[g])
        these_subs[g].grid(axis='y')
        sns.despine(ax=these_subs_2mic[g])
        these_subs_2mic[g].grid(axis='y')
        these_subs[g].set_xticks([])
    spots = [(offsets[i]+offsets[i+1])/2 for i in range(len(offsets)-1)]
    these_subs[g].set_xticks([-2200,-800]+spots)
    these_subs[g].set_xticklabels(['2-micron', 'chrMito']+chromos[:-2], fontsize=7) 
        
    fig.savefig(figout+well+'_allgens_depth.png', background='transparent', bbox_inches='tight', pad_inches=0.1)
    fig_2mic.savefig(figout+well+'_allgens_depth_2micron.png', background='transparent', bbox_inches='tight', pad_inches=0.1)
    pl.close("all")
    
# Part 2 plotting putative CNVs
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
figout2 = outdir + 'coverage_check_plots/'
all_cnvs = pd.read_csv(outdir+'_all_cnvs_to_check.csv')
all_cnvs['Window_ID'] = all_cnvs.apply(lambda r: r['CHROM']+'_'+str(int(r['CNV_start']))+'_'+str(int(r['CNV_end'])), axis=1)

def plot_window(row):
    well, chromo, s, e = row['Well'], row['CHROM'], row['CNV_start'], row['CNV_end']
    window_len = np.max([e-s, 4000])
    fig, these_subs = pl.subplots(6, 1, figsize=(13, 5), dpi=100, sharex=True, sharey=True)
    td = pd.read_csv(outdir+row['Well']+'_depth_summary.tsv', delimiter='\t')
    ttd = td[(td['CHROM']==chromo) & (td['window']*500 > s-window_len) & (td['window']*500 < e+window_len)]
    for g in range(len(seq_gens)):
        gen = seq_gens[g]
        if 'standardized_depth_'+str(gen) in ttd:
            these_subs[g].plot(ttd['window']*500, ttd['relative_depth_'+str(gen)], lw=1, c=cbs[0])
            these_subs[g].plot(ttd['window']*500, ttd['inferred_state_'+str(gen)], lw=1, c='k')
    sns.despine()
    these_subs[g].grid(axis='y')
    fig.suptitle(well + ' ' + row['Window_ID'])
    fig.savefig(figout2+row['Window_ID']+'_'+well+'.png', background='transparent', bbox_inches='tight', pad_inches=0.1)
    pl.close("all")

for j, row in all_cnvs[all_cnvs['Well'].apply(lambda w: w[1]==whichplate)].iterrows():
    plot_window(row)