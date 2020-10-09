import pandas as pd
import numpy as np
from scipy import stats as sci_stats
from glob import glob
from collections import defaultdict, Counter
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_gaf
from statsmodels.stats.multitest import fdrcorrection as benjamini_hochberg
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('outdir', help='directory where chromo dbs etc and future output will go')
args = parser.parse_args()
outdir = args.outdir

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

def in_telomere(row):
    if row['CHROM'] not in ['2-micron', 'chrMito']:
        left_edge = telo_lens['TEL'+str(romans[row['CHROM'][3:]]).zfill(2)+'L']
        right_edge = chromo_lens[row['CHROM']] - telo_lens['TEL'+str(romans[row['CHROM'][3:]]).zfill(2)+'R']
        if row['POS'] < left_edge or row['POS'] > right_edge:
            return True
    return False

## FUNCTIONS FOR GETTING THE # NONSYN OPPORTUNITIES FOR MULTIPLICITY

nt2codon = {
    'TTT': 'F', 'TTC': 'F',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'TAT': 'Y', 'TAC': 'Y',
    'TAA': '*', 'TAG': '*', 'TGA': '*',
    'TGT': 'C', 'TGC': 'C',
    'TGG': 'W',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def get_attrib(row, attrib):
    if row['type']=='gene':
        if attrib+'=' in row['attributes']:
            return row['attributes'].split(attrib+'=')[1].split(';')[0]
    return ''

def read_fasta(fasta_file):
    """
    Reads a fasta file and returns a dictionary with seqid keys and sequence values
    """
    fd = dict()
    with open(fasta_file, 'r') as infile:
        for line in infile:
            if '>' in line:
                current_key = line[1:].strip()
                fd[current_key] = ''
            else:
                fd[current_key] += line.strip()
    return fd

def reverse_transcribe(seq):
    """reverse transcribes a dna sequence (does not convert any non-atcg/ATCG characters)"""
    watson_crick = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join([watson_crick.setdefault(c, c) for c in seq[::-1]])

class SeqInfoGetter:
    
    def __init__(self, gff_file, fasta_file):
        gff_cols = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        self.gff = pd.read_csv(gff_file, delimiter='\t', skiprows=1, header=None, names=gff_cols)
        self.gff['ORF'] = self.gff.apply(lambda row: get_attrib(row, "ID"), axis=1)
        self.genes = self.gff[self.gff['ORF']!='']
        self.genes['Gene'] = self.genes.apply(lambda row: get_attrib(row, "gene"), axis=1)
        self.chromo_seqs = read_fasta(fasta_file)
    
    def get_nt_seq(self, element_name, element_type):
        td = self.genes[self.genes[element_type]==element_name]
        if len(td) != 1:
            print(len(td), 'hits, aborting.')
            return None
        else:
            row = td.iloc[0]
            cs = self.chromo_seqs[row['seqid']]
            if row['strand'] == '+':
                return cs[row['start']-1:row['end']]
            else:
                return reverse_transcribe(cs[row['start']-1:row['end']])
            
    def get_aa_seq(self, element_name, element_type):
        nt_s = self.get_nt_seq(element_name, element_type)
        if nt_s:
            aas = ''
            for i in range(len(nt_s)//3):
                aas += nt2codon[nt_s[i*3:(i+1)*3]]
            if len(nt_s) % 3 != 0:
                aas += '-leftover->' + nt_s[-1*(len(nt_s) % 3):]
            return aas
            
    def get_mutational_opps2(self, element_name, element_type, verbose=False):
        nt_s = self.get_nt_seq(element_name, element_type)
        if nt_s:
            if len(nt_s) % 3 != 0:
                if verbose:
                    print('Warning: seq len not a multiple of 3', element_name)
                    print(self.genes[self.genes[element_type]==element_name].iloc[0]['Gene'])
                    print(self.get_aa_seq(element_name, element_type))
           
            syn, nonsyn = 0, 0
            for i in range(len(nt_s)//3):
                codon_seq = nt_s[i*3:(i+1)*3]
                codes_for = nt2codon[codon_seq]
                for j in range(3):
                    for nt in 'ATCG':
                        if nt != codon_seq[j]:
                            if nt2codon[codon_seq[:j]+nt+codon_seq[j+1:]] == codes_for:
                                syn += 1
                            else:
                                nonsyn += 1
            return nonsyn / (syn+nonsyn)

def genes_affected(a):
    keep_anns = []
    for ann in str(a).split(','):
        if '|' in ann:
            sa = ann.split('|')
            for mut_type in sa[1].split('&'):
                if sa[1] not in ['SV', 'synonymous', 'noncoding']:
                    if 'Dubious' not in sa[5]:
                        keep_anns.append(sa[4])
    return ';'.join(keep_anns)

def genes_affected_syn(a):
    keep_anns = []
    for ann in str(a).split(','):
        if '|' in ann:
            sa = ann.split('|')
            for mut_type in sa[1].split('&'):
                if sa[1] == 'synonymous':
                    if 'Dubious' not in sa[5]:
                        keep_anns.append(sa[4 ])
    return ';'.join(keep_anns)

def simplify_mutation_type(m):
    mut_type_list = [mu_simplify[i] for i in str(m).split('&')]
    mut_types_in_order = ['SV','indel','nonsense','missense','synonymous','noncoding']
    mut_types_here = [mt for mt in mut_types_in_order if mt in mut_type_list]
    if len(mut_types_here) > 0:
        # take the first one (they are arranged in precedence order)
        return mut_types_here[0]
    else:
        print('ok')
        return 'NA'

def simple_ann(a):
    sa = a.split('|')
    # format: ALT|simple mutation type|Putative_impact|Gene/ORF name (to be used)|ORF name|proteinAA position|briefDescription
    if len(sa) > 13:
        aa_pos = sa[13]
    else:
        aa_pos = ''
    return '|'.join([str(i).replace(',', '_') for i in [sa[0], simplify_mutation_type(sa[1]), sa[2], o2g.get(sa[3], [sa[3]])[0], sa[3], aa_pos, o2g.get(sa[3], [None, sa[3]])[1]]])

def annotation_parser(row):
    # input is annotation column
    # goal is to remove annotations that are "boring" - intergenic or more than 100 bp upstream
    # simplifies annotations to the format: ALT|simple mutation type|Putative_impact|Gene/ORF name (to be used)|ORF name|proteinAA position|briefDescription
    # only preserves the most "important" annotation
    # SVs are a special case, and are not annotated specifically since they may affect many many genes:
    if pd.notnull(row['SVTYPE']):
        return row['ALT'] + '|SV|||||'
    annotations = row['ANN']
    anns = []
    for a in str(annotations).split(','):
        sa = a.split('|')
        if len(sa) > 1:
            effs = sa[1]
            # excluding intergenic or downstream annotations
            if effs not in ['intergenic_region', 'downstream_gene_variant']:
                if 'upstream' in effs and len(sa) > 14: # if upstream, must be within 100 bp of the feature to keep it
                    distance_to_feature = sa[14]
                    if distance_to_feature != '':
                        if int(distance_to_feature) < 100:
                            anns.append(simple_ann(a))
                else:
                    anns.append(simple_ann(a))
    return ','.join(anns)

def fixed_by(r, gen, fixed_freq_thresh, prefix=''):
    # find gen index
    g_start = seq_gens.index(gen)
    # if fixed at the previous generation, automatically fixed now
    if g_start > 0:
        if r[prefix + 'fixed_by_'+str(seq_gens[g_start-1])]:
            return True

    states = []
    for g in seq_gens[g_start:]:
        if 'G'+str(g)+'_alt_counts'  in r:
            ref, alt = r['G'+str(g)+'_ref_counts'], r['G'+str(g)+'_alt_counts']
            tot = ref + alt
            if tot >= 5:
                if alt/tot >= fixed_freq_thresh:
                    states.append('fixed')
                else:
                    states.append('not fixed')
            else:
                states.append('unknown')
    if states[0] == 'fixed' and 'not fixed' not in states:
        return True
    else:
        return False
    
def present_at(r, gen, freq_thresh):
    # enforce that if it is fixed it's present (which might not be true due to low counts otherwise)
    if r['fixed_by_'+str(gen)]:
        return True
    ref, alt = r['G'+str(gen)+'_ref_counts'], r['G'+str(gen)+'_alt_counts']
    tot = ref + alt
    if tot >= 5:
        if alt/tot >= freq_thresh:
            return True
        else:
            return False
    else:
        return False

def test_for_merge(test_row, mg_rows):
    # Tests that a mutation row is the same or different from other rows that are nearby (in the genome)
    # at each time point, using Fisher's exact test
    pvals = []
    for gen in seq_gens:
        if 'G'+str(gen)+'_alt_counts' in test_row:
            ref, alt = test_row['G'+str(gen)+'_ref_counts'], test_row['G'+str(gen)+'_alt_counts']
            mg_ref, mg_alt = 0, 0
            for r in mg_rows:
                mg_ref += r['G'+str(gen)+'_ref_counts']
                mg_alt += r['G'+str(gen)+'_alt_counts']
            pvals.append(sci_stats.fisher_exact([[alt, ref], [mg_alt, mg_ref]])[1])
    return len([i for i in pvals if i<0.01])

def get_mutation_type(ann):
    if len(str(ann).split('|'))>1:
        if str(ann).split('|')[1]=='SV':
            return 'SV'
        splitter = str(ann).split(',')
        for mt in mutation_types_in_consequence_order:
            for s in splitter:
                if len(s.split('|'))>1:
                    if s.split('|')[1] == mt:
                        return mt
    # if none of the listed types, call it noncoding
    return 'noncoding'
    
def get_mutation_impact(ann):
    if len(str(ann).split('|'))>1:
        if str(ann).split('|')[1]=='SV':
            return 'SV'
        splitter = str(ann).split(',')
        for mt in mutation_impacts_in_consequence_order:
            for s in splitter:
                if len(s.split('|'))>2:
                    if s.split('|')[2] == mt:
                        return mt
    # if none of the listed types, call it NA
    return 'NA'

def get_group_mutation_type(row, mg_counts, mg_to_type):
    if mg_counts[row['mutation_group']] == 1:
        return row['mutation_type']
    else:
        all_mut_types_in_group = mg_to_type[row['mutation_group']]
        for mt in mutation_types_in_consequence_order:
            if mt in all_mut_types_in_group:
                return mt
            
def get_group_mutation_impact(row, mg_counts, mg_to_impact):
    if mg_counts[row['mutation_group']] == 1:
        return row['mutation_impact']
    else:
        all_mut_impacts_in_group = mg_to_impact[row['mutation_group']]
        for mt in mutation_impacts_in_consequence_order:
            if mt in all_mut_impacts_in_group:
                return mt
            
### FUNCTIONS FOR PART 2 ###

def get_gene_hit_probs(data, pop_hits):
    # Given a G x N array of G genes and N populations, 
    # with 1's where a population has a gene hit else 0's,
    # returns a matrix of probability of fixing in each population
    # based on the number of mutatations in each population (pop_hits, len N array)
    # and the number of populations with a mutation in that gene
    gene_hits = data.sum(axis=1) + 0.1        # +0.1 pseudocount, pops w/ this gene hit, N_i
    total_hits = np.sum(pop_hits)             # total muts, M
    gene_hit_probs = gene_hits / total_hits   # N_i / M
    # TOO SLOW: prob_no_hits = np.array([[(1-i)**j for j in pop_hits] for i in gene_hit_probs])
    prob_no_hits = np.tile(1-gene_hit_probs, (len(pop_hits), 1)).T**np.tile(pop_hits, (len(gene_hits), 1)) 
    return 1-prob_no_hits

def split_probs(pop_hits, mhd, well_sets):
    # given sets of well indices, computes hit probabilities separately for each set, 
    # then stitches them back together in order
    well_to_probs = dict()
    full_probs = np.zeros_like(mhd*1.1) # multiplication to force it to be a float array
    for well_set in well_sets:
        full_probs[:,well_set] = get_gene_hit_probs(mhd[:,well_set], pop_hits[well_set])
    return full_probs

def draw_data(probs):
    draw = np.random.rand(probs.shape[0], probs.shape[1])
    return np.where(draw<probs, 1, 0)

def get_log_like(hits, probs):
    polarized_hits = (hits-0.5)*2 # changes hits from 0 and 1 to -1 and 1
    # if hit = 1, yeilds prob(hit)
    # if hit = 0, yeilds 1-prob(hit)
    log_like = np.log( (1 - hits) + (probs*polarized_hits) )
    return np.sum(log_like, axis=1)

def simulate_hits(n, probs_to_sim, pop_hits):
    ll_rec = []
    for i in range(n):
        sim = draw_data(probs_to_sim)
        pmtmp = dict()
        pmtmp['model1'] = get_gene_hit_probs(sim, pop_hits)
        pmtmp['model2'] = split_probs(pop_hits, sim, [[w for w in range(len(wells)) if well_to_strain[wells[w]]==s] for s in strains])
        pmtmp['model3'] = split_probs(pop_hits, sim, [[w for w in range(len(wells)) if wells[w][:2]==p] for p in plates])
        pmtmp['model4'] = split_probs(pop_hits, sim, [[w for w in range(len(wells)) if wells[w][:2]==p and well_to_strain[wells[w]]==s] for p in plates for s in strains])
        ll_rec.append([get_log_like(sim, pmtmp[model]) for model in pmtmp])
    return ll_rec

def decide_on_main_dependence(row):
    pcols = [i for i in row.keys() if 'ratio_p_corr' in i]
    if np.min(row[pcols]) < 0.05:
        if row['model4_LL_ratio_p_corrected'] < 0.05:
            if row['AIC_model4_v_2']<0 and row['AIC_model4_v_3']<0:
                return 'both'
        if row['model2_LL_ratio_p_corrected'] < 0.05 and row['model2_LL_ratio'] > row['model3_LL_ratio']:
            return 'strain'
        elif row['model3_LL_ratio_p_corrected'] < 0.05 and row['model2_LL_ratio'] < row['model3_LL_ratio']:
            return 'environment'
    return 'NA'

## BEGIN MAIN ##

plates = ['P1', 'P2', 'P3']
strains = ['a', 'alpha', 'diploid']
wi = pd.read_csv('../accessory_files/VLTE_by_well_info.csv')[['plate.well', 'contam', 'strain']]
wi['plate.well'] = wi['plate.well'].apply(lambda p: p[:2]+p[3:])
well_to_strain = {i[0]: i[1] for i in wi.as_matrix(['plate.well', 'strain'])}

seq_tool = SeqInfoGetter('../../Output/WGS/reference/w303_vlte.gff', '../../Output/WGS/reference/w303_vlte.fasta')

gene_info = pd.read_csv('../accessory_files/yeast_gene_annotations.tsv', delimiter='\t')
gene_info = gene_info[gene_info['featureType']=='ORF'].loc[gene_info['briefDescription'].apply(lambda bd: ('Putative protein' not in bd) and ('Dubious open reading frame' not in bd))]
gene_info['size'] = gene_info['end']-gene_info['start']

o2g = {i[0]:i[1:] for i in gene_info.as_matrix(['ORF', 'Gene_ORF', 'briefDescription', 'SGDID', 'chromosome'])}
orfs_w_sizes = {i[0]:i[1] for i in gene_info.as_matrix(['ORF', 'size'])}
gene_info_orfs_only = gene_info[gene_info['featureType']=='ORF'].loc[gene_info['briefDescription'].apply(lambda bd: ('Putative protein' not in bd) and ('Dubious open reading frame' not in bd))]
gene_info_orfs_only = gene_info_orfs_only[gene_info['ORF'].isin(list(seq_tool.genes['ORF']))]
gene_info_orfs_only['nonsyn_opportunity'] = gene_info_orfs_only.apply(lambda row: seq_tool.get_mutational_opps2(row['ORF'], "ORF")*row['size'], axis=1) # Yeilds % of nonsyn/all random mutations
mean_orf_size = np.mean(gene_info_orfs_only['size'])
mean_orf_nonsyn_opps = np.mean(gene_info_orfs_only['nonsyn_opportunity'])

mutation_types_in_consequence_order = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
mutation_impacts_in_consequence_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
mut_type_simple = {
    'missense': ['missense_variant'],
    'nonsense': ['stop_lost', 'stop_gained', 'start_lost'],
    'synonymous': ['synonymous_variant', 'stop_retained_variant', 'initiator_codon_variant'],
    'noncoding': ['upstream_gene_variant', 'splice_region_variant', 'intron_variant', 'splice_acceptor_variant', ''],
    'indel': ['conservative_inframe_insertion', 'conservative_inframe_deletion', 'disruptive_inframe_insertion', 'disruptive_inframe_deletion', 
           'frameshift_variant'],
    'SV': ['SV']
 }

mut_types = [i for i in mut_type_simple]
mu_simplify = dict()
for m in mut_type_simple:
    mu_simplify.update({i:m for i in mut_type_simple[m]})

all_wells = [f.split('/')[-1].split('_')[0] for f in glob(outdir + '/well_output/*filtered.tsv')]    
excluded_wells = ['P1B03'] # excluding a haploid that autodiploidized from haploid mutation total calculations
wells = [w for w in all_wells if w not in excluded_wells]
seq_gens = [70, 1410, 2640, 5150, 7530, 10150]


# dictionary to keep track of nonsynonymous mutations in different ORFs (across populations)
orf_hit_info = defaultdict(Counter)
orf_chromos = dict()
for well in all_wells:
    if well_to_strain[well] == 'diploid':
        thresh1, thresh2 = 0.1, 0.4
    else:
        thresh1, thresh2 = 0.1, 0.9
    muts = pd.read_csv(outdir+'/well_output/' + well + '_filtered.tsv', delimiter='\t').sort_values(['CHROM', 'POS'])
    muts['ANN_simpler'] = muts.apply(lambda r: annotation_parser(r), axis=1)
    muts['mutation_type'] = muts['ANN_simpler'].apply(lambda a: get_mutation_type(a))
    muts['mutation_impact'] = muts['ANN_simpler'].apply(lambda a: get_mutation_impact(a))
    muts['ORF_hit'] = muts['ANN_simpler'].apply(lambda a: genes_affected(a))
    muts['ORF_hit_synonymous'] = muts['ANN_simpler'].apply(lambda a: genes_affected_syn(a))
    fixed_rec, present_rec = [], []
    for g in seq_gens:
        muts['fixed_by_'+str(g)] = muts.apply(lambda row: fixed_by(row, g, thresh2), axis=1)
        if well_to_strain[well] == 'diploid':
            muts['LOH_fixed_by_'+str(g)] = muts.apply(lambda row: fixed_by(row, g, 0.9, prefix='LOH_'), axis=1)
        else:
            muts['LOH_fixed_by_'+str(g)] = ['NA']*len(muts)
        muts['present_at_'+str(g)] = muts.apply(lambda row: fixed_by(row, g, thresh1), axis=1)

    # Adding a column called "mutation_group" that tags together mutations that seem like they are part of the same event
    # Meaning they are within 25 bp of another mutation in the group and have allele counts consistent w being at the same frequency
    # at each timepoint
    muts['ID'] = muts.apply(lambda row: row['CHROM'] + '_' + str(row['POS']) + '_' + str(row['ALT']), axis=1)
    mutid_to_mutation_group = dict()
    mg_num = 1
    big_mat = []
    mutation_group_to_rows = dict()
    for row_ind, row in muts.iterrows():
        no_merge = True
        if pd.isnull(row['SVTYPE']): # Not going to merge any SVs
            for mg in mutation_group_to_rows:
                near_mg = False
                for mg_row in mutation_group_to_rows[mg]:
                    if (mg_row['CHROM'] == row['CHROM']) and (np.abs(mg_row['POS'] - row['POS'])<=25):
                        near_mg = True
                if near_mg:
                    if test_for_merge(row, mutation_group_to_rows[mg]) < 1:
                        mutation_group_to_rows[mg].append(row)
                        no_merge = False
                        mutid_to_mutation_group[row['ID']] = mg
                        break
        if no_merge:
            mutation_group_to_rows[mg_num] = [row]
            mutid_to_mutation_group[row['ID']] = mg_num
            mg_num += 1
    muts['mutation_group'] = muts['ID'].apply(lambda mid: mutid_to_mutation_group[mid])
    # this give the mutation type for the group if mutations are grouped by proximity (putatively one event)
    # we will make all of them in the group have that mutation type in this column
    mutation_group_counts = dict(muts['mutation_group'].value_counts())
    mg_to_mut_types = dict()
    mg_to_mut_impacts = dict()
    for mg in [i for i in mutation_group_counts if mutation_group_counts[i]>1]:
        mg_to_mut_types[mg] = list(muts[muts['mutation_group']==mg]['mutation_type'])
        mg_to_mut_impacts[mg] = list(muts[muts['mutation_group']==mg]['mutation_impact'])
    muts['group_mutation_type'] = muts.apply(lambda r: get_group_mutation_type(r, mutation_group_counts, mg_to_mut_types), axis=1)
    muts['group_mutation_impact'] = muts.apply(lambda r: get_group_mutation_impact(r, mutation_group_counts, mg_to_mut_impacts), axis=1)
    muts['in_telomere'] = muts.apply(lambda r: in_telomere(r), axis=1)
    muts.to_csv(outdir+'/processed_well_output/'+well+'_processed.tsv', index=False, sep='\t')
    
    ### ALL FURTHER ANALYSIS EXCLUDES 2-micron and telomeric mutations, and excludes the well that autodiploidized
    if well not in excluded_wells:
        ## Adding to ORF hit dictionary
        mgs_seen = set()
        # Exclude from analysis mutations in the 2-micron plasmid and telomeres, and SVs
        td = muts[pd.isnull(muts['SVTYPE'])& (~muts['in_telomere']) & (~(muts['CHROM']=='2-micron'))]
        for jnk, row in td[pd.notnull(td['ORF_hit'])].iterrows():
            if row['mutation_group'] not in mgs_seen:
                mgs_seen.add(row['mutation_group'])
                for orf in str(row['ORF_hit']).split(';'):
                    if orf != '':
                        orf_chromos[orf] = row['CHROM']
                        tmp_dict = orf_hit_info[orf]
                        tmp_dict[well+'_present'] += 1
                        if row['fixed_by_10150']:
                            tmp_dict['num_hits'] += 1
                            tmp_dict[well] += 1
                            if well_to_strain[well] == 'diploid':
                                tmp_dict['dip_hits'] += 1
                                if row['group_mutation_impact'] == 'HIGH':
                                    tmp_dict['high_impact'] += 1
                                    if row['LOH_fixed_by_10150']:
                                        tmp_dict['LOH'] += 1
                                        tmp_dict['LOH_high_impact'] += 1
                                    else:
                                        tmp_dict['no_LOH'] += 1
                                        tmp_dict['no_LOH_high_impact'] += 1
                                else:
                                    if row['LOH_fixed_by_10150']:
                                        tmp_dict['LOH'] += 1
                                    else:
                                        tmp_dict['no_LOH'] += 1
                            else:
                                tmp_dict['hap_hits'] += 1
                                if row['group_mutation_impact'] == 'HIGH':
                                    tmp_dict['haploid_high_impact'] += 1
                                    tmp_dict['high_impact'] += 1


    
#### PART 2 looking at multi hit genes (ORFs) ###

orf_cols = ['num_hits', 'dip_hits', 'hap_hits', 'high_impact', 'LOH', 'no_LOH', 'LOH_high_impact', 'no_LOH_high_impact', 'haploid_high_impact'] + wells + [w+'_present' for w in wells]
big_mat = [[orf, o2g.get(orf, [orf])[0], orf_chromos[orf]] + [orf_hit_info[orf][c] for c in orf_cols] for orf in orf_hit_info]
hit_df = pd.DataFrame(big_mat, columns=['ORF', 'Gene_ORF', 'CHROM']+orf_cols)
hit_df['pops_hit'] = hit_df.apply(lambda row: len([w for w in wells if row[w]>0]), axis=1)
hit_df['size'] = hit_df['ORF'].apply(lambda o: orfs_w_sizes.get(o, mean_orf_size))
def get_nonsyn_opp(row):
    if row['ORF'] in list(seq_tool.genes['ORF']):
        return seq_tool.get_mutational_opps2(row['ORF'], "ORF")*row['size']
    else:
        # if we don't have info in our gff, return nan
        return np.nan

hit_df['nonsyn_opportunity'] = hit_df.apply(lambda row: get_nonsyn_opp(row), axis=1) # Yeilds % of nonsyn/all random mutations
hit_df['multiplicity'] = mean_orf_nonsyn_opps * (hit_df['num_hits']/hit_df['nonsyn_opportunity'])

pop_to_genes_mat = [[w, ';'.join(list(hit_df[hit_df[w]>0]['Gene_ORF']))] for w in wells]
pd.DataFrame(pop_to_genes_mat, columns=['plate_well', 'genehits_fixed_10K']).to_csv('../../Output/Browser/pops_to_genes_fixed.tsv', index=False, sep='\t')

## LOOKING FOR WHETHER ORF HITS DEPEND ON ENVIRONMENT OR STRAIN BACKGROUND ##

# count total hits for each population
pop_hits = np.sum(hit_df[wells], axis=0)
# focus on multi-hit orfs
multi_hit_data = hit_df[hit_df['pops_hit']>=6]
orf_names = list(multi_hit_data['ORF'])
mh_dat = np.clip(np.array(multi_hit_data[wells]), 0, 1) # just presence/absence of a hit, not worrying about 2+ hits in any one population
prob_models = dict()
prob_models['model1'] = get_gene_hit_probs(mh_dat, pop_hits) # no dependence
prob_models['model2'] = split_probs(pop_hits, mh_dat, [[w for w in range(len(wells)) if well_to_strain[wells[w]]==s] for s in strains]) # prob of hit depends on strain
prob_models['model3'] = split_probs(pop_hits, mh_dat, [[w for w in range(len(wells)) if wells[w][:2]==p] for p in plates]) # prob of hit depends on env
prob_models['model4'] = split_probs(pop_hits, mh_dat, [[w for w in range(len(wells)) if wells[w][:2]==p and well_to_strain[wells[w]]==s] for p in plates for s in strains]) # prob of hit depends on both

# calculate log-likelihood of each model per orf
actual_ll = [get_log_like(mh_dat, prob_models[model]) for model in prob_models]
# simulate orf hits under model one, calculate log-likelihood of each model per orf for every sim
sim_log_likes = simulate_hits(10000, prob_models['model1'], pop_hits)
# compute log-likelihood ratios
gene_to_sim_ll_ratios = dict()
gene_to_data_ll_ratios = dict()
for g in range(len(orf_names)):
    tg = [[],[],[]] # log-likelihood ratios of model 2 / 1 , 3 / 1, and 4 / 1
    for rec in sim_log_likes:
        for i in range(3):
            tg[i].append(rec[i+1][g]-rec[0][g])
    gene_to_sim_ll_ratios[g] = tg
    gene_to_data_ll_ratios[g] = [actual_ll[i+1][g]-actual_ll[0][g] for i in range(3)]

# record log-likes per orf and calculate p-values based on simulated data
orfs_to_results = dict()
for g in range(len(orf_names)):
    tmp = gene_to_data_ll_ratios[g]
    for i in range(3):
        percentile = len([ll_ratio for ll_ratio in gene_to_sim_ll_ratios[g][i] if ll_ratio > gene_to_data_ll_ratios[g][i]])/len(gene_to_sim_ll_ratios[g][i])
        tmp.append(percentile)
    orfs_to_results[orf_names[g]] = tmp
    
for i in range(3):
    multi_hit_data['model'+str(i+2)+'_LL_ratio'] = multi_hit_data['ORF'].apply(lambda g: orfs_to_results.get(g, [np.nan]*6)[i])
    multi_hit_data['model'+str(i+2)+'_LL_ratio_p'] = multi_hit_data['ORF'].apply(lambda g: orfs_to_results.get(g, [np.nan]*6)[i+3])

all_pvals = list(multi_hit_data['model2_LL_ratio_p'])+list(multi_hit_data['model3_LL_ratio_p'])+list(multi_hit_data['model4_LL_ratio_p'])
corrected_sig_test = benjamini_hochberg(all_pvals, alpha=0.05)
for i in range(3):
    multi_hit_data['model'+str(i+2)+'_LL_ratio_p_corrected'] = corrected_sig_test[1][i*len(multi_hit_data):(i+1)*len(multi_hit_data)]

# AIC = 2k - 2LL where k is num parameters and LL is log likelihood. Model 4 has 9 parameters, vs. 3 in model 2 or 3
# AIC_4_vs_2 = 18 - 2LL4 - (6 - 2LL2) = 12 - 2(LL4-LL2)
# LL4-LL2 is the same as model4_LL_ratio-model2_LL_ratio (since both are just the LL minus LL1)
# AIC_4_vs_3 calculated the same way. If this AIC comparison is less than 0, model 4 has a lower AIC and is favored
for i in range(2):
    multi_hit_data['AIC_model4_v_'+str(i+2)] = multi_hit_data.apply(lambda r: 12 - 2*(r['model4_LL_ratio']-r['model'+str(i+2)+'_LL_ratio']), axis=1)
    
multi_hit_data['dependent_on'] = multi_hit_data.apply(lambda r: decide_on_main_dependence(r), axis=1)

hit_df.to_csv(outdir+'/gene_hit_data.tsv', index=False, sep='\t')
multi_hit_data.to_csv(outdir+'/multi_hit_genes.tsv', index=False, sep='\t')

## GO ANALYSIS

# Convert ORF names to SGDIDs for GO analysis
multi_hit_sgdids = list(gene_info[gene_info['ORF'].isin(orf_names)]['SGDID'])

obodag = GODag("../accessory_files/go-basic.obo")  # http://geneontology.org/ontology/go-basic.obo
goid_to_gene_list = defaultdict(list)
genename_2_id = dict()
with open('../accessory_files/gene_association.sgd', 'r') as infile:
    for line in infile:
        if line[0] != '!':
            s = line.split('\t')
            goid_to_gene_list[s[4]].append(s[1])
            genename_2_id[s[2]] = s[1]
id_2_genename = {genename_2_id[i]: i for i in genename_2_id}
# Only looking at "biological process" GO terms
geneid2gos_yeast = read_gaf('../accessory_files/gene_association.sgd', namespace='BP')
ids = [i for i in geneid2gos_yeast.keys()]
background_set = [genename_2_id[i]for i in genename_2_id]
goeaobj = GOEnrichmentStudy(
        background_set, # List of all genes in analysis
        geneid2gos_yeast, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method

goea_results_all = goeaobj.run_study(multi_hit_sgdids, keep_if=lambda x: x.p_uncorrected < 0.05)
go_results = sorted(goea_results_all, key=lambda r: r.p_fdr_bh)

cols = ['GO ID', 'GO term', 'pval_uncorrected', 'pval_benjamini_hochberg', 'num_hits', 'num_in_group', 'hits']
big_mat = []
for res in go_results:
    big_mat.append([res.GO, res.name, res.p_uncorrected, res.p_fdr_bh, res.ratio_in_study[0], res.ratio_in_pop[0], 
                    ';'.join([id_2_genename[i] for i in res.study_items])])
        
pd.DataFrame(big_mat, columns=cols).sort_values(by='pval_benjamini_hochberg').to_csv(outdir+'/GO_enrichments.tsv', sep='\t', index=False)
