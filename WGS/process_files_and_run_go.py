import pandas as pd
import numpy as np
from scipy import stats as sci_stats
from glob import glob
from collections import defaultdict
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_gaf

def genes_affected(a):
    keep_anns = []
    for ann in str(a).split(';'):
        if '|' in ann:
            anns = ann.split('|')
            for mut_type in anns[1].split('&'):
                if mu_simplify[mut_type] not in ['synonymous', 'noncoding']:
                    if 'Dubious' not in anns[5]:
                        keep_anns.append(anns[2])
    return ';'.join(keep_anns)

def genes_affected_syn(a):
    keep_anns = []
    for ann in str(a).split(';'):
        if '|' in ann:
            anns = ann.split('|')
            for mut_type in anns[1].split('&'):
                if mu_simplify[mut_type]=='synonymous':
                    if 'Dubious' not in anns[5]:
                        keep_anns.append(anns[2])
    return ';'.join(keep_anns)

def fixed_by(r, gen, freq_thresh):
    # find gen index
    g_start = seq_gens.index(gen)
    # if fixed at the previous generation, automatically fixed now
    if g_start > 0:
        if r['fixed_by_'+str(seq_gens[g_start-1])]:
            return True

    states = []
    for col in ['G'+str(g)+'_allele_counts' for g in seq_gens[g_start:]]:
        if col in r:
            ref, alt = int(r[col].split(',')[0]), int(r[col].split(',')[1])
            tot = ref + alt
            if tot >= 5:
                if alt/tot > freq_thresh:
                    states.append('fixed')
                else:
                    states.append('not fixed')
            else:
                states.append('unknown')
    if 'fixed' in states and 'not fixed' not in states:
        return True
    else:
        return False
    
def present_at(r, gen, freq_thresh):
    ac = r['G'+str(gen)+'_allele_counts']
    # enforce that if it is fixed it's present (which might not be true due to low counts otherwise)
    if r['fixed_by_'+str(gen)]:
        return True
    ref, alt = int(ac.split(',')[0]), int(ac.split(',')[1])
    tot = ref + alt
    if tot >= 5:
        if alt/tot > freq_thresh:
            return True
        else:
            return False
    else:
        return False


def test_for_merge(test_row, mg_rows, c2i):
    # Tests that a mutation row is the same or different from other rows that are nearby (in the genome)
    # at each time point, using Fisher's exact test
    pvals = []
    for gen in seq_gens:
        if 'G'+str(gen)+'_allele_counts' in c2i:
            ref, alt = [int(i) for i in test_row[c2i['G'+str(gen)+'_allele_counts']].split(',')]
            mg_ref, mg_alt = 0, 0
            for r in mg_rows:
                mg_ref += int(r[c2i['G'+str(gen)+'_allele_counts']].split(',')[0])
                mg_alt += int(r[c2i['G'+str(gen)+'_allele_counts']].split(',')[1])
            pvals.append(sci_stats.fisher_exact([[alt, ref], [mg_alt, mg_ref]])[1])
    return len([i for i in pvals if i<0.01])


wi = pd.read_csv('../accessory_files/VLTE_by_well_info.csv')[['plate.well', 'contam', 'strain']]
wi['plate.well'] = wi['plate.well'].apply(lambda p: p[:2]+p[3:])
well_to_strain = {i[0]: i[1] for i in wi.as_matrix(['plate.well', 'strain'])}

mut_type_simple = {
    'missense': ['missense_variant'],
    'nonsense': ['stop_lost', 'stop_gained', 'start_lost'],
    'synonymous': ['synonymous_variant', 'stop_retained_variant', 'initiator_codon_variant'],
    'noncoding': ['upstream_gene_variant', 'splice_region_variant', 'intron_variant', 'splice_acceptor_variant', ''],
    'indel': ['conservative_inframe_insertion', 'conservative_inframe_deletion', 'disruptive_inframe_insertion', 'disruptive_inframe_deletion', 
           'frameshift_variant']
 }
mut_types = [i for i in mut_type_simple]
mu_simplify = dict()
for m in mut_type_simple:
    mu_simplify.update({i:m for i in mut_type_simple[m]})

wells = [f.split('/')[-1].split('_')[0] for f in glob('../../Output/WGS/well_output/*filtered.tsv')]    
seq_gens = [70, 1410, 2640, 5150, 7530, 10150]
wells_to_10k_fixed_genes = dict()
wells_to_10k_fixed_genes_syn = dict()
for well in wells:
    if well_to_strain[well] == 'diploid':
        thresh1, thresh2 = 0.1, 0.4
    else:
        thresh1, thresh2 = 0.1, 0.95
    muts = pd.read_csv('../../Output/WGS/well_output/' + well + '_filtered.tsv', delimiter='\t').sort_values(['CHROM', 'POS'])
    muts['ORF_hit'] = muts['ANN'].apply(lambda a: genes_affected(a))
    muts['ORF_hit_synonymous'] = muts['ANN'].apply(lambda a: genes_affected_syn(a))
    fixed_rec, present_rec = [], []
    for g in seq_gens:
        muts['fixed_by_'+str(g)] = muts.apply(lambda row: fixed_by(row, g, thresh2), axis=1)
        muts['present_at_'+str(g)] = muts.apply(lambda row: fixed_by(row, g, thresh1), axis=1)
    wells_to_10k_fixed_genes[well], wells_to_10k_fixed_genes_syn[well] = [], []
    for entry in list(muts[muts['fixed_by_10150']]['ORF_hit']):
        if entry != '':
            wells_to_10k_fixed_genes[well] += entry.split(';')

    # Adding a column called "mutation_group" that tags together mutations that seem like they are part of the same event
    # Meaning they are within 50 bp of another mutation in the group and have allele counts consistent w being at the same frequency
    # at each timepoint
    cols = muts.columns
    c2i = {cols[i]: i for i in range(len(cols))}
    mg_num = 1
    big_mat = []
    mutation_group_to_rows = dict()
    for row in muts.as_matrix(cols):
        no_merge = True
        for mg in mutation_group_to_rows:
            near_mg = False
            for mg_row in mutation_group_to_rows[mg]:
                if (mg_row[c2i['CHROM']] == row[c2i['CHROM']]) and (np.abs(mg_row[c2i['POS']] - row[c2i['POS']])<50):
                    near_mg = True
            if near_mg:
                if test_for_merge(row, mutation_group_to_rows[mg], c2i) < 1:
                    mutation_group_to_rows[mg].append(row)
                    no_merge = False
                    big_mat.append(list(row) + [mg])
                    break
        if no_merge:
            mutation_group_to_rows[mg_num] = [row]
            big_mat.append(list(row) + [mg_num])
            mg_num += 1
    pd.DataFrame(big_mat, columns=list(cols)+['mutation_group']).to_csv('../../Output/WGS/processed_well_output/' + well + '_processed.tsv', index=False, sep='\t')

genes_to_pops = defaultdict(set)
for well in wells:
    for gene in wells_to_10k_fixed_genes[well]:
        genes_to_pops[gene].add(well)

multi_hit = [g for g in genes_to_pops if len(genes_to_pops[g])>=6]
gene_to_pop_mat = [[g, ';'.join(list(genes_to_pops[g]))] for g in genes_to_pops]
pd.DataFrame(gene_to_pop_mat, columns=['Gene', 'Populations_w_fixation_at_10K']).to_csv('../../Output/WGS/genes_to_pops_hit.tsv', index=False, sep='\t')

# Convert ORF names to SGDIDs for GO analysis
gene_info = pd.read_csv('../accessory_files/yeast_gene_annotations.tsv', delimiter='\t')
multi_hit_sgdids = list(gene_info[gene_info['ORF'].isin(multi_hit)]['SGDID'])

## GO Analysis
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
        
pd.DataFrame(big_mat, columns=cols).sort_values(by='pval_benjamini_hochberg').to_csv('../../Output/WGS/GO_enrichments.tsv', sep='\t')