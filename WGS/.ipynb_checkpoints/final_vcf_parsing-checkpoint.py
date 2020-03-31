from glob import glob
import pandas as pd
import numpy as np

def annotation_parser(a):
    # input is a single annotation
    sa = a.split('|')
    allele, effs, gene_hit = sa[0], sa[1], sa[3]
    extra_ann = list(orf2ann.setdefault(gene_hit, ['']*2))
    if effs not in ['intergenic_region', 'downstream_gene_variant']:
        if 'upstream' in effs:
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


orf2ann = {i[0]: [str(j) for j in i[1:]] for i in pd.read_csv('yeast_gene_annotations.csv').as_matrix(['ORF', 'Gene', 'briefDescription'])}

flist = glob('../../Output/WGS/chromo_vcfs/*ann.vcf')

d = pd.concat([pd.read_csv(fin, skiprows=62, delimiter='\t') for fin in flist])
d['ANN'] = d['INFO'].apply(lambda i: format_annotations(i))
d = d.rename(columns={'#CHROM': 'CHROM'})
orig_cols = cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
new_cols = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'ANN', 'FORMAT']
d[new_cols + [i for i in d if i not in orig_cols and i not in new_cols]].to_csv('../../Output/WGS/FINAL_FULL_VARIANTS.tsv', index=False, sep='\t')
print("done!")
