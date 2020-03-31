import pandas as pd

d1 = pd.read_csv('SampleSheet_Batch1.csv', skiprows=24)
d2 = pd.read_csv('SampleSheet_Batch2.csv')
d1['suffixes'] = ['_S' + str(i) for i in range(1, len(d1)+1)]
d2['suffixes'] = ['_S' + str(i) + '_L001' for i in range(1, len(d2)+1)]
d1['samp'] = d1['Sample_ID'].apply(lambda s: '_'.join(s.split('_')[4:]))
d2['samp'] = d2['Sample_ID'].apply(lambda s: '_'.join(s.split('_')[4:]))
d1['sample_vcf_combined_option'] = d1['samp'].apply(lambda s: '../../Output/WGS/combined_work/' + s + '.g.vcf')
d1['sample_vcf'] = d1['samp'].apply(lambda s: '../../Output/WGS/work/' + s + '.g.vcf')
d2['sample_vcf_combined_option'] = d2['samp'].apply(lambda s: '../../Output/WGS/work/' + s + '.g.vcf')
d2['sample_vcf'] = d2['samp'].apply(lambda s: '../../Output/WGS/work/' + s + '.g.vcf')
d = pd.concat([d1[['Sample_ID', 'suffixes', 'samp', 'sample_vcf', 'sample_vcf_combined_option']], d2[['Sample_ID', 'suffixes', 'samp', 'sample_vcf', 'sample_vcf_combined_option']]])
d = d[d['Sample_ID'].apply(lambda s: 'P4' not in s and 'G1410_P1C07' not in s)]
d['prefixes'] = d['Sample_ID'].apply(lambda s: '_'.join(s.split('_')[:4]))
d['gen'] = d['samp'].apply(lambda s: int(s.split('_')[0][1:]))
d = d.sort_values(by='gen')
d[['samp']].to_csv('Samples.txt', header=None, index=False)
d[['prefixes']].to_csv('Sample_Prefixes.txt', header=None, index=False)
d[['suffixes']].to_csv('Sample_Suffixes.txt', header=None, index=False)
d[['samp', 'sample_vcf']].to_csv('vcf_map.txt', header=None, index=False, sep='\t')
d[['samp', 'sample_vcf_combined_option']].to_csv('vcf_map_combined_option.txt', header=None, index=False, sep='\t')
