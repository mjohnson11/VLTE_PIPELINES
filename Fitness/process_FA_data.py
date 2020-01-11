# A script to process FACS data from 384-well plates and measure fitness using that

# importing
from matplotlib import pyplot as pl
import numpy as np
from scipy import stats as sci_stats
import pandas as pd
import csv
# this was installed like: pip install FlowCytometryTools
import FlowCytometryTools as fct
import seaborn as sns
from glob import glob
sns.set_style('white')

# These are controls we originally planned to use to measure ancestral fitness. We had a ton of issues with these, in part probably due to to mutations 
# present in our ancestral stocks that we wanted to be clonal. In any case, we only use the 2490A strain to standardize fitness, the rest we ignore, but
# we still flag them in this analysis
control_d = {
    'F06': 'v2N-1',
    'A10': 'v2N-2',
    'A12': 'va-1',
    'G12': 'valpha-1',
    'D04': '2490A'
}

# Gating on good cells - this excludes the smallest events, which sometimes are just some sort of junk in the media or machine
cell_gate = fct.PolyGate([(10000, 5500), (30000, 7000), (300000, 9000), (300000, 7000), (100000, 6000), (40000, 5500)], ['FSC-A', 'SSC-A'])
# Gates for reference for each plate / machine - the reference FACS results are different depending on the type of media and which machine we used, so the gates are different
# We are gating on FITC and PE here, even though the ref doesn't glow red - it just ended up working better than FITC and SSC in terms of separating the two clouds
p23_fortessa = fct.PolyGate([(6500, -1), (6500, 3500), (6000, 5000), (6500, 6300), (6700, 6800), (10000, 12000), (10000, -1)], ['FITC-A', 'PE-A'])
p1_fortessa = fct.PolyGate([(6200, -1), (5700, 5000), (6000, 5700), (8000, 7700), (10000, 11000), (10000, -1)], ['FITC-A', 'PE-A'])
p23_lsr2 = fct.PolyGate([(6200, -1), (6000, 6500), (6500, 7350), (8000, 9000), (10000, 11000), (10000, -1)], ['FITC-A', 'PE-A'])
p1_lsr2 = fct.PolyGate([(6200, -1), (5900, 6300), (6500, 7200), (8000, 8500), (10000, 11000), (10000, -1)], ['FITC-A', 'PE-A'])
pe_fitc_ref_gates = {'P2_R2': p23_fortessa, 'P3_R2': p23_fortessa, 'P2_R1': p23_lsr2, 'P3_R1': p23_lsr2, 'P1_R1': p1_fortessa, 'P1_R2': p1_lsr2}

# Defining the order of plates within each 384 well plate
plate_sets = {'VGFA1': [['Gen550', 'Gen1410'], ['Gen5150', 'Gen10150']],
              'VGFA2':[['Gen70', 'Gen2640'], ['Gen3630', 'Gen7530']]}
# Fitness assays were done in two batches - two sets of 4 generations
gens = {'VGFA1': [550, 1410, 5150, 10150], 'VGFA2': [70, 2640, 3630, 7530]}
all_gens = [70, 550, 1410, 2640, 3630, 5150, 7530, 10150]

# blanks to use to calculate the average percentage of reference in the ref gate
blanks_for_ref_perc = {
    'VGFA1_P1': [(5150, 'A11'), (5150, 'G03'), (5150, 'H03'), (10150, 'A11'), (10150, 'G03'), (10150, 'H03')],
    'VGFA1_P2': [(5150, 'A08'), (5150, 'D07'), (5150, 'F10'), (10150, 'A08'), (10150, 'D07'), (10150, 'F10')],
    'VGFA1_P3': [(5150, 'B02'), (5150, 'B03'), (5150, 'B04'), (10150, 'B02'), (10150, 'B03'), (10150, 'B04')],
    'VGFA2_P1': [(3630, 'A11'), (3630, 'G03'), (3630, 'H03'), (7530, 'A11'), (7530, 'G03'), (7530, 'H03')],
    'VGFA2_P2': [(3630, 'A08'), (3630, 'D07'), (3630, 'F10'), (7530, 'A08'), (7530, 'D07'), (7530, 'F10')],
    'VGFA2_P3': [(3630, 'B02'), (3630, 'B03'), (3630, 'B04'), (7530, 'B02'), (7530, 'B03'), (7530, 'B04')]
}

# Tells which replicate was done on which machine
rep_dict = {'P1': {'R2': 'LSR2', 'R1': 'Fortessa'}, 'P2': {'R1': 'LSR2', 'R2': 'Fortessa'}, 'P3': {'R1': 'LSR2', 'R2': 'Fortessa'}}

def well_to_plate_well(well, plates):
    # From 384 well plate to plate name and 96 well plate well name
    row = ord(well[0])-65
    col = int(well[1:])
    plate = plates[row % 2][(col-1) %2]
    plate_row = chr((row-2)//2+66)
    plate_col = (col-1)//2+1
    return plate + '_' + plate_row + str(plate_col).zfill(2)

def get_plate_data(dir_base, plate_set, plate_type='384', lsr2=False):
    # given a directory base like 'FA_timepoint_0/Specimen_001_', this will try to read in all the files corresponding to each well in a 96 or 384 well plate
    # and will return a dictionary like td['well_id'] = dataframe with facs info
    print('reading from', dir_base)
    td = dict()
    wells_missed = []
    if plate_type == '384':
        let_top, col_top = 16, 25
    elif plate_type == '96':
        let_top, col_top = 8, 13
    else:
        print('unrecognized plate type')
        return None
    for let in [chr(i+65) for i in range(let_top)]:
        for col in range(1, col_top):
            orig_well = let + str(col).zfill(2)
            if plate_type == '384':
                well = well_to_plate_well(orig_well, plate_set)
            else:
                well = orig_well
            flist = glob(dir_base + '*' + orig_well + '*.fcs')
            try:
                assert len(flist) == 1
                # Reading in file and immediately gating on good cells
                samp = fct.FCMeasurement(ID=well, datafile=flist[0]).transform('tlog', channels=['FITC-A', 'SSC-A', 'PE-A']).gate(cell_gate)
                td[well] = samp
            except AssertionError:
                wells_missed.append(well)
                td[well] = None
    if len(wells_missed) > 0:
        print('Missed files for', len(wells_missed), 'wells:', ' '.join(wells_missed))
    return td

def get_ref_counts(df):
    ref_counts = df.gate(use_gate).shape[0]
    total_counts = df.shape[0]
    if total_counts < 1000: # Excluding timepoints with less than 1000 reads
        freq = np.nan
        #print('Low counts for a sample...')
    else:
        freq = ref_counts/total_counts
    return ref_counts, total_counts-ref_counts, freq

def get_fit(row, tps, gen_string):
    ref_freqs = np.array([row[gen_string + '_Ref_Freq_T' + str(t)] for t in tps if pd.notnull(row[gen_string + '_Ref_Freq_T' + str(t)])])
    times = np.array([t for t in tps if pd.notnull(row[gen_string + '_Ref_Freq_T' + str(t)])])*10
    # excluding time intervals where ref or test is >95% in both timepoints
    use_tp_until = len(times)
    for t in range(1, len(times)):
        if (ref_freqs[t] > 0.95 and ref_freqs[t-1] > 0.95) or (ref_freqs[t] < 0.05 and ref_freqs[t-1] > 0.05):
            use_tp_until = t-1
            break
    if use_tp_until > 0:
        test_freqs = 1-ref_freqs
        return sci_stats.linregress(times[:use_tp_until+1], np.log(test_freqs[:use_tp_until+1]))[0] - sci_stats.linregress(times[:use_tp_until+1], np.log(ref_freqs[:use_tp_until+1]))[0]

def pick_s_use(row, gen_string):
    if pd.notnull(row[gen_string + '_s']):
        return row[gen_string + '_s']
    else:
        return row[gen_string + '_s_with_T0']
    
def process_files(assay_base, assay_name):
    ## READING FILES
    tps = [0,1,2,3]
    if 'VGFA2' in assay_base:
        dirs = [assay_base+'_D'+str(tp+2)+'_B2_T'+str(tp) for tp in tps]
    else:
        dirs = [assay_base+'_D'+str(tp+2)+'_T'+str(tp) for tp in tps]
    dir_d = {d.split('/')[-1]: d for d in dirs}
    dat_d = dict()
    for d in dir_d:
        if 'LSR2' in dir_d[d].split('/')[1]:
            dat_d[d] = get_plate_data(dir_d[d] + '/Specimen_001_', plate_sets[dir_d[d].split('/')[0]], lsr2=True)
        else:
            dat_d[d] = get_plate_data(dir_d[d] + '/Specimen_001_', plate_sets[dir_d[d].split('/')[0]])

    ## It is very clear from the data that we mixed up the plates for 1410 and 5150 when FACSing the last time point. I will fix that here:      
    if assay_name == 'VGFA1_P3_R1':
        td = dat_d['VGFA1_P3_R1_D5_T3']
        for well in [i.split('_')[1] for i in td if 'Gen5150' in i]:
            tmp = td['Gen1410_' + well]
            td['Gen1410_' + well] = td['Gen5150_' + well]
            td['Gen5150_' + well] = tmp

    ## PROCESSING / GETTING REF FREQUENCY DATA
    mat = []
    ref_freq_in_blanks = {d.split('/')[-1]: np.mean([get_ref_counts(dat_d[d.split('/')[-1]]['Gen'+str(g)+'_'+well])[-1] for g, well in blanks_for_ref_perc[assay_name[:8]]]) for d in dirs}
    print('Blank reference freqs:', ref_freq_in_blanks)
    for row in range(8):
        for col in range(12):
            well = chr(row+65) + str(col+1).zfill(2)
            tmp = [well]
            tmp_colnames = ['Well']
            for gen in gens[assay_name.split('_')[0]]:
                g = 'Gen' + str(gen) 
                for d in dirs:
                    dat_d_name = d.split('/')[-1]
                    result = get_ref_counts(dat_d[dat_d_name][g+'_'+well]) 
                    tmp += list(result) + [(1/ref_freq_in_blanks[dat_d_name])*result[2]]
                    tp = dat_d_name[-1]
                    tmp_colnames += [g+'_Ref_Counts_T'+tp, g+'_NonRef_Counts_T'+tp, g+'_Uncorrected_Ref_Freq_T'+tp, g+'_Ref_Freq_T'+tp]
            mat.append(tmp)
                
    td = pd.DataFrame(mat, columns=tmp_colnames)
    for gen in gens[assay_name.split('_')[0]]:
        g = 'Gen' + str(gen) 
        td[g+'_s'] = td.apply(lambda r: get_fit(r, tps, g), axis=1)
        td[g+'_s_with_T0'] = td.apply(lambda r: get_fit(r, tps, g, exclude_zero=False), axis=1)

    return dat_d, td

def get_strain(well, ex_d):
    if well in control_d:
        return control_d[well]
    else:
        if ex_d[well][0] == 'No':
            return ex_d[well][1]
        else:
            return 'BAD'
        
## PLOTTING
    
def plot_well_fitc_pe(sub, df, plot_gate):
    df.plot(['FITC-A', 'PE-A'], ax=sub, gates=[plot_gate])
    sub.set_xlim([0, 10000])
    sub.set_ylim([0, 10000])
    sub.set_xlabel('')
    sub.set_ylabel('')
    sub.set_xticks([])
    sub.set_yticks([])

def plot_well_report(well, outname=None):
    gen2exp = {i: 'VGFA1' for i in [550, 1410, 5150, 10150]}
    gen2exp.update({i: 'VGFA2' for i in [70, 2640, 3630, 7530]})
    use_row = fd.loc[fd['Well']==well].iloc[0]
    fig, subps = pl.subplots(10, 8, figsize=(24, 25))
    fig.suptitle(plate + ' ' + well, fontsize=20, y=0.92)
    pl.subplots_adjust(hspace=0.2, wspace=0.2)
    all_gens = [70, 550, 1410, 2640, 3630, 5150, 7530, 10150]
    for g in range(len(all_gens)):
        gen = all_gens[g]
        exp = gen2exp[gen]
        freq_subs = [subps[i][g] for i in range(2)]
        facs_subs = [subps[i][g] for i in range(2, 10)]
        for rep_num in range(2):
            tps = [0,1,2,3]
            freqs = np.array([use_row['Gen' + str(gen) + '_Ref_Freq_T' + str(t) + '_R' + str(rep_num+1)] for t in tps])
            freq_subs[rep_num].plot(np.array(tps)*10, 1-freqs)
            freq_subs[rep_num].annotate('s = ' + str(use_row['Gen' + str(gen) + '_s_with_T0_R'+str(rep_num+1)])[:6], xy=(0.1, 0.8), xycoords='axes fraction', fontsize=15)
            freq_subs[rep_num].set_ylim([0, 1])
            freq_subs[rep_num].set_xlim([0, 30])
            if rep_num == 0:
                freq_subs[rep_num].set_title('Gen ' + str(gen), fontsize=17)
            if g != 0:
                freq_subs[rep_num].set_yticks([])
            else:
                freq_subs[rep_num].set_yticks([0, 0.5, 1])
                freq_subs[rep_num].tick_params(which='major', axis='y', labelsize=14)
                freq_subs[rep_num].set_ylabel('R' + str(rep_num+1), rotation='horizontal', fontsize=17, labelpad=35)
            freq_subs[rep_num].set_xticks([])
            sns.despine(ax=freq_subs[rep_num])
            # plotting PE vs. FITC
            assay = exp + '_' + plate + '_R' + str(rep_num+1)
            tf = facs_data[assay]
            exp_to_b = {'VGFA1': '', 'VGFA2': '_B2'}
            for t in range(4):
                plot_well_fitc_pe(facs_subs[rep_num*4+t], tf[assay + '_D' + str(t+2) + exp_to_b[exp] + '_T' + str(t)]['Gen'+str(gen) + '_' + well], pe_fitc_ref_gates[plate + '_R' + str(rep_num+1)])
                if g == 0:
                    facs_subs[rep_num*4+t].set_ylabel('R' + str(rep_num+1) + ' T' + str(t), rotation='horizontal', fontsize=17, labelpad=50)
                sns.despine(ax=facs_subs[rep_num*4+t], bottom=True, left=True)
    if outname:
        fig.savefig(outname, background='transparent', bbox_inches='tight', pad_inches=0.1)
        pl.close('all')
          
## BEGIN MAIN ##
for plate in ['P1', 'P2', 'P3']:
    print(plate)
    outname = 'browser3/' + plate + '_freq_and_s_data.csv'
    outcorr = 'browser3/' + plate + '_s_corr.png'

    facs_data = dict()
    freq_data = dict()
    for v in ['VGFA1', 'VGFA2']:
        tmp_freq_data = dict()
        for rep in ['R1', 'R2']:
            machine = rep_dict[plate][rep]
            
            a_base = v+'/'+v+'_'+machine+'/'+ '_'.join([v, plate, rep])
            a_name = a_base.split('/')[-1]
            use_gate = pe_fitc_ref_gates[plate + '_' + rep]
            facs_data[a_name], tmp_freq_data[rep] = process_files(a_base, a_name)
        freq_data[v] = tmp_freq_data['R1'].merge(tmp_freq_data['R2'], on='Well', how='inner', suffixes=('_R1', '_R2'))
        # averaging s columns
        s_cols = ['s_with_T0', 's']
        for scol in s_cols:
            # Standardizing s by subtracting off the s value for the unlabeled strain 2490A for the 8 reps in that assay
            well_2490_s = []
            for g in gens[v]:
                freq_data[v]['Gen'+str(g)+'_'+scol] = np.nanmean(freq_data[v][['Gen'+str(g)+'_'+scol+'_R1', 'Gen'+str(g)+'_'+scol+'_R2']], axis=1)
                freq_data[v]['Gen'+str(g)+'_'+scol+'_range'] = np.abs(freq_data[v]['Gen'+str(g)+'_'+scol+'_R1'] - freq_data[v]['Gen'+str(g)+'_'+scol+'_R2'])
                well_2490_s += list(freq_data[v].loc[freq_data[v]['Well']=='D04']['Gen'+str(g)+'_'+scol])
            print(v, '2490A ', scol, ' vals:', well_2490_s, 'mean:', np.nanmean(well_2490_s))
            for g in gens[v]:
                freq_data[v]['Gen'+str(g)+'_'+scol + '_scaled'] = freq_data[v]['Gen'+str(g)+'_'+scol] - np.nanmean(well_2490_s)
                freq_data[v]['Gen'+str(g)+'_'+scol + '_scaling_stddev'] = [np.nanstd(well_2490_s)]*len(freq_data[v])
    fd = freq_data['VGFA1'].merge(freq_data['VGFA2'], on='Well', how='inner')
    
    ## Adding strain / well annotation
    bwi = pd.read_csv('VLTE_by_well_info.csv')
    exclude_dict = {i[0]: i[1:] for i in bwi[bwi['plate']==plate].as_matrix(['well', 'contam', 'strain'])}
    fd['strain'] = fd['Well'].apply(lambda w: get_strain(w, exclude_dict))
    fd.to_csv(outname, index=False)
    
    ## Plotting s correlations and s over time
    color_by_strain = {
        'a': '#FFB000',
        'diploid': 'k',
        'alpha': '#648FFF',
    }
    f, subps = pl.subplots(2, 2, figsize=(14, 10))
    subs = subps[0]
    traj_subs = subps[1]
    pl.subplots_adjust(hspace=0.3, wspace=0.3)
    td = fd[fd['strain']!='BAD']
    for p in range(2):
        for row in fd.as_matrix(['strain'] + ['Gen'+str(g)+'_'+s_cols[p] + '_scaled' for g in all_gens]):
            if row[0] in color_by_strain:
                traj_subs[p].plot(all_gens, row[1:], color=color_by_strain[row[0]])
    
        subs[p].plot([-0.2, 0.14], [-0.2, 0.14], linestyle='dashed', c='k', alpha=0.5)
        for gen in all_gens:
            subs[p].scatter(td['Gen' + str(gen) + '_' + s_cols[p] + '_R1'], td['Gen' + str(gen) + '_' + s_cols[p] + '_R2'], s=20, alpha=0.6, label='Gen'+str(gen))
        subs[p].set_title(plate, fontsize='20', y=0.8)
        subs[p].tick_params(axis='both', which='major', labelsize=16)
        subs[p].set_xlabel('R1_' + s_cols[p], fontsize=16)
        subs[p].set_ylabel('R2_' + s_cols[p], fontsize=16)
        subs[p].set_ylim([-0.20, 0.14])
        subs[p].set_xlim([-0.20, 0.14])
    subs[1].legend()
    sns.despine()
    f.savefig(outcorr, background='transparent', bbox_inches='tight', pad_inches=0.1)
        
    for row in range(8):
        for col in range(12):
            well = chr(row+65) + str(col+1).zfill(2)
            plot_well_report(well, outname='browser3/report_graphs/' + plate + '_' + well + '.png')