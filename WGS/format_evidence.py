import pandas as pd

def alignment(r):
    # uses the cigar to format the sequence to show in an alignment
    # adds - for deletions, puts insertions in parentheses (will be shown only by a red line that shows the insertion on hover)
    seq, cigar = r['SEQ'], r['CIGAR']
    new_seq = ''
    current_nums = ''
    current_spot = 0
    for character in cigar:
        if character in ['M', 'D', 'I', 'S']:
            if character == 'S':
                return new_seq
            elif character == 'M':
                new_seq += seq[current_spot:current_spot+int(current_nums)]
                current_spot += int(current_nums)
            elif character == 'D':
                new_seq += '-'*int(current_nums)
            else:
                new_seq += '(' + seq[current_spot:current_spot+int(current_nums)] + ')'
                current_spot += int(current_nums)
            current_nums = ''
        else:
            current_nums += character
    return new_seq


td = pd.read_csv('../../Output/Browser/evidence_sams/P1B02.tsv', delimiter='\t')
td = td[td['CIGAR'].apply(lambda c: 'S' in c or 'I' in c or 'D' in c)].iloc[:500]
td['a'] = td.apply(lambda row: alignment(row), axis=1)
td.to_csv('test.tsv', sep='\t')