import argparse
import copy
import os
import pandas as pd
import re
import subprocess
import sys
from Bio import SeqIO

#list with ambiguous nucleotides
ambig_nt = [
    'n', 'N',  # A, T, G, C
    'r', 'R',  # Purine (A or G)
    'y', 'Y',  # Pyrimidine (T or C)
    'k', 'K',  # Keto (G or T)
    'm', 'M',  # Amino (A or C)
    's', 'S',  # Strong interaction (3 H bonds) (G or C)
    'w', 'W',  # Weak interaction (2 H bonds) (A or T)
    'b', 'B',  # Not A (C or G or T)
    'd', 'D',  # Not C (A or G or T)
    'h', 'H',  # Not G (A or C or T)
    'v', 'V',  # Not T or U (A or C or G)
    'total'
]

def resolve_ambiguous(input_file, output_dir, window, path_to_blast,
                     evalue, word_size, max_ambiguous, max_ambiguous_row):
    '''
    Resolves ambiguous nucleotides in nucleotide sequences according to
    consensus in most related sequences to the region with ambiguous nt.
    Removes sequences which have more than user-specified amount of 
    ambiguous nucleotides total/in a row.

    Input:
        input_file - path to file with alignment of nt sequences in fasta-format
        path_blast - path to bin directory of BLAST
    '''
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        output_dir += "\\"
    else:
        output_dir += "/"

    #alignment with nucleotide sequences in fasta-format
    fasta_al = list(SeqIO.parse(open(input_file), "fasta"))

    #the length of window around the ambiguous nucleotide to cut from original sequence
    #window = 500

    print('------Finding sequences with ambiguous characters------')

    #list with SeqIO objects which sequences contain less ambiguous characters that specified threshold
    fasta_al_less_amb = []

    # list with SeqIO objects which sequences 
    # correspond to slice surrounding ambiguous character
    list_slices = []
    ambiguous_nucleotides = "RYWSKMBDHVN"
    pattern = r"[" + ambiguous_nucleotides + "]{" + str(max_ambiguous_row) + ",}"

    for rec in fasta_al.copy():
        # total number of ambiguous nucleotides in sequence
        amb_total = len(re.findall(r"[RYSWKMBDHVN]", str(rec.seq)))
        if amb_total == 0:
            # adds records with no ambiguous characters to the new alignment
            fasta_al_less_amb.append(rec)
        else:
            # checking whether the number of ambiguous characters exceed specified threshold
            rec_seq_len = len(re.sub("-","", str(rec.seq)))
            print('The number of ambiguous nucleotides in {}: {} ({}%)'.format(rec.id, amb_total, round(amb_total/rec_seq_len,2)))
            if amb_total > max_ambiguous:
                print(rec.name, 'exceeded threshold')
                continue
            elif (re.findall(pattern, str(rec.seq))):
                print(rec.name, 'has {} or more n-nt in a row'.format(max_ambiguous_row))
                continue
            
            else:
                print(rec.name)
                # if the number of ambiguous nucleotides doesn't exceed the threshold
                # add copy record to a new list
                fasta_al_less_amb.append(rec)
                # finds positions of ambiguous nucleotides in sequence
                starts = [m.start() for m in re.finditer(r"[RYSWKMBDHVN]", str(rec.seq))]
                # for each ambiguous nt creates a slice with length=window surrounding this nt
                i=0
                # list with start and end positions of slices
                slices = []
                while i < len(starts):
                    # starts of ambiguous nt in current slice
                    current_starts = []
                    current_starts.append(str(starts[i]+1))
                    # takes the start of sequence if amb nt is closer than window/2
                    #  to the beginning of seq
                    if starts[i] < window/2:
                        st = 0
                        e = window
                    else:
                            # takes the end of the sequence if amb nt is closer than window/2
                        # to the end of sequences
                        if len(rec.seq) - starts[i] < window/2:
                            st = len(rec.seq) - window -1
                            e = len(rec.seq)
                            
                        # takes the slice [starts[i]-window/2, starts[i]+window/2]
                        else:
                            st = int(starts[i]-window/2)
                            e = int(starts[i]+window/2)

                    # creates slice object
                    cur_slice_rec = copy.deepcopy(rec[st:e])
                    cur_slice_rec.description = ''
                
                    #appends sliced sequence to the list
                    list_slices.append(cur_slice_rec)
                
                    slices.append([st,e])

                    if i+1 < len(starts):
                        for j in range(i+1, len(starts),1):
                            if st+int(window/5)<starts[j] and starts[j]<e-int(window/5):
                                current_starts.append(str(starts[j]+1))
                                i += 1
                                if j == len(starts) - 1:
                                    i +=1
                                continue
                            else:
                                i = i + 1
                                break                    
                    else:
                        i += 1

                    if e != len(rec.seq):
                        cur_slice_rec.id = rec.id + "_" + ":".join([str(st+1)]+current_starts+[str(e+1)])
                    else:
                        cur_slice_rec.id = rec.id + "_" + ":".join([str(st+1)]+current_starts+[str(e)])

    # filename for fasta-file with slices
    file_name_slices = os.path.splitext(input_file)[0] + "_slices.fasta"
    # writes slices to fasta_file
    SeqIO.write(list_slices, file_name_slices, "fasta")

    file_name_less_amb = os.path.splitext(input_file)[0] + "_less_amb.fasta"
    with open(file_name_less_amb,'w') as file_less_amb:
        SeqIO.write(fasta_al_less_amb, file_less_amb, "fasta")
    file_less_amb.close()

    # commands for creating local database and blast slices against it
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        makeblast_command = '{}makeblastdb.exe -in {} -dbtype nucl -out {}local_db'.format(path_to_blast, file_name_less_amb, output_dir)
        print(makeblast_command)
        blastn_command = '{blast_path}blastn.exe -db {out_path}local_db -query {input} -outfmt 6 -out \
                            {out_path}blast.out -strand plus -evalue 1e-20 -word_size 7 -max_target_seqs 30'.format(blast_path = path_to_blast, \
                            input = file_name_slices, out_path = output_dir)
    else:
        makeblast_command = '{}makeblastdb -in {} -dbtype nucl -out {}local_db'.format(path_to_blast, file_name_less_amb, output_dir)
        blastn_command = '{blast_path}blastn -db {out_path}local_db -query {input} -outfmt 6 -out \
                            {out_path}blast.out -strand plus -evalue 1e-20 -word_size 7 -max_target_seqs 30'.format(blast_path = path_to_blast, \
                            input = file_name_slices, out_path = output_dir)
    subprocess.run(makeblast_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run(blastn_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    subprocess.call(makeblast_command, shell=True)

    print('------Creating blast database------')
    # blast against reference sequences
    with open(output_dir+'blast.err', 'w') as stderr_file:
        subprocess.call(blastn_command, shell=True, stderr=stderr_file)

    # dataframe with blast results

    blast_output = pd.read_csv(output_dir+'blast.out', sep='\t', header = None, \
                                names=['qseqid','sseqid','pident','length','mismatch',\
                                'gapopen','qstart','qend','sstart','send','evalue','bitscore'])
    # create dictionary with sequences instead of list
    fasta_al_less_amb = SeqIO.to_dict(fasta_al_less_amb)

    # flag indicates whether the sequence has been resolved
    flag = 0

    current_seq_id = ''
    print('------Resolving nucleotides------')
    for _, row in blast_output.iterrows():
        # Changes flag when meets new sequence with amb nt
        if row['qseqid'] != current_seq_id:
            current_seq_id = row['qseqid']
            print('---------')
            print(current_seq_id)
            # ID of sequence with amb nt according to original fasta
            current_seq_id_orig = "_".join(row['qseqid'].split('_')[:-1])
            # Start position of slice (enumeration from 0)
            start = int(row['qseqid'].split('_')[-1].split(':')[0]) - 1
            # End positions of slice
            end = int(row['qseqid'].split('_')[-1].split(':')[-1]) - 1
            # Ambiguous positions within slice
            left_amb_pos = [int(x) - 1 for x in row['qseqid'].split('_')[-1].split(':')[1:-1]]
            #print("left amb pos")
            #print(left_amb_pos)
            flag = 0
        # Skips the row if amb nts have been resolved
        if flag != 1:
            if row['sseqid'] == '_'.join(row['qseqid'].split('_')[:-1]):
                continue
            else:
                # Relative start of ambiguous character in window, counting from 0
                cur_start = start + int(row['qstart']) - 1
                rel_amb_pos = [x - cur_start for x in left_amb_pos]
                # Position corresponding to ambiguous character in reference sequence
                ref_pos = [int(row['sstart']) - 1 + x for x in rel_amb_pos]

                # Nucleotides in reference sequence in positions that are ambiguous
                ref_res_nuc = [fasta_al_less_amb[row['sseqid']].seq[x] for x in ref_pos]

                left_amb_pos_copy = left_amb_pos.copy()
                # Changes amb nt to the ones in the reference sequence
                for i in range(len(left_amb_pos)):
                    if ref_res_nuc[i] not in ambig_nt:
                        fasta_al_less_amb[current_seq_id_orig].seq = fasta_al_less_amb[current_seq_id_orig].seq[:left_amb_pos[i]] + ref_res_nuc[i] + fasta_al_less_amb[current_seq_id_orig].seq[left_amb_pos[i] + 1:]
                        left_amb_pos_copy.remove(left_amb_pos[i])
                    else:
                        continue
                left_amb_pos = left_amb_pos_copy[:]
                if len(left_amb_pos) == 0:
                    flag = 1
                else:
                    print(current_seq_id, 'reference has amb')

    # Filename for fasta-file with sequences resolved
    file_name_less_amb = os.path.splitext(input_file)[0] + "_less_amb.fasta"
    SeqIO.write(fasta_al_less_amb.values(), file_name_less_amb, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Path to file with alignment of nt sequences in fasta-format", required=True)
    parser.add_argument("-pout", "--path_out", type=str,
                        help="Path to directory for output alignment without amb nucleotides. The directory of input file by default.")
    parser.add_argument("-w", "--window", type=int, default=100,
                        help="Window size")
    parser.add_argument("-evalue", "--evalue", type=float, default=1e-20,
                        help="BLASTn E-value")
    parser.add_argument("-word_size", "--word_size", type=int, default=7,
                        help="BLASTn word size")
    parser.add_argument("-pb", "--path_blast", type=str,
                        help="Path to BLAST", required=True)
    parser.add_argument("-max_ambiguous", "--max_ambiguous", type=float, default=5,
                        help="Maximum number of ambiguous nucleotides allowed")
    parser.add_argument("-max_ambiguous_row", "--max_ambiguous_row", type=int, default=3,
                        help="Maximum number of consecutive ambiguous nucleotides allowed")

    args = parser.parse_args()

    args.input_file = os.path.realpath(args.input_file)
    if not args.path_out:
        args.path_out = os.path.split(args.input_file)[0]

    resolve_ambiguous(args.input_file, args.path_out, args.window, args.path_blast, args.evalue, args.word_size, args.max_ambiguous, args.max_ambiguous_row)

#path_to_blast = "J:\\Programs\\blast-2.11.0+\\bin\\"
