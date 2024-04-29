import argparse
from Bio import SeqIO
from matplotlib import pyplot as plt

'''
 Plots distribution of ambiguous nucleotides in sequences


Input:
input_file - file with sequences in fasta format

Output:

distribution olot for every type of ambiguous nucleotide
'$input_file'_res.fasta - file with sequences that have no ambiguous nucleotides

'''
def search_amb(input_file):


    #list with ambiguous nucleotides
    ambig_nt = ['n',  # A, T, G, C 
                'r', # Purine (A or G)
                'y', # Pyrimidine (T or C)
                'k', # Keto (G or T)
                'm', # Amino (A or C)
                's', # Strong interaction (3 H bonds) (G or C)
                'w', # Weak interaction (2 H bonds) (A or T)
                'b', # Not A (C or G or T)
                'd', # Not C (A or G or T)
                'h', # Not G (A or C or T)
                'v', # Not T or U (A or C or G)
                'total'
                ]
    # counts of each type of ambiguous nucleotide
    amb_counts = {}

    for nt in ambig_nt:
        amb_counts[nt] = []
    
    # checks presence of ambiguous nucleotides in sequences from input_file
    new_records = []
    with open(input_file) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        print(len(records))
        for record in records.copy():
            # total number of ambiguous nucleotides
            total_amb = 0
            for nt in amb_counts.keys():
                amb_counts[nt].append(record.seq.lower().count(nt))
                total_amb += record.seq.lower().count(nt)
            amb_counts['total'].append(total_amb)
            
            # filter sequences with ambiguous nucleotides
            if total_amb < 1:
                new_records.append(record)
            else:
                print(record.id + ' ' + str(total_amb))
    SeqIO.write(new_records, input_file.replace('.fasta', '_res.fasta'),format = "fasta")
    
    # plots distribution histograms
    for key in ambig_nt:
        if len(amb_counts[key]) != 0:
            plt.hist(amb_counts[key], list(range(max(amb_counts[key])+5)), log = True)
            plt.title('Distribution of {}'.format(key))
            plt.xlabel('Number of ambiguous characters in sequence')
            plt.ylabel('Number of sequences')
            plt.show()
    handle.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    args = parser.parse_args()

    search_amb(args.input_file)