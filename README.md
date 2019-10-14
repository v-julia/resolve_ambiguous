# resolve_ambiguous
Script that resolves ambigugous nucleotides in sequences from alignment in fasta-format using the rest sequences:
 - 	discards sequences with too many ambiguous characters (e.g. 0.2%), which may be a sign of poor sequence quality rather than natural heterogeneity; 
 - 	in the remaining sequences, identifies a region (e.g. 100 nt) with an ambiguous character and blasts it against the dataset;
 - 	identifis the frequency of different nucleotides at this position in the most closely related sequences;
 - 	replaces the ambiguous character with the most common nucleotide among the most closely related sequences.

