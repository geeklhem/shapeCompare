#OPTIONS
#STORE OPTIONS 


# --------------------------------------------------------------------------
# LOADING OPTIONS
# --------------------------------------------------------------------------

# Fasta file contain all RNA sequences you could feed to the software, 
# By default the first one is used.
# You can specify the sequence to use for each *.shape file by creating
# a "seq.txt" file inside containing a name matching one of the sequences names 
# of this file. 
FASTA_FILE = "../data/sequences/7SK.fasta"

# Same file but with "-" for sequence alignement. You can generate it with muscle. (http://drive5.com/muscle/ or http://www.ebi.ac.uk/Tools/msa/muscle/)
FASTA_ALG = "../data/sequences/7SK_alg.fasta"

# Minimal proportion of bases to match between the .shape and the .fasta file.
MATCH_THRESHOLD = 0.8
