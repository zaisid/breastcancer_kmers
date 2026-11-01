###########################################
##    Created by ZAINAB SIDDIQUI         ##
##    On 10/26/2026                      ##
##    Test of K-merization Packages      ##
###########################################

# imports
install.packages("kmer")
library(kmer)

# DATA LOADING

# test sequence, approx. ~20 kbp
data_dir <- "~/Downloads/Crocs/Data/synteny/regions/gator/LOC_to_ARHG.fa"

# preprocess to get only the nucleotide data
og_seq <- read_file(data_dir)
seq <- unlist(strsplit(og_seq, "\n"))[-1]
seq[1:4]

test_seq <- seq[1]
cat(test_seq)
# Output:
# GGAAAGCAGGACTTTATTGGCAGCAAGATCTCCCCCCCACCCCCGCATAAAAACTCCATT

#convert to character-by-character list
test_nuc <- unlist(strsplit(test_seq, split = ""))
print(test_nuc)

# k=2? --> dinucleotides are important for GC analysis or smthn
kcount(test_nuc, k=2)
# Output:
#    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
# S1  7  3  4  4  6 11  1  3  3  4  3  0  2  3  1  4


# REAL SAMPLE TEST

# load data
small_sample <- readLines("../breastcancer_kmers/data/SRR41124782_smallsample.txt")

# row vector for holding kmer sum for sample
k=2
kmers_sum <- matrix(rep(0,k^4),nrow=1) 

# loop over lines in txt file & sum
for( line in small_sample ) {
    cat("Sequence: ")
    print(substr(line,1,10))
    
    # convert to list of bps
    line_list <- unlist(strsplit(line, split = ""))
    
    kmer_counts <- kcount(line_list, k=2)
    kmers_sum <- kmers_sum + kmer_counts
}
print(kmers_sum)
# Output:
#     AA AC  AG AT CA CC CG CT GA GC GG GT TA TC  TG  TT
# S1 169 64 112 95 96 87 34 92 98 74 94 93 76 86 117 103

# ^^ Successfully got sum of kmers for the same sample across different sequences













