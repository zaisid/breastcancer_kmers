################################################
##    Created by ZAINAB SIDDIQUI              ##
##    On 11/16/2026                           ##
##    Building Full SRA Accession             ##
##       Loading Pipeline                     ##
################################################

# Converting/editing SRAaccession_test.R to loop over all samples

library(data.table)
library(Biostrings)

# ---------------------- Modularize -----------------------

test_id <- read.delim("../breastcancer_kmers/data/abc_ids.txt",
                      header = FALSE)[2,] # SRR41124800

data_dir <- file.path("..",       #../breastcancer_kmers/data/prefetched
                      "breastcancer_kmers",
                      "data",
                      "prefetched") 

load_SRA <- function(sra_id,data_dir,verbose=TRUE) {
    ##PARAMS:
    # sra_id   = (character) ID of SRA accession file to load in
    # data_dir = (character) file path (e.g., using file.path()) of directory
    #            to download SRA files to
    # verbose  = (boolean, default TRUE) controls whether system messages are shown
    # 
    ##RETURNS:
    # accession_id  = (character) ID of SRA accession file downloaded
    # accession_dir = (character) file path of directory with SRA files
    
    if (verbose) {
        message("Downloading SRA file...")}
    system(paste("prefetch", sra_id,"-O",data_dir))

    accession_id  <- list.files(data_dir) # ID that's being used to label the download
    accession_dir <- file.path(data_dir,accession_id)
    
    # .sra file
    fetched_file <- file.path(accession_dir,list.files(accession_dir))
    
    # extract with fastq-dump
    if (verbose){
        message("Extracting FASTQ data...")}
    system(paste("fasterq-dump --split-files", fetched_file,"-O",accession_dir))
    # delete R2
    if (verbose){
        message("Deleting reverse read file...")}
    system(paste("rm",file.path(accession_dir,
                                list.files(accession_dir,pattern="_2\\.fastq$")))) 
    
    if (verbose){
        message("Deleting SRA file...")} 
    system(paste("rm",fetched_file)) 
    
    cat("File",accession_id,"downloaded")
    
    return(accession_id,accession_dir)
}


clean_SRA <- function(accession_id,accession_dir,verbose=TRUE) {
    ##PARAMS:
    # accession_id  = (character) ID of SRA accession file downloaded
    # accession_dir = (character) file path of directory with SRA files
    # verbose       = (boolean, default TRUE) controls whether system messages are shown
    # 
    ##RETURNS:
    # cleaned_fastq_file = (character) file path of raw .txt of accession file
    
    uncleaned_fastq_file <- file.path(accession_dir, 
                                      list.files(accession_dir,
                                                 pattern = "_1\\.fastq$")
    )
    
    cleaned_fastq_file <- file.path(accession_dir,paste0(accession_id,"_seqonly.txt"))
    
    # take the 2nd row of every 4 to isolate sequences
    if (verbose) {
        message("Extracting sequences from FASTQ dump...")}
    system(paste("awk 'NR % 4 == 2'", uncleaned_fastq_file,">",cleaned_fastq_file))

    # delete original file
    if (verbose) {
        message("Deleting original fastq file...")}
    system(paste("rm",uncleaned_fastq_file))
    
    return(cleaned_fastq_file)
}

sample_SRA <- function(accession_id,
                       accession_dir,
                       cleaned_fastq_file,
                       sample_size=1e6) {
    ##PARAMS:
    # accession_id       = (character) ID of SRA accession file downloaded
    # accession_dir      = (character) file path of directory with SRA files
    # cleaned_fastq_file = (character) file path of raw .txt of accession file
    # sample_size        = (double) file path of raw .txt of accession file
    # 
    ##RETURNS:
    # sampled_fastq_file = (character) path of sampled .txt lines from SRA
    
    # random sampling of sequence-lines
    n = 1e6 # num lines to sample
    sampled_fastq_file <- file.path(accession_dir,
                                    paste0(accession_id,"_sample",
                                           gsub("\\+0", "", 
                                                format(n, scientific = TRUE)),
                                           ".txt"))
    
    lines <- fread(cleaned_fastq_file, header = FALSE)[[1]] 
    sampled_lines <- sample(lines,n) 
    rm(lines) 
    
    writeLines(sampled_lines, sampled_fastq_file)
    
    return(sampled_fastq_file)
}


kmerize_SRA <- function(sampled_fastq_file,k=2,chunk_size=50000,
                        verbose = TRUE) {
    ##PARAMS:
    # sampled_fastq_file = (character) path of sampled .txt lines from SRA
    # k                  = (double, default = 2) k-mer size
    # chunk_size         = (double, default = 50,000) number of sequences/lines
    #                      to be loaded in at a time; tune to memory capability
    # verbose            = (boolean, default TRUE) controls whether system 
    #                      messages are shown
    # 
    ##RETURNS:
    # kmer_freq = frequency vector of kmers in given file
    # kmer_sums = raw sums  vector of kmers in given file
    
    # Initialize aggregate counts
    kmer_sums <- numeric(4^k)
    names(kmer_sums) <- oligonucleotideFrequency(
        DNAString("ACGT"),
        width = k
    ) |> names()
    
    # Open sampled sequence file
    con <- file(sampled_fastq_file, open = "r")
    chunk_num <- 0
    
    repeat {
        # Read next chunk
        seq_chunk <- readLines(con, n = chunk_size)
        
        # Stop when file is exhausted
        if (length(seq_chunk) == 0) break
        chunk_num <- chunk_num + 1
        
        # Convert to DNA sequences
        dna_chunk <- DNAStringSet(seq_chunk)
        
        # Aggregate k-mer counts across chunk
        chunk_counts <- oligonucleotideFrequency(
                        dna_chunk,
                        width = k,
                        simplify.as = "collapsed")
        
        # Accumulate counts
        kmer_sums[names(chunk_counts)] <-
            kmer_sums[names(chunk_counts)] + chunk_counts
        
        # Progress indicator
        if (chunk_num %% 10 == 0 && verbose) {
            message("Processed ", chunk_num * chunk_size,
                    " sequences")}
        
        rm(seq_chunk, dna_chunk, chunk_counts)
    }
    close(con)
    
    # Normalize into frequencies 
    kmer_freq <- kmer_sums / sum(kmer_sums)
    
    return(kmer_freq,kmer_sums)
}











