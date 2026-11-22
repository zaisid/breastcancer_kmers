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

    #accession_id  <- list.files(data_dir) # ID that's being used to label the download
    accession_dir <- file.path(data_dir,sra_id)
    
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
    
    cat("File",sra_id,"downloaded")
    
    return(accession_dir)
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
        message("Deleting original FASTQ file...")}
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
    
    # Output file
    sampled_fastq_file <- file.path(
        accession_dir,
        paste0(accession_id, "_sample",
               as.character(sample_size),
               ".txt")
    )
    
    # Temporary files
    index_file <- tempfile()
    sorted_index_file <- tempfile()
    
    # Count total sequences without loading them into R
    line_count <- system2(
        "wc",
        args = c("-l", shQuote(cleaned_fastq_file)),
        stdout = TRUE)
    
    total_lines <- as.numeric(
        strsplit(trimws(line_count[1]), "\\s+")[[1]][1])
    
    # Validate sample size
    if (sample_size > total_lines) {
        stop("Sample size exceeds number of available sequences.")}
    
    # Randomly select line numbers
    status <- system2(
        "gshuf",
        args = c(
            "-i", paste0("1-", total_lines),
            "-n", as.character(as.integer(sample_size))),
        stdout = index_file
    )
    
    if (status != 0) stop("Random sampling failed.")
    
    # Sort selected line numbers
    status <- system2(
        "sort",
        args = c("-n", index_file),
        stdout = sorted_index_file
    )
    
    if (status != 0) stop("Sorting sampled indices failed.")
    
    # Extract only selected sequences
    awk_script <- paste(
        'BEGIN { i=1 }',
        'NR==FNR { selected[++n]=$1; next }',
        'n > 0 && FNR == selected[i] { print; i++ }',
        sep = " "
    )
    
    status <- system2(
        "awk",
        args = c(
            shQuote(awk_script),
            shQuote(sorted_index_file),
            shQuote(cleaned_fastq_file)
        ),
        stdout = sampled_fastq_file
    )
    
    # Cleanup
    unlink(c(index_file, sorted_index_file))
    
    if (status != 0) {
        stop("Sequence extraction failed.")}
    
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
            message("Processed ", format(chunk_num * chunk_size,
                                         scientific=FALSE),
                    " sequences")}
        
        rm(seq_chunk, dna_chunk, chunk_counts)
    }
    close(con)
    
    # Normalize into frequencies 
    kmer_freq <- kmer_sums / sum(kmer_sums)
    
    return(list(
        kmer_freq = kmer_freq,
        kmer_sums = kmer_sums
    ))
}



load_and_kmerize <- function(sra_id, data_dir, verbose = TRUE,
                             sample_size = 1e6, k = 2,
                             chunk_size = 50000, sums = FALSE) {
    # full pipeline, from SRA ID to kmer vector
    # outputs kmer_freq vector (also kmer_sums if sum is set to TRUE; default FALSE)
    
    # Helper for consistent error handling
    handle_error <- function(stage, error) {
        message("[ERROR] ", sra_id, " | ", stage, ": ",
                conditionMessage(error))
        return(NULL)
    }
    
    # -------------------- Check for existing sample ---------------------
    
    # same path construction as in sample_SRA()
    existing_sample_file <- file.path(
        data_dir, sra_id,
        paste0(sra_id, "_sample", as.character(sample_size), ".txt"))
    
    sample_exists <- FALSE
    if (file.exists(existing_sample_file)) {
        line_count <- system2(
            "wc",
            args = c("-l", shQuote(existing_sample_file)),
            stdout = TRUE)
        n_lines <- as.numeric(
            strsplit(trimws(line_count[1]), "\\s+")[[1]][1])
        sample_exists <- (n_lines == sample_size)
    }
    
    if (sample_exists) {
        
        if (verbose) {
            message("Sample file already present for ", sra_id,
                    "; skipping to k-merizing...")}
        sampled_fastq_file <- existing_sample_file
        
    } else {
        
        # -------------------- Downloading ---------------------
        
        if (verbose) {
            message("Downloading file: ", sra_id)}
        
        accession_dir <- tryCatch(
            load_SRA(sra_id, data_dir, verbose),
            error = function(e) handle_error("SRA download", e))
        
        if (is.null(accession_dir)) return(NULL)
        
        # -------------------- Cleaning ---------------------
        
        if (verbose) {
            message("Cleaning FASTQ data...")}
        
        cleaned_fastq_file <- tryCatch(
            clean_SRA(sra_id, accession_dir, verbose),
            error = function(e) handle_error("SRA cleaning", e)
        )
        
        if (is.null(cleaned_fastq_file)) return(NULL)
        
        # -------------------- Sampling ---------------------
        
        if (verbose) {
            message("Sampling ", format(sample_size, scientific = FALSE),
                    " sequences... ")}
        
        sampled_fastq_file <- tryCatch(
            sample_SRA(sra_id, accession_dir,
                       cleaned_fastq_file, sample_size),
            error = function(e) handle_error("SRA sampling", e))
        
        if (is.null(sampled_fastq_file)) return(NULL)
    }
    
    # -------------------- K-merizing ---------------------
    
    if (verbose) {
        message("K-merizing...")}
    
    kmer_results <- tryCatch(
        kmerize_SRA(sampled_fastq_file, k, chunk_size, verbose),
        error = function(e) handle_error("K-merization", e))
    
    if (is.null(kmer_results)) return(NULL)
    
    # -------------------- Output ---------------------
    
    if (sums) {
        return(kmer_results)}
    
    return(kmer_results$kmer_freq)
}

#system.time(load_and_kmerize(test_id,data_dir))
#    user  system  elapsed 
# 481.727  64.875  352.151 (old)
# 350.685  33.849  337.301 (new)



# --------------------- Handling List of IDs ---------------------

id_path <- "../breastcancer_kmers/data/abc_ids.txt"

kmerize_id_list <- function(id_path, class_label, verbose = TRUE,
                            sample_size = 1e6, k = 2,
                            chunk_size = 50000, sums = FALSE) {
    # PARAMS: id_path <- file path of deliminated id list
    #         class_label <- (character) class assigned to every id in this list
    # OUTPUT: results_df of kmer_freqs across each id in id_path file
    ids <- read.delim(id_path, header = FALSE)[[1]]
    
    # create results df
    kmers <- oligonucleotideFrequency(
        DNAString("ACGT"),
        width = k
    ) |> names() # get row names (ids)
    
    # Initialize empty results dataframe
    results_df <- data.frame(
        matrix(
            nrow = 0,
            ncol = length(kmers),
            dimnames = list(NULL, kmers)
        ),
        class = character(0),
        check.names = FALSE)    
    
    for (id in ids) {
        data_dir <- file.path(dirname(id_path),"prefetched")
        
        kmer_freq <- load_and_kmerize(id, data_dir, verbose, 
                                      sample_size, k, chunk_size, sums)
        
        if (!is.null(kmer_freq)) {
            results_df[id, kmers] <- as.list(kmer_freq)
            results_df[id, "class"] <- class_label}
    }
    return(results_df)
}


#abc_kmers <- kmerize_id_list(id_path,"ABC")
#abc_kmers

#abca_kmers <- kmerize_id_list("../breastcancer_kmers/data/abca_ids.txt","ABCA")
#abca_kmers

# save as .csv
#k2_s1e6_freqs <- rbind(abc_kmers,abca_kmers)
#write.csv(k2_s1e6_freqs, 
#          "../breastcancer_kmers/data/k2_s1e6_freqs.csv", 
#          row.names = TRUE)

