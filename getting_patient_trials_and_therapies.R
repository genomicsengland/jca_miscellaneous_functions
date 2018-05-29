require(tidyverse)

getting_patient_trials_and_therapies <- function(filepath){

  #these column names may need to be adapted to the input file
  col_names <- c("sample", "disease_type", "disease_subtype", "tumour_type", "clinical_sample_date_time", "coords", "gene_name", "ensembl_id", "allele_origin", "tier", "cdna_change", "protein_change")
  
  #filepath is to output of get_info....py scripts on LSF; colnames below may need to be edited based on LSF script output
  ncols <- max(count.fields(filepath, sep = "\t"))
  get_vardata <- function(filepath, ncols){
    column.names <- c(col_names, paste0("V_", seq_len(ncols - length(col_names))))
    
    #read in table
    vardata <- read.table(filepath, sep = "\t", col.names = column.names, fill = T, stringsAsFactors = F) %>%
      mutate(disease_type = gsub(" ", "_", disease_type), disease_subtype = gsub(" ", "_", disease_subtype))
    
    return(vardata)
  }
  
  #read in the reference table that Nirupa made that links the tumour types in the Genomoncology links to the GEL tumour types
  keytable <- read.table("/Users/johnambrose/Documents/trial_counting/trial.counting/genomoncology_gel_linkingtable.tsv", sep = "\t", header = T)
  
  #get rid of unneeded cols
  vardata <- get_vardata(filepath, ncols)
  vardata <- vardata %>% select(-c(coords, ensembl_id, cdna_change, protein_change))
  
  #convert everything to uppercase to make matching easier
  for(i in 1:ncol(vardata)){
    vardata[,i] <- toupper(vardata[,i])
  }
  for(i in 1:ncol(keytable)){
    keytable[,i] <- toupper(keytable[,i])
  }
  
  #column number for the start of the url columns
  urlstart <- which(names(vardata) == "V_1")
  urlend <- ncol(vardata)
  
  #column numbers to provide the words to match
  matchindices <- which(names(vardata) == "disease_type" | names(vardata) == "disease_subtype")
  
  #function to count number of occurences of certain words in given vector
  #filter bit splits out different flavours of the trial links (those beginning with Trial vs therapeutic)
  #basically matchpattern generates the start of the trial link with various options of the genomoncology tumour type
  # pattern is TRIAL/THERAPEUTIC (<tumour type>): then False or True for gene level or variant level match
  # based on the gel tumour type, then we iterate through each of those and get the matches, difficult doing them
  # all together as the matchpattern has brackets in it so get interpreted as regex stuff
  wordmatch <- function(filter, words, input, varlevel = "FALSE"){
    matchpattern <- paste0(filter, " (", keytable$Genomoncology[keytable$gel %in% words], "):", varlevel)
    matches <- c()
    for(i in matchpattern){
      matches <- c(matches, grep(i, input, fixed = T))
    }
    return(length(matches))
  }
  
  #do the counting on each row
  trial_counts_disease_type.gene <- apply(vardata, 1, function(x) wordmatch("TRIAL", x[matchindices], x[urlstart:urlend]))
  trial_counts_disease_type.var <- apply(vardata, 1, function(x) wordmatch("TRIAL", x[matchindices], x[urlstart:urlend], varlevel = "TRUE"))
  thera_counts_disease_type.var <- apply(vardata, 1, function(x) wordmatch("THERAPEUTIC", x[matchindices], x[urlstart:urlend], varlevel = "TRUE"))
  
  #Combine with original data
  vardata <- vardata[,-c(urlstart:urlend)]
  vardata$trial_counts_disease_type_gene <- trial_counts_disease_type.gene + trial_counts_disease_type.var
  
  return(vardata)
}