getting_patients_for_focus4 <- function(){
  #getting patients with suitable mutations for focus4 trial (BRAF, KRAS/NRAS, PIK3CA)
  
  #N.B. THE 7 SAMPLES THAT ARE MISSING FROM THE VARDATA COMPARED TO THE MUTATION BURDEN AND OTHER FILES ARE BELOW
  #THEY DO NOT SEEM TO HAVE ANY TIER1 VARIANTS WHICH IS WHY THEY ARE NOT OUTPUT BY ~/jca_scripts/get_info.v2.py:
  # vardata_name <- vardata[,1:5] %>% filter(ttype == "COLORECTAL") %>% mutate(name = paste0(sample,"_",ttype)) %>% select(name) %>% unique()
  # patients_and_ldp_codes_name <- patients_and_ldp_codes %>% filter(tissue == "COLORECTAL") %>% mutate(name = paste0(sample,"_",tissue)) %>% select(name) %>% unique()
  # patients_and_ldp_codes_name$name[which(!(patients_and_ldp_codes_name$name %in% vardata_name$name))] %>% matrix()
  # [1,] "LP3000240-DNA_G08_COLORECTAL"
  # [2,] "LP3000572-DNA_B01_COLORECTAL"
  # [3,] "LP3000608-DNA_A02_COLORECTAL"
  # [4,] "LP3000343-DNA_A01_COLORECTAL"
  # [5,] "LP3000208-DNA_A05_COLORECTAL"
  # [6,] "LP3000241-DNA_A03_COLORECTAL"
  # [7,] "LP3000241-DNA_F04_COLORECTAL"
  
  #set up vardata file inputs from various runs
  datafilename_alona <- "/Users/johnambrose/Documents/trial_counting/trial.counting/variants.txt"
  datafilename_alona_batch4 <- "/Users/johnambrose/Documents/trial_counting/trial.counting/variants.batch4.txt"
  datafilename_cip_api <- "/Users/johnambrose/Documents/trial_counting/trial.counting/variants.05022018.txt"
  
  ncols <- max(max(count.fields(datafilename_alona, sep = "\t")), max(count.fields(datafilename_alona_batch4, sep = "\t")), max(count.fields(datafilename_cip_api, sep = "\t")))
  get_vardata <- function(datafilename, ncols){
    column.names <- c("sample", "ttype", "subttype", "coords", "gene", "ensemble", "tier", "cDNAchange", "protchange", paste0("V_", seq_len(ncols - 9)))
    
    #read in table
    vardata <- read.table(datafilename, sep = "\t", col.names = column.names, fill = T, stringsAsFactors = F) %>%
      dplyr::mutate(ttype = gsub(" ", "_", ttype), subttype = gsub(" ", "_", subttype))
    
    return(vardata)
  }
  
  vardata_combined <- rbind(get_vardata(datafilename_alona, ncols), get_vardata(datafilename_alona_batch4, ncols), get_vardata(datafilename_cip_api, ncols))
  vardata <- vardata_combined
  
  #read in the reference table that Nirupa made that links the tumour types in the Genomoncology links to the GEL tumour types
  keytable <- read.table("/Users/johnambrose/Documents/trial_counting/trial.counting/genomoncology_gel_linkingtable.tsv", sep = "\t", header = T)
  
  #get rid of crap
  drops <- c("coords", "cDNAchange", "protchange")
  vardata <- vardata[,!names(vardata) %in% drops]
  
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
  matchindices <- which(names(vardata) == "ttype" | names(vardata) == "subttype")
  
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
  trialcountsttype.gene <- apply(vardata, 1, function(x) wordmatch("TRIAL", x[matchindices], x[urlstart:urlend]))
  trialcountsttype.var <- apply(vardata, 1, function(x) wordmatch("TRIAL", x[matchindices], x[urlstart:urlend], varlevel = "TRUE"))
  theracountsttype.var <- apply(vardata, 1, function(x) wordmatch("THERAPEUTIC", x[matchindices], x[urlstart:urlend], varlevel = "TRUE"))
  
  #Combine with original data
  vardata <- vardata[,-c(urlstart:urlend)]
  vardata$trial_counts_ttype_gene <- trialcountsttype.gene + trialcountsttype.var
  
  vardata_filtered <- vardata %>%
    filter(ttype == "COLORECTAL" & (gene == "BRAF" | gene == "KRAS" | gene == "NRAS" | gene == "PIK3CA") & trial_counts_ttype_gene > 0) %>% 
    select(-subttype, -ensemble, -tier)
  
  return(vardata_filtered)
}