getting_library_type <- function(con, list, study, germline_study){
  library(opencgaR)
  library(getPass)
  
  # con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")
  # con <- opencgaLogin(opencga = con, userid = "jambrose", passwd = getPass())
  
  pca_stats_row <- data.frame(sample_name = character(0),
                              library_type = character(0),
                              row.names = NULL)
  pca_stats <- pca_stats_row
  
  pca_stats <- do.call(rbind, lapply(list, function(name){
    file <- try(fileClient(OpencgaR = con, action = "search", params = list(study = study, name = paste0(name,".bam"), uri="~file:///genomes/by_date/", include = "name,uri,stats,attributes")), silent = TRUE)
    if(is.data.frame(file) && nrow(file) > 0){
      print(file$name)
      if(!is.null(file$stats$SAMPLE_METADATA$sample_library_type)){
        library_type <- file$stats$SAMPLE_METADATA$sample_library_type
      }else{
        library_type <- "NA"
      }
      pca_stats_row <- data.frame(sample_name = name,
                                  library_type = library_type,
                                  row.names = NULL)
    }else{
      message <- paste0("not_in_opencga")
      pca_stats_row <- data.frame(sample_name = name,
                                  library_type = message,
                                  row.names = NULL)
    }
    return(pca_stats_row)
  }))
  return(pca_stats)
}