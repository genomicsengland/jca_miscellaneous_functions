getting_disease_type <- function(samples){
  start_time <- Sys.time()
  print(start_time)
  sample_disease_type_row <- data.frame(sample_name = character(0),
                                        disease = character(0),
                                        disease_subtype = character(0),
                                        row.names = NULL)
  sample_disease_type <- sample_disease_type_row
  
  sample_disease_type <- do.call(rbind,lapply(samples$annotationSets, function(x){
    if(is.data.frame(x) && nrow(x) > 0){
      if(x %>% pluck("annotations") %>% pluck(1) %>% filter(name == "sampleId") %>% select(value) %>% nrow > 0){
        sample_name <- x %>% pluck("annotations") %>% pluck(1) %>% filter(name == "sampleId") %>% select(value)
      }else{
        sample_name <- "NA"
      }
      if(x %>% pluck("annotations") %>% pluck(1) %>% filter(name == "diseaseType") %>% select(value) %>% nrow > 0){
        disease <- x %>% pluck("annotations") %>% pluck(1) %>% filter(name == "diseaseType") %>% select(value)
      }else{
        disease <- "NA"
      }
      if(x %>% pluck("annotations") %>% pluck(1) %>% filter(name == "diseaseSubType") %>% select(value) %>% nrow > 0){
        disease_subtype <- x %>% pluck("annotations") %>% pluck(1) %>% filter(name == "diseaseSubType") %>% select(value)
      }else{
        disease_subtype <- "NA"
      }
      sample_disease_type_row <- data.frame(sample_name = sample_name,
                                            disease = disease,
                                            disease_subtype = disease_subtype,
                                            row.names = NULL)
      names(sample_disease_type_row) <- c("sample_name", "disease", "disease_subtype")
    }
    return(sample_disease_type_row)
  }))
  end_time <- Sys.time()
  print(end_time)
  print(end_time - start_time)
  return(sample_disease_type)
}