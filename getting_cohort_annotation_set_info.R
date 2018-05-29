getting_cohort_annotation_set_info <- function(con, study, cohorts, latest_cohorts_only = TRUE){
  
  library(janitor)
  library(getPass)
  library(dplyr)
  library(tidyr)
  library(opencgaR)
  
  # con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")
  # con <- opencgaLogin(opencga = con, userid = "jambrose", passwd = getPass())
  
  options(stringsAsFactors = FALSE)
  
  anno_rows <- data.frame(cohort_id = character(0),
                          cohort_name = character(0),
                          cohort_creation_datetime = character(0),
                          individual_id = character(0),
                          ready_for_analysis = character(0),
                          tumour_sample = character(0),
                          germline_sample = character(0),
                          row.names = NULL)
  
  anno_all <- anno_rows
  
  for(i in 1:nrow(cohorts)){
    print(paste0("cohort ", i, " of ", nrow(cohorts)))
    cohort <- cohorts[i,]
    cohort_id <- cohort$id
    cohort_name <- cohort$name
    cohort_creation_datetime <- cohort$creationDate
    cohort_annotation_sets <- cohort$annotationSets
    for(j in 1:length(cohort_annotation_sets)){
      annotations <- cohort_annotation_sets[[j]]$annotations
      if(!is.null(annotations)){
        individual_id <- unlist(annotations[[1]]["value"][(annotations[[1]]["name"]=="individualId")])
        ready_for_analysis <- unlist(annotations[[1]]["value"][(annotations[[1]]["name"]=="readyForAnalysis")])
        matchedSamples <- annotations[[1]]["value"][(annotations[[1]]["name"]=="matchedSamples")]
        #go through each row of the matched samples to capture FF+GL and FFPE+GL, etc.
        for(k in 1:nrow(matchedSamples[[1]]["tumourSampleId"])){
          tumour_sample <- matchedSamples[[1]]["tumourSampleId"][k,]
          germline_sample <- matchedSamples[[1]]["germlineSampleId"][k,]
          pair_name <- paste0(tumour_sample, "_", germline_sample)
        }
      }else{
        individual_id <- "NA"
        ready_for_analysis <- "NA"
        pair_name <- "NA"
      }
      anno_rows <- data.frame(cohort_id = cohort_id,
                              cohort_name = cohort_name,
                              cohort_creation_datetime = cohort_creation_datetime,
                              individual_id = individual_id,
                              ready_for_analysis = ready_for_analysis,
                              tumour_sample = tumour_sample,
                              germline_sample = germline_sample,
                              row.names = NULL)
    }
    anno_all <- rbind(anno_all, anno_rows)
  }
  
  if(latest_cohorts_only == TRUE){
    anno_all_unique <- anno_all %>% 
      group_by(individual_id) %>%
      filter(cohort_id == max(cohort_id)) %>%
      ungroup() %>% 
      unique()
  }else{
    anno_all_unique <- anno_all
  }
  
  return(anno_all_unique)
}