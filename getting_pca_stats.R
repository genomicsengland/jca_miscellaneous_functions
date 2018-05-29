is_null_check <- function(value){
  if(!is.null(value)){
    result <- value
  }else{
    result <- "NA"
  }
  return(result)
}

getting_pca_stats <- function(con, list, study, germline_study){
  library(opencgaR)
  library(getPass)
  
  # con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")
  # con <- opencgaLogin(opencga = con, userid = "jambrose", passwd = getPass())
  
  pca_stats_row <- data.frame(sample_name = character(0),
                              germline_sample = character(0),
                              ldp_code = character(0),
                              lab_sample_id = character(0),
                              library_type = character(0),
                              gc_drop = character(0),
                              at_drop = character(0),
                              coverage_homogeneity = character(0),
                              chimeric_percentage = character(0),
                              average_fragment_size = character(0),
                              mapping_rate = character(0),
                              log10_snvs = character(0),
                              log10_indels = character(0),
                              disease = character(0),
                              disease_subtype = character(0),
                              tumour_delivery_id = character(0),
                              germline_delivery_id = character(0),
                              bertha_version = character(0),
                              row.names = NULL)
  
  pca_stats <- pca_stats_row
  
  pca_stats <- do.call(rbind, lapply(list, function(name){
    file <- try(fileClient(OpencgaR = con, action = "search", params = list(study = study, name = paste0(name,".bam"), uri="~file:///genomes/by_date/", include = "name,uri,stats,attributes")), silent = TRUE)
    if(is.data.frame(file) && nrow(file) > 0){
      print(file$name)
      library_type <- is_null_check(file$stats$SAMPLE_METADATA$sample_library_type)
      gc_drop <- is_null_check(file$stats$AT_GC_DROP$gc_drop)
      at_drop <- is_null_check(file$stats$AT_GC_DROP$at_drop)
      coverage_homogeneity <- is_null_check(round(file$stats$WHOLE_GENOME_COVERAGE$coverageSummary[[1]][file$stats$WHOLE_GENOME_COVERAGE$coverageSummary[[1]]$scope=="autosomes",]$localRMSD, 2))
      chimeric_percentage <- is_null_check(file$stats$SAMTOOLS_STATS_CALCULATIONS_ALL$SAMTOOLS_PAIRS_ON_DIFFERENT_CHROMOSOMES_PERCENT)
      average_fragment_size <- is_null_check(file$stats$SAMTOOLS_STATS_ALL$SAMTOOLS_INSERT_SIZE_AVERAGE)
      mapping_rate <- is_null_check(file$stats$SAMTOOLS_STATS_CALCULATIONS_ALL$SAMTOOLS_READS_MAPPED_PERCENT)
      log10_snvs <- is_null_check(round(log(file$stats$ILLUMINA_SUMMARY_REPORT_CANCER$pair_stats$SOMATIC_SNVS, 10), digits = 2))
      log10_indels <- is_null_check(round(log(file$stats$ILLUMINA_SUMMARY_REPORT_CANCER$pair_stats$SOMATIC_INDELS, 10), digits = 2))
      snvs_per_mb <- is_null_check(round(file$stats$ILLUMINA_SUMMARY_REPORT_CANCER$pair_stats$SOMATIC_SNVS/(3077073773/1000000), digits = 2))

      #this is not always there: germline_sample <- file$stats$SAMPLE_METADATA$sample_matched_sample
      if(!is.null(file$stats$ILLUMINA_SUMMARY_REPORT_CANCER$pair_stats$NORMAL_SAMPLE_NAME)){
        germline_sample <- file$stats$ILLUMINA_SUMMARY_REPORT_CANCER$pair_stats$NORMAL_SAMPLE_NAME
      }else if(!is.null(file$stats$SAMPLE_METADATA$sample_matched_sample)){
        germline_sample <- file$stats$SAMPLE_METADATA$sample_matched_sample
      }else{
        germline_sample <- NULL
      }
      if(!is.null(germline_sample)){
        file_germline <- try(fileClient(OpencgaR = con, action = "search", params = list(study = germline_study, name = paste0(germline_sample,".bam"), uri="~file:///genomes/by_date/", include = "name,uri,stats,attributes")), silent = TRUE)
        germline_delivery_id <- unlist(file_germline$attributes$deliveryId)
        if(length(germline_delivery_id) == 1){
          #check there is only one deliveryId that is not the germline, i.e. the tumour id
          if(length(which(!(unlist(file$attributes$deliveryId) %in% unlist(file_germline$attributes$deliveryId)))) == 1){
            tumour_delivery_id <- unlist(file$attributes$deliveryId)[which(!(unlist(file$attributes$deliveryId) %in% unlist(file_germline$attributes$deliveryId)))]
          }else{
            stop(paste0("Error: more than one non-germline deliveryId for tumour sample ", file$name))
          }
        }else{
          stop(paste0("Error: more than one germline deliveryId for tumour sample ", file$name))
        }
      }else{
        print(paste0("Error finding germline_sample for ", file$name))
        germline_sample <- "NA"
        tumour_delivery_id <- "NA"
        germline_delivery_id <- "NA"
      }

      bertha_version <- file$attributes$processes %>%
        as.data.frame() %>%
        arrange(desc(processDate)) %>%
        select(version) %>%
        slice(1) %>%
        rename(bertha_version = version)

      sampleAnnotationSet <- sampleAnnotationsetClient(con, sample = name, action = "search", params = list(study = study))
      if(is.data.frame(sampleAnnotationSet) && nrow(sampleAnnotationSet) > 0){
        lab_sample_id <- unlist(sampleAnnotationSet$annotations[[1]]["value"][(sampleAnnotationSet$annotations[[1]]["name"]=="labSampleId")])
        disease <- unlist(sampleAnnotationSet$annotations[[1]]["value"][(sampleAnnotationSet$annotations[[1]]["name"]=="diseaseType")])
        disease_subtype <- unlist(sampleAnnotationSet$annotations[[1]]["value"][(sampleAnnotationSet$annotations[[1]]["name"]=="diseaseSubType")])
        ldp_code <- unlist(sampleAnnotationSet$annotations[[1]]["value"][(sampleAnnotationSet$annotations[[1]]["name"]=="LDPCode")])
      }else{
        lab_sample_id <- "no_annotation_set_in_opencga"
        disease <- "no_annotation_set_in_opencga"
        disease_subtype <- "no_annotation_set_in_opencga"
        ldp_code <- "no_annotation_set_in_opencga"
      }
      
      pca_stats_row <- data.frame(sample_name = name,
                                  germline_sample = germline_sample,
                                  ldp_code = ldp_code,
                                  lab_sample_id = lab_sample_id,
                                  library_type = library_type,
                                  gc_drop = gc_drop,
                                  at_drop = at_drop,
                                  coverage_homogeneity = coverage_homogeneity,
                                  chimeric_percentage = chimeric_percentage,
                                  average_fragment_size = average_fragment_size,
                                  mapping_rate = mapping_rate,
                                  log10_snvs = log10_snvs,
                                  log10_indels = log10_indels,
                                  disease = disease,
                                  disease_subtype = disease_subtype,
                                  tumour_delivery_id = tumour_delivery_id,
                                  germline_delivery_id = germline_delivery_id,
                                  bertha_version = bertha_version,
                                  row.names = NULL)
    }else{
      message <- paste0("not_in_opencga")
      pca_stats_row <- data.frame(sample_name = name,
                                  germline_sample = message,
                                  ldp_code = message,
                                  sample_id = message,
                                  library_type = message,
                                  gc_drop = message,
                                  at_drop = message,
                                  coverage_homogeneity = message,
                                  chimeric_percentage = message,
                                  average_fragment_size = message,
                                  mapping_rate = message,
                                  log10_snvs = message,
                                  log10_indels = message,
                                  disease = disease,
                                  disease_subtype = disease_subtype,
                                  tumour_delivery_id = message,
                                  germline_delivery_id = message,
                                  bertha_version = message,
                                  row.names = NULL)
    }
    return(pca_stats_row)
  }))
  return(pca_stats)
}