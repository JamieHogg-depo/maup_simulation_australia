##
# Combine results
##

.libPaths("r_lib")
library(tidyverse)
library(sf)

#' @param dir character vector specifying where the specific dated folders can be found
#' @param cur_date character vector in the format yearmonthdate
#' @param list_element character vector specifying what elements of `which_list` should be used in the bind_rows() function
#' @param which_list character vector specifying which (of the two) lists to work with
#' @param subset_name character vector specifying the name for the resulting .rds name

combine_results <- function(dir,
                            cur_date,
                            list_element, 
                            which_list,
                            subset_name, 
                            drop_geo = FALSE){
  
  r_loc <- paste0(dir, "/", cur_date, "/r")
  file_names <- list.files(r_loc)
  cur_files <- file_names[which(str_detect(file_names, which_list))]

  for(i in 1:length(cur_files)){
  assign(cur_files[i],
       readRDS(paste0(r_loc, "/", cur_files[i])))
  }
  what_to_eval <- paste0("list(", 
  paste(
  paste0(cur_files, "$", list_element), 
  collapse = ", "),
  ")" )
  if(drop_geo == FALSE){
    saveRDS(
      bind_rows(eval(parse(text=what_to_eval))), 
      file=paste0(dir, "/", cur_date, "/", which_list, "_", subset_name, ".rds")
      ) 
  } else {
    saveRDS(
      bind_rows(eval(parse(text=what_to_eval))) %>% 
        st_drop_geometry(), 
      file=paste0(dir, "/", cur_date, "/", which_list, "_", subset_name, ".rds")
    ) 
  }

  message("File ", 
          paste0(dir, "/", cur_date, "/", which_list, "_", subset_name, ".rds"), 
          " saved!")
}

combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "data", 
                which_list = "storeddatalist",
                subset_name = "sa1data")
				
combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "MI", 
                which_list = "storeddatalist",
                subset_name = "sa1MI")
				
combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "sa1model", 
                which_list = "storeddatalist",
                subset_name = "sa1models")

combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "sa2model", 
                which_list = "storeddatalist",
                subset_name = "sa2models")
				
combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "sa3model", 
                which_list = "storeddatalist",
                subset_name = "sa3models")
				
combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "sa4model", 
                which_list = "storeddatalist",
                subset_name = "sa4models")
				
combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "MI", 
                which_list = "simlist",
                subset_name = "MI")
				
combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "OD", 
                which_list = "simlist",
                subset_name = "OD")
				
combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "model", 
                which_list = "simlist",
                subset_name = "models")
				
combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "df_agg", 
                which_list = "simlist",
                subset_name = "df_agg")

combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "df_agg", 
                which_list = "simlist",
                subset_name = "df_agg_nogeo",
                drop_geo = TRUE)
				
combine_results(dir = "CDS_MAUP2/outputs",
                cur_date = cur_date,
                list_element = "data_zone_metrics", 
                which_list = "simlist",
                subset_name = "data_zone_metrics")


