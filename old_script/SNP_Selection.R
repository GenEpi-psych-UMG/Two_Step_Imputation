#######################
#
# R scripted algorithm for variants filtration lists creation 
#
# (c) M.Kamal Nasr 2024 mohammed.nasr@uni-greifswald.de
#
#######################

####Get SNP_Lists to be filtered out from the main Panel and extract it from the add on panel.
if(!require("data.table")) install.packages("data.table")
if(!require("dplyr")) install.packages("dplyr")

add_quotes_to_column <- function(df, column_name) {
  # Check if the specified column exists
  if (!column_name %in% colnames(df)) {
    stop(paste("Column", column_name, "not found in the data frame"))
  }
  
  # Add quotes to each element in the specified column
  df[[column_name]] <- paste0('"', df[[column_name]], '"')
  
  return(df)
}


arguments = commandArgs(trailingOnly = TRUE)

if (length(arguments) == 0) {
  stop("No input pathway provided.")
}

input_pathway <- arguments[1]
print(paste("Loading ", input_pathway, sep = ""))
all_data <- fread(input_pathway)

if (exists("all_data") && is.data.frame(all_data)) {
  print("The dataframe 'all_data' is loaded.")
} else {
  print("The dataframe 'all_data' is not loaded.")
}

for (i in unique(all_data$chr)) {
  
  input <- all_data[all_data$chr == i,]

  
  if (nrow(input) == 1) {
    #stop("Not enough panels")
    Filter_out_main_path <- input$Path[1]
      Filter_out_main_info <- fread(Filter_out_main_path)
      
      ##Select unique variants 
      Filter_out_main_info$POS <- sub(".*:(\\d+):.*", "\\1", Filter_out_main_info$SNPID)
      Filter_out_main_info$chr <- sub("(\\d+):.*", "\\1", Filter_out_main_info$SNPID)
      Filter_out_main_info <- Filter_out_main_info[as.numeric(Filter_out_main_info$R2) < 0.9,] ##Threshold for imputation quality
      Filter_out_main <- subset(Filter_out_main_info, select = c("chr", "POS"))
      write.table(Filter_out_main[,1:2], paste("Selected_filterout_main_chr", Filter_out_main$chr[1], ".txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

  }else {
    ##Read files
    if (nrow(input) == 2) {
      Main_panel_path <- input$Path[1]
      Main_panel_info <- fread(Main_panel_path)
      Addon_panel_path <- input$Path[2]
      Addon_panel_info <- fread(Addon_panel_path)
      
      ##Select unique variants 
      Addon_panel_info$POS <- sub(".*:(\\d+):.*", "\\1", Addon_panel_info$SNPID)
      Addon_panel_info$chr <- sub("(\\d+):.*", "\\1", Addon_panel_info$SNPID)
      #Addon_panel_info$Reference <- arguments[2,2]
      #Main_panel_info$Reference <- arguments[1,2]
      Addon_unique <- Addon_panel_info[!(Addon_panel_info$SNPID %in% Main_panel_info$SNPID),]
      Addon_unique <- Addon_unique[as.numeric(Addon_unique$R2) >= 0.9,] ##Threshold for imputation quality
      Addon_unique <- subset(Addon_unique, select = c("chr", "POS"))
      ###Define overlapping variants with higher R2 in the add-on panel
      Panel_merge <- merge(Main_panel_info, Addon_panel_info, by = c("SNPID", "ID"))
      
      ##Select filter-out SNPs in main SNP-list
      Filter_out_main <- Panel_merge[Panel_merge$R2.y > Panel_merge$R2.x]
      Filter_out_main <- Filter_out_main[as.numeric(Filter_out_main$R2.y) >= 0.9,] ##Threshold for imputation quality
      Filter_out_main$POS <- sub(".*:(\\d+):.*", "\\1", Filter_out_main$SNPID)
      Filter_out_main$chr <- sub("(\\d+):.*", "\\1", Filter_out_main$SNPID)
      Filter_out_main <- subset(Filter_out_main, select = c("chr", "POS"))
      Filter_out_main$reference <- input[1,2]
      write.table(Filter_out_main, paste("Selected_filterout_main_chr", Filter_out_main$chr[1], "_with_reference.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
      write.table(Filter_out_main[,1:2], paste("Selected_filterout_main_chr", Filter_out_main$chr[1], ".txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
      
      ##Select keep SNPs in the add-on SNP-list
      
      addon_keep <- data.frame(
        chr = Filter_out_main$chr,
        POS = Filter_out_main$POS)
      
      Unique_keep <- data.frame(
        chr = Addon_unique$chr,
        POS = Addon_unique$POS)
      
      
      
      Keep_list <- rbind(Unique_keep, addon_keep)
      Keep_list$reference <- input$Reference[2]
      write.table(Keep_list, paste("Selected_keep_list_main_chr", Keep_list$chr[1], "_with_reference.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
      write.table(Keep_list[,1:2], paste("Selected_keep_list_main_chr", Keep_list$chr[1],"_", input$Reference[2],".txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
      
    }else {
      ##Create list of all add-on panels
      #input <- add_quotes_to_column(input, "Path")
      dfs <- lapply(input$Path[2:nrow(input)], function(file_path) {
        # Debugging: print the file path being read
        print(paste("Reading file:", file_path))
        
        # Read the file using fread
        fread(file_path)
      })
      cohort_names <- input$Reference[2:nrow(input)]
      names(dfs) <- cohort_names
      dfs <- lapply(names(dfs), function(name) {
        df <- dfs[[name]]
        df$Reference <- name ##column contain the name of the reference panel
        return(df)
      })
      #names(dfs) <- tools::file_path_sans_ext(basename(input[2:nrow(input),2]))
      
      Combined_dfs <- do.call(rbind, dfs)
      
      
      ##Select the unique SNPIDs and make a datframe of its info
      duplicates <- Combined_dfs[duplicated(Combined_dfs$SNPID) | duplicated(Combined_dfs$SNPID, fromLast = TRUE), ]
      uniques <- Combined_dfs[!duplicated(Combined_dfs$SNPID) & !duplicated(Combined_dfs$SNPID, fromLast = TRUE), ]
      
      ######test for uniqueness and duplication 
      if (nrow(duplicates) + nrow(uniques) != nrow(Combined_dfs)) {
        stop("duplicates and uniques are not generated properly")
      }
      High_duplicates <- duplicates %>%
        group_by(SNPID) %>%
        arrange(desc(R2)) %>%
        slice(1)
      
      Main_panel_path <- input$Path[1]
      Main_panel_info <- fread(Main_panel_path)
      Addon_panel_info <- rbind(uniques, High_duplicates)
      Main_panel_info$Reference_main <- input[1,2]
      
      Addon_panel_info$POS <- sub(".*:(\\d+):.*", "\\1", Addon_panel_info$SNPID)
      Addon_panel_info$chr <- sub("(\\d+):.*", "\\1", Addon_panel_info$SNPID)
      Addon_unique <- Addon_panel_info[!(Addon_panel_info$SNPID %in% Main_panel_info$SNPID),]
      Addon_unique <- Addon_unique[as.numeric(Addon_unique$R2) >= 0.9,] ##Threshold for imputation quality
      Addon_unique <- subset(Addon_unique, select = c("chr", "POS", "Reference"))
      ###Define overlapping variants with higher R2 in the add-on panel
      Panel_merge <- merge(Main_panel_info, Addon_panel_info, by = c("SNPID", "ID"))
      
      ##Select filter-out SNPs in main SNP-list
      Filter_out_main <- Panel_merge[Panel_merge$R2.y > Panel_merge$R2.x]
      Filter_out_main <- Filter_out_main[as.numeric(Filter_out_main$R2.y) >= 0.9,] ##Threshold for imputation quality
      Filter_out_main <- subset(Filter_out_main, select = c("chr", "POS", "Reference","Reference_main"))
      
      write.table(Filter_out_main, paste("Selected_filterout_main_chr", Filter_out_main$chr[1], "_with_reference.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
      write.table(Filter_out_main[,1:2], paste("Selected_filterout_main_chr", Filter_out_main$chr[1], ".txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
      
      
      ##Select keep SNPs in the add-on SNP-list
      addon_keep <- data.frame(
        chr = Filter_out_main$chr,
        POS = Filter_out_main$POS,
        Reference = Filter_out_main$Reference)
      
      Unique_keep <- data.frame(
        chr = Addon_unique$chr,
        POS = Addon_unique$POS,
        Reference = Addon_unique$Reference)
      
      Keep_list <- rbind(Unique_keep, addon_keep)
      write.table(Keep_list, paste("Selected_keep_list_main_chr", Keep_list$chr[1], "_with_reference.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
      ###Write SNP lists for each add-on panel
      unique_references <- unique(Keep_list$Reference)
      
      for (ref in unique_references) {
        subset_df <- subset(Keep_list, Reference == ref)
        write.table(subset_df[,1:2], paste("Selected_keep_list_main_chr", Keep_list$chr[1],"_", ref ,".txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
        
      }
      
      
      
    }
    print(paste("SNPs lists for chr", Filter_out_main$chr[1], "are generated", sep = " "))
    
    rm(list = setdiff(ls(), c("arguments", "all_data")))
  }
}



