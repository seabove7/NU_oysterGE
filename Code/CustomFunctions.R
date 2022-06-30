

##################################################################
###################### plotPCA modifications ######################
##################################################################

### Updated 'plotPCA' (from DESeq2) to save all PCs for plasticity measures

plotPCA_allPCs <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object)) # Variance estimates for each row (column) in a matrix.
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))] # orders genes from largest -> smallest and pulls ntop (default 500) genes
  pca <- prcomp(t(assay(object)[select, ]), center = TRUE, scale. = FALSE) # performs PCA on ntop (default 500) from dataset
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(group = group, intgroup.df, name = colnames(object), pca$x)
  if (returnData) {
    attr(d, "percentVar") <- percentVar # add [1:2] if only want 1st two PC variances
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}


################################################################################
################### Sample plasticity based on PC distances ####################
################################################################################

## Code pulled from github: https://github.com/seabove7/RandomFun/tree/main/Plasticity_function

PCAplast <- function(pca, data, sample_ID = NA, num_pca = "all", control_col, control_lvl = "none", group = NA) {
  
  # rename the user input info
  pca_df <- pca
  data_df <- data
  control_name <- control_col
  control_lvl <- control_lvl
  group_col <- group
  
  # ifelse statement to pull PCs from dataframe or prcomp objects
  if(class(pca_df) == "prcomp"){
    pca_dist <- pca_df$x # grab PC distances from prcomp() object
  } else {
    pca_dist <- pca_df # grab PC distances from data.frame object
  }
  
  
  # check for correct number of PCAs provided
  if(class(num_pca) == "numeric") {
    if(num_pca < 2) { # will throw error if too few PCAs requested
      stop("please select more than 2 PCs to calculate distance")
    } 
  }
  
  if(class(num_pca) == "numeric") {
    if(num_pca > (pca_dist %>% dplyr::select(starts_with("PCA")) %>% ncol())) { # will throw error if too many PCAs requested
      stop(paste(num_pca, "PCs requested for calculation, but only", (pca_dist %>% dplyr::select(starts_with("PCA")) %>% ncol()), "PCs available. Select appropriate number of PCAs for calculation."))
    } 
  }
  
  
  # oder the dataframe to ensure correct pairing after calculating distance 
  if(sample_ID %in% colnames(data_df)){
    data_df <- data_df[order(data_df[[sample_ID]]),]
  } else {
    data_df <- data_df[order(row.names(data_df)),]
  }
  
  
  # combine the datasets
  dist_df <- cbind(data_df, pca_dist) 
  
  
  # if there is no control level, modify so function pulls all levels
  if(control_lvl == "none") {
    control_lvl = levels(dist_df[[control_name]])
  } else {
    control_lvl = control_lvl
  }
  
  
  # make dataframe of control grouping only
  if(!is.na(group_col)) { # calculate mean per grouping ID (if provided)
    mean_control <- dist_df %>%
      filter(dist_df[[control_name]] == list(control_lvl)[[1]]) %>% 
      rename_with(tolower) %>% # renames all pc's with lowercase 'PC' (just to differentiate from all sample PCs)
      dplyr::select(colnames((dist_df %>% rename_with(tolower))[tolower(group_col)]), starts_with("pc")) %>% 
      dplyr::select(-contains("pco")) # use this to remove 'pco2' column. may need to add custom ones of these
    
    # add the control PCA values to treatment samples per grouping
    dist_df2 <- left_join(dist_df, mean_control, by.x = group_col, by.y = tolower(group_col))
    
  } else { # calculate mean per control treatment
    mean_control <- dist_df %>%
      filter(dist_df[[control_name]] == list(control_lvl)[[1]]) %>% 
      rename_with(tolower) %>% # renames all pc's with lowercase 'PC' (just to differentiate from all sample PCs)
      dplyr::select(colnames((dist_df %>% rename_with(tolower))[tolower(control_name)]), starts_with("pc")) %>% # select just the PCs 
      #group_by() 
      dplyr::select(-contains("pco")) %>% # use this to remove 'pco2' column. may need to add custom ones of these
      group_by_at(vars(tolower(control_name))) %>% 
      summarise_if(is.numeric, mean) %>% 
      dplyr::select(-tolower(control_name))
    
    # add the control PCA values to all samples 
    dist_df2 <- merge(dist_df, mean_control, all = TRUE)
    
  }
  
  
  
  # again, reorder data
  if(sample_ID %in% colnames(data_df)){
    dist_df2 <- dist_df2[order(dist_df2[[sample_ID]]),]
  } else {
    rownames(dist_df2) <- rownames(data_df)
    dist_df2 <- dist_df2[order(row.names(dist_df2)),]
  }
  
  
  ### Calculate sample (PCA) distances from control (pca) using all PCAs
  # make dataframe to populate with pca distances
  full_calc_dist <- data.frame(control_name = dist_df2[control_name])
  
  if(num_pca == "all") {
    ## forloop that will calculate distances between control and sample for all PCs (n will be total number)
    for(n in 1:(dist_df %>% dplyr::select(starts_with("PC") & -contains("pco")) %>% ncol())){
      # makes the PCA column name for control (lowercase) and sample (uppercase)
      PC_col <- paste0("PC", n)
      pc_col <- paste0("pc", n)
      
      # pulls the PC column for control (lowercase) and sample (uppercase)
      PCx <- dist_df2[PC_col]
      pcx <- dist_df2[pc_col]
      
      pca_calc_dist <- data.frame((PCx - pcx)^2) # calculates the distance between 2 PCs
      full_calc_dist <- cbind(full_calc_dist, pca_calc_dist) # add that distance to running dataframe
    }
  } else {
    ## forloop that will calculate distances between control and sample for SPECIFIED # of PCs (n will be total number)
    for(n in 1:as.numeric(num_pca)){
      # makes the PC column name for control (lowercase) and sample (uppercase)
      PC_col <- paste0("PC", n)
      pc_col <- paste0("pc", n)
      
      # pulls the PC column for control (lowercase) and sample (uppercase)
      PCx <- dist_df2[PC_col]
      pcx <- dist_df2[pc_col]
      
      pca_calc_dist <- data.frame((PCx - pcx)^2) # calculates the distance between 2 PCs
      full_calc_dist <- cbind(full_calc_dist, pca_calc_dist) # add that distance to running dataframe
    }
  }
  
  
  ## final distance calculation (adds all PCA distances and takes squareroot)
  distance <- full_calc_dist %>% 
    mutate(dis_sum = rowSums(across(where(is.numeric)))) %>% 
    mutate(dist = sqrt(dis_sum)) %>% 
    dplyr::select(matches("dist"))
  
  
  ## combine the calculated distance with the metadata and remove controls for final dataframe
  dist_df <- data_df %>% 
    bind_cols(distance) %>% 
    filter(!is.na(dist)) 
  
  ## removes the control levels
  if(length(control_lvl) > 1){
    dist_df <- dist_df
  } else {
    dist_df <- dist_df %>%
      filter(dist_df[[control_name]] != control_lvl)  %>% 
      droplevels()
  }
  
}


