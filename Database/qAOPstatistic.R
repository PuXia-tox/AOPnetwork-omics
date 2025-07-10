qAOPstatistic <- function(cmke) {
  
  message("Parsing AOP network...")
  steatosis_aop <- convert_cytoscape_to_aop("/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/aopwiki-KER_2020_09_21.cyjs")
  nodes <- convert_aop_to_graph(steatosis_aop) %>%
    as("matrix") %>%
    colnames() %>%
    sapply(function(x) getAOPNodeName(steatosis_aop, x)) %>%
    gsub("\\(.*\\)", "", .) %>%
    gsub(",", "", .) %>%
    toupper() %>%
    cbind(matrix(.), names(.)) %>%
    .[, c(1, 3)]
  
  # Extract POD stats and MIE paths
  pod_stats <- sapply(cmke, function(x) x[[1]])
  mie_paths <- sapply(cmke, function(x) x[[2]])
  mie_nodes <- sapply(mie_paths, function(x) x[1])
  
  # Split POD strings into components
  pod_split <- strsplit(pod_stats, " ")
  pod_min    <- sapply(pod_split, `[`, 1)
  pod_count  <- sapply(pod_split, `[`, 2)
  pod_mean   <- sapply(pod_split, `[`, 3)
  pod_median <- sapply(pod_split, `[`, 4)
  
  # Combine KE and POD summaries
  summary_matrix <- cbind(mie_nodes, pod_min, pod_mean, pod_median)
  mie_summary <- data.frame(
    KE = names(table(mie_nodes)),
    min.aop = tapply(as.numeric(pod_min), mie_nodes, min),
    mean.aop = tapply(as.numeric(pod_mean), mie_nodes, min),
    median.aop = tapply(as.numeric(pod_median), mie_nodes, min),
    stringsAsFactors = FALSE
  )
  
  # Merge KE frequency and annotation info
  mie_freq <- data.frame(KE = names(table(mie_nodes)), num = as.integer(table(mie_nodes)), stringsAsFactors = FALSE)
  mie_annot <- nodes[nodes[, 2] %in% mie_freq$KE, , drop = FALSE]
  colnames(mie_annot) <- c("name", "KE")
  
  stat_table <- merge(mie_freq, mie_annot, by = "KE")
  stat_table <- merge(stat_table, mie_summary, by = "KE")
  stat_table <- stat_table[order(as.numeric(stat_table$median.aop)), ]
  
  # Generate all KE event frequency table
  all_KEs <- table(unlist(sapply(cmke, `[[`, 2)))
  all_KE_percent <- round(sort(all_KEs, decreasing = TRUE) / length(cmke) * 100, 2)
  
  # Map KE names
  KE_names <- sapply(names(all_KE_percent), function(k) {
    match_row <- which(nodes[, 2] == k)
    if (length(match_row)) nodes[match_row[1], 1] else NA
  })
  
  all_aop <- data.frame(
    nodes = names(all_KE_percent),
    name = KE_names,
    percent_cmke = as.numeric(all_KE_percent),
    matched = "n",
    mie = "n",
    ao = "n",
    stringsAsFactors = FALSE
  )
  
  # Flag matched nodes
  matched_nodes <- unique(unlist(sapply(cmke, function(x) x[[3]])))
  all_aop$matched[all_aop$nodes %in% matched_nodes] <- "yes"
  
  # Flag MIEs and AOs
  mie_set <- unique(sapply(mie_paths, function(x) x[1]))
  ao_set  <- unique(sapply(mie_paths, function(x) x[length(x)]))
  
  all_aop$mie[all_aop$nodes %in% mie_set] <- "yes"
  all_aop$ao[all_aop$nodes %in% ao_set]  <- "yes"
  
  # Estimate minimum POD per KE
  pod_vals <- sapply(cmke, function(z) strsplit(z[[1]], " ")[[1]][3])
  pod_numeric <- as.numeric(pod_vals)
  
  KE_pod <- sapply(all_aop$nodes, function(k) {
    matched_indices <- which(sapply(cmke, function(z) k %in% z[[3]]))
    if (length(matched_indices)) {
      min(pod_numeric[matched_indices], na.rm = TRUE)
    } else {
      NA
    }
  })
  
  all_aop$pod <- KE_pod
  
  return(all_aop)
}