end_time_infected <- function(table){
  I_list <- c()
  for (i in 1:nrow(table)){
    # if in i list and recover remove
    if(table$infector[i] %in% I_list && table$affected[i] == 0){
      I_list <- I_list[!I_list %in% table$infector[i]]
    }
    else{
      I_list <- c(I_list, table$affected[i])
    }
  }
  return(I_list)
}

# take in table from full sim and n which is number of infected for your subtable
sublist <- function(table, I_list, n){
  sub_Ilist <- sample(I_list, n)
  sub_Ilist <- unlist(sub_Ilist)
  
  sampled_indices <- c()
  
  for (i in sub_Ilist){
    sampled_indices <- c(sampled_indices, which(table$affected== i))
  }
  # Retrieve the heights for the sampled indices
  sampled_heights <- table$heights[sampled_indices]
  # print(sampled_heights)
  
  # Create a data frame with the sampled indices and their corresponding heights
  sampled_df <- data.frame(indices = sub_Ilist, heights = sampled_heights)
  
  # Sort the data frame by heights in descending order
  sampled_df <- sampled_df[order(-sampled_df$heights), ]
  
  # Extract the sorted indices
  sorted_sampled_indices <- sampled_df$indices
  
  return(sorted_sampled_indices)
}


# Function to get ancestors of a given id
MRCA <- function(table,id1, id2) {
  get_ancestors <- function(id) {
    ancestors <- c(id)
    index <- which(table$affected == id)
    
    while (length(index) > 0) {
      infector <- table$infector[[index]]
      ancestors <- c(ancestors, infector)
      index <- which(table$affected == infector)
    }
    return(ancestors)
  }
  
  # Get ancestors for both ids
  ancestors1 <- get_ancestors(id1)
  ancestors2 <- get_ancestors(id2)
  
  intersection <- intersect(ancestors1, ancestors2)
  
  ancestors1 <- setdiff(ancestors1, intersection)
  ancestors2 <- setdiff(ancestors2, intersection)
  
  A1 <- ancestors1[length(ancestors1)]
  A2 <- ancestors2[length(ancestors2)]
  
  ind1 <- 0
  ind2 <- 0
  
  mrca <- 0
  
  
  if(length(A1) == 0 && length(A2) == 0){
    mrca <- 0
  }
  
  else if(length(A1) > 0 && length(A2) == 0){
    mrca <- A1
  }
  
  else if(length(A1) == 0 && length(A2) > 0){
    mrca <- A2
  }
  
  else if(length(A1) > 0 && length(A2) > 0){
    ind1 <- which(table$affected == ancestors1[length(ancestors1)])
    ind2 <- which(table$affected == ancestors2[length(ancestors2)])
    
    if (!(length(ind1) == 0 && length(ind2) == 0)){
      if (ind1 < ind2){
        mrca <- A1
      } else{
        mrca <- A2
      }
    } 
  }
  
  if (length(mrca) == 0){
    mrca <- 0
  }
  
  # Return the intersection of both ancestors lists
  return(mrca)
}


# get matrix displaying time of infection via mrca from the ids. 
MRCA_MATRIX <- function(table, sample){
  n <- length(sample)
  # print(sample)
  mrca_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(sample, sample))
  
  mrca <- NA
  for (i in seq_along(sample)){
    for (j in seq_along(sample)){
      mrca <- MRCA(table, sample[i], sample[j])
      mrca_matrix[i, j] <- mrca
    }
  }
  # dont care about lower triangle set to 0s
  lower_tri_index <- lower.tri(mrca_matrix)
  mrca_matrix[lower_tri_index] <- 0

  # convert to times
  for (i in 1:nrow(mrca_matrix)){
    for (j in 1:nrow(mrca_matrix)){
      if (mrca_matrix[i, j] == 0){
        mrca_matrix[i, j] <- 0
      } 
      else{
        mrca_time <- table$time[which(table$affected == mrca_matrix[i, j])]
        mrca_matrix[i, j] <- mrca_time 
      }
    }
  }
  return(mrca_matrix)
}

subtable <- function(table, mrca_matrix, sample){
  p_0 <- as.numeric(table$affected[1])

  if (sample[1] != p_0){
    sampled_p_0 <- sample[1]
    sample[1] <- p_0
  }
  else {
    sampled_p_0 <- FALSE
  }
  
  drawn <- c()
  
  get_drawn_ancestors <- function(id) {
    ancestors <- c(id)
    index <- which(as.numeric(table$affected) == as.numeric(id))
    
    while (length(index) > 0 && !(table$infector[[index]] %in% drawn || table$infector[[index]] == p_0)) {
      infector <- table$infector[[index]]
      ancestors <- c(ancestors, infector)
      index <- which(table$affected == infector)
    }
    
    drawn <<- c(drawn, ancestors[length(ancestors)])
    return(ancestors)
  }
  
  for (i in 1:length(sample)){
    get_drawn_ancestors(sample[i])
  }
  
  if (drawn[1] != p_0){
    drawn[1] <- p_0
  }
  
  subtable <- table
  subtable <- subtable[subtable$affected %in% drawn | subtable$affected == 0, ]
  
  # correct ids
  if (sampled_p_0){
    sample[1] <- sampled_p_0
  }
  sample_w_0 <- c(0, sample)
  drawn_w_0 <- c(0, drawn)
  
  id_replacement <- setNames(sample_w_0, drawn_w_0)
  
  subtable$affected <- id_replacement[as.character(subtable$affected)]
  subtable$infector <- id_replacement[as.character(subtable$infector)]
  
  # correct times with matrix
  for (i in 2:nrow(subtable)){
    row <- subtable$infector[[i]]
    col <- subtable$affected[[i]]
    
    row_index <- which(sample == row)
    col_index <- which(sample == col)

    subtable$time[i] <- mrca_matrix[row_index, col_index]
  }
  
  # ensure who infects is correct
  for (i in 3:nrow(subtable)){
    time <- subtable$time[i]
    inf <- subtable$infector[i]
    while (inf != sample[1] && time < subtable$time[which(subtable$affected == inf)]){
      inf <- subtable$infector[which(subtable$affected == inf)]
    }
    subtable$infector[i] <- inf
  }
  subtable <- subtable[order(subtable$time), ]
  
  # correct counts
  for (i in 1:nrow(subtable)){
    subtable$S[i] <- nrow(subtable) - i
    subtable$I[i] <- i
    subtable$R[i] <- 0
  }
  return (subtable)
}


subtree <- function(table, n){
  # get list of those infected at end time
  I_list <- end_time_infected(table)
  if (length(I_list) < n){
    stop("Error: not enough infected individuals remain alive at sample time, 
         try a different simulation or a smaller sample size.")
  }
  
  # get table as only infection events not including p_0 infection
  filtered_table <- table[2:nrow(table), ]
  filtered_table <- filtered_table[filtered_table$affected != 0, ]
  small_table <- rbind(table[1, ], filtered_table)
  # get sublist of sampled infected at sample time
  sublist <- sublist(small_table, I_list, n)
  table$affected[1] <- 0
  # get MRCA information and store in matrix
  mat <- MRCA_MATRIX(small_table, sublist)
  
  # get subtable for plotting
  subtable <- subtable(small_table, mat, sublist)
  
  return(subtable)
}
