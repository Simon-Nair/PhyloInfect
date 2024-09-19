add_height <- function(table){
  table$heights <- NA
  table$heights <- rep(NA, nrow(table))
  table$heights[1] <- 1 
  
  prepped_table <- table
  return(prepped_table)
}

genPhyloInfect <- function(table, end_time=10, recent=table$affected[1], time_of_infection = 0, counter = 0){

  # get indices of events involving recent id
  indices <- rev(which(as.numeric(unlist(table$infector)) == as.numeric(unlist(recent)))) 
  
  segment_count <- 0  # Initialize segment count for this id
  
  # if recovery
  if (length(indices) > 0 && table$affected[indices[1]] == 0){ 
    recovery_time <- table$time[indices[1]] # get recovery time
    # draw line
    segments(time_of_infection, height, recovery_time, height, lwd=1)
    segments(time_of_infection, height, time_of_infection,  height + 0.016 * counter, lwd=1,)
    text(recovery_time + gap, height + 0.003, labels = recent, cex = 0.55)
    indices <- indices[-1] #remove recovery
    segment_count <- segment_count + 1
  }
  # if patient does not recover
  else {
    # draw lines till end time since no recovery
    segments(time_of_infection, height, end_time, height, lwd=1)
    segments(time_of_infection, height, time_of_infection,  height + 0.016 * counter, lwd=1,)
    text(end_time + gap, height + 0.003, labels = recent, cex = 0.55)
    segment_count <- segment_count + 1
  }
  
  for (i in indices){
    if (gap == 0.3){
      gap <<-  0.1
    }else{
      gap <<-  0.3
    }
    height <<- height - 0.016 # decrement height for next line to be drawn
    prepped_table$heights[i] <<- height
    segment_count <- segment_count + genPhyloInfect(table, end_time, table$affected[i], table$time[i], segment_count)
  }
  return(segment_count)
}


phyloInfect <- function(table, title="", end_time=10){
  height <<- 1
  gap <<- 0.3
  ticks <- seq(0, end_time, by = 0.01)
  label_ticks <- seq(0, end_time, by = 1)
  labels <- as.character(label_ticks)
  plot(1, type = "n", xlim = c(0, end_time + 0.1), ylim = c(0, 1), main=title, xlab = "Time", ylab = "", axes = FALSE)
  axis(1, at = label_ticks, labels = label_ticks)
  
  
  prepped_table <<- add_height(table)
  
  genPhyloInfect(prepped_table, end_time)
}