split.taxon <- function(x,
                        num_classes = 6,
                        class_method = c("quantile","equal","natural","stddev")[1],
                        rij_taxon,
                        name_feat = "sp_gen_dim",
                        return_class_midpoint = F) {
  
    # Limit x values to PUs where the species is present
  x.sp <- x[which(rij_taxon==1)]
  
  # Calculate breaks between classes
  interval <- ( max(x.sp) - min(x.sp) ) / num_classes # used only for equal interval
  breaks <- switch(class_method,
      "quantile" =  c(min(x.sp), quantile_breaks(num_classes,data.frame(x.sp=x.sp)), max(x.sp)),
      "equal" = c(min(x.sp), min(x.sp)+interval*c(1:(num_classes-1)), max(x.sp)),
      "natural" =  c(min(x.sp), natural_breaks(num_classes,data.frame(x.sp=x.sp)), max(x.sp)),
      "stddev" = c(min(x.sp), stddev_breaks(data.frame(x.sp=x.sp)), max(x.sp))
  )
  
  # If breaks are not all unique:
  if (length(breaks) != length(unique(breaks))) {
      breaks <- unique(breaks)
      warning(paste("There are non-unique breaks for",class_method,num_classes,"classes. Using",length(breaks)-1,"classes instead"))
      # stop(paste("There are non-unique breaks for",class_method,num_classes,"classes. Using",length(breaks)-1,"classes instead"))
      num_classes <- length(breaks) - 1
      }
  
  # Classify x values into num_classes
  x_fac <- cut(x, breaks, include.lowest=T)
  
  # Define x_classes matrix containing the distribution of the taxon split by columns
  x_classes <- matrix(NA,nrow=length(x),ncol=num_classes)
  # Loop on columns: fill the x_classes with taxon presence/absence if classified into that factor level
  for (i in 1 : num_classes) {
    x_classes[,i] <- ifelse(x_fac==levels(x_fac)[i],rij_taxon,0)
  }
  # PUs that were not classified into any factor level (bcs outside of species range) are assigned 0
  x_classes[which(is.na(x_classes),arr.ind=T)] <- 0
  # Assign names to the classes
  colnames(x_classes) <- paste0(name_feat,"_",c(1:num_classes))
  
  # return
  if (return_class_midpoint) list(midpoints = Midx(breaks), st_matrix = x_classes) else x_classes
}


split.taxon.multi <- function(x,
                        num_classes = 6,
                        rij_taxon,
                        name_feat = "sp_gen_dim",
                        return_class_midpoint = F) {
    
    nvar <- ncol(x)
    x$name <- rownames(x)
    
    # Limit x values to PUs where the species is present
    x.sp <- x[which(rij_taxon==1),]
    
    a <- kmeans(x.sp[,1:nvar], num_classes)
    x.sp$cluster <- a$cluster
    aa <- data.frame(name=x.sp$name,cluster=factor(a$cluster))
    
    x %>% left_join(aa,by="name") %>% pull(cluster) -> x_fac
    
    # Define x_classes matrix containing the distribution of the taxon split by columns
    x_classes <- matrix(NA,nrow=nrow(x),ncol=num_classes)
    # Loop on columns: fill the x_classes with taxon presence/absence if classified into that factor level
    for (i in 1 : num_classes) {
        x_classes[,i] <- ifelse(x_fac==levels(x_fac)[i],rij_taxon,0)
    }
    # PUs that were not classified into any factor level (bcs outside of species range) are assigned 0
    x_classes[which(is.na(x_classes),arr.ind=T)] <- 0
    # Assign names to the classes
    colnames(x_classes) <- paste0(name_feat,"_",c(1:num_classes))
    
    # return
    if (return_class_midpoint) list(midpoints = a$centers, st_matrix = x_classes) else x_classes
}
