
## TODO: export and document these functions (but leave them as internal)

sampledf <- function(df, size, replace=FALSE, prob=NULL, group, status){
  
  if (!missing(group) && !missing(status))
    elig.ids <- df$ids[df$group %in% group & df$status %in% status]
  
  if (missing(group) && !missing(status))
    elig.ids <- df$ids[df$status %in% status]
  
  if (!missing(group) && missing(status))
    elig.ids <- df$ids[df$group %in% group]
  
  if (missing(group) && missing(status))
    elig.ids <- df$ids
  
  if (length(elig.ids) > 1) {
    ids <- sample(elig.ids, size, replace, prob)
  } else {
    if (size > 0) {
      ids <- elig.ids
    } else {
      ids <- NULL
    }
  }
  
  return(ids)
}


get.prevICM <- function(df, odf, type, nsteps, at, set.df=TRUE) {
  
  groups <- length(unique(df$group))
  status <- df$status
  
  if (set.df == TRUE) {
    if (groups == 1) {
      if (at == 1) {
        odf <- data.frame(s.num.g1 = sum(status == 0, na.rm=TRUE))
        odf$i.num.g1 <- sum(status == 1, na.rm=TRUE)
        odf$num.g1 <- odf$s.num.g1+odf$i.num.g1
        if (type == 'SIR') {
          odf$r.num.g1 <- sum(status == 2, na.rm=TRUE)
          odf$num.g1 <- odf$num.g1 + odf$r.num.g1
        }
      } else {
        odf$s.num.g1[at] <- sum(status == 0, na.rm=TRUE)
        odf$i.num.g1[at] <- sum(status == 1, na.rm=TRUE)
        odf$num.g1[at] <- odf$s.num.g1[at] + odf$i.num.g1[at]
        if (type == 'SIR') {
          odf$r.num.g1[at] <- sum(status == 2, na.rm=TRUE)
          odf$num.g1[at] <- odf$num.g1[at] + odf$r.num.g1[at]
        }
      }
    } else {
      ids.g1 <- df$ids[df$group == 1]
      ids.g2 <- df$ids[df$group == 2]
      if (at == 1) {
        odf <- data.frame(s.num.g1 = sum(status[ids.g1] == 0, na.rm=TRUE))
        odf$i.num.g1 <- sum(status[ids.g1] == 1, na.rm=TRUE)
        odf$s.num.g2 <- sum(status[ids.g2] == 0, na.rm=TRUE)
        odf$i.num.g2 <- sum(status[ids.g2] == 1, na.rm=TRUE)
        odf$s.num <- odf$s.num.g1 + odf$s.num.g2
        odf$i.num <- odf$i.num.g1 + odf$i.num.g2
        odf$num.g1 <- odf$s.num.g1 + odf$i.num.g1
        odf$num.g2 <- odf$s.num.g2 + odf$i.num.g2
        if (type == 'SIR') {
          odf$r.num.g1 <- sum(status[ids.g1] == 2, na.rm=TRUE)
          odf$r.num.g2 <- sum(status[ids.g2] == 2, na.rm=TRUE)
          odf$r.num <- odf$r.num.g1 + odf$r.num.g2
          odf$num.g1 <- odf$num.g1 + odf$r.num.g1
          odf$num.g2 <- odf$num.g2 + odf$r.num.g2
        }
      } else {
        odf$s.num.g1[at] <- sum(status[ids.g1] == 0, na.rm=TRUE)
        odf$i.num.g1[at] <- sum(status[ids.g1] == 1, na.rm=TRUE)
        odf$s.num.g2[at] <- sum(status[ids.g2] == 0, na.rm=TRUE)
        odf$i.num.g2[at] <- sum(status[ids.g2] == 1, na.rm=TRUE)
        odf$s.num[at] <- odf$s.num.g1[at] + odf$s.num.g2[at]
        odf$i.num[at] <- odf$i.num.g1[at] + odf$i.num.g2[at]
        odf$num.g1[at] <- odf$s.num.g1[at] + odf$i.num.g1[at]
        odf$num.g2[at] <- odf$s.num.g2[at] + odf$i.num.g2[at]
        if (type == 'SIR') {
          odf$r.num.g1[at] <- sum(status[ids.g1] == 2, na.rm=TRUE)
          odf$r.num.g2[at] <- sum(status[ids.g2] == 2, na.rm=TRUE)
          odf$r.num[at] <- odf$r.num.g1[at] + odf$r.num.g2[at]
          odf$num.g1[at] <- odf$num.g1[at] + odf$r.num.g1[at]
          odf$num.g2[at] <- odf$num.g2[at] + odf$r.num.g2[at]
        }
      }
    }
    if (at == 1) {
      empt.odf <- data.frame(matrix(rep(NA, ncol(odf)*(nsteps-1)), nrow=nsteps-1))
      names(empt.odf) <- names(odf)
      odf <- rbind(odf, empt.odf)
    }
    return(odf)
  }
  
}

ssample <- function(x, size, replace = FALSE, prob = NULL) {
  
  if (length(x) > 1) 
    return(sample(x, size, replace, prob))
  
  if (length(x) == 1 && size > 0)
    return(x)
  
  if (length(x) == 1 && size == 0)
    return(NULL)
}