entropySupCal <- function(maplabels = maplabels, u = u){

    entS <- 0
    for(i in  1:ncol(u)){
        tab <- table(maplabels, u[,i])
        #Leave at zero if one of clusters has 0 observations
        if(min(rowSums(tab)) > 0) {
            tab <- tab*log(tab/rowSums(tab))
            tab <- sum(ifelse(tab == "NaN", 0, tab))
        } else {
            tab <- 0
        }
        entS <- entS + tab
        }

   return(entS)
}

