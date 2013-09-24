entropyMarCal <- function(maplabels = maplabels, u = u){
    entM <- 0
    for(i in  1:ncol(u)){
        tab <- table(maplabels, u[,i])
        #Leave at zero if one of clusters has 0 observations
        if(colSums(tab)[2] > 0) {
            divid <- colSums(tab)
            pen <- tab[,2]*log(tab[,2]/divid[2])
            pen <- sum(ifelse(pen == "NaN", 0, pen))
        } else {
            pen <- 0
        }
        entM <- entM + pen
        }
   return(entM)
}

