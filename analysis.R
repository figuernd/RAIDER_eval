sortSeedBySen <- function(F, c) {
    F2 = subset(F, chr==c & tool == "raider")
    return(F2[order(F2$tpr), "seed"])
}