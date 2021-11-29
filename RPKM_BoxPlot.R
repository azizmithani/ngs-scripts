data <- as.matrix(read.table("RPKM_Per_Gene_Afg.txt", sep = "\t", header=T))
RPKM <- matrix(as.numeric(data[,-1]),nrow=nrow(data), byrow=FALSE)
colnames(RPKM) <- colnames(data)[-1]
rownames(RPKM) <- data[,1]

# normalised RPKM
nRPKM <- RPKM
nRPKM[,3:ncol(RPKM)] <- nRPKM[,3:ncol(RPKM)]/apply(RPKM[,3:ncol(RPKM)],1,sum)

# median
mRPKM <- apply(RPKM[,3:ncol(RPKM)],1,median)


#################################################################

# absolute deviation
adRPKM <- RPKM
adRPKM[,3:ncol(RPKM)] <-abs(mRPKM - adRPKM[,3:ncol(RPKM)])

# median absolute deviation
madRPKM <- apply(RPKM[,3:ncol(adRPKM)],1,median)
# sort RPKM decreasingly according to median values
RPKM <- RPKM[sort(madRPKM, decreasing = T, index.return = T)$ix,]
madRPKM <- sort(madRPKM, decreasing = T)

#################################################################



#################################################################
# remove values with zero median
zero <- which (mRPKM == 0)
RPKM <- RPKM[-zero,]
mRPKM <- mRPKM[-zero]

# sort RPKM according to median values
RPKM <- RPKM[sort(mRPKM, decreasing = F, index.return = T)$ix,]
mRPKM <- sort(mRPKM, decreasing = F)

#################################################################


# sort RPKM decreasingly according to median values
RPKM <- RPKM[sort(mRPKM, decreasing = T, index.return = T)$ix,]
mRPKM <- sort(mRPKM, decreasing = T)

# Unigenes that are variable (at least n times than the median or at most 1/n times the median)
n <- 2.0
RPKM.Var <- RPKM[apply(RPKM[,3:ncol(RPKM)] >= mRPKM*n | RPKM[,3:ncol(RPKM)] <= mRPKM/n, 1, sum) > 0,]

# Plotting
RPKM.Unigene <- c()
RPKM.RPKM <- c()
RPKM.Unigene.Var <- c()
RPKM.RPKM.Var <- c()
RPKM.Range<- 1:100
for (i in 3:ncol(RPKM)) {
	RPKM.Unigene <- c(RPKM.Unigene, paste(format(RPKM.Range), rownames(RPKM)[RPKM.Range]))
	RPKM.Unigene.Var <- c(RPKM.Unigene.Var, paste(format(RPKM.Range), rownames(RPKM.Var)[RPKM.Range]))
	RPKM.RPKM <- c(RPKM.RPKM, RPKM[RPKM.Range,i])
	RPKM.RPKM.Var <- c(RPKM.RPKM.Var, RPKM.Var[RPKM.Range,i])
}
RPKM1 <- data.frame(Unigene=RPKM.Unigene, RPKM=log2(RPKM.RPKM))
RPKM1.Var <- data.frame(Unigene=RPKM.Unigene.Var, RPKM=log2(RPKM.RPKM.Var))
pdf(file="RPKM.pdf", width=10, height=7)
par(mar=c(7,4,3,1))
# All Unigenes
boxplot(log2(t(RPKM[RPKM.Range,3:8])), las=2, cex.axis=0.6, main=paste("Top", RPKM.Range[length(RPKM.Range)], "unigenes"), ylab="Log2(RPKM)")
for (i in 1:6) {
	stripchart(RPKM ~ Unigene, data=RPKM1[(i-1)*length(RPKM.Range) + RPKM.Range,], vertical = TRUE, las =2, cex.axis=0.6, cex = 0.6, 	add=TRUE, col="RED")
}
# Variable unigenes only
boxplot(log2(t(RPKM.Var[RPKM.Range,3:8])), las=2, cex.axis=0.6, main=paste("Top", RPKM.Range[length(RPKM.Range)], "variable unigenes"), ylab="Log2(RPKM)")
stripchart(RPKM ~ Unigene, data=RPKM1.Var, vertical = TRUE, las =2, cex.axis=0.6, cex = 0.6, add=TRUE, col="RED")
dev.off()


# test - each line with different colour
pdf(file="RPKM1.pdf", width=10, height=7)
par(mar=c(7,4,3,1))
# All Unigenes
boxplot(log2(t(RPKM[RPKM.Range,3:8])), las=2, cex.axis=0.6, main=paste("Top", RPKM.Range[length(RPKM.Range)], "unigenes"), ylab="Log2(RPKM)")
for (i in 1:6) {
	stripchart(RPKM ~ Unigene, data=RPKM1[(i-1)*length(RPKM.Range) + RPKM.Range,], vertical = TRUE, las =2, cex.axis=0.6, cex = 0.6, add=TRUE, col=i)
}
# Variable unigenes only
boxplot(log2(t(RPKM.Var[RPKM.Range,3:8])), las=2, cex.axis=0.6, main=paste("Top", RPKM.Range[length(RPKM.Range)], "variable unigenes"), ylab="Log2(RPKM)")
stripchart(RPKM ~ Unigene, data=RPKM1.Var, vertical = TRUE, las =2, cex.axis=0.6, cex = 0.6, add=TRUE, col="RED")
dev.off()
