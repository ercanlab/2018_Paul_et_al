#############
#Loess span validation
#############


#plot on all trnas
#significant intra and inter
#
#condensin peaks
#TfIIC bidning peaks
#cen as well
if (!require("fANCOVA")) {
  install.packages("fANCOVA")
  library(fANCOVA)
}
biocLite('fANCOVA')
poliiivscondensin<-cbind(smc4_chip_trna_for_circos,poliii_trna_score_circos[match(smc4_chip_trna_for_circos[,2],poliii_trna_score_circos[,2]),3])
ind<-which(is.na(poliiivscondensin[,4])=='TRUE')


pdf('~/Box Sync/Lab/Data_Analysis/tRNA/polIII_at_trna_correlation_with_smc4.pdf')
par(mar=c( 5 , 5 , 2 , 2 ) )
plot(poliiivscondensin[,4],poliiivscondensin[,3], ylab='SMC4 ChIP', xlab='PolIII ChIP', pch='.')
abline(lm( as.numeric(poliiivscondensin[,3]) ~ as.numeric(poliiivscondensin[,4])),col='red')
text(50,1.5,lab='R2 = 0.0001436')
cor(cbind(as.numeric(poliiivscondensin[-ind,3]),as.numeric(poliiivscondensin[-ind,4])), method='pearson')
ml<-lm( as.numeric(poliiivscondensin[,4]) ~ as.numeric(poliiivscondensin[,3]))
dev.off()
summary(ml)



#Lets optimise loess model

df<-data.frame(cbind(as.numeric(poliiivscondensin[,4]),as.numeric(poliiivscondensin[,3])),stringsAsFactors = F)
df<-df[-which(is.na(df[,1])=='TRUE'),]

#Set up useful variables to handle cross-validation loop.
span.seq <- seq(from = 0.15, to = 0.95, by = 0.05) #explores range of spans
k <- 10 #number of folds
set.seed(1) # replicate results
folds <- sample(x=(1:k), size = nrow(df), replace = TRUE)
cv.error.mtrx <- matrix(rep(x = NA, times = k * length(span.seq)), 
                        nrow = length(span.seq), ncol = k)


#Run a nested for loop iterating over each span possibility in span.seq, and each fold in folds
for(i in 1:length(span.seq)) {
  for(j in 1:k) {
    loess.fit <- loess(df[folds != j, 2] ~ df[folds != j, 1], span = span.seq[i])
    preds <- predict(loess.fit, df[folds == j,1 ])
    cv.error.mtrx[i, j] <- mean((df[folds == j,1] - preds)^2, na.rm = TRUE)
    # some predictions result in `NA` because of the `x` ranges in each fold
  }
}


#Calculate average cross-validation mean square error from each of the 10 folds:
cv.errors <- rowMeans(cv.error.mtrx,na.rm = T)  

#Find which span resulted in the lowest MSE.
best.span.i <- which.min(cv.errors)
best.span.i<-15
span.seq[best.span.i]

#Plot your results.
plot(x = span.seq, y = cv.errors, type = "l", main = "CV Plot")
points(x = span.seq, y = cv.errors, pch = 20, cex = 0.75, col = "blue")
points(x = span.seq[best.span.i], y = cv.errors[best.span.i], pch = 20, cex = 1, col = "red")

best.loess.fit <- loess(df[,2]~df[,1], span = span.seq[best.span.i])

x.seq <- seq(from = min(df[,1]), to = max(df[,1]), length = 100)

plot(x = df[,1], y = df[,2], main = "Best Span Plot")
lines(x = x.seq, y = predict(best.loess.fit,x.seq), col = "red", lwd = 2)
