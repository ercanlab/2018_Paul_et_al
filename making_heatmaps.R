######################################
# Heatmaps from interaction matrices #
######################################

#Date
#18.04.26


#Load in required packages 
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(RColorBrewer)
}
if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(reshape2)
}
if (!require("DescTools")) {
  install.packages("DescTools", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(DescTools)
}
if (!require("viridis")) {
  install.packages("viridis", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(viridis)
}

#Prevent scientific notation
options(scipen=999)

#Read in the fends tables - Do this for MRP23-MRP27
MRP23<-data.matrix(read.table("~/Box Sync/Lab/Data/Hi-C/ts_rep1/MRP23.gz", stringsAsFactors = F, sep=" "))
MRP24<-data.matrix(read.table("~/Box Sync/Lab/Data/Hi-C/ts_rep1/MRP24.gz", stringsAsFactors = F, sep=" "))
MRP25<-data.matrix(read.table("~/Box Sync/Lab/Data/Hi-C/ts_rep1/MRP25.gz", stringsAsFactors = F, sep=" "))
MRP26<-data.matrix(read.table("~/Box Sync/Lab/Data/Hi-C/ts_rep1/MRP26.gz", stringsAsFactors = F, sep=" "))
print('all read in')

#Define bins and chromosome sizes
chr.size<-c(203893,794508,342718,1490682,602514,284456,1067526,544538,435585,719294,687260,1008248,908607,812465, 1054033,921188)

res<-10000

chr_bins<-list("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13"
               , "chr14", "chr15", "chr16")

bin_per_chr<-floor(chr.size/res+1)

incrementalbins<-matrix(0,length(chr_bins),1)
for (i in 1:(length(chr_bins)-1)){
  incrementalbins[i+1,]<-incrementalbins[i,]+bin_per_chr[i]
}
incrementalbins<-rbind(incrementalbins,"1196")
incrementalbins<-as.numeric(incrementalbins)

print('bins are set')

#Normalise the matrices
(MRP23/sum(MRP23))*100->norm_MRP23
(MRP24/sum(MRP24))*100->norm_MRP24
(MRP25/sum(MRP25))*100->norm_MRP25
(MRP26/sum(MRP26))*100->norm_MRP26

print('matrices normalised')

#Break down into each chromosome
for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP23_chr",h), norm_MRP23[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP24_chr",h), norm_MRP24[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP25_chr",h), norm_MRP25[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP26_chr",h), norm_MRP26[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

print('per chr matrix')

#Give me the matrix for chr 14 through 16 
assign("norm_MRP23_chr14_16", norm_MRP23[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP24_chr14_16", norm_MRP24[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP25_chr14_16", norm_MRP25[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP26_chr14_16", norm_MRP26[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])

print('per 3 chr matrix')

#Set up parameters for heatmaps
size=1196^2

my_palette <- colorRampPalette(c("white", "yellow","orange", "red", "black" ))(n = 119)

# (optional) defines the color breaks manually for a "skewed" color transition
cols = c(seq(0,(0.000049*1428025)/size,length=30),
         seq((0.00005*1428025)/size,(0.00009*1428025)/size,length=30),
         seq((0.0001*1428025)/size,(0.00049*1428025)/size,length=20),
         seq((0.0005*1428025)/size,(0.0009*1428025)/size,length=10),
         seq((0.001*1428025)/size,(0.00149*1428025)/size,length=10),             
         seq((0.0015*1428025)/size,(0.002*1428025)/size,length=20))  

#This sets up the edges of the chr to draw as a black line
incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

#Draw out the heatmaps

png("~/Desktop/MRP23_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 

heatmap.2(norm_MRP23_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)

dev.off()
png("~/Desktop/MRP24_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP24_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

png("~/Desktop/MRP25_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP25_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins2), v=incrementalbins2)
)

dev.off()

png("~/Desktop/MRP26_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP26_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

# Need to set up edges for genome-wide matrix and chr 12 matrix
incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

assign("norm_MRP23_chr12", norm_MRP23[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP25_chr12", norm_MRP25[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])

incrementalbins5<-incrementalbins
incrementalbins5[1]<-1
incrementalbins5<-incrementalbins5-0.5
incrementalbins5[4]<-incrementalbins5[4]+1


incrementalbins6<-incrementalbins5
incrementalbins6[1]<-1
incrementalbins6[2]<-(incrementalbins5[17]-incrementalbins5[16])+0.5
incrementalbins6[3]<-(incrementalbins5[17]-incrementalbins5[15])+0.5
incrementalbins6[4]<-(incrementalbins5[17]-incrementalbins5[14])+0.5
incrementalbins6[5]<-(incrementalbins5[17]-incrementalbins5[13])+0.5
incrementalbins6[6]<-(incrementalbins5[17]-incrementalbins5[12])+0.5
incrementalbins6[7]<-(incrementalbins5[17]-incrementalbins5[11])+0.5
incrementalbins6[8]<-(incrementalbins5[17]-incrementalbins5[10])+0.5
incrementalbins6[9]<-(incrementalbins5[17]-incrementalbins5[9])+0.5
incrementalbins6[10]<-(incrementalbins5[17]-incrementalbins5[8])+0.5
incrementalbins6[11]<-(incrementalbins5[17]-incrementalbins5[7])+0.5
incrementalbins6[12]<-(incrementalbins5[17]-incrementalbins5[6])+0.5
incrementalbins6[13]<-(incrementalbins5[17]-incrementalbins5[5])+0.5
incrementalbins6[14]<-(incrementalbins5[17]-incrementalbins5[4])+0.5
incrementalbins6[15]<-(incrementalbins5[17]-incrementalbins5[3])+0.5
incrementalbins6[16]<-(incrementalbins5[17]-incrementalbins5[2])+0.5
incrementalbins6[17]<-(incrementalbins5[17]-incrementalbins5[1])+0.5

incrementalbins5[1]<-incrementalbins5[1]+0.3
incrementalbins6[1]<-incrementalbins6[1]+0.25
incrementalbins5[16]<-incrementalbins5[16]-0.2
incrementalbins6[16]<-incrementalbins6[16]-0.3


heatmap.2(norm_MRP23,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins6), v=incrementalbins5)
)


heatmap.2(norm_MRP25,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins6), v=incrementalbins5)
)

#Just ChrXII

heatmap.2(norm_MRP23_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
)

heatmap.2(norm_MRP25_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
)

#Give me the matrix for chr 14 through 16 

assign("norm_MRP23_chr12", norm_MRP23[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP24_chr12", norm_MRP24[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP25_chr12", norm_MRP25[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP26_chr12", norm_MRP26[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])

print('per 3 chr matrix')

mean_MRP23<-mean(as.matrix(norm_MRP23))
mean_MRP24<-mean(as.matrix(norm_MRP24))
mean_MRP25<-mean(as.matrix(norm_MRP25))
mean_MRP26<-mean(as.matrix(norm_MRP26))

sd_MRP23<-sd(as.matrix(norm_MRP23))
sd_MRP24<-sd(as.matrix(norm_MRP24))
sd_MRP25<-sd(as.matrix(norm_MRP25))
sd_MRP26<-sd(as.matrix(norm_MRP26))

zscore_MRP23<-(norm_MRP23-mean_MRP23)/sd_MRP23
zscore_MRP24<-(norm_MRP24-mean_MRP24)/sd_MRP24
zscore_MRP25<-(norm_MRP25-mean_MRP25)/sd_MRP25
zscore_MRP26<-(norm_MRP26-mean_MRP26)/sd_MRP26

(zscore_MRP25-zscore_MRP23)->zscore_subtract_MRP23_25
(zscore_MRP24-zscore_MRP23)->zscore_subtract_MRP23_24
(zscore_MRP25-zscore_MRP24)->zscore_subtract_MRP24_25

assign("zscore_subtract_MRP23_25_chr14_16", zscore_subtract_MRP23_25[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP23_24_chr14_16", zscore_subtract_MRP23_24[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP24_25_chr14_16", zscore_subtract_MRP23_24[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])

# (optional) defines the color breaks manually for a "skewed" color transition
range=(0.000015*4428025)/size

col_breaks = c(seq((-5),-4,length=40),
               seq(((-3.99)),(-(3)),length=40),
               seq(-(2.99),-(2),length=40),
               seq(-(1.99),-(1),length=40),
               seq(-(0.99),0,length=40),
               seq(0.01,1,length=40),
               seq(1.01,2,length=40),
               seq(2.01,3,length=40),
               seq(3.01,4,length=40),
               seq(4.01,5,length=40))

incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

# creates a own color palette from red to green
my_palette <- colorRampPalette(c(magma(15)[4],'darkorchid3','white','palegreen3','darkgreen'))(n = 399)

png("~/Desktop/zscore_subtract_MRP23_25_chr14_16.png",# create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP23_25_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()
png("~/Desktop/zscore_subtract_MRP24_25_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP24_25_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)

dev.off()
png("~/Desktop/zscore_subtract_MRP23_24_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 

heatmap.2(zscore_subtract_MRP23_24_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)


dev.off()




#Lets do it it with just ChrXII

assign("zscore_subtract_MRP23_25_chr12", zscore_subtract_MRP23_25[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("zscore_subtract_MRP23_24_chr12", zscore_subtract_MRP23_24[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("zscore_subtract_MRP24_25_chr12", zscore_subtract_MRP23_24[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])


# creates a own color palette from red to green
#my_palette <- colorRampPalette(c("darkblue","blue4", "blue","cadetblue1", "white","darkgoldenrod1", "firebrick1", "red4", "darkred"))(n = 399)
# (optional) defines the color breaks manually for a "skewed" color transition
range=(0.000015*4428025)/size

col_breaks = c(seq((-5),-4,length=40),
               seq(((-3.99)),(-(3)),length=40),
               seq(-(2.99),-(2),length=40),
               seq(-(1.99),-(1),length=40),
               seq(-(0.99),0,length=40),
               seq(0.01,1,length=40),
               seq(1.01,2,length=40),
               seq(2.01,3,length=40),
               seq(3.01,4,length=40),
               seq(4.01,5,length=40))

incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3


png("~/Desktop/zscore_subtract_MRP23_25_chr12.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP23_25_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=incrementalbins_12, v=incrementalbins_12)
)
dev.off()


png("~/Desktop/zscore_subtract_MRP23_24_chr12.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP23_24_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=incrementalbins_12, v=incrementalbins_12)
)
dev.off()

png("~/Desktop/zscore_subtract_MRP24_25_chr12.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP24_25_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=incrementalbins_12, v=incrementalbins_12)
)
dev.off()



#Start again


MRP7<-data.matrix(read.table("~/Box Sync/Lab/Data/3C-seq/Mirny_pipeline/MRP7/heatmap.gz", stringsAsFactors = F, sep=" "))
MRP8<-data.matrix(read.table("~/Box Sync/Lab/Data/3C-seq/Mirny_pipeline/MRP8/heatmap.gz", stringsAsFactors = F, sep=" "))
MRP9<-data.matrix(read.table("~/Box Sync/Lab/Data/3C-seq/Mirny_pipeline/MRP9/heatmap.gz", stringsAsFactors = F, sep=" "))
MRP10<-data.matrix(read.table("~/Box Sync/Lab/Data/3C-seq/Mirny_pipeline/MRP10/heatmap.gz", stringsAsFactors = F, sep=" "))
print('all read in')

#Define bins and chromosome sizes
chr.size<-c(203893,794508,342718,1490682,602514,284456,1067526,544538,435585,719294,687260,1008248,908607,812465, 1054033,921188)

res<-10000

chr_bins<-list("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13"
               , "chr14", "chr15", "chr16")

bin_per_chr<-floor(chr.size/res+1)

incrementalbins<-matrix(0,length(chr_bins),1)
for (i in 1:(length(chr_bins)-1)){
  incrementalbins[i+1,]<-incrementalbins[i,]+bin_per_chr[i]
}
incrementalbins<-rbind(incrementalbins,"1196")
incrementalbins<-as.numeric(incrementalbins)

print('bins are set')

#Normalise the matrices
(MRP7/sum(MRP7))*100->norm_MRP7
(MRP8/sum(MRP8))*100->norm_MRP8
(MRP9/sum(MRP9))*100->norm_MRP9
(MRP10/sum(MRP10))*100->norm_MRP10

print('matrices normalised')

#Break down into each chromosome
for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP7_chr",h), norm_MRP7[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP8_chr",h), norm_MRP8[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP9_chr",h), norm_MRP9[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP10_chr",h), norm_MRP10[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

print('per chr matrix')

#Give me the matrix for chr 14 through 16 
assign("norm_MRP7_chr14_16", norm_MRP7[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP8_chr14_16", norm_MRP8[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP9_chr14_16", norm_MRP9[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP10_chr14_16", norm_MRP10[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])

print('per 3 chr matrix')

#Set up parameters for heatmaps
size=1196^2

my_palette <- colorRampPalette(c("white", "yellow","orange", "red", "black" ))(n = 119)

# (optional) defines the color breaks manually for a "skewed" color transition
cols = c(seq(0,(0.000049*1428025)/size,length=30),
         seq((0.00005*1428025)/size,(0.00009*1428025)/size,length=30),
         seq((0.0001*1428025)/size,(0.00049*1428025)/size,length=20),
         seq((0.0005*1428025)/size,(0.0009*1428025)/size,length=10),
         seq((0.001*1428025)/size,(0.00149*1428025)/size,length=10),             
         seq((0.0015*1428025)/size,(0.002*1428025)/size,length=20))  

#This sets up the edges of the chr to draw as a black line
incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

#Draw out the heatmaps

png("~/Desktop/MRP7_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 

heatmap.2(norm_MRP7_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)

dev.off()
png("~/Desktop/MRP8_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP8_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

png("~/Desktop/MRP9_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP9_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins2), v=incrementalbins2)
)

dev.off()

png("~/Desktop/MRP10_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP10_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

# Need to set up edges for genome-wide matrix and chr 12 matrix
incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

assign("norm_MRP7_chr12", norm_MRP7[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP9_chr12", norm_MRP9[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])

incrementalbins5<-incrementalbins
incrementalbins5[1]<-1
incrementalbins5<-incrementalbins5-0.5
incrementalbins5[4]<-incrementalbins5[4]+1


incrementalbins6<-incrementalbins5
incrementalbins6[1]<-1
incrementalbins6[2]<-(incrementalbins5[17]-incrementalbins5[16])+0.5
incrementalbins6[3]<-(incrementalbins5[17]-incrementalbins5[15])+0.5
incrementalbins6[4]<-(incrementalbins5[17]-incrementalbins5[14])+0.5
incrementalbins6[5]<-(incrementalbins5[17]-incrementalbins5[13])+0.5
incrementalbins6[6]<-(incrementalbins5[17]-incrementalbins5[12])+0.5
incrementalbins6[7]<-(incrementalbins5[17]-incrementalbins5[11])+0.5
incrementalbins6[8]<-(incrementalbins5[17]-incrementalbins5[10])+0.5
incrementalbins6[9]<-(incrementalbins5[17]-incrementalbins5[9])+0.5
incrementalbins6[10]<-(incrementalbins5[17]-incrementalbins5[8])+0.5
incrementalbins6[11]<-(incrementalbins5[17]-incrementalbins5[7])+0.5
incrementalbins6[12]<-(incrementalbins5[17]-incrementalbins5[6])+0.5
incrementalbins6[13]<-(incrementalbins5[17]-incrementalbins5[5])+0.5
incrementalbins6[14]<-(incrementalbins5[17]-incrementalbins5[4])+0.5
incrementalbins6[15]<-(incrementalbins5[17]-incrementalbins5[3])+0.5
incrementalbins6[16]<-(incrementalbins5[17]-incrementalbins5[2])+0.5
incrementalbins6[17]<-(incrementalbins5[17]-incrementalbins5[1])+0.5

incrementalbins5[1]<-incrementalbins5[1]+0.3
incrementalbins6[1]<-incrementalbins6[1]+0.25
incrementalbins5[16]<-incrementalbins5[16]-0.2
incrementalbins6[16]<-incrementalbins6[16]-0.3


heatmap.2(norm_MRP7,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins6), v=incrementalbins5)
)


heatmap.2(norm_MRP9,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins6), v=incrementalbins5)
)

#Just ChrXII

heatmap.2(norm_MRP7_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
)

heatmap.2(norm_MRP10_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
)

#Give me the matrix for chr 14 through 16 

assign("norm_MRP7_chr12", norm_MRP7[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP8_chr12", norm_MRP8[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP9_chr12", norm_MRP9[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP10_chr12", norm_MRP10[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])

print('per 3 chr matrix')

mean_MRP7<-mean(as.matrix(norm_MRP7))
mean_MRP8<-mean(as.matrix(norm_MRP8))
mean_MRP9<-mean(as.matrix(norm_MRP9))
mean_MRP10<-mean(as.matrix(norm_MRP10))

sd_MRP7<-sd(as.matrix(norm_MRP7))
sd_MRP8<-sd(as.matrix(norm_MRP8))
sd_MRP9<-sd(as.matrix(norm_MRP9))
sd_MRP10<-sd(as.matrix(norm_MRP10))

zscore_MRP7<-(norm_MRP7-mean_MRP7)/sd_MRP7
zscore_MRP8<-(norm_MRP8-mean_MRP8)/sd_MRP8
zscore_MRP9<-(norm_MRP9-mean_MRP9)/sd_MRP9
zscore_MRP10<-(norm_MRP10-mean_MRP10)/sd_MRP10

(zscore_MRP9-zscore_MRP7)->zscore_subtract_MRP7_9
(zscore_MRP8-zscore_MRP7)->zscore_subtract_MRP7_8
(zscore_MRP9-zscore_MRP8)->zscore_subtract_MRP8_9
(zscore_MRP10-zscore_MRP7)->zscore_subtract_MRP10_7

assign("zscore_subtract_MRP7_9_chr14_16", zscore_subtract_MRP7_9[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP7_8_chr14_16", zscore_subtract_MRP7_8[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP8_9_chr14_16", zscore_subtract_MRP8_9[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP10_7_chr14_16", zscore_subtract_MRP10_7[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])


# (optional) defines the color breaks manually for a "skewed" color transition
range=(0.000015*4428025)/size

col_breaks = c(seq((-5),-4,length=40),
               seq(((-3.99)),(-(3)),length=40),
               seq(-(2.99),-(2),length=40),
               seq(-(1.99),-(1),length=40),
               seq(-(0.99),0,length=40),
               seq(0.01,1,length=40),
               seq(1.01,2,length=40),
               seq(2.01,3,length=40),
               seq(3.01,4,length=40),
               seq(4.01,5,length=40))

incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

# creates a own color palette from red to green
my_palette <- colorRampPalette(c(magma(15)[4],'darkorchid3','white','palegreen3','darkgreen'))(n = 399)

png("~/Desktop/zscore_subtract_MRP7_10_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP10_7_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

#Lets do it it with just ChrXII

assign("zscore_subtract_MRP10_7_chr12", zscore_subtract_MRP10_7[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])


# creates a own color palette from red to green
#my_palette <- colorRampPalette(c("darkblue","blue4", "blue","cadetblue1", "white","darkgoldenrod1", "firebrick1", "red4", "darkred"))(n = 399)
# (optional) defines the color breaks manually for a "skewed" color transition
range=(0.000015*4428025)/size

col_breaks = c(seq((-5),-4,length=40),
               seq(((-3.99)),(-(3)),length=40),
               seq(-(2.99),-(2),length=40),
               seq(-(1.99),-(1),length=40),
               seq(-(0.99),0,length=40),
               seq(0.01,1,length=40),
               seq(1.01,2,length=40),
               seq(2.01,3,length=40),
               seq(3.01,4,length=40),
               seq(4.01,5,length=40))

incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3


png("~/Desktop/zscore_subtract_MRP7_10_chr12.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP10_7_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          #add.expr = abline(h=incrementalbins_12, v=incrementalbins_12)
)
dev.off()

######################################
# Heatmaps from interaction matrices #
######################################

#Date
#18.04.26


#Load in required packages 
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(RColorBrewer)
}
if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(reshape2)
}
if (!require("DescTools")) {
  install.packages("DescTools", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(DescTools)
}
if (!require("viridis")) {
  install.packages("viridis", dependencies = TRUE, repos='http://cran.us.r-project.org')
  library(viridis)
}

#Prevent scientific notation
options(scipen=999)

#Read in the fends tables - Do this for MRP23-MRP27
MRP23<-data.matrix(read.table("~/Box Sync/Lab/Data/Hi-C/ts_rep1/MRP23.gz", stringsAsFactors = F, sep=" "))
MRP24<-data.matrix(read.table("~/Box Sync/Lab/Data/Hi-C/ts_rep1/MRP24.gz", stringsAsFactors = F, sep=" "))
MRP25<-data.matrix(read.table("~/Box Sync/Lab/Data/Hi-C/ts_rep1/MRP25.gz", stringsAsFactors = F, sep=" "))
MRP26<-data.matrix(read.table("~/Box Sync/Lab/Data/Hi-C/ts_rep1/MRP26.gz", stringsAsFactors = F, sep=" "))
print('all read in')

#Define bins and chromosome sizes
chr.size<-c(203893,794508,342718,1490682,602514,284456,1067526,544538,435585,719294,687260,1008248,908607,812465, 1054033,921188)

res<-10000

chr_bins<-list("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13"
               , "chr14", "chr15", "chr16")

bin_per_chr<-floor(chr.size/res+1)

incrementalbins<-matrix(0,length(chr_bins),1)
for (i in 1:(length(chr_bins)-1)){
  incrementalbins[i+1,]<-incrementalbins[i,]+bin_per_chr[i]
}
incrementalbins<-rbind(incrementalbins,"1196")
incrementalbins<-as.numeric(incrementalbins)

print('bins are set')

#Normalise the matrices
(MRP23/sum(MRP23))*100->norm_MRP23
(MRP24/sum(MRP24))*100->norm_MRP24
(MRP25/sum(MRP25))*100->norm_MRP25
(MRP26/sum(MRP26))*100->norm_MRP26

print('matrices normalised')

#Break down into each chromosome
for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP23_chr",h), norm_MRP23[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP24_chr",h), norm_MRP24[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP25_chr",h), norm_MRP25[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP26_chr",h), norm_MRP26[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

print('per chr matrix')

#Give me the matrix for chr 14 through 16 
assign("norm_MRP23_chr14_16", norm_MRP23[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP24_chr14_16", norm_MRP24[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP25_chr14_16", norm_MRP25[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP26_chr14_16", norm_MRP26[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])

print('per 3 chr matrix')

#Set up parameters for heatmaps
size=1196^2

my_palette <- colorRampPalette(c("white", "yellow","orange", "red", "black" ))(n = 119)

# (optional) defines the color breaks manually for a "skewed" color transition
cols = c(seq(0,(0.000049*1428025)/size,length=30),
         seq((0.00005*1428025)/size,(0.00009*1428025)/size,length=30),
         seq((0.0001*1428025)/size,(0.00049*1428025)/size,length=20),
         seq((0.0005*1428025)/size,(0.0009*1428025)/size,length=10),
         seq((0.001*1428025)/size,(0.00149*1428025)/size,length=10),             
         seq((0.0015*1428025)/size,(0.002*1428025)/size,length=20))  

#This sets up the edges of the chr to draw as a black line
incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

#Draw out the heatmaps

png("~/Desktop/MRP23_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 

heatmap.2(norm_MRP23_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)

dev.off()
png("~/Desktop/MRP24_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP24_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

png("~/Desktop/MRP25_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP25_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins2), v=incrementalbins2)
)

dev.off()

png("~/Desktop/MRP26_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP26_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

# Need to set up edges for genome-wide matrix and chr 12 matrix
incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

assign("norm_MRP23_chr12", norm_MRP23[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP25_chr12", norm_MRP25[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])

incrementalbins5<-incrementalbins
incrementalbins5[1]<-1
incrementalbins5<-incrementalbins5-0.5
incrementalbins5[4]<-incrementalbins5[4]+1


incrementalbins6<-incrementalbins5
incrementalbins6[1]<-1
incrementalbins6[2]<-(incrementalbins5[17]-incrementalbins5[16])+0.5
incrementalbins6[3]<-(incrementalbins5[17]-incrementalbins5[15])+0.5
incrementalbins6[4]<-(incrementalbins5[17]-incrementalbins5[14])+0.5
incrementalbins6[5]<-(incrementalbins5[17]-incrementalbins5[13])+0.5
incrementalbins6[6]<-(incrementalbins5[17]-incrementalbins5[12])+0.5
incrementalbins6[7]<-(incrementalbins5[17]-incrementalbins5[11])+0.5
incrementalbins6[8]<-(incrementalbins5[17]-incrementalbins5[10])+0.5
incrementalbins6[9]<-(incrementalbins5[17]-incrementalbins5[9])+0.5
incrementalbins6[10]<-(incrementalbins5[17]-incrementalbins5[8])+0.5
incrementalbins6[11]<-(incrementalbins5[17]-incrementalbins5[7])+0.5
incrementalbins6[12]<-(incrementalbins5[17]-incrementalbins5[6])+0.5
incrementalbins6[13]<-(incrementalbins5[17]-incrementalbins5[5])+0.5
incrementalbins6[14]<-(incrementalbins5[17]-incrementalbins5[4])+0.5
incrementalbins6[15]<-(incrementalbins5[17]-incrementalbins5[3])+0.5
incrementalbins6[16]<-(incrementalbins5[17]-incrementalbins5[2])+0.5
incrementalbins6[17]<-(incrementalbins5[17]-incrementalbins5[1])+0.5

incrementalbins5[1]<-incrementalbins5[1]+0.3
incrementalbins6[1]<-incrementalbins6[1]+0.25
incrementalbins5[16]<-incrementalbins5[16]-0.2
incrementalbins6[16]<-incrementalbins6[16]-0.3


heatmap.2(norm_MRP23,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins6), v=incrementalbins5)
)


heatmap.2(norm_MRP25,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins6), v=incrementalbins5)
)

#Just ChrXII

heatmap.2(norm_MRP23_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
)

heatmap.2(norm_MRP25_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
)

#Give me the matrix for chr 14 through 16 

assign("norm_MRP23_chr12", norm_MRP23[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP24_chr12", norm_MRP24[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP25_chr12", norm_MRP25[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP26_chr12", norm_MRP26[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])

print('per 3 chr matrix')

mean_MRP23<-mean(as.matrix(norm_MRP23))
mean_MRP24<-mean(as.matrix(norm_MRP24))
mean_MRP25<-mean(as.matrix(norm_MRP25))
mean_MRP26<-mean(as.matrix(norm_MRP26))

sd_MRP23<-sd(as.matrix(norm_MRP23))
sd_MRP24<-sd(as.matrix(norm_MRP24))
sd_MRP25<-sd(as.matrix(norm_MRP25))
sd_MRP26<-sd(as.matrix(norm_MRP26))

zscore_MRP23<-(norm_MRP23-mean_MRP23)/sd_MRP23
zscore_MRP24<-(norm_MRP24-mean_MRP24)/sd_MRP24
zscore_MRP25<-(norm_MRP25-mean_MRP25)/sd_MRP25
zscore_MRP26<-(norm_MRP26-mean_MRP26)/sd_MRP26

(zscore_MRP25-zscore_MRP23)->zscore_subtract_MRP23_25
(zscore_MRP24-zscore_MRP23)->zscore_subtract_MRP23_24
(zscore_MRP25-zscore_MRP24)->zscore_subtract_MRP24_25

assign("zscore_subtract_MRP23_25_chr14_16", zscore_subtract_MRP23_25[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP23_24_chr14_16", zscore_subtract_MRP23_24[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP24_25_chr14_16", zscore_subtract_MRP23_24[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])

# (optional) defines the color breaks manually for a "skewed" color transition
range=(0.000015*4428025)/size

col_breaks = c(seq((-5),-4,length=40),
               seq(((-3.99)),(-(3)),length=40),
               seq(-(2.99),-(2),length=40),
               seq(-(1.99),-(1),length=40),
               seq(-(0.99),0,length=40),
               seq(0.01,1,length=40),
               seq(1.01,2,length=40),
               seq(2.01,3,length=40),
               seq(3.01,4,length=40),
               seq(4.01,5,length=40))

incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

# creates a own color palette from red to green
my_palette <- colorRampPalette(c(magma(15)[4],'darkorchid3','white','palegreen3','darkgreen'))(n = 399)

png("~/Desktop/zscore_subtract_MRP23_25_chr14_16.png",# create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP23_25_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()
png("~/Desktop/zscore_subtract_MRP24_25_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP24_25_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)

dev.off()
png("~/Desktop/zscore_subtract_MRP23_24_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 

heatmap.2(zscore_subtract_MRP23_24_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)


dev.off()




#Lets do it it with just ChrXII

assign("zscore_subtract_MRP23_25_chr12", zscore_subtract_MRP23_25[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("zscore_subtract_MRP23_24_chr12", zscore_subtract_MRP23_24[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("zscore_subtract_MRP24_25_chr12", zscore_subtract_MRP23_24[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])


# creates a own color palette from red to green
#my_palette <- colorRampPalette(c("darkblue","blue4", "blue","cadetblue1", "white","darkgoldenrod1", "firebrick1", "red4", "darkred"))(n = 399)
# (optional) defines the color breaks manually for a "skewed" color transition
range=(0.000015*4428025)/size

col_breaks = c(seq((-5),-4,length=40),
               seq(((-3.99)),(-(3)),length=40),
               seq(-(2.99),-(2),length=40),
               seq(-(1.99),-(1),length=40),
               seq(-(0.99),0,length=40),
               seq(0.01,1,length=40),
               seq(1.01,2,length=40),
               seq(2.01,3,length=40),
               seq(3.01,4,length=40),
               seq(4.01,5,length=40))

incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3


png("~/Desktop/zscore_subtract_MRP23_25_chr12.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP23_25_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=incrementalbins_12, v=incrementalbins_12)
)
dev.off()


png("~/Desktop/zscore_subtract_MRP23_24_chr12.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP23_24_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=incrementalbins_12, v=incrementalbins_12)
)
dev.off()

png("~/Desktop/zscore_subtract_MRP24_25_chr12.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP24_25_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=incrementalbins_12, v=incrementalbins_12)
)
dev.off()



#Start again


MRP7<-data.matrix(read.table("~/Box Sync/Lab/Data/3C-seq/Mirny_pipeline/MRP7/heatmap.gz", stringsAsFactors = F, sep=" "))
MRP8<-data.matrix(read.table("~/Box Sync/Lab/Data/3C-seq/Mirny_pipeline/MRP8/heatmap.gz", stringsAsFactors = F, sep=" "))
MRP9<-data.matrix(read.table("~/Box Sync/Lab/Data/3C-seq/Mirny_pipeline/MRP9/heatmap.gz", stringsAsFactors = F, sep=" "))
MRP10<-data.matrix(read.table("~/Box Sync/Lab/Data/3C-seq/Mirny_pipeline/MRP10/heatmap.gz", stringsAsFactors = F, sep=" "))
print('all read in')

#Define bins and chromosome sizes
chr.size<-c(203893,794508,342718,1490682,602514,284456,1067526,544538,435585,719294,687260,1008248,908607,812465, 1054033,921188)

res<-10000

chr_bins<-list("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13"
               , "chr14", "chr15", "chr16")

bin_per_chr<-floor(chr.size/res+1)

incrementalbins<-matrix(0,length(chr_bins),1)
for (i in 1:(length(chr_bins)-1)){
  incrementalbins[i+1,]<-incrementalbins[i,]+bin_per_chr[i]
}
incrementalbins<-rbind(incrementalbins,"1196")
incrementalbins<-as.numeric(incrementalbins)

print('bins are set')

#Normalise the matrices
(MRP7/sum(MRP7))*100->norm_MRP7
(MRP8/sum(MRP8))*100->norm_MRP8
(MRP9/sum(MRP9))*100->norm_MRP9
(MRP10/sum(MRP10))*100->norm_MRP10

print('matrices normalised')

#Break down into each chromosome
for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP7_chr",h), norm_MRP7[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP8_chr",h), norm_MRP8[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP9_chr",h), norm_MRP9[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("norm_MRP10_chr",h), norm_MRP10[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
}

print('per chr matrix')

#Give me the matrix for chr 14 through 16 
assign("norm_MRP7_chr14_16", norm_MRP7[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP8_chr14_16", norm_MRP8[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP9_chr14_16", norm_MRP9[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("norm_MRP10_chr14_16", norm_MRP10[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])

print('per 3 chr matrix')

#Set up parameters for heatmaps
size=1196^2

my_palette <- colorRampPalette(c("white", "yellow","orange", "red", "black" ))(n = 119)

# (optional) defines the color breaks manually for a "skewed" color transition
cols = c(seq(0,(0.000049*1428025)/size,length=30),
         seq((0.00005*1428025)/size,(0.00009*1428025)/size,length=30),
         seq((0.0001*1428025)/size,(0.00049*1428025)/size,length=20),
         seq((0.0005*1428025)/size,(0.0009*1428025)/size,length=10),
         seq((0.001*1428025)/size,(0.00149*1428025)/size,length=10),             
         seq((0.0015*1428025)/size,(0.002*1428025)/size,length=20))  

#This sets up the edges of the chr to draw as a black line
incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

#Draw out the heatmaps

png("~/Desktop/MRP7_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 

heatmap.2(norm_MRP7_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)

dev.off()
png("~/Desktop/MRP8_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP8_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

png("~/Desktop/MRP9_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP9_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins2), v=incrementalbins2)
)

dev.off()

png("~/Desktop/MRP10_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)

heatmap.2(norm_MRP10_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

# Need to set up edges for genome-wide matrix and chr 12 matrix
incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

assign("norm_MRP7_chr12", norm_MRP7[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP9_chr12", norm_MRP9[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])

incrementalbins5<-incrementalbins
incrementalbins5[1]<-1
incrementalbins5<-incrementalbins5-0.5
incrementalbins5[4]<-incrementalbins5[4]+1


incrementalbins6<-incrementalbins5
incrementalbins6[1]<-1
incrementalbins6[2]<-(incrementalbins5[17]-incrementalbins5[16])+0.5
incrementalbins6[3]<-(incrementalbins5[17]-incrementalbins5[15])+0.5
incrementalbins6[4]<-(incrementalbins5[17]-incrementalbins5[14])+0.5
incrementalbins6[5]<-(incrementalbins5[17]-incrementalbins5[13])+0.5
incrementalbins6[6]<-(incrementalbins5[17]-incrementalbins5[12])+0.5
incrementalbins6[7]<-(incrementalbins5[17]-incrementalbins5[11])+0.5
incrementalbins6[8]<-(incrementalbins5[17]-incrementalbins5[10])+0.5
incrementalbins6[9]<-(incrementalbins5[17]-incrementalbins5[9])+0.5
incrementalbins6[10]<-(incrementalbins5[17]-incrementalbins5[8])+0.5
incrementalbins6[11]<-(incrementalbins5[17]-incrementalbins5[7])+0.5
incrementalbins6[12]<-(incrementalbins5[17]-incrementalbins5[6])+0.5
incrementalbins6[13]<-(incrementalbins5[17]-incrementalbins5[5])+0.5
incrementalbins6[14]<-(incrementalbins5[17]-incrementalbins5[4])+0.5
incrementalbins6[15]<-(incrementalbins5[17]-incrementalbins5[3])+0.5
incrementalbins6[16]<-(incrementalbins5[17]-incrementalbins5[2])+0.5
incrementalbins6[17]<-(incrementalbins5[17]-incrementalbins5[1])+0.5

incrementalbins5[1]<-incrementalbins5[1]+0.3
incrementalbins6[1]<-incrementalbins6[1]+0.25
incrementalbins5[16]<-incrementalbins5[16]-0.2
incrementalbins6[16]<-incrementalbins6[16]-0.3


heatmap.2(norm_MRP7,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins6), v=incrementalbins5)
)


heatmap.2(norm_MRP9,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins6), v=incrementalbins5)
)

#Just ChrXII

heatmap.2(norm_MRP7_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
)

heatmap.2(norm_MRP10_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=cols,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
)

#Give me the matrix for chr 14 through 16 

assign("norm_MRP7_chr12", norm_MRP7[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP8_chr12", norm_MRP8[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP9_chr12", norm_MRP9[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])
assign("norm_MRP10_chr12", norm_MRP10[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])

print('per 3 chr matrix')

mean_MRP7<-mean(as.matrix(norm_MRP7))
mean_MRP8<-mean(as.matrix(norm_MRP8))
mean_MRP9<-mean(as.matrix(norm_MRP9))
mean_MRP10<-mean(as.matrix(norm_MRP10))

sd_MRP7<-sd(as.matrix(norm_MRP7))
sd_MRP8<-sd(as.matrix(norm_MRP8))
sd_MRP9<-sd(as.matrix(norm_MRP9))
sd_MRP10<-sd(as.matrix(norm_MRP10))

zscore_MRP7<-(norm_MRP7-mean_MRP7)/sd_MRP7
zscore_MRP8<-(norm_MRP8-mean_MRP8)/sd_MRP8
zscore_MRP9<-(norm_MRP9-mean_MRP9)/sd_MRP9
zscore_MRP10<-(norm_MRP10-mean_MRP10)/sd_MRP10

(zscore_MRP9-zscore_MRP7)->zscore_subtract_MRP7_9
(zscore_MRP8-zscore_MRP7)->zscore_subtract_MRP7_8
(zscore_MRP9-zscore_MRP8)->zscore_subtract_MRP8_9
(zscore_MRP10-zscore_MRP7)->zscore_subtract_MRP10_7

assign("zscore_subtract_MRP7_9_chr14_16", zscore_subtract_MRP7_9[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP7_8_chr14_16", zscore_subtract_MRP7_8[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP8_9_chr14_16", zscore_subtract_MRP8_9[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])
assign("zscore_subtract_MRP10_7_chr14_16", zscore_subtract_MRP10_7[(incrementalbins[14]+1):incrementalbins[16+1],(incrementalbins[14]+1):incrementalbins[16+1]])


# (optional) defines the color breaks manually for a "skewed" color transition
range=(0.000015*4428025)/size

col_breaks = c(seq((-5),-4,length=40),
               seq(((-3.99)),(-(3)),length=40),
               seq(-(2.99),-(2),length=40),
               seq(-(1.99),-(1),length=40),
               seq(-(0.99),0,length=40),
               seq(0.01,1,length=40),
               seq(1.01,2,length=40),
               seq(2.01,3,length=40),
               seq(3.01,4,length=40),
               seq(4.01,5,length=40))

incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3

# creates a own color palette from red to green
my_palette <- colorRampPalette(c(magma(15)[4],'darkorchid3','white','palegreen3','darkgreen'))(n = 399)

png("~/Desktop/zscore_subtract_MRP7_10_chr14_16.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP10_7_chr14_16,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          add.expr = abline(h=(incrementalbins3), v=incrementalbins2)
)
dev.off()

#Lets do it it with just ChrXII

assign("zscore_subtract_MRP10_7_chr12", zscore_subtract_MRP10_7[(incrementalbins[12]+1):incrementalbins[12+1],(incrementalbins[12]+1):incrementalbins[12+1]])


# creates a own color palette from red to green
#my_palette <- colorRampPalette(c("darkblue","blue4", "blue","cadetblue1", "white","darkgoldenrod1", "firebrick1", "red4", "darkred"))(n = 399)
# (optional) defines the color breaks manually for a "skewed" color transition
range=(0.000015*4428025)/size

col_breaks = c(seq((-5),-4,length=40),
               seq(((-3.99)),(-(3)),length=40),
               seq(-(2.99),-(2),length=40),
               seq(-(1.99),-(1),length=40),
               seq(-(0.99),0,length=40),
               seq(0.01,1,length=40),
               seq(1.01,2,length=40),
               seq(2.01,3,length=40),
               seq(3.01,4,length=40),
               seq(4.01,5,length=40))

incrementalbins2<-(incrementalbins[14:17]-915)
incrementalbins2[1]<-1
incrementalbins2<-incrementalbins2-0.5
incrementalbins2[4]<-incrementalbins2[4]+1

incrementalbins3<-incrementalbins2
incrementalbins3[2]<-(incrementalbins2[4]-incrementalbins2[3])+0.5
incrementalbins3[3]<-(incrementalbins2[4]-incrementalbins2[2])+0.5

incrementalbins2[1]<-incrementalbins2[1]+0.3
incrementalbins3[1]<-incrementalbins3[1]+0.25
incrementalbins2[4]<-incrementalbins2[4]-0.2
incrementalbins3[4]<-incrementalbins3[4]-0.3


png("~/Desktop/zscore_subtract_MRP7_10_chr12.png",   # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2(zscore_subtract_MRP10_7_chr12,
          main = "TITLE", # heat map title
          notecol="black",      # change font color of cell labels to black
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,
          dendrogram="none",
          labRow = FALSE,
          labCol = FALSE,
          colsep=FALSE,
          rowsep=FALSE,
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          trace="none",
          key=TRUE,
          key.xlab="% of interactions",
          key.ylab=NULL,
          density.info="none",
          key.ytickfun = NULL,
          key.title = NA,
          #add.expr = abline(h=incrementalbins_12, v=incrementalbins_12)
)
dev.off()

