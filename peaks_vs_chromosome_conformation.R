#######################################################
# How does condensin inform local conformation change #
#######################################################

#Do the regions that change most correlate with peaks of condnesin binding.

#Load plyr for
if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}

options(stringsAsFactors = FALSE)

#Read in the chip data and assign peaks to the same bins used in conformation analysis 
smc4_chip_peaks<-read.table('~/Box Sync/Lab/Paper/chip_stuff/ChIP_data-selected/AH6408K-031517-SK1Yue-B3W3-MACS2/AH6408K-031517-SK1Yue-PM_B3W3_MACS2_peaks.narrowPeak', stringsAsFactors = F, header = F)

smc4_chip_peak_bins<-data.frame(cbind(smc4_chip_peaks[,1],ceiling(smc4_chip_peaks[,2]/10000)),stringsAsFactors = F)
colnames(smc4_chip_peak_bins)<-c('chr','bin')

smc4_chip_peak_bins_counts<-smc4_chip_peak_bins %>% group_by(chr, bin) %>% tally()

#Lets read in the conformation matrices

wtrapa<-read.table('~/Box Sync/Lab/Data_Analysis/Heatmaps/MRP7-10/MRP7/heatmap', stringsAsFactors = F)
conrapa<-read.table('~/Box Sync/Lab/Data_Analysis/Heatmaps/MRP7-10/MRP10/heatmap', stringsAsFactors = F)

wtrapa_norm<-data.matrix(wtrapa/sum(wtrapa))
conrapa_norm<-data.matrix(conrapa/sum(conrapa))

norm_log2FC<-log2(conrapa_norm/wtrapa_norm)
norm_log2FC[which(norm_log2FC=='Inf'|norm_log2FC=='-Inf'|norm_log2FC=='NaN')]<-0

#Lets get the chromosome sizes
chr.size<-c(203893,794508,342718,1490682,602514,284456,1067526,544538,435585,719294,687260,1008248,908607,812465, 1054033,921188)
res<-10000
chr_bins_number<-list("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13"
               , "chr14", "chr15", "chr16")
chr_bins_roman<-c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII"
            , "chrXIV", "chrXV", "chrXVI")
bin_per_chr<-floor(chr.size/res+1)
incrementalbins<-matrix(0,length(chr_bins),1)
for (i in 1:(length(chr_bins)-1)){
  incrementalbins[i+1,]<-incrementalbins[i,]+bin_per_chr[i]
}
incrementalbins<-rbind(incrementalbins, 1196)

#Take each peak bin and convert into indices that apply to genome wide matrix
smc4_chip_peak_bins_counts<-cbind(data.frame(smc4_chip_peak_bins_counts, stringsAsFactors = F),rep('NA',nrow(smc4_chip_peak_bins_counts),stringsAsFactors = F),rep('NA',nrow(smc4_chip_peak_bins_counts),stringsAsFactors = F))
colnames(smc4_chip_peak_bins_counts)<-c('chr','bin','peaks','bin_increm')
for (i in 1:16){
ind<-which(smc4_chip_peak_bins_counts[,1]==chr_bins_roman[i])
smc4_chip_peak_bins_counts[ind,4]<-as.numeric(smc4_chip_peak_bins_counts[ind,2])+incrementalbins[i]
smc4_chip_peak_bins_counts[ind,5]<-i
print(i)
}

#I need a list of every interaction
temp1<-rep('NA',7)
for (i in 1:nrow(smc4_chip_peak_bins_counts)){
  left<-smc4_chip_peak_bins_counts[-i,c(5,4,3)]
  right<-smc4_chip_peak_bins_counts[i,c(5,4,3)]
  temp2<-rep('NA',nrow(left))
  for (j in 1:nrow(left)){
    temp2[j]<-norm_log2FC[as.numeric(left[j,2]),as.numeric(right[2])]
  }
  left<-cbind(left,temp2)
  temp1<-rbind(temp1,(cbind(right,left)))
}
allvsall<-temp1[-1,]

#Plot out
#plot(allvsall[,3],allvsall[,6], xlab='ChIP peaks 1', ylab='ChIP peaks 2')

#Plot with color gradients
rbPal <- colorRampPalette(c('blue','white','red'))
layout(matrix(1:1))
par(mar = c(7, 6, 4,2) + 0.2)
Col <- rbPal(90)[(as.numeric(cut((as.numeric(allvsall[,7])),breaks = 81)))+9]

#plot(jitter(as.numeric(allvsall[,3]), amount=0.3),jitter(as.numeric(allvsall[,6]), amount=0.3), col=Col, pch=4, xlab='# ChIP Peaks', ylab='# ChIP Peaks')

#pdf('~/Box Sync/Lab/Data_Analysis/ChIP_conformation_correlation/peaksvspeak.pdf')
par(mar = c(7, 6, 8,8) + 0.2)
plot(jitter(as.numeric(allvsall[,3]), amount=0.31),jitter(as.numeric(allvsall[,6]), amount=0.31), col=Col, pch=4, xlab='# ChIP Peaks', ylab='# ChIP Peaks', main='LogFC between bins with differnt amounts of SMC4 peaks')
legend_image <- as.raster(rev(matrix(rbPal(90), ncol=1)))
rasterImage(legend_image, 16, 3, 15,11, xpd=T)
text(x=16.5, y = seq(3.470588,10.529412,l=9), labels = seq(-4,4,by=1), xpd=T)
text(x=16, y= 12.5, labels='Log2FC\n(Brn1-FRB/WT)', xpd=T)
#dev.off()


# 3D plot
source("https://bioconductor.org/biocLite.R")
biocLite('scatterplot3d')

scatterplot3d(allvsall[,3],allvsall[,7],allvsall[,6], point.col=Col, xlab='ChIP peaks 1', zlab='ChIP peaks 2', ylab='Interaction Score' ,pch=19)

biocLite('rgl')
library(rgl)
biocLite('evd')
library(evd)
biocLite('car')
library(car)
scatter3d(allvsall[,3],allvsall[,7],allvsall[,6], color=Col, xlab='ChIP peaks 1', zlab='ChIP peaks 2', ylab='Interaction Score' ,pch=19)


#Try again but I want to do it agin with only intra chromosomal
temp1<-rep('NA',7)
for (i in 1:nrow(smc4_chip_peak_bins_counts)){
  left_A<-smc4_chip_peak_bins_counts[-i,c(5,4,3)]
  right<-smc4_chip_peak_bins_counts[i,c(5,4,3)]
  left_B<-left_A[which(as.numeric(left_A[,1])==unlist(right[1])),]
  temp2<-rep('NA',nrow(left_B))
  for (j in 1:nrow(left_B)){
    temp2[j]<-norm_log2FC[as.numeric(left_B[j,2]),as.numeric(right[2])]
  }
  left_B<-cbind(left_B,temp2)
  temp1<-rbind(temp1,(cbind(right,left_B)))
}
allvsall_intra_only<-temp1[-1,]


#Plot with color gradients
rbPal <- colorRampPalette(c('blue','white','red'))
layout(matrix(1:1))
par(mar = c(7, 6, 4,2) + 0.2)
Col <- rbPal(90)[(as.numeric(cut((as.numeric(allvsall_intra_only[,7])),breaks = 77)))+9]

#pdf('~/Box Sync/Lab/Data_Analysis/ChIP_conformation_correlation/peaksvspeak_intra_only.pdf')
par(mar = c(7, 6, 8,8) + 0.2)
plot(jitter(as.numeric(allvsall_intra_only[,3]), amount=0.31),jitter(as.numeric(allvsall_intra_only[,6]), amount=0.31), col=Col, pch=4, xlab='# ChIP Peaks', ylab='# ChIP Peaks', main='LogFC between bins with differnt amounts of SMC4 peaks')
legend_image <- as.raster(rev(matrix(rbPal(90), ncol=1)))
rasterImage(legend_image, 16, 3, 15,11, xpd=T)
text(x=16.5, y = seq(3.470588,10.529412,l=9), labels = seq(-4,4,by=1), xpd=T)
text(x=16, y= 12.5, labels='Log2FC\n(Brn1-FRB/WT)', xpd=T)
#dev.off()



#Lets look at intra scores
#First extract all intra scores
temp1<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  temp2<-c('NA')
  for (j in 1:length(chr_ind)){
    my_scores<-as.numeric(norm_log2FC[as.numeric(smc4_chip_peak_bins_counts[j,4]),chr_range])
    temp2<-c(temp2,my_scores)
    print(paste0(j,'A'))}
  temp2<-temp2[-1]
  temp1<-c(temp1,temp2)
  print(paste0(i,'B'))
}
all_intra_peakvspeakorother_scores<-temp1[-1]

all_intra_peakvspeakorother_scores<-(as.numeric(all_intra_peakvspeakorother_scores))


#Next just get scores from between peaks

coords<-cbind(as.numeric(allvsall_intra_only[,2]),as.numeric(allvsall_intra_only[,5]))
all_intra_peakvspeak_scores<-norm_log2FC[coords]


#Next get scores betwen other regions
temp1<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  onewaysubset<-norm_log2FC[chr_range,chr_range][-as.numeric(smc4_chip_peak_bins_counts[chr_ind,2]),]
  twowaysubset<-onewaysubset[-as.numeric(smc4_chip_peak_bins_counts[chr_ind,2])]
  temp1<-c(temp1,as.vector(twowaysubset))
}  
all_intra_othervsother_scores<-temp1[-1]
all_intra_othervsother_scores<-as.numeric(all_intra_othervsother_scores)  



head(all_intra_peakvspeak_scores)
boxplot(all_intra_othervsother_scores,all_intra_peakvspeakorother_scores, outline=F)

boxplot(all_intra_othervsother_scores[-which(all_intra_othervsother_scores==0)],all_intra_peakvspeakorother_scores[-which(all_intra_peakvspeakorother_scores==0)],all_intra_peakvspeak_scores[-which(all_intra_peakvspeak_scores==0)], outline=F)

boxplot(all_intra_peakvspeak_scores,all_intra_peakvspeakorother_scores,all_intra_othervsother_scores, outline=F)

#Subset to within 100Kb
temp1<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  left<-as.numeric(smc4_chip_peak_bins_counts[chr_ind,2])-10
  right<-as.numeric(smc4_chip_peak_bins_counts[chr_ind,2])+10
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  chr_range2<-bin_per_chr[i]
  left[which(left<1)]<-1
  right[which(right>max(chr_range2))]<-max(chr_range2)
  temp2<-c('NA')
  for (j in 1:length(left)){
    my_scores<-norm_log2FC[chr_range,chr_range][as.numeric(smc4_chip_peak_bins_counts[chr_ind[i][j],2]),left[j]:right[j]]
    temp2<-c(temp2,my_scores)
    print(paste0(j,'A'))}
  temp2<-temp2[-1]
  temp1<-c(temp1,temp2)
  print(paste0(i,'B'))
}
all_100kb_intra_peakvspeakorother_scores<-temp1[-1]
all_100kb_intra_peakvspeakorother_scores<-as.numeric(all_100kb_intra_peakvspeakorother_scores)

temp1<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  left<-as.numeric(smc4_chip_peak_bins_counts[chr_ind,2])-10
  right<-as.numeric(smc4_chip_peak_bins_counts[chr_ind,2])+10
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  chr_range2<-bin_per_chr[i]
  left[which(left<1)]<-1
  right[which(right>max(chr_range2))]<-max(chr_range2)
  temp2<-c('NA')
  for (j in 1:length(left)){
    
    #to do this options. maybe exlcuse just the lines across. 
    #then call out the border pieces
    
    onewaysubset<-as.vector(norm_log2FC[chr_range,chr_range][-as.numeric(smc4_chip_peak_bins_counts[chr_ind[i][j],2]),])
    twowaysubset<-as.vector(norm_log2FC[chr_range,chr_range][as.numeric(smc4_chip_peak_bins_counts[chr_ind[i][j],2]),-(left[j]:right[j])])
    
    temp2<-c(temp2,onewaysubset,twowaysubset)
    print(paste0(j,'A'))}
  temp2<-temp2[-1]
  temp1<-c(temp1,temp2)
  print(paste0(i,'B'))
}
all_100kb_intra_othervsother_scores<-temp1[-1]
all_100kb_intra_othervsother_scores<-as.numeric(all_100kb_intra_othervsother_scores)


#peaks to peaks less then 100kb
coords<-cbind(as.numeric(allvsall_intra_only[,2]),as.numeric(allvsall_intra_only[,5]))
coords[which(abs(coords[,1]-coords[,2])<11),]->close_coords
all_100kb_intra_peakvspeak_scores<-norm_log2FC[close_coords]

pdf('~/Box Sync/Lab/Data_Analysis/ChIP_conformation_correlation/condensin_peak_and_nopeak_changes.pdf')
par(mar = c(8, 6, 2,2) + 0.2)
boxplot(all_intra_peakvspeak_scores,all_intra_peakvspeakorother_scores,all_intra_othervsother_scores,all_100kb_intra_peakvspeak_scores,all_100kb_intra_peakvspeakorother_scores,all_100kb_intra_othervsother_scores, outline=F, xaxt='n', ylab='log2FC(Brn1-FRB/WT)',ylim=c(-2.5,2.5))
abline(h=0, col='gray')
axis(1,at=1:6,lab=c('Peak to peak','Peak to intra','No peak to\n no peak','Peak to peak \n<100kb','Peak to intra \n< 00kb','No peak to \nno peak \n<100kb'),las=2)
dev.off()

wilcox.test(all_100kb_intra_othervsother_scores,all_100kb_intra_peakvspeak_scores)
wilcox.test(all_intra_othervsother_scores,all_intra_peakvspeak_scores)
wilcox.test(all_intra_peakvspeakorother_scores,all_intra_peakvspeak_scores)
wilcox.test(all_intra_peakvspeakorother_scores,all_intra_othervsother_scores)

boxplot(all_intra_peakvspeak_scores[-which(all_intra_peakvspeak_scores==0)],all_intra_peakvspeakorother_scores[-which(all_intra_peakvspeakorother_scores==0)],all_intra_othervsother_scores[-which(all_intra_othervsother_scores==0)],all_100kb_intra_peakvspeak_scores[-which(all_100kb_intra_peakvspeak_scores==0)],all_100kb_intra_peakvspeakorother_scores[-which(all_100kb_intra_peakvspeakorother_scores==0)],all_100kb_intra_othervsother_scores[-which(all_100kb_intra_othervsother_scores==0)], outline=F)



#Redo all but exclude XII
smc4_chip_peak_bins_counts_no12<-smc4_chip_peak_bins_counts[-which(smc4_chip_peak_bins_counts[,5]==12),]

#I need a list of every interaction
temp1<-rep('NA',7)
for (i in 1:nrow(smc4_chip_peak_bins_counts_no12)){
  left<-smc4_chip_peak_bins_counts_no12[-i,c(5,4,3)]
  right<-smc4_chip_peak_bins_counts_no12[i,c(5,4,3)]
  temp2<-rep('NA',nrow(left))
  for (j in 1:nrow(left)){
    temp2[j]<-norm_log2FC[as.numeric(left[j,2]),as.numeric(right[2])]
  }
  left<-cbind(left,temp2)
  temp1<-rbind(temp1,(cbind(right,left)))
}
allvsall<-temp1[-1,]

#Plot out
#plot(allvsall[,3],allvsall[,6], xlab='ChIP peaks 1', ylab='ChIP peaks 2')

#Plot with color gradients
rbPal <- colorRampPalette(c('blue','white','red'))
layout(matrix(1:1))
par(mar = c(7, 6, 4,2) + 0.2)
Col <- rbPal(90)[(as.numeric(cut((as.numeric(allvsall[,7])),breaks = 81)))+9]

#plot(jitter(as.numeric(allvsall[,3]), amount=0.3),jitter(as.numeric(allvsall[,6]), amount=0.3), col=Col, pch=4, xlab='# ChIP Peaks', ylab='# ChIP Peaks')

pdf('~/Box Sync/Lab/Data_Analysis/ChIP_conformation_correlation/peaksvspeak_no12.pdf')
par(mar = c(7, 6, 8,8) + 0.2)
plot(jitter(as.numeric(allvsall[,3]), amount=0.31),jitter(as.numeric(allvsall[,6]), amount=0.31), col=Col, pch=4, xlab='# ChIP Peaks', ylab='# ChIP Peaks', main='LogFC between bins with different \n amounts of SMC4 peaks (no C12)')
legend_image <- as.raster(rev(matrix(rbPal(90), ncol=1)))
rasterImage(legend_image, 16, 3, 15,11, xpd=T)
text(x=16.5, y = seq(3.470588,10.529412,l=9), labels = seq(-4,4,by=1), xpd=T)
text(x=16, y= 12.5, labels='Log2FC\n(Brn1-FRB/WT)', xpd=T)
dev.off()





#Try again but I want to do it agin with only intra chromosomal
temp1<-rep('NA',7)
for (i in 1:nrow(smc4_chip_peak_bins_counts_no12)){
  left_A<-smc4_chip_peak_bins_counts_no12[-i,c(5,4,3)]
  right<-smc4_chip_peak_bins_counts_no12[i,c(5,4,3)]
  left_B<-left_A[which(as.numeric(left_A[,1])==unlist(right[1])),]
  temp2<-rep('NA',nrow(left_B))
  for (j in 1:nrow(left_B)){
    temp2[j]<-norm_log2FC[as.numeric(left_B[j,2]),as.numeric(right[2])]
  }
  left_B<-cbind(left_B,temp2)
  temp1<-rbind(temp1,(cbind(right,left_B)))
}
allvsall_intra_only<-temp1[-1,]


#Plot with color gradients
rbPal <- colorRampPalette(c('blue','white','red'))
layout(matrix(1:1))
par(mar = c(7, 6, 4,2) + 0.2)
Col <- rbPal(90)[(as.numeric(cut((as.numeric(allvsall_intra_only[,7])),breaks = 77)))+9]

pdf('~/Box Sync/Lab/Data_Analysis/ChIP_conformation_correlation/peaksvspeak_intra_only_no12.pdf')
par(mar = c(7, 6, 8,8) + 0.2)
plot(jitter(as.numeric(allvsall_intra_only[,3]), amount=0.31),jitter(as.numeric(allvsall_intra_only[,6]), amount=0.31), col=Col, pch=4, xlab='# ChIP Peaks', ylab='# ChIP Peaks', main='LogFC between bins with different \n amounts of SMC4 peaks (no C12)')
legend_image <- as.raster(rev(matrix(rbPal(90), ncol=1)))
rasterImage(legend_image, 16, 3, 15,11, xpd=T)
text(x=16.5, y = seq(3.470588,10.529412,l=9), labels = seq(-4,4,by=1), xpd=T)
text(x=16, y= 12.5, labels='Log2FC\n(Brn1-FRB/WT)', xpd=T)
dev.off()




#####HERE!!!!







#include the inter inter interacitons. 



#Lets look at intra scores
#First extract all intra scores peakvspeak
temp1<-c('NA')
temp2<-c('NA')
temp3<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  my_scores<-as.numeric(wtrapa_norm[as.numeric(smc4_chip_peak_bins_counts[chr_ind,4]),as.numeric(smc4_chip_peak_bins_counts[chr_ind,4])])
  temp1<-c(temp1,my_scores)
  my_scores<-as.numeric(conrapa_norm[as.numeric(smc4_chip_peak_bins_counts[chr_ind,4]),as.numeric(smc4_chip_peak_bins_counts[chr_ind,4])])
  temp2<-c(temp2,my_scores)
  my_scores<-as.numeric(norm_log2FC[as.numeric(smc4_chip_peak_bins_counts[chr_ind,4]),as.numeric(smc4_chip_peak_bins_counts[chr_ind,4])])
  temp3<-c(temp3,my_scores)
}
all_intra_peakvspeak_scores_wtrapa_norm<-temp1[-1]
all_intra_peakvspeak_scores_wtrapa_norm<-(as.numeric(all_intra_peakvspeak_scores_wtrapa_norm))
all_intra_peakvspeak_scores_conrapa_norm<-temp2[-1]
all_intra_peakvspeak_scores_conrapa_norm<-(as.numeric(all_intra_peakvspeak_scores_conrapa_norm))
all_intra_peakvspeak_scores_log2FC_norm<-temp3[-1]
all_intra_peakvspeak_scores_log2FC_norm<-(as.numeric(all_intra_peakvspeak_scores_log2FC_norm))


#Next look at intra no peaks chromosomal
temp1<-c('NA','NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  no_peaks<-cbind(chr_bins_roman[i],chr_range[is.na(match(chr_range,as.numeric(smc4_chip_peak_bins_counts[chr_ind,4])))])
  temp1<-rbind(temp1,no_peaks)
}
no_peaks_increment_bins<-temp1[-1,]


temp1<-c('NA')
temp2<-c('NA')
temp3<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_ind_no_peak<-which(no_peaks_increment_bins[,1]==chr_bins_roman[i])
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  my_scores<-as.numeric(wtrapa_norm[as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2]),as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])])
  temp1<-c(temp1,my_scores)
  my_scores<-as.numeric(conrapa_norm[as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2]),as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])])
  temp2<-c(temp2,my_scores)
  my_scores<-as.numeric(norm_log2FC[as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2]),as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])])
  temp3<-c(temp3,my_scores)
}
all_intra_nopeakvsnopeak_scores_wtrapa_norm<-temp1[-1]
all_intra_nopeakvsnopeak_scores_wtrapa_norm<-(as.numeric(all_intra_nopeakvsnopeak_scores_wtrapa_norm))
all_intra_nopeakvsnopeak_scores_conrapa_norm<-temp2[-1]
all_intra_nopeakvsnopeak_scores_conrapa_norm<-(as.numeric(all_intra_nopeakvsnopeak_scores_conrapa_norm))
all_intra_nopeakvsnopeak_scores_log2FC_norm<-temp3[-1]
all_intra_nopeakvsnopeak_scores_log2FC_norm<-(as.numeric(all_intra_nopeakvsnopeak_scores_log2FC_norm))

#Peak vs no peak intra
temp1<-c('NA')
temp2<-c('NA')
temp3<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_ind_no_peak<-which(no_peaks_increment_bins[,1]==chr_bins_roman[i])
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  my_scores<-as.numeric(wtrapa_norm[as.numeric(smc4_chip_peak_bins_counts[chr_ind,4]),as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])])
  temp1<-c(temp1,my_scores)
  my_scores<-as.numeric(conrapa_norm[as.numeric(smc4_chip_peak_bins_counts[chr_ind,4]),as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])])
  temp2<-c(temp2,my_scores)
  my_scores<-as.numeric(norm_log2FC[as.numeric(smc4_chip_peak_bins_counts[chr_ind,4]),as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])])
  temp3<-c(temp3,my_scores)
}
all_intra_peakvsnopeak_scores_wtrapa_norm<-temp1[-1]
all_intra_peakvsnopeak_scores_wtrapa_norm<-(as.numeric(all_intra_peakvsnopeak_scores_wtrapa_norm))
all_intra_peakvsnopeak_scores_conrapa_norm<-temp2[-1]
all_intra_peakvsnopeak_scores_conrapa_norm<-(as.numeric(all_intra_peakvsnopeak_scores_conrapa_norm))
all_intra_peakvsnopeak_scores_log2FC_norm<-temp3[-1]
all_intra_peakvsnopeak_scores_log2FC_norm<-(as.numeric(all_intra_peakvsnopeak_scores_log2FC_norm))


#First extract all inter peakvsnopeak
temp1<-c('NA')
temp2<-c('NA')
temp3<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  chr_ind_no_peak<-which(no_peaks_increment_bins[,1]==chr_bins_roman[i])
  for (j in 1:length(chr_ind)){  
    my_scores<-as.numeric(wtrapa_norm[as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4]),as.numeric(no_peaks_increment_bins[is.na(match(as.numeric(no_peaks_increment_bins[,2]),chr_range)),2])])
    temp1<-c(temp1,my_scores)
    my_scores<-as.numeric(conrapa_norm[as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4]),as.numeric(no_peaks_increment_bins[is.na(match(as.numeric(no_peaks_increment_bins[,2]),chr_range)),2])])
    temp2<-c(temp2,my_scores)
    my_scores<-as.numeric(norm_log2FC[as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4]),as.numeric(no_peaks_increment_bins[is.na(match(as.numeric(no_peaks_increment_bins[,2]),chr_range)),2])])
    temp3<-c(temp3,my_scores)
  }
  print(i)}
all_inter_peakvsnopeak_scores_wtrapa_norm<-temp1[-1]
all_inter_peakvsnopeak_scores_wtrapa_norm<-(as.numeric(all_inter_peakvsnopeak_scores_wtrapa_norm))
all_inter_peakvsnopeak_scores_conrapa_norm<-temp2[-1]
all_inter_peakvsnopeak_scores_conrapa_norm<-(as.numeric(all_inter_peakvsnopeak_scores_conrapa_norm))
all_inter_peakvsnopeak_scores_log2FC_norm<-temp3[-1]
all_inter_peakvsnopeak_scores_log2FC_norm<-(as.numeric(all_inter_peakvsnopeak_scores_log2FC_norm))

#Peak vs peak inter
temp1<-c('NA')
temp2<-c('NA')
temp3<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  for (j in 1:length(chr_ind)){  
    my_scores<-as.numeric(wtrapa_norm[as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4]),as.numeric(smc4_chip_peak_bins_counts[,4][is.na(match(smc4_chip_peak_bins_counts[,4],chr_range))])])
    temp1<-c(temp1,my_scores)
    my_scores<-as.numeric(conrapa_norm[as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4]),as.numeric(smc4_chip_peak_bins_counts[,4][is.na(match(smc4_chip_peak_bins_counts[,4],chr_range))])])
    temp2<-c(temp2,my_scores)
    my_scores<-as.numeric(norm_log2FC[as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4]),as.numeric(smc4_chip_peak_bins_counts[,4][is.na(match(smc4_chip_peak_bins_counts[,4],chr_range))])])
    temp3<-c(temp3,my_scores)
  }
  print(i)}
all_inter_peakvspeak_scores_wtrapa_norm<-temp1[-1]
all_inter_peakvspeak_scores_wtrapa_norm<-(as.numeric(all_inter_peakvspeak_scores_wtrapa_norm))
all_inter_peakvspeak_scores_conrapa_norm<-temp2[-1]
all_inter_peakvspeak_scores_conrapa_norm<-(as.numeric(all_inter_peakvspeak_scores_conrapa_norm))
all_inter_peakvspeak_scores_log2FC_norm<-temp3[-1]
all_inter_peakvspeak_scores_log2FC_norm<-(as.numeric(all_inter_peakvspeak_scores_log2FC_norm))

#Nopeak vs nopeak inter
temp1<-c('NA')
temp2<-c('NA')
temp3<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  chr_ind_no_peak<-which(no_peaks_increment_bins[,1]==chr_bins_roman[i])
  for (j in 1:length(chr_ind_no_peak)){  
    my_scores<-as.numeric(wtrapa_norm[as.numeric(no_peaks_increment_bins[chr_ind_no_peak[j],2]),as.numeric(no_peaks_increment_bins[is.na(match(as.numeric(no_peaks_increment_bins[,2]),chr_range)),2])])
    temp1<-c(temp1,my_scores)
    my_scores<-as.numeric(conrapa_norm[as.numeric(no_peaks_increment_bins[chr_ind_no_peak[j],2]),as.numeric(no_peaks_increment_bins[is.na(match(as.numeric(no_peaks_increment_bins[,2]),chr_range)),2])])
    temp2<-c(temp2,my_scores)
    my_scores<-as.numeric(norm_log2FC[as.numeric(no_peaks_increment_bins[chr_ind_no_peak[j],2]),as.numeric(no_peaks_increment_bins[is.na(match(as.numeric(no_peaks_increment_bins[,2]),chr_range)),2])])
    temp3<-c(temp3,my_scores)
  }
  print(i)}
all_inter_nopeakvsnopeak_scores_wtrapa_norm<-temp1[-1]
all_inter_nopeakvsnopeak_scores_wtrapa_norm<-(as.numeric(all_inter_nopeakvsnopeak_scores_wtrapa_norm))
all_inter_nopeakvsnopeak_scores_conrapa_norm<-temp2[-1]
all_inter_nopeakvsnopeak_scores_conrapa_norm<-(as.numeric(all_inter_nopeakvsnopeak_scores_conrapa_norm))
all_inter_nopeakvsnopeak_scores_log2FC_norm<-temp3[-1]
all_inter_nopeakvsnopeak_scores_log2FC_norm<-(as.numeric(all_inter_nopeakvsnopeak_scores_log2FC_norm))


#First extract all intra <100kb peakvsnopeak 
temp1<-c('NA')
temp2<-c('NA')
temp3<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  chr_ind_no_peak<-which(no_peaks_increment_bins[,1]==chr_bins_roman[i])
  for (j in 1:length(chr_ind)){  
    range100kb<-c(as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4])-10, as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4])+10)
    if(range100kb[1]<(incrementalbins[i]+1)){range100kb[1]<-(incrementalbins[i]+1)}
    if(range100kb[2]>(incrementalbins[i+1])){range100kb[2]<-(incrementalbins[i+1])}
    left<-as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4])
    right<-as.numeric(no_peaks_increment_bins[chr_ind_no_peak,][which(as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])>range100kb[1] & as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])<range100kb[2]),2])
    my_scores<-as.numeric(wtrapa_norm[left,right])
    temp1<-c(temp1,my_scores)
    my_scores<-as.numeric(conrapa_norm[left,right])
    temp2<-c(temp2,my_scores)
    my_scores<-as.numeric(norm_log2FC[left,right])
    temp3<-c(temp3,my_scores)
  }
  print(i)}
intra_100kb_peakvsnopeak_scores_wtrapa_norm<-temp1[-1]
intra_100kb_peakvsnopeak_scores_wtrapa_norm<-(as.numeric(intra_100kb_peakvsnopeak_scores_wtrapa_norm))
intra_100kb_peakvsnopeak_scores_conrapa_norm<-temp2[-1]
intra_100kb_peakvsnopeak_scores_conrapa_norm<-(as.numeric(intra_100kb_peakvsnopeak_scores_conrapa_norm))
intra_100kb_peakvsnopeak_scores_log2FC_norm<-temp3[-1]
intra_100kb_peakvsnopeak_scores_log2FC_norm<-(as.numeric(intra_100kb_peakvsnopeak_scores_log2FC_norm))


#First extract all intra <100kb peakvspeak
temp1<-c('NA')
temp2<-c('NA')
temp3<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  chr_ind_no_peak<-which(no_peaks_increment_bins[,1]==chr_bins_roman[i])
  for (j in 1:length(chr_ind)){  
    range100kb<-c(as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4])-10, as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4])+10)
    if(range100kb[1]<(incrementalbins[i]+1)){range100kb[1]<-(incrementalbins[i]+1)}
    if(range100kb[2]>(incrementalbins[i+1])){range100kb[2]<-(incrementalbins[i+1])}
    left<-as.numeric(smc4_chip_peak_bins_counts[chr_ind[j],4])
    right<-as.numeric(smc4_chip_peak_bins_counts[chr_ind[-j],][which(as.numeric(smc4_chip_peak_bins_counts[chr_ind[-j],4])>range100kb[1] & as.numeric(smc4_chip_peak_bins_counts[chr_ind[-j],4])<range100kb[2]),4])
    my_scores<-as.numeric(wtrapa_norm[left,right])
    temp1<-c(temp1,my_scores)
    my_scores<-as.numeric(conrapa_norm[left,right])
    temp2<-c(temp2,my_scores)
    my_scores<-as.numeric(norm_log2FC[left,right])
    temp3<-c(temp3,my_scores)
  }
  print(i)}
intra_100kb_peakvspeak_scores_wtrapa_norm<-temp1[-1]
intra_100kb_peakvspeak_scores_wtrapa_norm<-(as.numeric(intra_100kb_peakvspeak_scores_wtrapa_norm))
intra_100kb_peakvspeak_scores_conrapa_norm<-temp2[-1]
intra_100kb_peakvspeak_scores_conrapa_norm<-(as.numeric(intra_100kb_peakvspeak_scores_conrapa_norm))
intra_100kb_peakvspeak_scores_log2FC_norm<-temp3[-1]
intra_100kb_peakvspeak_scores_log2FC_norm<-(as.numeric(intra_100kb_peakvspeak_scores_log2FC_norm))

#First extract all intra <100kb nopeakvsnopeak  #####EDIT
temp1<-c('NA')
temp2<-c('NA')
temp3<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  chr_ind_no_peak<-which(no_peaks_increment_bins[,1]==chr_bins_roman[i])
  for (j in 1:length(chr_ind_no_peak)){  
    range100kb<-c(as.numeric(no_peaks_increment_bins[chr_ind_no_peak[j],2])-10, as.numeric(no_peaks_increment_bins[chr_ind_no_peak[j],2])+10)
    if(range100kb[1]<(incrementalbins[i]+1)){range100kb[1]<-(incrementalbins[i]+1)}
    if(range100kb[2]>(incrementalbins[i+1])){range100kb[2]<-(incrementalbins[i+1])}
    left<-as.numeric(no_peaks_increment_bins[chr_ind_no_peak[j],2])
    right<-as.numeric(no_peaks_increment_bins[chr_ind_no_peak,][which(as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])>range100kb[1] & as.numeric(no_peaks_increment_bins[chr_ind_no_peak,2])<range100kb[2]),2])
    my_scores<-as.numeric(wtrapa_norm[left,right])
    temp1<-c(temp1,my_scores)
    my_scores<-as.numeric(conrapa_norm[left,right])
    temp2<-c(temp2,my_scores)
    my_scores<-as.numeric(norm_log2FC[left,right])
    temp3<-c(temp3,my_scores)
  }
  print(i)}
intra_100kb_nopeakvsnopeak_scores_wtrapa_norm<-temp1[-1]
intra_100kb_nopeakvsnopeak_scores_wtrapa_norm<-(as.numeric(intra_100kb_nopeakvsnopeak_scores_wtrapa_norm))
intra_100kb_nopeakvsnopeak_scores_conrapa_norm<-temp2[-1]
intra_100kb_nopeakvsnopeak_scores_conrapa_norm<-(as.numeric(intra_100kb_nopeakvsnopeak_scores_conrapa_norm))
intra_100kb_nopeakvsnopeak_scores_log2FC_norm<-temp3[-1]
intra_100kb_nopeakvsnopeak_scores_log2FC_norm<-(as.numeric(intra_100kb_nopeakvsnopeak_scores_log2FC_norm))



pdf('~/Box Sync/Lab/Data_Analysis/ChIP_conformation_correlation/condensin_peak_and_nopeak_changes_no12.pdf')
par(mar = c(8, 6, 4,2) + 0.2)
boxplot(all_intra_peakvspeak_scores,all_intra_peakvspeakorother_scores,all_intra_othervsother_scores,all_100kb_intra_peakvspeak_scores,all_100kb_intra_peakvspeakorother_scores,all_100kb_intra_othervsother_scores, outline=F, xaxt='n', ylab='log2FC(Brn1-FRB/WT)', main='interactions changes between regions that \n do or do not contain condensin peaks',ylim=c(-2.5,2.5))
abline(h=0, col='gray')
axis(1,at=1:6,lab=c('Peak to peak','Peak to intra','No peak to\n no peak','Peak to peak \n<100kb','Peak to intra \n< 00kb','No peak to \nno peak \n<100kb'),las=2)
dev.off()












temp1<-c('NA')
for (i in (1:16)[-12]){
  chr_ind<-which(smc4_chip_peak_bins_counts_no12[,5]==i)
  left<-as.numeric(smc4_chip_peak_bins_counts_no12[chr_ind,2])-10
  right<-as.numeric(smc4_chip_peak_bins_counts_no12[chr_ind,2])+10
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  chr_range2<-bin_per_chr[i]
  left[which(left<1)]<-1
  right[which(right>max(chr_range2))]<-max(chr_range2)
  temp2<-c('NA')
  for (j in 1:length(left)){
    my_scores<-norm_log2FC[chr_range,chr_range][as.numeric(smc4_chip_peak_bins_counts_no12[chr_ind[i][j],2]),left[j]:right[j]]
    temp2<-c(temp2,my_scores)
    print(paste0(j,'A'))}
  temp2<-temp2[-1]
  temp1<-c(temp1,temp2)
  print(paste0(i,'B'))
}
all_100kb_intra_peakvspeakorother_scores<-temp1[-1]
all_100kb_intra_peakvspeakorother_scores<-as.numeric(all_100kb_intra_peakvspeakorother_scores)












#Extract all inter peakvsnopeak
temp1<-c('NA')
temp2<-c('NA')
for (i in 1:16){
  chr_ind<-which(smc4_chip_peak_bins_counts[,5]==i)
  chr_range<-(incrementalbins[i]+1):(incrementalbins[i+1])
  for (j in 1:length(chr_range)){  
    my_scores<-as.numeric(wtrapa_norm[as.numeric(chr_range[j]),seq(1,1196,by=1)[is.na(match(1:1196,chr_range))]])
    temp1<-c(temp1,my_scores)
    my_scores<-as.numeric(conrapa_norm[as.numeric(chr_range[j]),seq(1,1196,by=1)[is.na(match(1:1196,chr_range))]])
    temp2<-c(temp2,my_scores)
  }
  print(i)}
all_inter_peakvsnopeak_scores_wtrapa_norm<-temp1[-1]
all_inter_peakvsnopeak_scores_wtrapa_norm<-(as.numeric(all_inter_peakvsnopeak_scores_wtrapa_norm))
all_inter_peakvsnopeak_scores_conrapa_norm<-temp2[-1]
all_inter_peakvsnopeak_scores_conrapa_norm<-(as.numeric(all_inter_peakvsnopeak_scores_conrapa_norm))


pdf('~/Box Sync/Lab/Data_Analysis/ChIP_conformation_correlation/all_compariosn_logfc.pdf')
boxplot(
  intra_100kb_peakvspeak_scores_log2FC_norm,
  intra_100kb_peakvsnopeak_scores_log2FC_norm,
  intra_100kb_nopeakvsnopeak_scores_log2FC_norm,
  all_intra_peakvspeak_scores_log2FC_norm,
  all_intra_peakvsnopeak_scores_log2FC_norm,
  all_intra_nopeakvsnopeak_scores_log2FC_norm,
  all_inter_peakvspeak_scores_log2FC_norm,
  all_inter_peakvsnopeak_scores_log2FC_norm,
  all_inter_nopeakvsnopeak_scores_log2FC_norm,
 outline=F
)
dev.off()

all_intra_peakvspeak_scores_log2FC_norm_nozeros<-all_intra_peakvspeak_scores_log2FC_norm[-which(all_intra_peakvspeak_scores_log2FC_norm==0)]
all_intra_peakvsnopeak_scores_log2FC_norm_nozeros<-all_intra_peakvsnopeak_scores_log2FC_norm[-which(all_intra_peakvsnopeak_scores_log2FC_norm==0)]
all_intra_nopeakvsnopeak_scores_log2FC_norm_nozeros<-all_intra_nopeakvsnopeak_scores_log2FC_norm[-which(all_intra_nopeakvsnopeak_scores_log2FC_norm==0)]
all_inter_peakvspeak_scores_log2FC_norm_nozeros<-all_inter_peakvspeak_scores_log2FC_norm[-which(all_inter_peakvspeak_scores_log2FC_norm==0)]
all_inter_peakvsnopeak_scores_log2FC_norm_nozeros<-all_inter_peakvsnopeak_scores_log2FC_norm[-which(all_inter_peakvsnopeak_scores_log2FC_norm==0)]
all_inter_nopeakvsnopeak_scores_log2FC_norm_nozeros<-all_inter_nopeakvsnopeak_scores_log2FC_norm[-which(all_inter_nopeakvsnopeak_scores_log2FC_norm==0)]
intra_100kb_peakvspeak_scores_log2FC_norm_nozeros<-intra_100kb_peakvspeak_scores_log2FC_norm[-which(intra_100kb_peakvspeak_scores_log2FC_norm==0)]
intra_100kb_peakvsnopeak_scores_log2FC_norm_nozeros<-intra_100kb_peakvsnopeak_scores_log2FC_norm[-which(intra_100kb_peakvsnopeak_scores_log2FC_norm==0)]
intra_100kb_nopeakvsnopeak_scores_log2FC_norm_nozeros<-intra_100kb_nopeakvsnopeak_scores_log2FC_norm[-which(intra_100kb_nopeakvsnopeak_scores_log2FC_norm==0)]

boxplot(
  intra_100kb_peakvspeak_scores_log2FC_norm_nozeros,
  intra_100kb_peakvsnopeak_scores_log2FC_norm_nozeros,
  intra_100kb_nopeakvsnopeak_scores_log2FC_norm_nozeros,
  all_intra_peakvspeak_scores_log2FC_norm_nozeros,
  all_intra_peakvsnopeak_scores_log2FC_norm_nozeros,
  all_intra_nopeakvsnopeak_scores_log2FC_norm_nozeros,
  all_inter_peakvspeak_scores_log2FC_norm_nozeros,
  all_inter_peakvsnopeak_scores_log2FC_norm_nozeros,
  all_inter_nopeakvsnopeak_scores_log2FC_norm_nozeros,
 outline=F
)


pdf('~/Box Sync/Lab/Data_Analysis/ChIP_conformation_correlation/all_compariosn_depletions.pdf')
boxplot(
  intra_100kb_peakvspeak_scores_wtrapa_norm,
  intra_100kb_peakvspeak_scores_conrapa_norm,
  intra_100kb_peakvsnopeak_scores_wtrapa_norm,
  intra_100kb_peakvsnopeak_scores_conrapa_norm,
  intra_100kb_nopeakvsnopeak_scores_wtrapa_norm,
  intra_100kb_nopeakvsnopeak_scores_conrapa_norm,
all_intra_peakvspeak_scores_wtrapa_norm,
all_intra_peakvspeak_scores_conrapa_norm,
all_intra_peakvsnopeak_scores_wtrapa_norm,
all_intra_peakvsnopeak_scores_conrapa_norm,
all_intra_nopeakvsnopeak_scores_wtrapa_norm,
all_intra_nopeakvsnopeak_scores_conrapa_norm,
all_inter_peakvspeak_scores_wtrapa_norm,
all_inter_peakvspeak_scores_conrapa_norm,
all_inter_peakvsnopeak_scores_wtrapa_norm,
all_inter_peakvsnopeak_scores_conrapa_norm,
all_inter_nopeakvsnopeak_scores_wtrapa_norm,
all_inter_nopeakvsnopeak_scores_conrapa_norm,
 outline=F
)
dev.off()

boxplot( intra_100kb_peakvspeak_scores_wtrapa_norm,
         intra_100kb_peakvspeak_scores_conrapa_norm,
         intra_100kb_peakvsnopeak_scores_wtrapa_norm,
         intra_100kb_peakvsnopeak_scores_conrapa_norm,
         intra_100kb_nopeakvsnopeak_scores_wtrapa_norm,
         intra_100kb_nopeakvsnopeak_scores_conrapa_norm, outline=F
)
boxplot(
  all_intra_peakvspeak_scores_wtrapa_norm,
  all_intra_peakvspeak_scores_conrapa_norm,
  all_intra_peakvsnopeak_scores_wtrapa_norm,
  all_intra_peakvsnopeak_scores_conrapa_norm,
  all_intra_nopeakvsnopeak_scores_wtrapa_norm,
  all_intra_nopeakvsnopeak_scores_conrapa_norm,outline=F
)
boxplot(all_inter_peakvspeak_scores_wtrapa_norm,
  all_inter_peakvspeak_scores_conrapa_norm,
  all_inter_peakvsnopeak_scores_wtrapa_norm,
  all_inter_peakvsnopeak_scores_conrapa_norm,
  all_inter_nopeakvsnopeak_scores_wtrapa_norm,
  all_inter_nopeakvsnopeak_scores_conrapa_norm,outline=F
  )


df<-c(intra_100kb_peakvspeak_scores_log2FC_norm,
       intra_100kb_peakvsnopeak_scores_log2FC_norm,
       intra_100kb_nopeakvsnopeak_scores_log2FC_norm,
       all_intra_peakvspeak_scores_log2FC_norm,
       all_intra_peakvsnopeak_scores_log2FC_norm,
       all_intra_nopeakvsnopeak_scores_log2FC_norm,
       all_inter_peakvspeak_scores_log2FC_norm,
       all_inter_peakvsnopeak_scores_log2FC_norm,
       all_inter_nopeakvsnopeak_scores_log2FC_norm)

key<-c(rep(1,length(intra_100kb_peakvspeak_scores_log2FC_norm)),
       rep(2,length(intra_100kb_peakvsnopeak_scores_log2FC_norm)),
       rep(3,length(intra_100kb_nopeakvsnopeak_scores_log2FC_norm)),
rep(4,length(all_intra_peakvspeak_scores_log2FC_norm)),
rep(5,length(all_intra_peakvsnopeak_scores_log2FC_norm)),
rep(6,length(all_intra_nopeakvsnopeak_scores_log2FC_norm)),
rep(7,length(all_inter_peakvspeak_scores_log2FC_norm)),
rep(8,length(all_inter_peakvsnopeak_scores_log2FC_norm)),
rep(9,length(all_inter_nopeakvsnopeak_scores_log2FC_norm)))

pairwise.wilcox.test(df,key,p.adjust.method ="bonferroni")





df<-c(intra_100kb_peakvspeak_scores_log2FC_norm_nozeros,
      intra_100kb_peakvsnopeak_scores_log2FC_norm_nozeros,
      intra_100kb_nopeakvsnopeak_scores_log2FC_norm_nozeros,
      all_intra_peakvspeak_scores_log2FC_norm_nozeros,
      all_intra_peakvsnopeak_scores_log2FC_norm_nozeros,
      all_intra_nopeakvsnopeak_scores_log2FC_norm_nozeros,
      all_inter_peakvspeak_scores_log2FC_norm_nozeros,
      all_inter_peakvsnopeak_scores_log2FC_norm_nozeros,
      all_inter_nopeakvsnopeak_scores_log2FC_norm_nozeros)

key<-c(rep(1,length(intra_100kb_peakvspeak_scores_log2FC_norm_nozeros)),
       rep(2,length(intra_100kb_peakvsnopeak_scores_log2FC_norm_nozeros)),
       rep(3,length(intra_100kb_nopeakvsnopeak_scores_log2FC_norm_nozeros)),
       rep(4,length(all_intra_peakvspeak_scores_log2FC_norm_nozeros)),
       rep(5,length(all_intra_peakvsnopeak_scores_log2FC_norm_nozeros)),
       rep(6,length(all_intra_nopeakvsnopeak_scores_log2FC_norm_nozeros)),
       rep(7,length(all_inter_peakvspeak_scores_log2FC_norm_nozeros)),
       rep(8,length(all_inter_peakvsnopeak_scores_log2FC_norm_nozeros)),
       rep(9,length(all_inter_nopeakvsnopeak_scores_log2FC_norm_nozeros)))

g.sizes<-c(length(intra_100kb_peakvspeak_scores_log2FC_norm_nozeros),
          length(intra_100kb_peakvsnopeak_scores_log2FC_norm_nozeros),
           length(intra_100kb_nopeakvsnopeak_scores_log2FC_norm_nozeros),
          length(all_intra_peakvspeak_scores_log2FC_norm_nozeros),
           length(all_intra_peakvsnopeak_scores_log2FC_norm_nozeros),
           length(all_intra_nopeakvsnopeak_scores_log2FC_norm_nozeros),
        length(all_inter_peakvspeak_scores_log2FC_norm_nozeros),
       length(all_inter_peakvsnopeak_scores_log2FC_norm_nozeros),
           length(all_inter_nopeakvsnopeak_scores_log2FC_norm_nozeros))

pairwise.wilcox.test(df,key,p.adjust.method ="bonferroni")

