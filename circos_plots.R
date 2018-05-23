#Circos plots

#04.30.18

#https://jokergoo.github.io/circlize_book/book/index.html


#Define the chromosome length
chr.size<-c(203893,794508,342718,1490682,602514,284456,1067526,544538,435585,719294,687260,1008248,908607,812465, 1054033,921188)
chr_bins<-c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13"
               , "chr14", "chr15", "chr16")

#Load circlize to create circos plots
if (!require("circlize")) {
  install.packages("circlize")
  library(circlize)
}

#Initilaize the circos plotting tracks
circos.par("track.height" = 0.15)
circos.initialize(factors = chr_bins, xlim=cbind(rep(0,16),chr.size))
circos.track(factors = chr_bins, ylim = cbind(rep(0,16),rep(1,16)))

if (!require("RCircos")) {
  install.packages("RCircos")
  library(RCircos)
}

chr.exclude <- NULL
tracks.inside <- 5
tracks.outside <- 0
cyto.info <- data.frame(cbind(chr_bins,rep(0,16),chr.size),rep(0,16),rep('gneg',16),stringsAsFactors = F)
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.cyto <- RCircos.Get.Plot.Ideogram()
rcircos.position <- RCircos.Get.Plot.Positions()
RCircos.List.Plot.Parameters()

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$radius.len <- 20
rcircos.params$base.per.unit <- 500
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.List.Plot.Parameters()
RCircos.Set.Plot.Area()

out.file <- "~/Desktop/RCircosDemoHumanGenome.pdf"
pdf(file=out.file, height=8, width=8, compress=TRUE)
RCircos.Set.Plot.Area()

par(mai=c(0.25, 0.25, 0.25, 0.25))
plot.new()
plot.window(c(-2.5,2.5), c(-2.5, 2.5))
RCircos.Chromosome.Ideogram.Plot()
dev.off()

##
#Didnt work too well. Try differnt package
###

source("https://bioconductor.org/biocLite.R")
biocLite('OmicCircos')
library(OmicCircos)

options(stringsAsFactors = FALSE) 

#########Trial with other data sets
data(UCSC.hg19.chr)
head(UCSC.hg19.chr)
#express
data( TCGA.BC.gene.exp.2k.60 )
#fusion
data( TCGA.BC.fus ) 

#sim data
seg.num<-16
ind.num<-20
seg.po<-ceiling(chr.size/10000)
link.num<-10
link.pg.num<-4

sim.out<-sim.circos(seg=seg.num , po=seg.po , ind=ind.num , link=link.num , link.pg=link.pg.num)
#seg.frame gives chr bins. Need 10 kb bins
seg.f<-sim.out$seg.frame
#seg.mapping is a ex[pression data. Call it from the circos with col.v to get specific column
seg.v<-sim.out$seg.mapping
#seg.link gives links between bins
link.v<-sim.out$seg.link
link.pg.v<-sim.out$seg.link.pg

seg.name <-paste0( "chr", 1 : seg.num , sep="" )
db<-segAnglePo(seg.frame , seg=chr_bins)

colors <- rainbow(seg.num, alpha=0.5)
par(mar=c( 2 , 2 , 2 , 2 ) )
plot( c(1, 800), c(1, 800), type="n" , axes=FALSE, xlab="", ylab="", main="")
circos(R=400, cir=db, type="chr", col=colors, print.chr.lab=TRUE, W=4, scale=TRUE)# mapping=seg.v, col.v=3, type="l" ->maybe add this and change it to centrlomere lcaltions

circos(R=360, cir=db , W=40, mapping=seg.v, col.v=3, type="l", B=TRUE, col=colors[1], lwd=2, scale=TRUE)

circos(R=350, cir=db , W=40, mapping=link.v , type="link" , lwd=2, col=colors[c(1 , 7)])
names(sim.out)

#Lets do this with my samples

#seg.frame sets up the segmentation 
init<-rep(0,3)
for (i in 1:16){
init<-rbind(init,cbind(chr_bins[i],0:(ceiling(chr.size/10000)[i]-1),1:ceiling(chr.size/10000)[i]))
}
seg.frame<-init[-1,]
seg.frame<-cbind(seg.frame, rep('NA',nrow(seg.frame)), rep('NA',nrow(seg.frame)))
colnames(seg.frame)<-c('seg.name','seg.Start', 'seg.End','the.v','NO')
seg.frame<-(data.frame(seg.frame, stringsAsFactors = F))
db<-segAnglePo(seg.frame , seg=chr_bins)

#seg.v is scores for different data sets I want for each. 
#seg.name seg.po  name1 name2 ......
#chr01 1 1.344 1.2344

#smc4 binned score dist
smc4_chip_trna<-read.table('~/Desktop/tRNA_score', stringsAsFactors = F,row.names=NULL,header=T)
smc4_chip_trna<-data.frame((smc4_chip_trna)[-1,],stringsAsFactors = F)

smc4_chip_binned<-read.table('~/Desktop/binned_score', stringsAsFactors = F)
smc4_chip_binned<-data.frame(t(smc4_chip_binned)[-1,],stringsAsFactors = F)

#convert the roman numerals into digits.
chr_roman_to_chr_numeric<-function(input){
roman<-strsplit(input,'chr')[[1]][2]
number<-as.numeric(as.roman(roman))
if(number<10){
output<-paste0('chr0',number)
}else{
output<-paste0('chr',number)
}
return(output)
}

smc4_chip_trna[,1]<-apply(data.frame(smc4_chip_trna[,1],stringsAsFactors = F),1,chr_roman_to_chr_numeric)
#transform into circos format
smc4_chip_trna_for_circos<-data.frame(cbind(smc4_chip_trna[,1],as.numeric(smc4_chip_trna[,2])/10000,as.numeric(smc4_chip_trna[,4])))

#Get my centromere information
cen_info<-read.table("~/Box Sync/Lab/Data_Analysis/Genomes_and_Subsets/SK1/140303_chr_cen_INFO.txt", stringsAsFactors = F, sep="\t", header=T)
##transform into circos format
cen_info_circos<-data.frame(cbind(cen_info[,1],as.numeric(cen_info[,8])/10000,rep(1,16)),stringsAsFactors = F)
cen_info_circos<-data.frame(cbind(cen_info[,1],(as.numeric(cen_info[,4])/10000)-1,(as.numeric(cen_info[,5])/10000)+1,1:16),stringsAsFactors = F)
cen_info_circos[8,3]<-10

#Colour scale for trNA binding strength
rbPal1 <- colorRampPalette(c('white','red','darkred','black'))
Col_smc4 <- rbPal1(90)[(as.numeric(cut((as.numeric(smc4_chip_trna_for_circos[,3])),breaks = 90)))]


#poliii
#Lets read in poliii chip data and adapt it to make it relevant

#read in the file
poliii_chip<-read.csv('~/Box Sync/Lab/Data/Pol_III_Supplementary_Data_.csv',header = TRUE,as.is= TRUE, stringsAsFactors = F)

#remove the i from the names with an i
strip_i<-function(matrix){
  strsplit(matrix[1],'i')[[1]][2]} 

poliii_names<-apply(poliii_chip,1,strip_i)
poliii_chip[which(is.na(poliii_names)==TRUE),1]->poliii_names[which(is.na(poliii_names)==TRUE)]

#which ones are tRNAs
poliii_tRNA_index<-which(substr(poliii_names,1,1)=='t')

#which ones have tRNAs in the description
poliii_near_genes<-poliii_chip[grep('\\(' , poliii_chip[,3]),3]

strip_tRNA<-function(mat){
  space_split<-strsplit(as.character(mat[3])," ")[[1]]
  if(length(space_split)<2){
    output<-mat[3]
  }else{
    output<-space_split[grep('\\(',space_split)]
  }
  return(output)
}

poliii_nearby_trna<-apply(poliii_chip,1,strip_tRNA)


#Which ones that arent tRNAs have adjacent trnas

poliii_nearby_trna_index<-which(!grepl('(0)',poliii_nearby_trna)=='TRUE')

index_for_trna<-setdiff(poliii_nearby_trna_index,poliii_tRNA_index)


poliii_names_2<-poliii_names

poliii_names_2[index_for_trna][-which(poliii_nearby_trna[index_for_trna]=='')]<-poliii_nearby_trna[index_for_trna][-which(poliii_nearby_trna[index_for_trna]=='')]

cbind(poliii_names_2,poliii_chip[,2])->poliii_trna_score

strip_tRNA<-function(mat){
  strsplit(as.character(mat[9]),'=')[[1]][3]}

tRNA_gff<-cbind(tRNA_gff,apply(tRNA_gff,1,strip_tRNA))

unlist(poliii_names_2)->poliii_names_2

poliii_trna_score<-(cbind(poliii_trna_score,(tRNA_gff[match(poliii_trna_score[,1], tRNA_gff[,14]),c(1,4,12)])))
poliii_trna_score<-poliii_trna_score[which(as.numeric(poliii_trna_score[,5])>0),]


poliii_trna_score_circos<-data.frame(cbind(poliii_trna_score[,3],as.numeric(poliii_trna_score[,4])/10000,poliii_trna_score[,2]),stringsAsFactors = F)

poliii_trna_score_circos[,1]<-apply(data.frame(poliii_trna_score_circos[,1],stringsAsFactors = F),2,chr_roman_to_chr_numeric)

rbPal2 <- colorRampPalette(c('white','blue','darkblue','black'))
Col_poliii<- rbPal2(90)[(as.numeric(cut((as.numeric(poliii_trna_score_circos[,3])),breaks = 90)))]

smc4_chip_trna_for_circos<-cbind(smc4_chip_trna_for_circos,1)
poliii_trna_score_circos<-cbind(poliii_trna_score_circos,1)

smc4_chip_trna_for_circos2<-cbind(smc4_chip_trna_for_circos[,c(1,2)], (as.numeric(smc4_chip_trna_for_circos[,2])+2), 1:nrow(smc4_chip_trna_for_circos))
poliii_trna_score_circos2<-cbind(poliii_trna_score_circos[,c(1,2)], (as.numeric(poliii_trna_score_circos[,2])+2), 1:nrow(poliii_trna_score_circos))

#trna links
intra_len<-nrow(my_intras)
inter_len<-nrow(my_inters)
trna_links<-rbind(my_inters[,c(1,2,3,7,8,9)],my_intras[,c(1,2,3,7,8,9)])

trna_links[,1]<-apply(data.matrix(trna_links[,1]),1,chr_roman_to_chr_numeric)
trna_links[,4]<-apply(data.matrix(trna_links[,4]),1,chr_roman_to_chr_numeric)

trna_links[,2]<-ceiling(trna_links[,2]/10000)
trna_links[,5]<-ceiling(trna_links[,5]/10000)

for (i in 1:(intra_len+inter_len)){
trna_links[i,3]<-paste0('n',i)
trna_links[i,6]<-paste0('n',i)}

#link.col=c(rep('red',inter_len),rep('blue',intra_len))
link.col.intra=c(rep('blue',intra_len))
link.col.inter=c(rep('red',inter_len))

#Lets draw this circos plot out.
pdf('~/Box Sync/Lab/Data_Analysis/tRNA/circos_plot.pdf')
colors <- rainbow(seg.num, alpha=0.5)
par(mar=c( 1 , 1 , 8 , 10 ),bg = "white")
plot.new()
plot(c(1, 800), c(1, 800), type="n" , axes=FALSE, xlab="", ylab="", main="")
circos(R=280, cir=db, type="chr", col='darkgray', print.chr.lab=TRUE, W=5, scale=TRUE)# mapping=seg.v, col.v=3, type="l" ->maybe add this and change it to centrlomere lcaltions
circos(R=277,cir=db,type="arc2",mapping=cen_info_circos,col='#00A089',col.v=3,W=5,cex=2,cutoff=0,lwd=10)
circos(R=227,cir=db,type="arc2",mapping=smc4_chip_trna_for_circos2,col=Col_smc4,col.v=4,W=20,scale=F,cex=2,cutoff=0,lwd=10, B=TRUE)
circos(R=187,cir=db,type="arc2",mapping=poliii_trna_score_circos2,col=Col_poliii,col.v=4,W=20,scale=F,cex=2,cutoff=0,lwd=10, B=TRUE)
circos(R=150, cir=db , W=10, mapping=trna_links[1:inter_len,] , type="link" , lwd=1, col='black')
circos(R=150, cir=db , W=1, mapping=trna_links[(inter_len+1):nrow(trna_links),] , type="link" , lwd=1, col='#42B549')


legend_image <- as.raster(rev(matrix(rbPal2(90), ncol=1)))
rasterImage(legend_image, 750, 100,790,300, xpd=T)
text(x=820, y = seq(100,272.4137,l=4), labels = seq(0,1.5,by=0.5), xpd=T)
text(x=950, y = 200, labels = 'PolIII subunit\n enrichment', xpd=T)

legend_image <- as.raster(rev(matrix(rbPal1(90), ncol=1)))
rasterImage(legend_image, 750, 400,790,600, xpd=T)
text(x=820, y = seq(400,600,l=4), labels = seq(0,60,by=20), xpd=T)
text(x=950, y = 500, labels = 'SMC4-Pk9 \n enrichment', xpd=T)

dev.off()

