#Files read in
binned_overlap<-read.table('/scratch/mrp420/binned_overlap', stringsAsFactors = F, header = F)
tRNA_overlap<-read.table('/scratch/mrp420/tRNA_overlap', stringsAsFactors = F, header = F)
smc4_chip<-read.table('/scratch/mrp420/AH6408K-031517-SK1Yue-PM_B3W3_MACS2_FE.bdg', stringsAsFactors = F)

print('Files read in')

chr.names<-list("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", 
"chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")

#set fucntions

subset_1<-function(unique_locations, overlap, chip){
  ind<-which(overlap[,4]==unique_locations[1]&overlap[,5]==unique_locations[2])
  ind<-matrix(ind,length(ind),1)
  ind2<-apply(ind,1,subset_2,chip=smc4_chip, overlap=overlap)
  return(c(unique_locations,mean(as.numeric(ind2))))
  }

subset_2<-function(ind,chip, overlap){
  ind2<-which(chip[,1]==overlap[ind,1]&chip[,2]==overlap[ind,2])
  return(chip[ind2,4])
}

#Apply fucntion to trnas
#locations<-unique(tRNA_overlap[,4:6])
#tRNA_score<-t(apply(locations,1,subset_1,chip=smc4_chip, overlap=tRNA_overlap))
#print('Applied to trnas')
#write.table(tRNA_score,'/scratch/mrp420/tRNA_score',quote=F,row.names=F)

chr_strip<-function(){
  chr_num<-strsplit(as.character(mat),'score_')[[1]][2]
  return(chr_num)
}

print('File outputted')
locations<-unique(binned_overlap[,4:6])
key<-sample(1:16)
for(i in 1:length(key)){
j<-key[i]
ind<-grep('chr',list.files('/scratch/mrp420/'))
if(length(grep('chr',list.files('/scratch/mrp420/')))!=0){
seen_before<-data.frame(list.files('/scratch/mrp420/')[grep('chr',list.files('/scratch/mrp420/'))], ,stringsAsFactors=F)
chr_matched<-apply(seen_before,1,chr_strip)
}else{
chr_matched<-'NO MATCH'}
if(length(which(chr_matched==chr.names[[j]]))==0){
chr_locations<-which(locations[,1]==chr.names[[j]])
chr_chip<-which(smc4_chip[,1]==chr.names[[j]])
chr_overlap<-which(binned_overlap[,4]==chr.names[[j]])
bin_score<-t(apply(locations[chr_locations,],1,subset_1,chip=smc4_chip[chr_chip,], overlap=binned_overlap[chr_overlap,]))
print('Applied to all bins')
write.table(bin_score,paste0('/scratch/mrp420/binned_score_',chr.names[[j]]),quote=F,row.names=F)
print(chr.names[[j]])
print('File outputted')
}else{
print(chr.names[[j]])
print('Already done')
}}
