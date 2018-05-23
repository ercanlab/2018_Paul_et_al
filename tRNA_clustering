###################
# tRNA clustering #
###################

#Load plyr for
if (!require("plyr")) {
  install.packages("plyr")
  library(plyr)
}

#Load viridis color palette for graphing
if (!require("viridis")) {
  install.packages("viridis")
  library(viridis)
}

#Read in genomewide interaction files
#Interaction matrix files were generated using the mirny lab balancing algorithm through the script .......(insert name here). 
#Decay niormalised matrices were generated the decay_correction.R script and consists of matrices that have interactions normalised relative to decay distance
wtrapa<-read.table('~/Box Sync/Lab/Data_Analysis/Heatmaps/MRP7-10/MRP7/heatmap')
conrapa<-read.table('~/Box Sync/Lab/Data_Analysis/Heatmaps/MRP7-10/MRP10/heatmap')

wtrapa_decaynorm<-read.csv('~/Box Sync/Lab/Data_Analysis/Heatmaps/Decay_normalised/wtrapa_decay_normalised.csv')
conrapa_decaynorm<-read.csv('~/Box Sync/Lab/Data_Analysis/Heatmaps/Decay_normalised/conrapa_decay_normalised.csv')

wtrapa<-data.matrix(wtrapa)
conrapa<-data.matrix(conrapa)
wtrapa_decaynorm<-data.matrix(wtrapa_decaynorm)
conrapa_decaynorm<-data.matrix(conrapa_decaynorm)


#Normalise to total interaction counts per dataset
wtrapa_norm<-data.matrix(wtrapa/sum(wtrapa))
conrapa_norm<-data.matrix(conrapa/sum(conrapa))
wtrapa_decaynorm_norm<-data.matrix(wtrapa_decaynorm/sum(wtrapa_decaynorm))
conrapa_decaynorm_norm<-data.matrix(conrapa_decaynorm/sum(conrapa_decaynorm))

#read in trna positions and assign to bins
gff<-read.table("~/Box Sync/Lab/Data_Analysis/Genomes_and_Subsets/SK1/SK1_annotation_modified_v2.gff", stringsAsFactors = F, sep="\t", quote="")
SC_gff<-read.table("~/Box Sync/Lab/Data_Analysis/Genomes_and_Subsets/S288C/saccharomyces_cerevisiae.gff", stringsAsFactors = F, sep="\t", quote="")
SP_gff<-read.table("~/Box Sync/Lab/Data_Analysis/Genomes_and_Subsets/s_pombe/schizosaccharomyces_pombe.chr.gff3.txt", stringsAsFactors = F, sep="\t", quote="")

tRNA_gff<-gff[grep(")", gff[,9]),]
tRNA_gff<-tRNA_gff[-grep("ALPHA", tRNA_gff[,9]),]

#Break down the tRNA into subtypes
tRNA_type<-function(input){
  as.character(strsplit(strsplit(strsplit(input[9],')')[[1]][1],'=')[[1]][2],"\\(")[[1]])
}
tRNA_types<-apply(tRNA_gff,1,tRNA_type)
tRNA_gff<-(cbind(tRNA_gff,t(tRNA_types)))

#Extract chromosome starts. This is done by assessing how many bins in each chromosome. Then incrmentally assessing at which cordinate the chromomse starts in the interaction matrix
chr.size<-c(203893,794508,342718,1490682,602514,284456,1067526,544538,435585,719294,687260,1008248,908607,812465, 1054033,921188)
res<-10000
chr_bins<-list("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13"
               , "chr14", "chr15", "chr16")
bin_per_chr<-floor(chr.size/res+1)
incrementalbins<-matrix(0,length(chr_bins),1)
for (i in 1:(length(chr_bins)-1)){
  incrementalbins[i+1,]<-incrementalbins[i,]+bin_per_chr[i]
}
incrementalbins<-rbind(incrementalbins, 1196)

#For randomisation. Need range for each chromosome to move around in. This is because we are usinga window larger then the bin. 
all_bins<-as.numeric(as.matrix(rbind(incrementalbins,'1196')))[-1]
chr_ranges<-cbind(incrementalbins+1,all_bins)
chr_ranges[,1]+2->chr_ranges[,1]
chr_ranges[,2]-2->chr_ranges[,2]
chr_ranges_trna<-chr_ranges
chr_ranges_trna[,1]<-chr_ranges_trna[,1]-2
chr_ranges_trna[,2]<-chr_ranges_trna[,2]+2
chr_ranges_trna<-chr_ranges_trna[-17,]
chr_ranges<-chr_ranges[-17,]

#read in trna positions and assign to bins
chr_bins<-c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII"
            , "chrXIV", "chrXV", "chrXVI")
trna_pos<-matrix(,,)
for (i in 1:16){
  pos<-which(tRNA_gff[,1]==chr_bins[i])
  tRNA_gff[pos,4]
  bin_loci<-ceiling(tRNA_gff[pos,4]/10000)+incrementalbins[i]
  trna_pos<-c(trna_pos,bin_loci)
}

#Get gene names for regular gff
gene_name<-function(input){
  as.character(strsplit(input[9],'=')[[1]][3])
}
gene_names<-data.frame(apply(gff,1,gene_name), stringsAsFactors = F)
gff<-cbind(gff,gene_names)


#ETCs
etc_names<-c('ADE8','ARG8','BCK1','RAD2','YMR201C','TFC6','WTM2','RPB5')#'ZOD1')

etc_gff<-data.frame(gff[match(etc_names,unlist(gene_names)),],stringsAsFactors = F)
etc_gff[9,]<-c('chrXIII','NA','NA',84694,84751,'.','+','.','NA','ZOD1','ZOD1')

etc_pos<-matrix(,,)
for (i in 1:16){
  pos<-which(etc_gff[,1]==chr_bins[i])
  etc_gff[pos,4]
  bin_loci<-ceiling(as.numeric(etc_gff[pos,4])/10000)+incrementalbins[i]
  etc_pos<-c(etc_pos,bin_loci)
}
etc_pos<-etc_pos[-1]

#all bins have been ID'd. Need to remove first item in array as blank, then minimise to unique bins
trna_pos<-trna_pos[-1]
trna_pos_unique<-unique(trna_pos)

tRNA_gff<-cbind(tRNA_gff,trna_pos)
tRNA_gff<-cbind(tRNA_gff,vector(mode="logical",nrow(tRNA_gff)))

#Want to break this down into specifc subsets of tRNAs. Which tRNAs belong to each codon. 
#Identify and cycle through each subgroup of tRNA and then assign the appropariate gneome wide bin
tRNA_groups<-matrix(0,1,2)
tRNA_aacids<-unique(tRNA_gff[,10])
for (i in 1:length(tRNA_aacids)){
  index<-which(tRNA_gff[,10]==tRNA_aacids[i])
  pos<-matrix(0,length(index),2)
  pos[,1]<-paste0(as.character(tRNA_aacids[i]),"_all")
  pos[,2]<-tRNA_gff[index,12]
  tRNA_groups<-rbind(tRNA_groups,pos)
  tRNA_codon<-unique(tRNA_gff[index,11])
  for (j in 1:length(tRNA_codon)){
    index<-which(tRNA_gff[,11]==tRNA_codon[j])
    pos<-matrix(0,length(index),2)
    pos[,1]<-paste0(as.character(tRNA_aacids[i]),as.character(tRNA_codon[j]))
    pos[,2]<-tRNA_gff[index,12]
    tRNA_groups<-rbind(tRNA_groups,pos)
  }
}

#Assign them out to the matrix
tRNA_groups<-as.data.frame(tRNA_groups[-1,], stringsAsFactors = F)
tRNA_all<-matrix("t_all",length(unique(tRNA_groups[,2])),2)
tRNA_all[,2]<-unique(tRNA_groups[,2])
tRNA_groups<-rbind(tRNA_groups,tRNA_all[order(as.numeric(tRNA_all[,2])),])

check1<-cbind('etc',etc_pos)
colnames(check1)<-c('','')
tRNA_groups<-rbind(tRNA_groups, check1)


tRNA_groups<-rbind(tRNA_groups,cbind('t_etc',c(tRNA_groups[which(tRNA_groups[,1]=="t_all"),2],check1[,2])))


#where are centromeres and telomeres. trying to assess if tethered already
cen_info<-read.table("~/Box Sync/Lab/Data_Analysis/Genomes_and_Subsets/SK1/140303_chr_cen_INFO.txt", stringsAsFactors = F, sep="\t")
ce_pos<-as.numeric(cen_info[-1,8])
ce_pos<-(ce_pos/10000)
ce_pos<-ceiling(ce_pos)
cen_in_matrix<-ce_pos+incrementalbins[1:16]

#Take my tRNA file and see if it overlaps with centromere or telomere. 
for (i in 1:16){
  index1<-which((tRNA_gff[,12]<=cen_in_matrix[i]+1) & (tRNA_gff[,12]>=(cen_in_matrix[i]-1)))
  index2<-which((tRNA_gff[,12]<=chr_ranges[i,1]) & (tRNA_gff[,12]>=(chr_ranges[i,1]-2)))
  index3<-which((tRNA_gff[,12]<=chr_ranges[i,2]+2) & (tRNA_gff[,12]>=(chr_ranges[i,2])))
  index<-c(index1,index2,index3)
  tRNA_gff[index,13]<-'TRUE'
}

#Subset to tRNAs without tethers i.e centromere. 
subset<-which(tRNA_gff[,13]!='TRUE')
tRNA_gff_notether<-tRNA_gff[subset,]

#Take my no tether group
tRNA_groups_notether<-matrix(0,1,2)
tRNA_aacids<-unique(tRNA_gff_notether[,10])
for (i in 1:length(tRNA_aacids)){
  index<-which(tRNA_gff_notether[,10]==tRNA_aacids[i])
  pos<-matrix(0,length(index),2)
  pos[,1]<-paste0(as.character(tRNA_aacids[i]),"_all")
  pos[,2]<-tRNA_gff_notether[index,12]
  tRNA_groups_notether<-rbind(tRNA_groups_notether,pos)
  tRNA_codon<-unique(tRNA_gff_notether[index,11])
  for (j in 1:length(tRNA_codon)){
    index<-which(tRNA_gff_notether[,11]==tRNA_codon[j])
    pos<-matrix(0,length(index),2)
    pos[,1]<-paste0(as.character(tRNA_aacids[i]),as.character(tRNA_codon[j]))
    pos[,2]<-tRNA_gff_notether[index,12]
    tRNA_groups_notether<-rbind(tRNA_groups_notether,pos)
  }
}

tRNA_groups_notether<-as.data.frame(tRNA_groups_notether[-1,], stringsAsFactors = F)
tRNA_all<-matrix("t_all",length(unique(tRNA_groups_notether[,2])),2)
tRNA_all[,2]<-unique(tRNA_groups_notether[,2])
tRNA_groups_notether<-rbind(tRNA_groups_notether,tRNA_all[order(as.numeric(tRNA_all[,2])),])

#chr with trnas (no trnas on 9)
chrs<-1:16
chrs<-chrs[-9]

#Isolate the top 1%
wtrapa_norm_cutoff<-(round(length(wtrapa_norm)/100))
wtrapa_norm_oneperc<-wtrapa_norm[order(wtrapa_norm, decreasing=TRUE)][cutoff]
wtrapa_norm_coords<-which(wtrapa_norm>oneperc, arr.ind=TRUE)

index<-as.matrix(as.vector(unique(tRNA_groups[,1])))


tRNAsintop<-function(coords) {
my_matches<-matrix(0,1,2)
my_final_mat<-matrix(0,1,6)
#need to define potential coordinates of interaction for each tRNA group
for (i in 1:length(index)){
  #For each group, where are its members in the genome
  index2<-as.numeric(as.vector(tRNA_groups[which(tRNA_groups[,1]==index[i]),2]))
  my_mat<-matrix(0,1,2)  
  #What can the tRNA group interact with
  for (j in 1:length(index2)){
    my_mat2<-matrix(0,(length(index2)-1),2)    
    my_mat2[,1]<-index2[j]
    my_mat2[,2]<-index2[-j]
    my_mat<-rbind(my_mat, my_mat2)
  }
  my_mat<-my_mat[-1,]  
  #how many interactions are possible
  possible<-nrow(my_mat)
  #Put the coordinates and the potential interactions in a single alrge dataframe. 
  m3 <- rbind(my_mat, coords)
  #Everything that has a match  must do so as they are present in the 1% and in my list of potential interactions
  matches<-nrow(m3[duplicated(m3), , drop = FALSE])
  match_coords<-m3[duplicated(m3), , drop = FALSE]
  
  
  my_temp_mat<-matrix(0,1,6)
  #Output numbers
  my_temp_mat[1,1]<-possible
  my_temp_mat[1,2]<-matches
  my_temp_mat[1,3]<-index[i]
  #Proportional amount possible realtive to whole genome
  my_temp_mat[1,4]<-possible/1430416
  #Proportioanl amount of matches we see in the top 1% 
  my_temp_mat[1,5]<-matches/14304
  #Enrichment in top 1%
  my_temp_mat[1,6]<-as.numeric(my_temp_mat[1,5])/as.numeric(my_temp_mat[1,4])
  my_final_mat<-rbind(my_final_mat,my_temp_mat)
  my_matches<-rbind(my_matches,match_coords )
  print(i)
}
my_matches_final<-unique(my_matches[-1,])
colnames(my_final_mat)<-c('Genome wide','In top 1%','tRNA group','Genome wide proportion', '1% proportion', 'Enrichment')
output<-list(my_final_mat[-1,],my_matches_final)
return(output)
}

wtrapa_norm_tRNA_topone<-tRNAsintop(wtrapa_norm_coords)

#Isolate the top 1% - con rapa
conrapa_norm_cutoff<-(round(length(conrapa_norm)/100))
conrapa_norm_oneperc<-conrapa_norm[order(conrapa_norm, decreasing=TRUE)][cutoff]
conrapa_norm_coords<-which(conrapa_norm>oneperc, arr.ind=TRUE)
conrapa_norm_tRNA_topone<-tRNAsintop(conrapa_norm_coords)

#Isolate the top 1% - wt rapa decay norm
wtrapa_decaynorm_norm_cutoff<-(round(length(wtrapa_decaynorm_norm)/100))
wtrapa_decaynorm_norm_oneperc<-wtrapa_decaynorm_norm[order(wtrapa_decaynorm_norm, decreasing=TRUE)][cutoff]
wtrapa_decaynorm_norm_coords<-which(wtrapa_decaynorm_norm>oneperc, arr.ind=TRUE)
wtrapa_decaynorm_norm_tRNA_topone<-tRNAsintop(wtrapa_decaynorm_norm_coords)

#Isolate the top 1% - con rapa decay norm
conrapa_decaynorm_norm_cutoff<-(round(length(conrapa_decaynorm_norm)/100))
conrapa_decaynorm_norm_oneperc<-conrapa_decaynorm_norm[order(conrapa_decaynorm_norm, decreasing=TRUE)][cutoff]
conrapa_decaynorm_norm_coords<-which(conrapa_decaynorm_norm>oneperc, arr.ind=TRUE)
conrapa_decaynorm_norm_tRNA_topone<-tRNAsintop(conrapa_decaynorm_norm_coords)








#Extract chromosome starts. This is done by assessing how many bins in each chromosome. Then incrmentally assessing at which cordinate the chromomse starts in the interaction matrix
chr.size<-c(203893,794508,342718,1490682,602514,284456,1067526,544538,435585,719294,687260,1008248,908607,812465, 1054033,921188)
res<-10000
chr_bins<-list("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13"
               , "chr14", "chr15", "chr16")
bin_per_chr<-floor(chr.size/res+1)
incrementalbins<-matrix(0,length(chr_bins),1)
for (i in 1:(length(chr_bins)-1)){
  incrementalbins[i+1,]<-incrementalbins[i,]+bin_per_chr[i]
}
incrementalbins<-rbind(incrementalbins, 1196)

#Are the trNA interactions inter or intrachromosmsal that we are losing?
inter_intra<-matrix(0,nrow(conrapa_decaynorm_norm_tRNA_topone[[2]]),1)
for (i in 1:nrow(conrapa_decaynorm_norm_tRNA_topone[[2]])){
left<-max(which(incrementalbins<conrapa_decaynorm_norm_tRNA_topone[[2]][i,1]))
right<-max(which(incrementalbins<conrapa_decaynorm_norm_tRNA_topone[[2]][i,2]))
if (left==right){inter_intra[i,]<-'intra'
}else{
  inter_intra[i,]<-'inter'
  }
}

conrapa_decaynorm_norm_tRNA_topone[[2]]<-cbind(conrapa_decaynorm_norm_tRNA_topone[[2]],inter_intra)


inter_intra<-matrix(0,nrow(wtrapa_decaynorm_norm_tRNA_topone[[2]]),1)
for (i in 1:nrow(wtrapa_decaynorm_norm_tRNA_topone[[2]])){
  left<-max(which(incrementalbins<wtrapa_decaynorm_norm_tRNA_topone[[2]][i,1]))
  right<-max(which(incrementalbins<wtrapa_decaynorm_norm_tRNA_topone[[2]][i,2]))
  if (left==right){inter_intra[i,]<-'intra'
  }else{
    inter_intra[i,]<-'inter'
  }
}

wtrapa_decaynorm_norm_tRNA_topone[[2]]<-cbind(wtrapa_decaynorm_norm_tRNA_topone[[2]],inter_intra)


inter_intra<-matrix(0,nrow(conrapa_norm_tRNA_topone[[2]]),1)
for (i in 1:nrow(conrapa_norm_tRNA_topone[[2]])){
  left<-max(which(incrementalbins<conrapa_norm_tRNA_topone[[2]][i,1]))
  right<-max(which(incrementalbins<conrapa_norm_tRNA_topone[[2]][i,2]))
  if (left==right){inter_intra[i,]<-'intra'
  }else{
    inter_intra[i,]<-'inter'
  }
}

conrapa_norm_tRNA_topone[[2]]<-cbind(conrapa_norm_tRNA_topone[[2]],inter_intra)


inter_intra<-matrix(0,nrow(wtrapa_norm_tRNA_topone[[2]]),1)
for (i in 1:nrow(wtrapa_norm_tRNA_topone[[2]])){
  left<-max(which(incrementalbins<wtrapa_norm_tRNA_topone[[2]][i,1]))
  right<-max(which(incrementalbins<wtrapa_norm_tRNA_topone[[2]][i,2]))
  if (left==right){inter_intra[i,]<-'intra'
  }else{
    inter_intra[i,]<-'inter'
  }
}

wtrapa_norm_tRNA_topone[[2]]<-cbind(wtrapa_norm_tRNA_topone[[2]],inter_intra)

#What category of change am I finding in my tRNAs

#Subset the data to interactions that are found in both datasets
#Creates an overlapping matrix for each point and sees if they are present in both. 
#Then apply a which staemtn to see if conssitencies derive from same rows. 
#The indices generated can then be used to index the tRNA list positively to extract tRNAs that maintain interaction, or negatively to find those that lose intraction
maintained_ind<-unique(which( outer( data.frame(wtrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[,1], data.frame(conrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[,1], "==") & 
       outer(data.frame(wtrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[,2], data.frame(conrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[,2], "=="), 
       arr.ind=TRUE ))

lost<-data.frame(wtrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[! 1:NROW(wtrapa_norm_tRNA_topone[[2]]) %in% maintained_ind, ]

maintained<-data.frame(wtrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[1:NROW(wtrapa_norm_tRNA_topone[[2]]) %in% maintained_ind, ]

#Do the same agin but start with the condensin perturbed tRNA dataset. Now that it is reversed, finding those that change between datsets, give us a gained interactions dataset.
gained_ind<-unique(which( outer( data.frame(conrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[,1], data.frame(wtrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[,1], "==") & 
                            outer(data.frame(conrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[,2], data.frame(wtrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[,2], "=="), 
                          arr.ind=TRUE ))

gained<-data.frame(wtrapa_norm_tRNA_topone[[2]], stringsAsFactors = F)[! 1:NROW(conrapa_norm_tRNA_topone[[2]]) %in% gained_ind, ]

#What are these tRNAS? I want a GFF file for each group set. 

gained_details<-cbind(tRNA_gff[match(as.numeric(gained[,1]),tRNA_gff[,12]),c(1,4,5,10,11,12)],tRNA_gff[match(as.numeric(gained[,2]),tRNA_gff[,12]),c(1,4,5,10,11,12)],gained[,3])
maintained_details<-cbind(tRNA_gff[match(as.numeric(maintained[,1]),tRNA_gff[,12]),c(1,4,5,10,11,12)],tRNA_gff[match(as.numeric(maintained[,2]),tRNA_gff[,12]),c(1,4,5,10,11,12)],maintained[,3])
lost_details<-cbind(tRNA_gff[match(as.numeric(lost[,1]),tRNA_gff[,12]),c(1,4,5,10,11,12)],tRNA_gff[match(as.numeric(lost[,2]),tRNA_gff[,12]),c(1,4,5,10,11,12)],lost[,3])


l_er<-length(which(lost_details[,13]=='inter'))
l_ra<-length(which(lost_details[,13]=='intra'))
              


gained_list<-c(gained_details[,6], gained_details[,12])
maintained_list<-c(maintained_details[,6], maintained_details[,12])
lost_list<-c(lost_details[,6], lost_details[,12])

main_gain_pre<-unique(maintained_list[match(gained_list,maintained_list)])
main_gain<-unique(main_gain_pre[-which(is.na(main_gain_pre)==TRUE)])

main_lost_pre<-unique(maintained_list[match(lost_list,maintained_list)])
main_lost<-unique(main_lost_pre[-which(is.na(main_lost_pre)==TRUE)])

gain_lost_pre<-unique(gained_list[match(lost_list,gained_list)])
gain_lost<-unique(gain_lost_pre[-which(is.na(gain_lost_pre)==TRUE)])

main_gain_lost_pre<-unique(maintained_list[match(gain_lost,maintained_list)])
main_gain_lost<-unique(main_gain_lost_pre[-which(is.na(main_gain_lost_pre)==TRUE)])



lost_only_pre<-unique(match(main_gain,lost_list))
lost_only_pre<-lost_only_pre[-which(is.na(lost_only_pre)==TRUE)]
lost_only<-unique(lost_list[-lost_only_pre])

maintained_only_pre<-unique(match(gain_lost,maintained_list))
maintained_only_pre<-maintained_only_pre[-which(is.na(maintained_only_pre)==TRUE)]
maintained_only<-unique(maintained_list[-maintained_only_pre])

gained_only_pre<-unique(match(main_lost,gained_list))
gained_only_pre<-gained_only_pre[-which(is.na(gained_only_pre)==TRUE)]
gained_only<-unique(gained_list[-gained_only_pre])


#Lets read in poliii chip data and adapt it to make it relvant

#read in the file
poliii_chip<-read.csv('~/Google_Drive/Lab/Data/Pol_III_Supplementary_Data_.csv',header = TRUE,as.is= TRUE, stringsAsFactors = F)

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

#poliii_trna_score<-(cbind(poliii_trna_score,matrix('NA',nrow(poliii_trna_score),1)))
poliii_trna_score<-(cbind(poliii_trna_score,(tRNA_gff[match(poliii_trna_score[,1], tRNA_gff[,14]),12])))

poliii_trna_score<-cbind(poliii_trna_score,(tRNA_gff[match(poliii_trna_score[,1], tRNA_gff[,9]),]))


head(poliii_trna_score)

pdf('~/Google_Drive/Lab/Data_Analysis/tRNA/polIII_at_trna.pdf')
par(mfrow=c(1,1))
par(mar=c(8.1,4.1,4.1,2.1))
boxplot(unlist(poliii_trna_score[match(lost_only,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(main_lost,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(maintained_only,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(main_gain,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(gained_only,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(gain_lost,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(main_gain_lost,poliii_trna_score[,3]),3]), ylab='Median ChIP enrichment',xaxt='n', main='Median enrichment of PolIII sub-units in groups of tRNAs \n that change clustering in response to condensin depletion', outline=F)
axis(1,at=c(1:7), label=c('Lost','Main + Lost','Main','Main + Gain','Gain','Gain + Lost', 'Main + Gain + Lost'), las=2)
dev.off()

all_clust_trnas<-unique(c(lost_only,gained_only,maintained_only,main_gain,main_gain_lost,main_lost,gain_lost))
all_trnas<-tRNA_gff[,12]

pdf('~/Google_Drive/Lab/Data_Analysis/tRNA/polIII_at_trna_subset_lost_only.pdf')
par(mfrow=c(1,1))
par(mar=c(10.1,4.1,4.1,2.1))
boxplot(unlist(poliii_trna_score[match(all_trnas,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(all_clust_trnas,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(lost_only,poliii_trna_score[,3]),3]), ylab='Median ChIP enrichment',xaxt='n', main='Median enrichment of PolIII sub-units in groups of tRNAs \n that change clustering in response to condensin depletion', outline=F)
axis(1,at=c(1:3), label=c('All tRNAS','Clustered tRNAs','Condensin-dependent \n clustered tRNAs'), las=2)
dev.off()



a<-unlist(poliii_trna_score[match(all_trnas,poliii_trna_score[,3]),3])
b<-unlist(poliii_trna_score[match(all_clust_trnas,poliii_trna_score[,3]),3])
c<-unlist(poliii_trna_score[match(lost_only,poliii_trna_score[,3]),3])

pdf('~/Google_Drive/Lab/Data_Analysis/tRNA/polIII_at_trna_subset_lost_only_jitterplot.pdf')
plot(jitter(rep(0,length(a)),amount=0.2),a,
     xlim=c(-0.5,2.5), ylim=range(0,1200),
     xaxt='n',xlab='', ylab='Median ChIP enrichment', main='Median enrichment of PolIII sub-units in groups of tRNAs \n that change clustering in response to condensin depletion', col=1, pch=8)
points(jitter(rep(1,length(b)), amount=0.2), b, col=2, pch=8)
points(jitter(rep(2,length(c)), amount=0.2), c, col=3, pch=8)
axis(1,at=c(0:2), label=c('All tRNAS','Clustered tRNAs','Condensin-dependent \n clustered tRNAs'), las=2)
dev.off()


a_median<-median(a)
b_median<-median(b)
c_median<-median(c)

pdf('~/Google_Drive/Lab/Data_Analysis/tRNA/polIII_at_trna_subset_lost_only_medians.pdf')
plot(c(a_median,b_median,c_median),xlim=c(-0.5,3.5), ylim=range(0,1200))
dev.off()




wilcox.test(a,b)
wilcox.test(b,c)
wilcox.test(a,c)


capture.output(summary(poliii_trna_score_indexed_trnas),file='~/Desktop/poliii_trna.gff')

lapply(poliii_trna_score_indexed_trnas, function(x) write.table( data.frame(x), '~/Desktop/poliii_trna.gff'  , append= T, sep=',' ))

#What is overlap between SMC4-chip signal and PolIII
smc4_chip<-read.table('~/Google_Drive/Lab/Paper/chip_stuff/ChIP_data-selected/AH6408K-031517-SK1Yue-B3W3-MACS2/AH6408K-031517-SK1Yue-PM_B3W3_MACS2_FE.bdg', stringsAsFactors = F, header = F)

poliii_trna_score_indexed<-cbind(poliii_trna_score,(tRNA_gff[match(poliii_trna_score[,3],tRNA_gff[,12]),c(1,4,5)]))
poliii_trna_score_indexed_trnas<-poliii_trna_score_indexed[which(poliii_trna_score_indexed[,3]>0),]

x <- vapply(output$Title, length, 1L)          ## How many items per list element
output <- output[rep(rownames(output), x), ]   ## Expand the data frame
output$Title <- unlist(output$Title, use.names = FALSE)  ## Replace with raw values


if (!require("bedr")) {
  install.packages("bedr")
  library(bedr)
}

poliii_trna.gff<-read.table('/scratch/mrp420/poliii_trna.gff',stringsAsFactors = F)
b.sorted <- bedr(engine = "bedtools", input = list(i = poliii_trna.gff[,c(4:6)]), method = "sort", params = "", check.chr = FALSE)

poliii_trna.gff


a.sorted <- bedr(engine = "bedtools", input = list(i = smc4_chip[,1:3]), method = "sort", params = "", check.chr = FALSE)
b.sorted <- bedr(engine = "bedtools", input = list(i = poliii_trna_score_indexed_trnas[,c(4:6)]), method = "sort", params = "", check.chr = FALSE)
overlap<-bedr(engine = "bedtools",method='intersect',input = list(a =a.sorted,  b=b.sorted), params = "-loj", check.chr = FALSE)


nooverlap<-overlap[which(overlap[,4]=='.'),]
overlap<-overlap[-which(overlap[,4]=='.'),]


#Take the duplicates out. then run down and reassign smc4 to the matrix. then can do the comaprison.  

overlap_with_chip<-cbind(overlap[1:2,],smc4_chip[match(overlap[1:2,1:3],smc4_chip[,1:3]),4],)




smc4_col<-matrix(0,nrow(poliii_trna_score_indexed_trnas),1)




for (i in 1:nrow(poliii_trna_score_indexed_trnas)){
  tRNA_GFF<-
    
    for loop this 
  find it matchs
  then get mean(poliii_trna_score_indexed_trnas)
  
  

head(poliii_trna_score)

pdf('~/Google_Drive/Lab/Data_Analysis/tRNA/polIII_at_trna.pdf')
par(mfrow=c(1,1))
par(mar=c(8.1,4.1,4.1,2.1))
boxplot(unlist(poliii_trna_score[match(lost_only,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(main_lost,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(maintained_only,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(main_gain,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(gained_only,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(gain_lost,poliii_trna_score[,3]),3]),unlist(poliii_trna_score[match(main_gain_lost,poliii_trna_score[,3]),3]), ylab='Median ChIP enrichment',xaxt='n', main='Median enrichment of PolIII sub-units in groups of tRNAs \n that change clustering in response to condensin depletion', outline=F)
axis(1,at=c(1:7), label=c('Lost','Main + Lost','Main','Main + Gain','Gain','Gain + Lost', 'Main + Gain + Lost'), las=2)
dev.off()

poliii_trna_score_indexed<-cbind(poliii_trna_score,(tRNA_gff[match(poliii_trna_score[,3],tRNA_gff[,12]),c(1,4,5)]))


#Forkhead origins
fh_ori<-read.csv('~/Google_Drive/Lab/Data/forkheadandorigins.csv', stringsAsFactors = F)

name_joiner<-function(name_inputs){
  paste0(name_inputs[1],'_',name_inputs[2])
}

final_names<-apply(fh_ori[,1:2], 1, name_joiner)

fh_ori_bed<-cbind(fh_ori[,3:5],final_names)

write.table(fh_ori_bed,'~/Google_Drive/Lab/Data/forkheadandorigins_saccer3.bed',sep='\t',quote=FALSE,row.names=F)

#Frokehead was converted to sk1 coordinates

#Load viridis color palette for graphing
if (!require("bedr")) {
  install.packages("bedr")
  library(viridis)
}

fh_ori_bed_sk1<-read.table('~/Google_Drive/Lab/Data/forkheadandorigins_sk1.bed', stringsAsFactors = F)


convert_roman<-function(input){
  paste0('chr',as.roman(strsplit(input[1],'chr')[[1]][2]))
}

fh_ori_bed_sk1[,1]<-apply(fh_ori_bed_sk1,1,convert_roman)
unique(fh_ori_bed_sk1[,1])->chrs

a.sorted <- bedr(engine = "bedtools", input = list(i = fh_ori_bed_sk1[,1:3]), method = "sort", params = "", check.chr = FALSE)
b.sorted <- bedr(engine = "bedtools", input = list(i = tRNA_gff[,c(1,4,5)]), method = "sort", params = "", check.chr = FALSE)
overlap<-bedr(engine = "bedtools",method='intersect',input = list(a =a.sorted,  b=b.sorted), params = "-loj", check.chr = FALSE)
nooverlap<-overlap[which(overlap[,4]=='.'),]
overlap<-overlap[-which(overlap[,4]=='.'),]

fh_ori_bed_sk1_nooverlap_1kb<-nooverlap[,1:3]
row.names(fh_ori_bed_sk1_nooverlap_1kb)<-c()
fh_ori_bed_sk1_nooverlap_1kb[,2]<-nooverlap[,2]-1000
fh_ori_bed_sk1_nooverlap_1kb[,3]<-nooverlap[,3]+1000
fh_ori_bed_sk1_nooverlap_1kb[which(fh_ori_bed_sk1_nooverlap_1kb[,2]<0),2]<-0
fh_ori_bed_sk1_nooverlap_1kb[which(fh_ori_bed_sk1_nooverlap_1kb[,3]<0),3]<-0

a.sorted <- bedr(engine = "bedtools", input = list(i = fh_ori_bed_sk1_nooverlap_1kb), method = "sort", params = "", check.chr = FALSE)
overlap_1kb<-bedr(engine = "bedtools",method='intersect',input = list(a =a.sorted,  b=b.sorted), params = "-loj", check.chr = FALSE)
nooverlap_1kb<-overlap_1kb[which(overlap_1kb[,4]=='.'),]
overlap_1kb<-overlap_1kb[-which(overlap_1kb[,4]=='.'),]


fh_ori_bed_sk1_nooverlap_10kb<-nooverlap[,1:3]
row.names(fh_ori_bed_sk1_nooverlap_10kb)<-c()
fh_ori_bed_sk1_nooverlap_10kb[,2]<-nooverlap[,2]-10000
fh_ori_bed_sk1_nooverlap_10kb[,3]<-nooverlap[,3]+10000
fh_ori_bed_sk1_nooverlap_10kb[which(fh_ori_bed_sk1_nooverlap_10kb[,2]<0),2]<-0
fh_ori_bed_sk1_nooverlap_10kb[which(fh_ori_bed_sk1_nooverlap_10kb[,3]<0),3]<-0

a.sorted <- bedr(engine = "bedtools", input = list(i = fh_ori_bed_sk1_nooverlap_10kb), method = "sort", params = "", check.chr = FALSE)
overlap_10kb<-bedr(engine = "bedtools",method='intersect',input = list(a =a.sorted,  b=b.sorted), params = "-loj", check.chr = FALSE)
nooverlap_10kb<-overlap_10kb[which(overlap_10kb[,4]=='.'),]
overlap_10kb<-overlap_10kb[-which(overlap_10kb[,4]=='.'),]

#tRNAs
bins_overlap<-tRNA_gff[match(overlap[,5],tRNA_gff[,c(4)]),12]
overlap_counts<-c(length(intersect(unique(lost_only),bins_overlap)),length(intersect(unique(main_lost),bins_overlap)),length(intersect(unique(maintained_only),bins_overlap)),length(intersect(unique(main_gain),bins_overlap)),length(intersect(unique(gained_only),bins_overlap)),length(intersect(unique(gain_lost),bins_overlap)),length(intersect(unique(main_gain_lost),bins_overlap)))
overlap_counts<-(overlap_counts/sum(overlap_counts))

bins_overlap_1kb<-tRNA_gff[match(overlap_1kb[,5],tRNA_gff[,c(4)]),12]
overlap_1kb_counts<-c(length(intersect(unique(lost_only),bins_overlap_1kb)),length(intersect(unique(main_lost),bins_overlap_1kb)),length(intersect(unique(maintained_only),bins_overlap_1kb)),length(intersect(unique(main_gain),bins_overlap_1kb)),length(intersect(unique(gained_only),bins_overlap_1kb)),length(intersect(unique(gain_lost),bins_overlap_1kb)),length(intersect(unique(main_gain_lost),bins_overlap_1kb)))
overlap_1kb_counts<-(overlap_1kb_counts/sum(overlap_1kb_counts))

bins_overlap_10kb<-tRNA_gff[match(overlap_10kb[,5],tRNA_gff[,c(4)]),12]
overlap_10kb_counts<-c(length(intersect(unique(lost_only),bins_overlap_10kb)),length(intersect(unique(main_lost),bins_overlap_10kb)),length(intersect(unique(maintained_only),bins_overlap_10kb)),length(intersect(unique(main_gain),bins_overlap_10kb)),length(intersect(unique(gained_only),bins_overlap_10kb)),length(intersect(unique(gain_lost),bins_overlap_10kb)),length(intersect(unique(main_gain_lost),bins_overlap_10kb)))
overlap_10kb_counts<-(overlap_10kb_counts/sum(overlap_10kb_counts))

bins_nooverlap_10kb<-tRNA_gff[-match(overlap_10kb[,5],tRNA_gff[,c(4)]),12]
nooverlap_10kb_counts<-c(length(intersect(unique(lost_only),bins_nooverlap_10kb)),length(intersect(unique(main_lost),bins_nooverlap_10kb)),length(intersect(unique(maintained_only),bins_nooverlap_10kb)),length(intersect(unique(main_gain),bins_nooverlap_10kb)),length(intersect(unique(gained_only),bins_nooverlap_10kb)),length(intersect(unique(gain_lost),bins_nooverlap_10kb)),length(intersect(unique(main_gain_lost),bins_nooverlap_10kb)))
nooverlap_10kb_counts<-(nooverlap_10kb_counts/sum(nooverlap_10kb_counts))



if (!require("wesanderson")) {
  install.packages("wesanderson")
  library(wesanderson)
}


col1<-c(wes_palette("GrandBudapest"),wes_palette("GrandBudapest2"))
pdf('~/Google_Drive/Lab/Data_Analysis/tRNA/origin_trna_adjacency.pdf')
par(mar=c(5,5,5,15),xpd = T)
barplot(data.matrix(cbind(overlap_counts,overlap_1kb_counts,overlap_10kb_counts,nooverlap_10kb_counts)), xaxt='n', ylab='% of tRNAs', main='Proportion of tRNAs that change their clustering, \n relative to their adjacency to origins', xlab='Distance between tRNA and adjacent origin', col=col1)
legend("topright",col = col1,cex = 0.8,lwd = 8, lty = 1,legend= c('Lost','Lost/Maintained', 'Maintained','Maintained/Gained','Gained','Lost/Gained','Lost/Maintained/Gained' ),inset =c(-0.8,0))
axis(1, at=c(0.7,1.9,3.1,4.3), label=c('Overlap', '<1kb', '<10kb', '>10kb'))
dev.off()



#Forkhead 

ARS_name<-function(table){
  as.numeric(strsplit(table[4],'[_:-]')[[1]][5])
}
trna_loc<-apply(fh_ori_bed_sk1,1,ARS_name)

fh_ori_bed_sk1<-cbind(fh_ori_bed_sk1,fh_ori[match(trna_loc,fh_ori[,4]),c(6,7,8,9)])

head(fh_ori_bed_sk1)

#The change in activity of each origin in each of the mutant strains (relative to WT) is indicated as −1 (decreased), 0 (unchanged), or 1 (increased). Names and coordinates were taken from OriDB (Nieduszynski et al., 2007).

fh_neg<-which(apply(fh_ori_bed_sk1[7:10],1,sum)<0)
fh_zero<-which(apply(fh_ori_bed_sk1[7:10],1,sum)==0)
fh_pos<-which(apply(fh_ori_bed_sk1[7:10],1,sum)>0)


fh_neg_overlap<-overlap[match(fh_ori_bed_sk1[fh_neg,2],overlap[,2]),]
fh_neg_overlap<-fh_neg_overlap[which(fh_neg_overlap[,2]>0),]
fh_neg_bins_overlap<-tRNA_gff[match(fh_neg_overlap[,5],tRNA_gff[,c(4)]),12]
fh_neg_overlap_counts<-c(length(intersect(unique(lost_only),fh_neg_bins_overlap)),length(intersect(unique(main_lost),fh_neg_bins_overlap)),length(intersect(unique(maintained_only),fh_neg_bins_overlap)),length(intersect(unique(main_gain),fh_neg_bins_overlap)),length(intersect(unique(gained_only),fh_neg_bins_overlap)),length(intersect(unique(gain_lost),fh_neg_bins_overlap)),length(intersect(unique(main_gain_lost),fh_neg_bins_overlap)))
fh_neg_overlap_counts<-(fh_neg_overlap_counts/sum(fh_neg_overlap_counts))


fh_neg_overlap_1kb<-overlap_1kb[match(fh_ori_bed_sk1[fh_neg,2],overlap_1kb[,2]+1000),]
fh_neg_overlap_1kb<-fh_neg_overlap_1kb[which(fh_neg_overlap_1kb[,2]>0),]
fh_neg_bins_overlap_1kb<-tRNA_gff[match(fh_neg_overlap_1kb[,5],tRNA_gff[,c(4)]),12]
fh_neg_overlap_1kb_counts<-c(length(intersect(unique(lost_only),fh_neg_bins_overlap_1kb)),length(intersect(unique(main_lost),fh_neg_bins_overlap_1kb)),length(intersect(unique(maintained_only),fh_neg_bins_overlap_1kb)),length(intersect(unique(main_gain),fh_neg_bins_overlap_1kb)),length(intersect(unique(gained_only),fh_neg_bins_overlap_1kb)),length(intersect(unique(gain_lost),fh_neg_bins_overlap_1kb)),length(intersect(unique(main_gain_lost),fh_neg_bins_overlap_1kb)))
fh_neg_overlap_1kb_counts<-(fh_neg_overlap_1kb_counts/sum(fh_neg_overlap_1kb_counts))

fh_neg_overlap_10kb<-overlap_10kb[match(fh_ori_bed_sk1[fh_neg,2],overlap_10kb[,2]+10000),]
fh_neg_overlap_10kb<-fh_neg_overlap_10kb[which(fh_neg_overlap_10kb[,2]>0),]
fh_neg_bins_overlap_10kb<-tRNA_gff[match(fh_neg_overlap_10kb[,5],tRNA_gff[,c(4)]),12]
fh_neg_overlap_10kb_counts<-c(length(intersect(unique(lost_only),fh_neg_bins_overlap_10kb)),length(intersect(unique(main_lost),fh_neg_bins_overlap_10kb)),length(intersect(unique(maintained_only),fh_neg_bins_overlap_10kb)),length(intersect(unique(main_gain),fh_neg_bins_overlap_10kb)),length(intersect(unique(gained_only),fh_neg_bins_overlap_10kb)),length(intersect(unique(gain_lost),fh_neg_bins_overlap_10kb)),length(intersect(unique(main_gain_lost),fh_neg_bins_overlap_10kb)))
fh_neg_overlap_10kb_counts<-(fh_neg_overlap_10kb_counts/sum(fh_neg_overlap_10kb_counts))

fh_neg_nooverlap_10kb<-overlap_10kb[match(fh_ori_bed_sk1[fh_neg,2],overlap_10kb[,2]+10000),]
fh_neg_nooverlap_10kb<-fh_neg_nooverlap_10kb[which(fh_neg_nooverlap_10kb[,2]>0),]
fh_neg_bins_nooverlap_10kb<-tRNA_gff[match(fh_neg_nooverlap_10kb[,5],tRNA_gff[,c(4)]),12]
fh_neg_nooverlap_10kb_counts<-c(length(intersect(unique(lost_only),fh_neg_bins_nooverlap_10kb)),length(intersect(unique(main_lost),fh_neg_bins_nooverlap_10kb)),length(intersect(unique(maintained_only),fh_neg_bins_nooverlap_10kb)),length(intersect(unique(main_gain),fh_neg_bins_nooverlap_10kb)),length(intersect(unique(gained_only),fh_neg_bins_nooverlap_10kb)),length(intersect(unique(gain_lost),fh_neg_bins_nooverlap_10kb)),length(intersect(unique(main_gain_lost),fh_neg_bins_nooverlap_10kb)))
fh_neg_nooverlap_10kb_counts<-(fh_neg_nooverlap_10kb_counts/sum(fh_neg_nooverlap_10kb_counts))

barplot(data.matrix(cbind(fh_neg_overlap_counts,fh_neg_overlap_1kb_counts,fh_neg_overlap_10kb_counts,fh_neg_nooverlap_10kb_counts)))

pdf('~/Google_Drive/Lab/Data_Analysis/tRNA/origin_trna_adjacency_forkhead_promotes.pdf')
par(mar=c(5,5,5,15),xpd = T)
barplot(data.matrix(cbind(fh_neg_overlap_counts[c(1,3,5)],fh_neg_overlap_1kb_counts[c(1,3,5)],fh_neg_overlap_10kb_counts[c(1,3,5)],fh_neg_nooverlap_10kb_counts[c(1,3,5)])), xaxt='n', ylab='% of tRNAs', main='Proportion of tRNAs that change their clustering, \n relative to their adjacency to Forkhead promoted origins', xlab='Distance between tRNA and adjacent origin that Forkhead promotes', col=col1)
legend("topright",col = col1,cex = 0.8,lwd = 8, lty = 1,legend= c('Lost','Lost/Maintained', 'Maintained','Maintained/Gained','Gained','Lost/Gained','Lost/Maintained/Gained' ),inset =c(-0.8,0))
axis(1, at=c(0.7,1.9,3.1,4.3), label=c('Overlap', '<1kb', '<10kb', '>10kb'))
dev.off()
barplot(data.matrix(cbind(fh_neg_overlap_counts[c(1,3,5)]*(1/sum(fh_neg_overlap_counts[c(1,3,5)])),fh_neg_overlap_1kb_counts[c(1,3,5)]*(1/sum(fh_neg_overlap_1kb_counts[c(1,3,5)])),fh_neg_overlap_10kb_counts[c(1,3,5)]*(1/sum(fh_neg_overlap_10kb_counts[c(1,3,5)])),fh_neg_nooverlap_10kb_counts[c(1,3,5)]*(1/sum(fh_neg_nooverlap_10kb_counts[c(1,3,5)])))))
barplot(data.matrix(cbind(fh_pos_overlap_counts[c(1,3,5)]*(1/sum(fh_pos_overlap_counts[c(1,3,5)])),fh_pos_overlap_1kb_counts[c(1,3,5)]*(1/sum(fh_pos_overlap_1kb_counts[c(1,3,5)])),fh_pos_overlap_10kb_counts[c(1,3,5)]*(1/sum(fh_pos_overlap_10kb_counts[c(1,3,5)])),fh_pos_nooverlap_10kb_counts[c(1,3,5)]*(1/sum(fh_pos_nooverlap_10kb_counts[c(1,3,5)])))))

fh_pos_overlap<-overlap[match(fh_ori_bed_sk1[fh_pos,2],overlap[,2]),]
fh_pos_overlap<-fh_pos_overlap[which(fh_pos_overlap[,2]>0),]
fh_pos_bins_overlap<-tRNA_gff[match(fh_pos_overlap[,5],tRNA_gff[,c(4)]),12]
fh_pos_overlap_counts<-c(length(intersect(unique(lost_only),fh_pos_bins_overlap)),length(intersect(unique(main_lost),fh_pos_bins_overlap)),length(intersect(unique(maintained_only),fh_pos_bins_overlap)),length(intersect(unique(main_gain),fh_pos_bins_overlap)),length(intersect(unique(gained_only),fh_pos_bins_overlap)),length(intersect(unique(gain_lost),fh_pos_bins_overlap)),length(intersect(unique(main_gain_lost),fh_pos_bins_overlap)))
fh_pos_overlap_counts<-(fh_pos_overlap_counts/sum(fh_pos_overlap_counts))


fh_pos_overlap_1kb<-overlap_1kb[match(fh_ori_bed_sk1[fh_pos,2],overlap_1kb[,2]+1000),]
fh_pos_overlap_1kb<-fh_pos_overlap_1kb[which(fh_pos_overlap_1kb[,2]>0),]
fh_pos_bins_overlap_1kb<-tRNA_gff[match(fh_pos_overlap_1kb[,5],tRNA_gff[,c(4)]),12]
fh_pos_overlap_1kb_counts<-c(length(intersect(unique(lost_only),fh_pos_bins_overlap_1kb)),length(intersect(unique(main_lost),fh_pos_bins_overlap_1kb)),length(intersect(unique(maintained_only),fh_pos_bins_overlap_1kb)),length(intersect(unique(main_gain),fh_pos_bins_overlap_1kb)),length(intersect(unique(gained_only),fh_pos_bins_overlap_1kb)),length(intersect(unique(gain_lost),fh_pos_bins_overlap_1kb)),length(intersect(unique(main_gain_lost),fh_pos_bins_overlap_1kb)))
fh_pos_overlap_1kb_counts<-(fh_pos_overlap_1kb_counts/sum(fh_pos_overlap_1kb_counts))

fh_pos_overlap_10kb<-overlap_10kb[match(fh_ori_bed_sk1[fh_pos,2],overlap_10kb[,2]+10000),]
fh_pos_overlap_10kb<-fh_pos_overlap_10kb[which(fh_pos_overlap_10kb[,2]>0),]
fh_pos_bins_overlap_10kb<-tRNA_gff[match(fh_pos_overlap_10kb[,5],tRNA_gff[,c(4)]),12]
fh_pos_overlap_10kb_counts<-c(length(intersect(unique(lost_only),fh_pos_bins_overlap_10kb)),length(intersect(unique(main_lost),fh_pos_bins_overlap_10kb)),length(intersect(unique(maintained_only),fh_pos_bins_overlap_10kb)),length(intersect(unique(main_gain),fh_pos_bins_overlap_10kb)),length(intersect(unique(gained_only),fh_pos_bins_overlap_10kb)),length(intersect(unique(gain_lost),fh_pos_bins_overlap_10kb)),length(intersect(unique(main_gain_lost),fh_pos_bins_overlap_10kb)))
fh_pos_overlap_10kb_counts<-(fh_pos_overlap_10kb_counts/sum(fh_pos_overlap_10kb_counts))

fh_pos_nooverlap_10kb<-overlap_10kb[match(fh_ori_bed_sk1[fh_pos,2],overlap_10kb[,2]+10000),]
fh_pos_nooverlap_10kb<-fh_pos_nooverlap_10kb[which(fh_pos_nooverlap_10kb[,2]>0),]
fh_pos_bins_nooverlap_10kb<-tRNA_gff[match(fh_pos_nooverlap_10kb[,5],tRNA_gff[,c(4)]),12]
fh_pos_nooverlap_10kb_counts<-c(length(intersect(unique(lost_only),fh_pos_bins_nooverlap_10kb)),length(intersect(unique(main_lost),fh_pos_bins_nooverlap_10kb)),length(intersect(unique(maintained_only),fh_pos_bins_nooverlap_10kb)),length(intersect(unique(main_gain),fh_pos_bins_nooverlap_10kb)),length(intersect(unique(gained_only),fh_pos_bins_nooverlap_10kb)),length(intersect(unique(gain_lost),fh_pos_bins_nooverlap_10kb)),length(intersect(unique(main_gain_lost),fh_pos_bins_nooverlap_10kb)))
fh_pos_nooverlap_10kb_counts<-(fh_pos_nooverlap_10kb_counts/sum(fh_pos_nooverlap_10kb_counts))

pdf('~/Google_Drive/Lab/Data_Analysis/tRNA/origin_trna_adjacency_forkhead_repress.pdf')
par(mar=c(5,5,5,15),xpd = T)
barplot(data.matrix(cbind(fh_pos_overlap_counts,fh_pos_overlap_1kb_counts,fh_pos_overlap_10kb_counts,fh_pos_nooverlap_10kb_counts)), xaxt='n', ylab='% of tRNAs', main='Proportion of tRNAs that change their clustering, \n relative to their adjacency to Forkhead repressed origins', xlab='Distance between tRNA and adjacent origin that Forkhead represses', col=col1)
legend("topright",col = col1,cex = 0.8,lwd = 8, lty = 1,legend= c('Lost','Lost/Maintained', 'Maintained','Maintained/Gained','Gained','Lost/Gained','Lost/Maintained/Gained' ),inset =c(-0.8,0))
axis(1, at=c(0.7,1.9,3.1,4.3), label=c('Overlap', '<1kb', '<10kb', '>10kb'))
dev.off()



fh_zero_overlap<-overlap[match(fh_ori_bed_sk1[fh_zero,2],overlap[,2]),]
fh_zero_overlap<-fh_zero_overlap[which(fh_zero_overlap[,2]>0),]
fh_zero_bins_overlap<-tRNA_gff[match(fh_zero_overlap[,5],tRNA_gff[,c(4)]),12]
fh_zero_overlap_counts<-c(length(intersect(unique(lost_only),fh_zero_bins_overlap)),length(intersect(unique(main_lost),fh_zero_bins_overlap)),length(intersect(unique(maintained_only),fh_zero_bins_overlap)),length(intersect(unique(main_gain),fh_zero_bins_overlap)),length(intersect(unique(gained_only),fh_zero_bins_overlap)),length(intersect(unique(gain_lost),fh_zero_bins_overlap)),length(intersect(unique(main_gain_lost),fh_zero_bins_overlap)))
fh_zero_overlap_counts<-(fh_zero_overlap_counts/sum(fh_zero_overlap_counts))


fh_zero_overlap_1kb<-overlap_1kb[match(fh_ori_bed_sk1[fh_zero,2],overlap_1kb[,2]+1000),]
fh_zero_overlap_1kb<-fh_zero_overlap_1kb[which(fh_zero_overlap_1kb[,2]>0),]
fh_zero_bins_overlap_1kb<-tRNA_gff[match(fh_zero_overlap_1kb[,5],tRNA_gff[,c(4)]),12]
fh_zero_overlap_1kb_counts<-c(length(intersect(unique(lost_only),fh_zero_bins_overlap_1kb)),length(intersect(unique(main_lost),fh_zero_bins_overlap_1kb)),length(intersect(unique(maintained_only),fh_zero_bins_overlap_1kb)),length(intersect(unique(main_gain),fh_zero_bins_overlap_1kb)),length(intersect(unique(gained_only),fh_zero_bins_overlap_1kb)),length(intersect(unique(gain_lost),fh_zero_bins_overlap_1kb)),length(intersect(unique(main_gain_lost),fh_zero_bins_overlap_1kb)))
fh_zero_overlap_1kb_counts<-(fh_zero_overlap_1kb_counts/sum(fh_zero_overlap_1kb_counts))

fh_zero_overlap_10kb<-overlap_10kb[match(fh_ori_bed_sk1[fh_zero,2],overlap_10kb[,2]+10000),]
fh_zero_overlap_10kb<-fh_zero_overlap_10kb[which(fh_zero_overlap_10kb[,2]>0),]
fh_zero_bins_overlap_10kb<-tRNA_gff[match(fh_zero_overlap_10kb[,5],tRNA_gff[,c(4)]),12]
fh_zero_overlap_10kb_counts<-c(length(intersect(unique(lost_only),fh_zero_bins_overlap_10kb)),length(intersect(unique(main_lost),fh_zero_bins_overlap_10kb)),length(intersect(unique(maintained_only),fh_zero_bins_overlap_10kb)),length(intersect(unique(main_gain),fh_zero_bins_overlap_10kb)),length(intersect(unique(gained_only),fh_zero_bins_overlap_10kb)),length(intersect(unique(gain_lost),fh_zero_bins_overlap_10kb)),length(intersect(unique(main_gain_lost),fh_zero_bins_overlap_10kb)))
fh_zero_overlap_10kb_counts<-(fh_zero_overlap_10kb_counts/sum(fh_zero_overlap_10kb_counts))

fh_zero_nooverlap_10kb<-overlap_10kb[match(fh_ori_bed_sk1[fh_zero,2],overlap_10kb[,2]+10000),]
fh_zero_nooverlap_10kb<-fh_zero_nooverlap_10kb[which(fh_zero_nooverlap_10kb[,2]>0),]
fh_zero_bins_nooverlap_10kb<-tRNA_gff[match(fh_zero_nooverlap_10kb[,5],tRNA_gff[,c(4)]),12]
fh_zero_nooverlap_10kb_counts<-c(length(intersect(unique(lost_only),fh_zero_bins_nooverlap_10kb)),length(intersect(unique(main_lost),fh_zero_bins_nooverlap_10kb)),length(intersect(unique(maintained_only),fh_zero_bins_nooverlap_10kb)),length(intersect(unique(main_gain),fh_zero_bins_nooverlap_10kb)),length(intersect(unique(gained_only),fh_zero_bins_nooverlap_10kb)),length(intersect(unique(gain_lost),fh_zero_bins_nooverlap_10kb)),length(intersect(unique(main_gain_lost),fh_zero_bins_nooverlap_10kb)))
fh_zero_nooverlap_10kb_counts<-(fh_zero_nooverlap_10kb_counts/sum(fh_zero_nooverlap_10kb_counts))


pdf('~/Google_Drive/Lab/Data_Analysis/tRNA/origin_trna_adjacency_forkhead_independent.pdf')
par(mar=c(5,5,5,15),xpd = T)
barplot(data.matrix(cbind(fh_zero_overlap_counts,fh_zero_overlap_1kb_counts,fh_zero_overlap_10kb_counts,fh_zero_nooverlap_10kb_counts)), xaxt='n', ylab='% of tRNAs', main='Proportion of tRNAs that change their clustering, \n relative to their adjacency to Forkhead independent origins', xlab='Distance between tRNA and adjacent origin that are Forkhead independent', col=col1)
legend("topright",col = col1,cex = 0.8,lwd = 8, lty = 1,legend= c('Lost','Lost/Maintained', 'Maintained','Maintained/Gained','Gained','Lost/Gained','Lost/Maintained/Gained' ),inset =c(-0.8,0))
axis(1, at=c(0.7,1.9,3.1,4.3), label=c('Overlap', '<1kb', '<10kb', '>10kb'))
dev.off()



WHICH == THE SAME AND ALSO








next up need to cross refenrence trnas with my list
see how many in each group. do stacked bar chart

overlap with trna groups to old trna groups. 
look at codnesin enrichment ordering. 
find which ne with introns. 

polii epxression data with subsetted data. 




see spread between groups. seee if bias depedning on how distances. 
repeat for forkhead. bound etc.

need to check how many are 


































left<-which(match(coords_MRP7[,1],as.numeric(unique(tRNA_gff[,12])))!='NA')
right<-which(match(coords_MRP7[,2],as.numeric(unique(tRNA_gff[,12])))!='NA')
trna_MRP7_coord<-coords_MRP7[left[which(match(left,right)!='NA')],]

left<-which(match(coords_MRP10[,1],as.numeric(unique(tRNA_gff[,12])))!='NA')
right<-which(match(coords_MRP10[,2],as.numeric(unique(tRNA_gff[,12])))!='NA')
trna_MRP10_coord<-coords_MRP10[left[which(match(left,right)!='NA')],]

my_index<-0
for (i in 1:nrow(tRNA_coords_MRP7)){
  out<-which(tRNA_coords_MRP10[,1]==tRNA_coords_MRP7[i,1]&tRNA_coords_MRP10[,2]==tRNA_coords_MRP7[i,2])
  if(length(out)>0){
    my_index<-c(my_index,out)
  }
}
my_index<-my_index[-1]
maintained<-tRNA_coords_MRP10[my_index,]
gained<-tRNA_coords_MRP10[-my_index,]

my_index<-0
for (i in 1:nrow(tRNA_coords_MRP10)){
  out<-which(tRNA_coords_MRP7[,1]==tRNA_coords_MRP10[i,1]&tRNA_coords_MRP7[,2]==tRNA_coords_MRP10[i,2])
  if(length(out)>0){
    my_index<-c(my_index,out)
  }
}
my_index<-my_index[-1]
lost<-tRNA_coords_MRP7[-my_index,]

tRNA_gff<-as.data.frame(tRNA_gff,stringsAsFactors=FALSE)


tRNA_gff_minimised<-cbind(as.character(tRNA_gff[,1]),as.numeric(tRNA_gff[,4]),as.character(tRNA_gff[,10]),as.numeric(tRNA_gff[,11]),as.numeric(tRNA_gff[,12]))
tRNA_gff_minimised<-as.data.frame(tRNA_gff_minimised, stringsAsFactors = F)

trna_gff_maintained<-matrix(0,1,10)
colnames(trna_gff_maintained)<-1:10

for(i in 1:nrow(maintained)){
  ind1<-which(tRNA_gff_minimised[,5]==maintained[i,1])
  ind2<-which(tRNA_gff_minimised[,5]==maintained[i,2])
  #l<-max(c(length(ind1),length(ind2)))
  trna_gff_temp<-cbind((tRNA_gff_minimised[ind1[1],c(1,2,3,4,5)]),(tRNA_gff_minimised[ind2[1],c(1,2,3,4,5)]))
  colnames(trna_gff_temp)<-1:10
  trna_gff_maintained<-rbind(trna_gff_maintained,trna_gff_temp)
print(i)
  }

trna_gff_lost<-matrix(0,1,10)
colnames(trna_gff_lost)<-1:10

for(i in 1:nrow(lost)){
  ind1<-which(tRNA_gff_minimised[,5]==lost[i,1])
  ind2<-which(tRNA_gff_minimised[,5]==lost[i,2])
  #l<-max(c(length(ind1),length(ind2)))
  trna_gff_temp<-cbind((tRNA_gff_minimised[ind1[1],c(1,2,3,4,5)]),(tRNA_gff_minimised[ind2[1],c(1,2,3,4,5)]))
  colnames(trna_gff_temp)<-1:10
  trna_gff_lost<-rbind(trna_gff_lost,trna_gff_temp)
  print(i)
}


trna_gff_gained<-matrix(0,1,10)
colnames(trna_gff_gained)<-1:10

for(i in 1:nrow(gained)){
  ind1<-which(tRNA_gff_minimised[,5]==gained[i,1])
  ind2<-which(tRNA_gff_minimised[,5]==gained[i,2])
  #l<-max(c(length(ind1),length(ind2)))
  trna_gff_temp<-cbind((tRNA_gff_minimised[ind1[1],c(1,2,3,4,5)]),(tRNA_gff_minimised[ind2[1],c(1,2,3,4,5)]))
  colnames(trna_gff_temp)<-1:10
  trna_gff_gained<-rbind(trna_gff_gained,trna_gff_temp)
  print(i)
}


trna_gff_gained_nas<-trna_gff_gained[complete.cases(trna_gff_gained),]
trna_gff_lost_nas<-trna_gff_lost[complete.cases(trna_gff_lost),]
trna_gff_maintained_nas<-trna_gff_maintained[complete.cases(trna_gff_maintained),]

intra_gained<-which((trna_gff_gained_nas[,1]==trna_gff_gained_nas[,6])=='TRUE')
inter_gained<-which((trna_gff_gained_nas[,1]==trna_gff_gained_nas[,6])!='TRUE')

intra_lost<-which((trna_gff_lost_nas[,1]==trna_gff_lost_nas[,6])=='TRUE')
inter_lost<-which((trna_gff_lost_nas[,1]==trna_gff_lost_nas[,6])!='TRUE')

intra_maintained<-which((trna_gff_maintained_nas[,1]==trna_gff_maintained_nas[,6])=='TRUE')
inter_maintained<-which((trna_gff_maintained_nas[,1]==trna_gff_maintained_nas[,6])!='TRUE')


intra_gained_prop<-length(intra_gained)/(length(intra_gained)+length(inter_gained))
inter_gained_prop<-length(inter_gained)/(length(intra_gained)+length(inter_gained))
intra_lost_prop<-length(intra_lost)/(length(intra_lost)+length(inter_lost))
inter_lost_prop<-length(inter_lost)/(length(intra_lost)+length(inter_lost))
intra_maintained_prop<-length(intra_maintained)/(length(intra_maintained)+length(inter_maintained))
inter_maintained_prop<-length(inter_maintained)/(length(intra_maintained)+length(inter_maintained))

data<-rbind(c(intra_gained_prop,intra_lost_prop,intra_maintained_prop),c(inter_gained_prop,inter_lost_prop,inter_maintained_prop))
rownames(data)<-c('intra','inter')
colnames(data)<-c('gained','lost','maintained')

barplot(data, col=colors()[c(89,12)] , border="white", space=0.04, font.axis=2, xlab="group")














#####lost gained intra inter. 

plot(0,0,ylim=c(1,17), xlim=c(1,1500000), yaxt="n", ylab="Chr number", xlab="Chr position")
axis(2,at=(1.25:16.25), labels=(16:1), las=1)

intra_gained<-trna_gff_gained_nas[intra_gained[-1],]
inter_gained<-trna_gff_gained_nas[inter_gained,]
intra_gained_fin<-cbind(matrix('intra_gained',nrow(intra_gained),1),intra_gained[,2],intra_gained[,1])
inter_gained_fin<-cbind(matrix('inter_gained',nrow(inter_gained),1),inter_gained[,2],inter_gained[,1])
unique(intra_gained_fin)->intra_gained_fin
unique(inter_gained_fin)->inter_gained_fin

intra_maintained<-trna_gff_maintained_nas[intra_maintained[-1],]
inter_maintained<-trna_gff_maintained_nas[inter_maintained,]
intra_maintained_fin<-cbind(matrix('intra_maintained',nrow(intra_maintained),1),intra_maintained[,2],intra_maintained[,1])
inter_maintained_fin<-cbind(matrix('inter_maintained',nrow(inter_maintained),1),inter_maintained[,2],inter_maintained[,1])
unique(intra_maintained_fin)->intra_maintained_fin
unique(inter_maintained_fin)->inter_maintained_fin

intra_lost<-trna_gff_lost_nas[intra_lost[-1],]
inter_lost<-trna_gff_lost_nas[inter_lost,]
intra_lost_fin<-cbind(matrix('intra_lost',nrow(intra_lost),1),intra_lost[,2],intra_lost[,1])
inter_lost_fin<-cbind(matrix('inter_lost',nrow(inter_lost),1),inter_lost[,2],inter_lost[,1])
unique(intra_lost_fin)->intra_lost_fin
unique(inter_lost_fin)->inter_lost_fin


select_few_deets<-rbind(intra_gained_fin,inter_gained_fin,intra_maintained_fin,inter_maintained_fin,intra_lost_fin,inter_lost_fin)


for (i in 1:16){
  rect(xleft=(1), ybottom=(17.5-i), xright=(chr.size[i]), ytop=(17.5-(i+0.5)))
}

chrs<-c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
my_trna_types1<-c("inter_lost","intra_gained","inter_gained","intra_maintained","inter_maintained","intra_lost")

my_cols <- brewer.pal(9, "Set1")

for (j in 1:length(my_trna_types1)){
  subset<-select_few_deets[which(select_few_deets[,1]==my_trna_types1[j]),]
  for (i in 1:nrow(subset)){
    index<-which(chrs==subset[i,3])
    segments(x0=(as.numeric(as.matrix(subset[i,2]))),y0=(17.5-index),y1=(17.5-(index+0.5)), lwd=4, col=my_cols[j])   
  }}


#subset

cen_gff<-gff[grep("CEN", gff[,9]),]
for(i in 1:16){
  points(x=cen_gff[i,4],y=(17.5-(i+0.25)), col='black')   
}

trna_nums<-matrix(0,1,1)
for (i in 1:length(my_select_few)){
  index<-which(tRNA_groups[,1]==my_select_few[i])
  trna_nums<-rbind(trna_nums,length(index))
}
trna_nums<-as.character(trna_nums[-1])

my_select_few<-unique(select_few_deets[,1])
trna_nums<-c(nrow(intra_gained_fin),nrow(inter_gained_fin),nrow(intra_maintained_fin),nrow(inter_maintained_fin),nrow(intra_lost_fin),nrow(inter_lost_fin))
my_cols2<-c("#377EB8" ,"#4DAF4A", "#984EA3" ,"#FF7F00", "#FFFF33", "#E41A1C")
pos=(9.25:1.25)
for (i in 1:length(my_select_few)){
  text(x=1200000,y=pos[i], my_select_few[i])
  segments(x0=1350000,x1=1400000,y0=pos[i], lwd=4, col=my_cols2[i]) 
  #text(x=1450000,y=pos[i], trna_nums[i])
}

("inter_lost","intra_gained","inter_gained","intra_maintained","inter_maintained","intra_lost")


intra_prob<-c(sum(trna_nums[c(1,3,5)]))/(c(sum(trna_nums[c(1,3,5)])+sum(trna_nums[c(2,4,6)])))
inter_prob<-c(sum(trna_nums[c(2,4,6)]))/(c(sum(trna_nums[c(1,3,5)])+sum(trna_nums[c(2,4,6)])))

res<-chisq.test(trna_nums[1:2], p=c(intra_prob,inter_prob))

