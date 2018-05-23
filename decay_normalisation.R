###################################
# Correction of decay by distance #
###################################

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
wtrapa<-read.table('~/Google_Drive/Lab/Data_Analysis/Heatmaps/MRP7-10/wtrapa')
conrapa<-read.table('~/Google_Drive/Lab/Data_Analysis/Heatmaps/MRP7-10/conrapa')

wtrapa<-data.matrix(wtrapa)
conrapa<-data.matrix(conrapa)

#Normalise to total interaction counts per dataset
wtrapa_norm<-data.matrix(wtrapa/sum(wtrapa))
conrapa_norm<-data.matrix(conrapa/sum(conrapa))

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

#Lets cure for decay curve. Break down into single chromosomes. Then get the scores at each distance. 
for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("pre_norm_wtrapa_norm_chr",h), wtrapa_norm[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
  chr_matrix<-get(paste0("pre_norm_wtrapa_norm_chr",h))
  my_dims<-nrow(chr_matrix)
  chr_decay_matrix<-matrix(0,150,3)
  for (i in 1:my_dims){
    for (j in 1:my_dims){
      if (abs(i-j)>0){
        sing_dist<-abs(i-j)
        chr_decay_matrix[sing_dist,2]<-chr_decay_matrix[sing_dist,2]+chr_matrix[i,j]
        chr_decay_matrix[sing_dist,3]<-chr_decay_matrix[sing_dist,3]+1
      }
      if (h==1){full_decay_matrix<-chr_decay_matrix
      }else{
        full_decay_matrix<-full_decay_matrix+chr_decay_matrix
      }}}}

#What is the mean score per point
decay<-full_decay_matrix[,2]/full_decay_matrix[,3]

#Normalise the data to the decay curve. Divide each point on norm  by the distance score for decay. 
for (h in 1:(length(incrementalbins)-1)){
  chr_matrix<-get(paste0("pre_norm_wtrapa_norm_chr",h))
  my_dims<-nrow(chr_matrix)
  for (i in 1:my_dims){
    for (j in 1:my_dims){
      if (abs(i-j)>0){
        chr_matrix[i,j]<-chr_matrix[i,j]/decay[abs(i-j)]
      }}}
  ind<-which(chr_matrix=='NaN')
  chr_matrix[ind]<-0
  assign(paste0("norm_wtrapa_norm_chr",h),chr_matrix)
}

#Place this back into the context of the whole genome matrix. Need to create matrix with just inter, to nomrlaise just inter-chromosomal 
wtrapa_norm->wtrapa_norm_decay_corrected
for (h in 1:(length(incrementalbins)-1)){
  wtrapa_norm_decay_corrected[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]]<-0
}

#What about the ean interchromomsoal interactions. Normalise for this next. 
mean_inter<-sum(wtrapa_norm_decay_corrected)/length(which(wtrapa_norm_decay_corrected>0))
wtrapa_norm_decay_corrected_final<-wtrapa_norm_decay_corrected
my_dims<-nrow(wtrapa_norm_decay_corrected_final)
for (i in 1:my_dims){
  for (j in 1:my_dims){
    if (abs(i-j)>0){
      wtrapa_norm_decay_corrected_final[i,j]<-wtrapa_norm_decay_corrected[i,j]/mean_inter
    }}}

# Put the normalised intra chromosmal matrices back in
for (h in 1:(length(incrementalbins)-1)){
  wtrapa_norm_decay_corrected_final[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]]<-get(paste0("norm_wtrapa_norm_chr",h))
}


#Lets cure for decay curve. Break down into single chromosomes. Then get the scores at each distance. 
for (h in 1:(length(incrementalbins)-1)){
  assign(paste0("pre_norm_conrapa_norm_chr",h), conrapa_norm[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]])
  chr_matrix<-get(paste0("pre_norm_conrapa_norm_chr",h))
  my_dims<-nrow(chr_matrix)
  chr_decay_matrix<-matrix(0,150,3)
  for (i in 1:my_dims){
    for (j in 1:my_dims){
      if (abs(i-j)>0){
        sing_dist<-abs(i-j)
        chr_decay_matrix[sing_dist,2]<-chr_decay_matrix[sing_dist,2]+chr_matrix[i,j]
        chr_decay_matrix[sing_dist,3]<-chr_decay_matrix[sing_dist,3]+1
      }
      if (h==1){full_decay_matrix<-chr_decay_matrix
      }else{
        full_decay_matrix<-full_decay_matrix+chr_decay_matrix
      }}}}

#What is the mean score per point
decay<-full_decay_matrix[,2]/full_decay_matrix[,3]

#Normalise the data to the decay curve. Divide each point on norm  by the distance score for decay. 
for (h in 1:(length(incrementalbins)-1)){
  chr_matrix<-get(paste0("pre_norm_conrapa_norm_chr",h))
  my_dims<-nrow(chr_matrix)
  for (i in 1:my_dims){
    for (j in 1:my_dims){
      if (abs(i-j)>0){
        chr_matrix[i,j]<-chr_matrix[i,j]/decay[abs(i-j)]
      }}}
  ind<-which(chr_matrix=='NaN')
  chr_matrix[ind]<-0
  assign(paste0("norm_conrapa_norm_chr",h),chr_matrix)
}

#Place this back into the context of the whole genome matrix. Need to create matrix with just inter, to nomrlaise just inter-chromosomal 
conrapa_norm->conrapa_norm_decay_corrected
for (h in 1:(length(incrementalbins)-1)){
  conrapa_norm_decay_corrected[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]]<-0
}

#What about the ean interchromomsoal interactions. Normalise for this next. 
mean_inter<-sum(conrapa_norm_decay_corrected)/length(which(conrapa_norm_decay_corrected>0))
conrapa_norm_decay_corrected_final<-conrapa_norm_decay_corrected
my_dims<-nrow(conrapa_norm_decay_corrected_final)
for (i in 1:my_dims){
  for (j in 1:my_dims){
    if (abs(i-j)>0){
      conrapa_norm_decay_corrected_final[i,j]<-conrapa_norm_decay_corrected[i,j]/mean_inter
    }}}

# Put the normalised intra chromosmal matrices back in
for (h in 1:(length(incrementalbins)-1)){
  conrapa_norm_decay_corrected_final[(incrementalbins[h]+1):incrementalbins[h+1],(incrementalbins[h]+1):incrementalbins[h+1]]<-get(paste0("norm_conrapa_norm_chr",h))
}

#Export data
write.csv(conrapa_norm_decay_corrected_final, '~/Box Sync/Lab/Data_Analysis/Heatmaps/Decay_normalised/conrapa_decay_normalised.csv', row.names = F)
write.csv(wtrapa_norm_decay_corrected_final, '~/Box Sync/Lab/Data_Analysis/Heatmaps/Decay_normalised/wtrapa_decay_normalised.csv', row.names = F)

