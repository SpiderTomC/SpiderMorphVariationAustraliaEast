#initial file read#
AllSpiders = read.csv( "your file path", stringsAsFactors=F)
allspiders= AllSpiders[,5:29]
allspiderslog=log(allspiders)
allspidersdist1=dist(scale(allspiderslog))
allspiderspcoa1= cmdscale(allspidersdist1,k=10, eig=T)
allspiderspcoa1$points[,1]=allspiderspcoa1$points[,1]*-1
#size adjusted matrix
adjspider= matrix(nrow = 553,ncol = 25)
for (i in 1:nrow(allspiders)) {
  rowratio<-mean(allspiders[,2],na.rm=T)/allspiders[i,2]
  for (j in 1:ncol(allspiders)){
    adjspider[i,j]<-allspiders[i,j]*rowratio
  }
}
#axis extreme comparisons
axis1max= allspiders[which(allspiderspcoa1$points[,1]>(mean(allspiderspcoa1$points[,1])+sd(allspiderspcoa1$points[,1]))),]
axis1min= allspiders[which(allspiderspcoa1$points[,1]<(mean(allspiderspcoa1$points[,1])-sd(allspiderspcoa1$points[,1]))),]
axis1maxmean=geometric.mean(axis1max)
axis1minmean=geometric.mean(axis1min)
axis2max= allspiders[which(allspiderspcoa1$points[,2]>(mean(allspiderspcoa1$points[,2])+sd(allspiderspcoa1$points[,2]))),]
axis2min= allspiders[which(allspiderspcoa1$points[,2]<(mean(allspiderspcoa1$points[,2])-sd(allspiderspcoa1$points[,2]))),]
axis2maxmean=geometric.mean(axis2max)
axis2minmean=geometric.mean(axis2min)
axis3max= allspiders[which(allspiderspcoa1$points[,3]>(mean(allspiderspcoa1$points[,3])+sd(allspiderspcoa1$points[,3]))),]
axis3min= allspiders[which(allspiderspcoa1$points[,3]<(mean(allspiderspcoa1$points[,3])-sd(allspiderspcoa1$points[,3]))),]
axis3maxmean=geometric.mean(axis3max)
axis3minmean=geometric.mean(axis3min)
axis4max= allspiders[which(allspiderspcoa1$points[,4]>(mean(allspiderspcoa1$points[,4])+sd(allspiderspcoa1$points[,4]))),]
axis4min= allspiders[which(allspiderspcoa1$points[,4]<(mean(allspiderspcoa1$points[,4])-sd(allspiderspcoa1$points[,4]))),]
axis4maxmean=geometric.mean(axis4max)
axis4minmean=geometric.mean(axis4min)
#with adjusted matrix
axis1maxadj= adjspider[which(adjspiderpcoa$points[,1]>(mean(adjspiderpcoa$points[,1])+sd(adjspiderpcoa$points[,1]))),]
axis1minadj= adjspider[which(adjspiderpcoa$points[,1]<(mean(adjspiderpcoa$points[,1])-sd(adjspiderpcoa$points[,1]))),]
axis1maxmeanadj=geometric.mean(axis1maxadj)
axis1minmeanadj=geometric.mean(axis1minadj)
axis2maxadj= adjspider[which(adjspiderpcoa$points[,2]>(mean(adjspiderpcoa$points[,2])+sd(adjspiderpcoa$points[,2]))),]
axis2minadj= adjspider[which(adjspiderpcoa$points[,2]<(mean(adjspiderpcoa$points[,2])-sd(adjspiderpcoa$points[,2]))),]
axis2maxmeanadj=geometric.mean(axis2maxadj)
axis2minmeanadj=geometric.mean(axis2minadj)
axis3maxadj= adjspider[which(adjspiderpcoa$points[,3]>(mean(adjspiderpcoa$points[,3])+sd(adjspiderpcoa$points[,3]))),]
axis3minadj= adjspider[which(adjspiderpcoa$points[,3]<(mean(adjspiderpcoa$points[,3])-sd(adjspiderpcoa$points[,3]))),]
axis3maxmeanadj=geometric.mean(axis3maxadj)
axis3minmeanadj=geometric.mean(axis3minadj)
#generation of site species lists#
squarescoord=cbind(squareslat[which(Cellcounts>220)],squareslong[which(Cellcounts>220)])
squarescoord=squarescoord[order(squarescoord[,1]),]
sitelats= c(-38.5,-37.75,-37.75,-37.5,-37,-37,-34.75,-34.5,-34.5,-34,-33.75,-33.5,-33.5,-33.5,-33.25,-33.25,-33,
            -32.25,-31.5,-31.5,-30.75,-28.25,-28,-27.5,-27.5,-27.5,-27.25,-27.25,-27.25,-26.75,-25.25,-24.25,-20.25,-17.25,
            -17,-16.75,-16.75,-16.5,-16)
sitelongs= c(143.5, 145.25,147, 145.75,149.25,149.5,146,150.75,151,151.25,151.25,150.5,151.25,151.5,151.5,151.75,151.75,
             152.75,152.75,159.25,152.25,153.25,153.25,153,153.25,153.5,151.25,153,153.25,153,150.25,151.25,148.75,146,
             145.75,145.75,146,145.5,145.25)
sitelatclim=c(-38.5,-37.75,-37.75,-37.5,-37,-37,-34.75,-34.5,-34.5,-34.004,-33.75,-33.5,-33.5,-33.5,-33.25,-33.29,-33,
              -32.25,-31.5,-31.55,-30.75,-28.25,-28,-27.5,-27.5,-27.54,-27.25,-27.25,-27.25,-26.75,-25.25,-24.25,-20.28,-17.25,
              -17,-16.75,-16.88,-16.5,-16)
sitelongclim=c(143.5, 145.25,147, 145.75,149.25,149.5,146,150.75,150.88,151.22,151.25,150.5,151.25,151.4,151.5,151.54,151.72,
               152.54,152.75,159.08,152.25,153.25,153.25,153,153.25,153.27,151.25,153,153.09,153,150.25,151.25,148.71,145.84,
               145.75,145.63,145.92,145.46,145.25)
sitenames= c("Great Otway National Park","Warrandyte","Avon Wilderness","Marysville","Rosemeath","South East Forest","Coleambally",
             "Avon/East Kangaloon","Wollongong","Botany Bay South","Sydney North Shore","Blue Mountains","Brisbane Water National Park",
             "Bouddi National Park","Central Coast","Tuggerah","Awabakal","Wallis Lake","Bago Bluff National Park","Lord Howe Island","Cunnawurra/Carrai","Limpinwood",
             "Lamington National Park","Southwest Brisbane","Southeast Brisbane","Raby Bay","Dalby/Kumbarilla","Northwest Brisbane","Moreton Bay",
             "Dularcha National Park","Precipice National Park","Wietalaba State Forest","Conway National Park","Gadgarra/Wooroonooran",
             "Atherton Tablelands","Kuranda State Forest","Green Island/Trinity Forest","Mowbray National Park","Daintree")
eastspidlist= list()
for (i in 1:39) {
  eastspidsite <- EastSpiders[which(EastSpiders$Latitude<= sitelats[i] &EastSpiders$Latitude>sitelats[i]-.25 & EastSpiders$Longitude<=sitelongs[i]& EastSpiders$Longitude>sitelongs[i]-.25),]
  eastnames= rep(sitenames[i],nrow(eastspidsite))
  eastspidsite=cbind(eastspidsite,eastnames)
  eastspidlist[[i]]<-eastspidsite
  
}
library(data.table)
allsitesoccu <-dplyr::bind_rows(eastspidlist)
setDT(allsitesoccu)
eastspidpres= dcast(allsitesoccu, allsitesoccu$Species~allsitesoccu$eastnames,length)
eastspidpres<-as.data.frame(eastspidpres)
specieslist=list()
for (i in 2:40) {
  speciessite<- eastspidpres[which(eastspidpres[,i]>=1),]
  speciessite<-(speciessite[,1])
  specieslist[[i]]<-(speciessite)
}
specieslist<-specieslist[2:40]
names(specieslist)[1:39]<-sort(sitenames)
#taxonomic groupings
allspiderfams<- list()
for (i in unique(AllSpiders$Family)){
  allfams<- which(AllSpiders[,1]==i)
  allspiderfams[[i]]<- allfams
}
mygalo= c(allspiderfams$Actinopodidae,allspiderfams$Atracidae,allspiderfams$Barychelidae,allspiderfams$Dipluridae,allspiderfams$Hexathelidae,allspiderfams$Idiopidae,allspiderfams$Nemesiidae,allspiderfams$Microstigmatidae)
synsperm= c(allspiderfams$Dysderidae,allspiderfams$Filistatidae,allspiderfams$Orsolobidae,allspiderfams$Pholcidae,allspiderfams$Segestriidae,allspiderfams$Scytodidae,allspiderfams$Gradungulidae,allspiderfams$Archaeidae)
araneio= c(allspiderfams$Nicodamidae,allspiderfams$Theridiidae,allspiderfams$Anapidae,allspiderfams$Araneidae,allspiderfams$Malkaridae,allspiderfams$Mimetidae,allspiderfams$Arkyidae,allspiderfams$Tetragnathidae,allspiderfams$Cyatholipidae,allspiderfams$Physoglenidae, allspiderfams$Linyphiidae)
entely= c(allspiderfams$Deinopidae,allspiderfams$Hersilidae,allspiderfams$Uloboridae)
marro= c(allspiderfams$Zodariidae, allspiderfams$Amaurobiidae,allspiderfams$Desidae,allspiderfams$Stiphidiidae,allspiderfams$Sparassidae)
calamis= c(allspiderfams$Zoropsidae,allspiderfams$Ctenidae,allspiderfams$Oxyopidae,allspiderfams$Pisauridae,allspiderfams$Lycosidae,allspiderfams$Thomisidae)
dionycha= c(allspiderfams$Prodidomidae,allspiderfams$Gallieniellidae,allspiderfams$Trochanteriidae,allspiderfams$Clubionidae,allspiderfams$Phrurolithidae,allspiderfams$Gnaphosidae,allspiderfams$Lamponidae,allspiderfams$Corinnidae,allspiderfams$Selenopidae,allspiderfams$Miturgidae,allspiderfams$Cheiracanthiidae,allspiderfams$Salticidae)
araneomorph= c(synsperm, araneio, entely, marro, calamis, dionycha)
phylogroups= list(mygalo,synsperm,araneio,entely,marro,calamis,dionycha)

#simplified phylogeny#
phylogenerator=matrix(nrow=7, ncol=7, dimnames = list(c("Mygalomorphae","Synspermiata","Araneoidea","Other Entelygynae","Marronoid clade","Oval calamistrum clade","Dionycha"),c(1,2,3,4,5,6,7)))
phylogenerator[1,]=c(1,0,0,0,0,0,0)
phylogenerator[2,]=c(1,1,0,0,0,0,0)
phylogenerator[3,]=c(1,1,1,0,0,0,0)
phylogenerator[4,]=c(1,1,1,1,0,0,0)
phylogenerator[5,]=c(1,1,1,1,1,0,0)
phylogenerator[6,]=c(1,1,1,1,1,1,0)
phylogenerator[7,]=c(1,1,1,1,1,1,1)
plot(hclust(vegdist(phylogenerator),method = "single"),cex=2.5)

png('Fig 01 Simplified phylogeny.png', width=800, height=1000)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(hclust(vegdist(phylogenerator),method = "single"),cex=2.2, axes=F, xlab='',ylab='',sub='', main='')
dev.off()

mygalopch= rep (16,64)
mygalopch[1:2]= rep (1,2)
mygalopch[3:12]=rep (15,10)
mygalopch[13:27]=rep (6,15)
mygalopch[28:34]=rep (17,7)
mygalopch[36:37]=rep (3,2)
mygalopch[38:45]=rep (4,8)
mygalopch[64]=rep(7,1)

synspermpch= rep(16,22)
synspermpch[1]=rep(1,1)
synspermpch[2]=rep(4,1)
synspermpch[3:6]=rep(15,4)
synspermpch[7:11]=rep(6,5)
synspermpch[12]=rep(17,1)
synspermpch[13]=rep(3,1)
synspermpch[14:20]=rep(12,7)

araneiopch=rep(16,134)
araneiopch[1:9]=rep(1,9)
araneiopch[10:27]=rep(4,18)
araneiopch[28:36]=rep(15,9)
araneiopch[37:99]=rep(6,63)
araneiopch[100:121]=rep(17,22)
araneiopch[122:130]=rep(3,9)

entelypch=rep(16,10)
entelypch[1:4]=rep(1,4)
entelypch[5:8]=rep(15,4)

marropch= rep(4,81)
marropch[1:22]=rep (16,22)
marropch[23:28]=rep (15,6)
marropch[29:39]=rep(6,11)
marropch[40:53]=rep(17,14)


calamispch=rep(16,98)
calamispch[1:7]=rep(1,7)
calamispch[8]=rep(15,1)
calamispch[9:15]=rep(6,7)
calamispch[16:21]=rep(17,6)
calamispch[22:65]=rep(4,44)

dionychapch=rep(16,144)
dionychapch[1:13]=rep(1,13)
dionychapch[14:20]=rep(15,7)
dionychapch[21:53]=rep(6,33)
dionychapch[54:65]=rep(17,12)
dionychapch[66:75]=rep(4,10)

#hunting groupings
ambush= which(AllSpiders[,3]=="Ambush")
ambushsilk=which(AllSpiders[,3]=="Ambush silk")
burrowers=which(AllSpiders[,3]=="Burrow")
hunters=which(AllSpiders[,3]=="Ground hunter")
irregularweb=which(AllSpiders[,3]=="Irregular web")
mimic=which(AllSpiders[,3]=="Mimicry")
orbweaver=which(AllSpiders[,3]=="Orb-web")
retreats=which(AllSpiders[,3]=="Retreat")
sheetwebs=which(AllSpiders[,3]=='Sheet web')
ambushers= c(ambush,ambushsilk)
sedentary=c(burrowers,retreats)
webusers=c(orbweaver,irregularweb,sheetwebs,ambushsilk)
nonweb=c(burrowers,ambush,mimic,hunters,retreats)



ambushpch=rep(16, 69)
ambushpch[1:61]=rep(17,61)
sedentpch=rep(16,82)
sedentpch[1:49]=rep(17,49)
#phylogenetic centroids#
spidpcoacoord=allspiderspcoa1$points[,1:2]
phylocentroids<-array()
for (j in 1:length(phylogroups)){
  phylosample= sample(allspiderspcoa1$points[phylogroups[[j]],1],length(phylogroups[[j]]))
  phylocentroids[j]<- mean(phylosample)
  
}
phylocentroids2<-array()
for (j in 1:length(phylogroups)){
  
  phylosample= sample(allspiderspcoa1$points[phylogroups[[j]],2],length(phylogroups[[j]]))
  phylocentroids2[j]<- mean(phylosample)
}


phylocentroids3<-array()
for (j in 1:length(phylogroups)){
  phylosample= sample(allspiderspcoa1$points[phylogroups[[j]],3],length(phylogroups[[j]]))
  phylocentroids3[j]<- mean(phylosample)
}


phylocentmat= matrix(cbind(phylocentroids,phylocentroids2, phylocentroids3), nrow=7, ncol = 3)
#distance calculations#
totcentroids<- array()
for(j in 1:3) {
  totsample= sample(allspiderspcoa1$points[,j],553)
  totcentroids[j]<- mean(totsample)
}
alldistance<-array()
for (i in 1:nrow(allspiders)){
  alldistance[i]<- dist(rbind(allspiderspcoa1$points[i,1:3],totcentroids),method="euclidean")^2
}
sum(alldistance)
#total variance explained for all three primary axes
compcentroids<- array()
for(j in 1:10) {
  compsample= sample(allspiderspcoa1$points[,j],553)
  compcentroids[j]<- mean(compsample)
}
compdistance<-array()
for (i in 1:nrow(allspiders)){
  compdistance[i]<- dist(rbind(allspiderspcoa1$points[i,1:10],compcentroids),method="euclidean")^2
}
sum(compdistance)
compaxisdistance<-array()
allaxisdistances<- array()
for (j in 1:10) {
  for (i in 1:nrow(allspiders)){
    compaxisdistance[i]<- dist(rbind(allspiderspcoa1$points[i,j],compcentroids[j]),method="euclidean")^2
  }
  allaxisdistances[j]<-sum(compaxisdistance)
}
#with size-adjusted
compcentroidsadj<- array()
for(j in 1:10) {
  compsample= sample(adjspiderpcoa$points[,j],553)
  compcentroidsadj[j]<- mean(compsample)
}
compdistanceadj<-array()
for (i in 1:nrow(allspiders)){
  compdistanceadj[i]<- dist(rbind(adjspiderpcoa$points[i,1:10],compcentroidsadj),method="euclidean")^2
}
sum(compdistanceadj)
compaxisdistanceadj<-array()
allaxisdistancesadj<- array()
for (j in 1:10) {
  for (i in 1:nrow(allspiders)){
    compaxisdistanceadj[i]<- dist(rbind(adjspiderpcoa$points[i,j],compcentroidsadj[j]),method="euclidean")^2
  }
  allaxisdistancesadj[j]<-sum(compaxisdistanceadj)
}
#variance explained actual
axisvariance<- array()
for (i in 1:10){
  axisvariance[i]<-(1-(allaxisdistances[i])/sum(compdistance))*100
}
#adjusted
axisvarianceadj<- array()
for (i in 1:10){
  axisvarianceadj[i]<-(1- sum(allaxisdistancesadj[i])/sum(compdistanceadj))*100
}
#phylogenetic distances
phylodistances <- array()
for (i in 1:7){
  group= unlist(phylogroups[[i]])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(allspiderspcoa1$points[group[j],1:3],phylocentmat[i,]))^2
  }
  phylodistances[i]<-sum(groupdist)
}
(1- sum(phylodistances)/sum(alldistance))*100
huntdistances <- array()
for (i in 1:9){
  group= which(AllSpiders[,3]==huntgroups[i])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(allspiderspcoa1$points[group[j],1:3],huntcentmat[i,]))^2
  }
  huntdistances[i]<-sum(groupdist)
}
(1- sum(huntdistances)/sum(alldistance))*100
#by axis variance
#axis 1
axisdistance<-array()
for (i in 1:nrow(allspiders)){
  axisdistance[i]<- dist(rbind(allspiderspcoa1$points[i,1],totcentroids[1]),method="euclidean")^2
}
sum(axisdistance)
phylodistaxis1 <- array()
for (i in 1:7){
  group= unlist(phylogroups[[i]])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(allspiderspcoa1$points[group[j],1],phylocentmat[i,1]))^2
  }
  phylodistaxis1[i]<-sum(groupdist)
}
(1- sum(phylodistaxis1)/sum(axisdistance))*100
huntdistaxis1 <- array()
for (i in 1:9){
  group= which(AllSpiders[,3]==huntgroups[i])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(allspiderspcoa1$points[group[j],1],huntcentmat[i,1]))^2
  }
  huntdistaxis1[i]<-sum(groupdist)
}
(1- sum(huntdistaxis1)/sum(axisdistance))*100
#axis 2
axisdistance2<-array()
for (i in 1:nrow(allspiders)){
  axisdistance2[i]<- dist(rbind(allspiderspcoa1$points[i,2],totcentroids[2]),method="euclidean")^2
}
sum(axisdistance2)
phylodistaxis2 <- array()
for (i in 1:7){
  group= unlist(phylogroups[[i]])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(allspiderspcoa1$points[group[j],2],phylocentmat[i,2]))^2
  }
  phylodistaxis2[i]<-sum(groupdist)
}
(1- sum(phylodistaxis2)/sum(axisdistance2))*100
huntdistaxis2 <- array()
for (i in 1:9){
  group= which(AllSpiders[,3]==huntgroups[i])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(allspiderspcoa1$points[group[j],2],huntcentmat[i,2]))^2
  }
  huntdistaxis2[i]<-sum(groupdist)
}
(1- sum(huntdistaxis2)/sum(axisdistance2))*100

#axis 3
axisdistance3<-array()
for (i in 1:nrow(allspiders)){
  axisdistance3[i]<- dist(rbind(allspiderspcoa1$points[i,3],totcentroids[3]),method="euclidean")^2
}
sum(axisdistance3)
phylodistaxis3 <- array()
for (i in 1:7){
  group= unlist(phylogroups[[i]])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(allspiderspcoa1$points[group[j],3],phylocentmat[i,3]))^2
  }
  phylodistaxis3[i]<-sum(groupdist)
}
(1- sum(phylodistaxis3)/sum(axisdistance3))*100
huntdistaxis3 <- array()
for (i in 1:9){
  group= which(AllSpiders[,3]==huntgroups[i])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(allspiderspcoa1$points[group[j],3],huntcentmat[i,3]))^2
  }
  huntdistaxis3[i]<-sum(groupdist)
}
(1- sum(huntdistaxis3)/sum(axisdistance3))*100
#with adjusted matrix
#phylogenetic centroids#
phylocentroidsadj<-array()
for (j in 1:length(phylogroups)){
  phylosample= sample(adjspiderpcoa$points[phylogroups[[j]],1],length(phylogroups[[j]]))
  phylocentroidsadj[j]<- mean(phylosample)
  
}
phylocentroidsadj2<-array()
for (j in 1:length(phylogroups)){
  
  phylosample= sample(adjspiderpcoa$points[phylogroups[[j]],2],length(phylogroups[[j]]))
  phylocentroidsadj2[j]<- mean(phylosample)
}


phylocentroidsadj3<-array()
for (j in 1:length(phylogroups)){
  phylosample= sample(adjspiderpcoa$points[phylogroups[[j]],3],length(phylogroups[[j]]))
  phylocentroidsadj3[j]<- mean(phylosample)
}


phylocentmatadj= matrix(cbind(phylocentroidsadj,phylocentroidsadj2, phylocentroidsadj3), nrow=7, ncol = 3)

#hunting centroids
huntcentroidsadj<-array()
for (j in 1:length(huntgroups)){
  huntsample= sample(adjspiderpcoa$points[which(AllSpiders[,3]==huntgroups[j]),1],nrow(allspiders[which(AllSpiders[,3]==huntgroups[j]),]))
  huntcentroidsadj[j]<- mean(huntsample)
}

huntcentroidsadj2<-array()
for (j in 1:length(huntgroups)){
  huntsample= sample(adjspiderpcoa$points[which(AllSpiders[,3]==huntgroups[j]),2],nrow(allspiders[which(AllSpiders[,3]==huntgroups[j]),]))
  huntcentroidsadj2[j]<- mean(huntsample)
}

huntcentroidsadj3<-array()
for (j in 1:length(huntgroups)){
  huntsample= sample(adjspiderpcoa$points[which(AllSpiders[,3]==huntgroups[j]),3],nrow(allspiders[which(AllSpiders[,3]==huntgroups[j]),]))
  huntcentroidsadj3[j]<- mean(huntsample)
}

huntcentmatadj= matrix(cbind(huntcentroidsadj,huntcentroidsadj2, huntcentroidsadj3), nrow=9, ncol = 3)

#distance calculations#
totcentroidsadj<- array()
for(j in 1:3) {
  totsample= sample(adjspiderpcoa$points[,j],553)
  totcentroidsadj[j]<- mean(totsample)
}
alldistanceadj<-array()
for (i in 1:nrow(adjspider)){
  alldistanceadj[i]<- dist(rbind(adjspiderpcoa$points[i,1:3],totcentroidsadj),method="euclidean")^2
}
sum(alldistanceadj)
phylodistancesadj <- array()
for (i in 1:7){
  group= unlist(phylogroups[[i]])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(adjspiderpcoa$points[group[j],1:3],phylocentmatadj[i,]))^2
  }
  phylodistancesadj[i]<-sum(groupdist)
}
(1- sum(phylodistancesadj)/sum(alldistanceadj))*100
huntdistancesadj <- array()
for (i in 1:9){
  group= which(AllSpiders[,3]==huntgroups[i])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(adjspiderpcoa$points[group[j],1:3],huntcentmatadj[i,]))^2
  }
  huntdistancesadj[i]<-sum(groupdist)
}
(1- sum(huntdistancesadj)/sum(alldistanceadj))*100
#total var explained by axis
allspiderspcoa1$eig[1]/sum(abs(allspiderspcoa1$eig))
allspiderspcoa1$eig[2]/sum(abs(allspiderspcoa1$eig))
allspiderspcoa1$eig[3]/sum(abs(allspiderspcoa1$eig))
#adjusted
adjspiderpcoa$eig[1]/sum(abs(adjspiderpcoa$eig))
adjspiderpcoa$eig[2]/sum(abs(adjspiderpcoa$eig))
adjspiderpcoa$eig[3]/sum(abs(adjspiderpcoa$eig))
#one by one axis variance
#axis 1
axisdistanceadj<-array()
for (i in 1:nrow(allspiders)){
  axisdistanceadj[i]<- dist(rbind(adjspiderpcoa$points[i,1],totcentroidsadj[1]),method="euclidean")^2
}
sum(axisdistanceadj)
phylodistaxis1adj <- array()
for (i in 1:7){
  group= unlist(phylogroups[[i]])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(adjspiderpcoa$points[group[j],1],phylocentmatadj[i,1]))^2
  }
  phylodistaxis1adj[i]<-sum(groupdist)
}
(1- sum(phylodistaxis1adj)/sum(axisdistanceadj))*100
huntdistaxis1adj <- array()
for (i in 1:9){
  group= which(AllSpiders[,3]==huntgroups[i])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(adjspiderpcoa$points[group[j],1],huntcentmatadj[i,1]))^2
  }
  huntdistaxis1adj[i]<-sum(groupdist)
}
(1- sum(huntdistaxis1adj)/sum(axisdistanceadj))*100
#axis 2
axisdistance2adj<-array()
for (i in 1:nrow(allspiders)){
  axisdistance2adj[i]<- dist(rbind(adjspiderpcoa$points[i,2],totcentroidsadj[2]),method="euclidean")^2
}
sum(axisdistance2adj)
phylodistaxis2adj <- array()
for (i in 1:7){
  group= unlist(phylogroups[[i]])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(adjspiderpcoa$points[group[j],2],phylocentmatadj[i,2]))^2
  }
  phylodistaxis2adj[i]<-sum(groupdist)
}
(1- sum(phylodistaxis2adj)/sum(axisdistance2adj))*100
huntdistaxis2adj <- array()
for (i in 1:9){
  group= which(AllSpiders[,3]==huntgroups[i])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(adjspiderpcoa$points[group[j],2],huntcentmatadj[i,2]))^2
  }
  huntdistaxis2adj[i]<-sum(groupdist)
}
(1- sum(huntdistaxis2adj)/sum(axisdistance2adj))*100

#axis 3
axisdistance3adj<-array()
for (i in 1:nrow(allspiders)){
  axisdistance3adj[i]<- dist(rbind(adjspiderpcoa$points[i,3],totcentroidsadj[3]),method="euclidean")^2
}
sum(axisdistance3adj)
phylodistaxis3adj <- array()
for (i in 1:7){
  group= unlist(phylogroups[[i]])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(adjspiderpcoa$points[group[j],3],phylocentmatadj[i,3]))^2
  }
  phylodistaxis3adj[i]<-sum(groupdist)
}
(1- sum(phylodistaxis3adj)/sum(axisdistance3adj))*100
huntdistaxis3adj <- array()
for (i in 1:9){
  group= which(AllSpiders[,3]==huntgroups[i])
  groupdist<- array()
  for (j in 1:length(group)) {
    groupdist[j]<-dist(rbind(adjspiderpcoa$points[group[j],3],huntcentmatadj[i,3]))^2
  }
  huntdistaxis3adj[i]<-sum(groupdist)
}
(1- sum(huntdistaxis3adj)/sum(axisdistance3adj))*100

#wilcox tests for each axis#
allspidphylo <- rep ("Dionycha",553)
allspidphylo[mygalo]= rep("Mygalomorph",64)
allspidphylo[synsperm]= rep("Synspermiata",22)
allspidphylo[entely]= rep("Entelygynes",10)
allspidphylo[araneio]= rep("Araneiodea",134)
allspidphylo[marro]= rep("Marronoid",81)
allspidphylo[calamis]= rep("Oval Calamistrium",98)
allspidgroups=unique(allspidphylo)
allres <- matrix(NA, nrow= 7, ncol = 7, dimnames = list(unique(allspidphylo),unique(allspidphylo)))

for (i in 1:ncol(allres)){
  for (j in 1:nrow(allres)){
    x<- wilcox.test(allspiderspcoa1$points[,1][allspidphylo == allspidgroups[i]], allspiderspcoa1$points[,1][allspidphylo == allspidgroups[j]])
    allres[i,j] <- x$p.value
  }
}
allres2 <- matrix(NA, nrow= 7, ncol = 7, dimnames = list(unique(allspidphylo),unique(allspidphylo)))
for (i in 1:ncol(allres2)){
  for (j in 1:nrow(allres2)){
    x<- wilcox.test(allspiderspcoa1$points[,2][allspidphylo == allspidgroups[i]], allspiderspcoa1$points[,2][allspidphylo == allspidgroups[j]])
    allres2[i,j] <- x$p.value
  }
}

allres3 <- matrix(NA, nrow= 7, ncol = 7, dimnames = list(unique(allspidphylo),unique(allspidphylo)))

for (i in 1:ncol(allres3)){
  for (j in 1:nrow(allres3)){
    x<- wilcox.test(allspiderspcoa1$points[,3][allspidphylo == allspidgroups[i]], allspiderspcoa1$points[,3][allspidphylo == allspidgroups[j]])
    allres3[i,j] <- x$p.value
  }
}

#hunting method tests#
allspidhunts=unique(AllSpiders[,3])
huntres <- matrix(NA, nrow= 9, ncol = 9, dimnames = list(unique(AllSpiders[,3]),unique(AllSpiders[,3])))

for (i in 1:ncol(huntres)){
  for (j in 1:nrow(huntres)){
    x<- wilcox.test(allspiderspcoa1$points[,1][AllSpiders[,3] == allspidhunts[i]], allspiderspcoa1$points[,1][AllSpiders[,3] == allspidhunts[j]])
    huntres[i,j] <- x$p.value
  }
}
huntres2 <- matrix(NA, nrow= 9, ncol = 9, dimnames = list(unique(AllSpiders[,3]),unique(AllSpiders[,3])))
for (i in 1:ncol(huntres2)){
  for (j in 1:nrow(huntres2)){
    x<- wilcox.test(allspiderspcoa1$points[,2][AllSpiders[,3] == allspidhunts[i]], allspiderspcoa1$points[,2][AllSpiders[,3] == allspidhunts[j]])
    huntres2[i,j] <- x$p.value
  }
}

huntres3 <- matrix(NA, nrow= 9, ncol = 9, dimnames = list(unique(AllSpiders[,3]),unique(AllSpiders[,3])))

for (i in 1:ncol(huntres3)){
  for (j in 1:nrow(huntres3)){
    x<- wilcox.test(allspiderspcoa1$points[,3][AllSpiders[,3] == allspidhunts[i]], allspiderspcoa1$points[,3][AllSpiders[,3] == allspidhunts[j]])
    huntres3[i,j] <- x$p.value
  }
}
#cluster analyses
library(cluster)
spidfit <- kmeans(dist(scale(log(allspiders))), 5)
spidfit1 <- kmeans(allspiderspcoa1$points, 5) #5 cluster solution
pamspidfit <- pam(dist(scale(log(allspiders))), 5)
pamspidfit1 <- pam(allspiderspcoa1$points, 5) #5 cluster solution
adjspidfit <- pam(dist(scale(log(adjspider))), 5)
adjspidfit1 <- kmeans(dist(scale(log(adjspider))), 5) #5 cluster solution
#chi-squared tests
#with regular kmeans
spidcluster1<- allspiderspcoa1$points[which(spidfit$cluster==1),1]
spidcluster<-matrix()
phylocluster=matrix()
phylocount=matrix()
chisquaredspid <-matrix(nrow = 7, ncol=5,dimnames = list(c("Mygalomorphae","Synspermiata","Araneiodea","Other Entelygynae","Marronoid","Oval Calamistrum","Dionycha"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")))
for (i in 1:5){
  spidcluster= allspiders[which(spidfit$cluster==i),]
  for (j in 1:length(phylogroups)){
    phylocluster<- allspiders[phylogroups[[j]],]
    phylocount<- spidcluster[which(rownames(phylocluster)%in%rownames(spidcluster)),]
    chisquaredspid[j,i]<-nrow(phylocount)
  }
}
spidcluster<-matrix()
huntcluster=matrix()
huntcount=matrix()
chisquaredspid2 <-matrix(nrow = 9, ncol=5,dimnames = list(c("Ambush","Ambush silk","Burrow","Ground hunter","Irregular web","Mimicry","Orb web","Retreat","Sheet web"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")))
for (i in 1:4){
  spidcluster= allspiders[which(spidfit$cluster==i),]
  for (j in 1:length(huntgroups)){
    huntcluster<- allspiders[AllSpiders[,3]==huntgroups[j],]
    huntcount<- spidcluster[which(rownames(huntcluster)%in%rownames(spidcluster)),]
    chisquaredspid2[j,i]<-nrow(huntcount)
  }
}
#with pam
spidcluster<-matrix()
phylocluster=matrix()
phylocount=matrix()
chisquaredspid3 <-matrix(nrow = 7, ncol=5,dimnames = list(c("Mygalomorphae","Synspermiata","Araneiodea","Other Entelygynae","Marronoid","Oval Calamistrum","Dionycha"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")))
for (i in 1:5){
  spidcluster= allspiders[which(pamspidfit$clustering==i),]
  for (j in 1:length(phylogroups)){
    phylocluster<- allspiders[phylogroups[[j]],]
    phylocount<- spidcluster[which(rownames(phylocluster)%in%rownames(spidcluster)),]
    chisquaredspid3[j,i]<-nrow(phylocount)
  }
}
spidcluster<-matrix()
huntcluster=matrix()
huntcount=matrix()
chisquaredspid4 <-matrix(nrow = 9, ncol=5,dimnames = list(c("Ambush","Ambush silk","Burrow","Ground hunter","Irregular web","Mimicry","Orb web","Retreat","Sheet web"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")))
for (i in 1:5){
  spidcluster= allspiders[which(pamspidfit$clustering==i),]
  for (j in 1:length(huntgroups)){
    huntcluster<- allspiders[AllSpiders[,3]==huntgroups[j],]
    huntcount<- spidcluster[which(rownames(huntcluster)%in%rownames(spidcluster)),]
    chisquaredspid4[j,i]<-nrow(huntcount)
  }
}

spidcluster<-matrix()
phylocluster=matrix()
phylocount=matrix()
chisquaredspid3b <-matrix(nrow = 7, ncol=5,dimnames = list(c("Mygalomorphae","Synspermiata","Araneiodea","Other Entelygynae","Marronoid","Oval Calamistrum","Dionycha"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")))
for (i in 1:5){
  spidcluster= allspiders[which(pamspidfit1$clustering==i),]
  for (j in 1:length(phylogroups)){
    phylocluster<- allspiders[phylogroups[[j]],]
    phylocount<- spidcluster[which(rownames(phylocluster)%in%rownames(spidcluster)),]
    chisquaredspid3b[j,i]<-nrow(phylocount)
  }
}
#using size adjusted matrix
spidcluster<-matrix()
phylocluster=matrix()
phylocount=matrix()
chisquaredspid5 <-matrix(nrow = 7, ncol=5,dimnames = list(c("Mygalomorphae","Synspermiata","Araneiodea","Other Entelygynae","Marronoid","Oval Calamistrum","Dionycha"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")))
for (i in 1:5){
  spidcluster= adjspider[which(adjspidfit$clustering==i),]
  for (j in 1:length(phylogroups)){
    phylocluster<- adjspider[phylogroups[[j]],]
    phylocount<- rownames(phylocluster)%in%rownames(spidcluster)
    
    chisquaredspid5[j,i]<-length(phylocount[which(phylocount==TRUE)])
    
  }
}

spidcluster<-matrix()
huntcluster=matrix()
huntcount=matrix()
chisquaredspid6 <-matrix(nrow = 9, ncol=5,dimnames = list(c("Ambush","Ambush silk","Burrow","Ground hunter","Irregular web","Mimicry","Orb web","Retreat","Sheet web"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4",'Cluster 5')))
for (i in 1:5){
  spidcluster= adjspider[which(adjspidfit$clustering==i),]
  for (j in 1:length(huntgroups)){
    huntcluster<- adjspider[AllSpiders[,3]==huntgroups[j],]
    huntcount<- rownames(huntcluster)%in%rownames(spidcluster)
    chisquaredspid6[j,i]<-length(huntcount[which(huntcount==TRUE)])
  }
}

spidcluster<-matrix()
phylocluster=matrix()
phylocount=matrix()
chisquaredspid5b <-matrix(nrow = 7, ncol=5,dimnames = list(c("Mygalomorphae","Synspermiata","Araneiodea","Other Entelygynae","Marronoid","Oval Calamistrum","Dionycha"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")))
for (i in 1:5){
  spidcluster= adjspider[which(adjspidfit1$clustering==i),]
  for (j in 1:length(phylogroups)){
    phylocluster<- adjspider[phylogroups[[j]],]
    phylocount<- rownames(phylocluster)%in%rownames(spidcluster)
    
    chisquaredspid5b[j,i]<-length(phylocount[which(phylocount==TRUE)])
    
  }
}

spidcluster<-matrix()
huntcluster=matrix()
huntcount=matrix()
chisquaredspid6b <-matrix(nrow = 9, ncol=5,dimnames = list(c("Ambush","Ambush silk","Burrow","Ground hunter","Irregular web","Mimicry","Orb web","Retreat","Sheet web"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4",'Cluster 5')))
for (i in 1:5){
  spidcluster= adjspider[which(adjspidfit1$clustering==i),]
  for (j in 1:length(huntgroups)){
    huntcluster<- adjspider[AllSpiders[,3]==huntgroups[j],]
    huntcount<- rownames(huntcluster)%in%rownames(spidcluster)
    chisquaredspid6b[j,i]<-length(huntcount[which(huntcount==TRUE)])
  }
}
#chisquared tables as percentages
chisquaredcent <-matrix(nrow = 7, ncol=5,dimnames = list(c("Mygalomorphae","Synspermiata","Araneiodea","Other Entelygynae","Marronoid","Oval Calamistrum","Dionycha"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")))
for (i in 1:5){
  chisquaredcent[,i]<- chisquaredspid[,i]
  for (j in 1:length(phylogroups)){
    chisquaredcent[j,i]<-chisquaredspid[j,i]/length(phylogroups[[j]])*100
  }
}

chisquaredcent5 <-matrix(nrow = 7, ncol=4,dimnames = list(c("Mygalomorphae","Synspermiata","Araneiodea","Other Entelygynae","Marronoid","Oval Calamistrum","Dionycha"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")))
for (i in 1:5){
  chisquaredcent5[,i]<- chisquaredspid5[,i]
  for (j in 1:length(phylogroups)){
    chisquaredcent5[j,i]<-chisquaredspid5[j,i]/length(phylogroups[[j]])*100
  }
}

chisquaredcent6 <-matrix(nrow = 9, ncol=4,dimnames = list(c("Ambush","Ambush silk","Burrow","Ground hunter","Irregular web","Mimicry","Orb web","Retreat","Sheet web"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")))
for (i in 1:5){
  chisquaredcent6[,i]<- chisquaredspid6[,i]
  for (j in 1:length(huntgroups)){
    chisquaredcent6[j,i]<-chisquaredspid6[j,i]/nrow(AllSpiders[which(AllSpiders[,3]==huntgroups[j]),])*100
  }
}
chisquaredcent4 <-matrix(nrow = 9, ncol=4,dimnames = list(c("Ambush","Ambush silk","Burrow","Ground hunter","Irregular web","Mimicry","Orb web","Retreat","Sheet web"),c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")))
for (i in 1:4){
  chisquaredcent4[,i]<- chisquaredspid4[,i]
  for (j in 1:length(huntgroups)){
    chisquaredcent4[j,i]<-chisquaredspid4[j,i]/nrow(AllSpiders[which(AllSpiders[,3]==huntgroups[j]),])*100
  }
}

adjfammeans<-matrix(nrow=length(unique(AllSpiders$Family)),ncol=25,dimnames = list(unique(AllSpiders$Family,meancol)))
fammeanadj=matrix()
adjmeanarray=array()
for (i in unique(AllSpiders$Family)){
  fammeanadj<- adjspider[which(AllSpiders[,1]==i),]
  if (is.vector(fammeanadj)==T){
    fammeanadj<-t(as.matrix(fammeanadj))
  }
  for (j in 1:25){
    adjmeanarray[j]<-geometric.mean(as.numeric(as.matrix(log(fammeanadj[,j]))),na.rm=T)
  }
  adjfammeans[i,]<-adjmeanarray
  
}
chisq.test(chisquaredspid, simulate.p.value = TRUE, B=1000)
chisq.test(chisquaredspid2,simulate.p.value = TRUE, B=1000)
chisq.test(chisquaredspid3,simulate.p.value = TRUE, B=1000)
chisq.test(chisquaredspid4,simulate.p.value = TRUE, B=1000)
chisq.test(chisquaredspid5,simulate.p.value = TRUE, B=1000)
chisq.test(chisquaredspid6,simulate.p.value = TRUE, B=1000)
chisq.test(chisquaredspid5b,simulate.p.value = TRUE, B=1000)
chisq.test(chisquaredspid6b,simulate.p.value = TRUE, B=1000)
#figures
png('Fig 01 Phylo Morphospace.png', width=800, height=800)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,2], cex=0,xlab="Relative size", ylab="Body/Leg length ratio")
points(allspiderspcoa1$points[mygalo,1],allspiderspcoa1$points[mygalo,2], cex=0.9,pch=mygalopch)
points(allspiderspcoa1$points[synsperm,1],allspiderspcoa1$points[synsperm,2], cex=0.9,pch=synspermpch,col="#D55E00")
points(allspiderspcoa1$points[araneio,1],allspiderspcoa1$points[araneio,2], cex=0.9,pch=araneiopch,col="#0072b2")
points(allspiderspcoa1$points[entely,1],allspiderspcoa1$points[entely,2], cex=0.9,col="#009e73",pch=entelypch)
points(allspiderspcoa1$points[marro,1],allspiderspcoa1$points[marro,2], cex=0.9,pch=marropch,col="#e79f00")
points(allspiderspcoa1$points[calamis,1],allspiderspcoa1$points[calamis,2], cex=0.9,pch=calamispch,col="#cc79a7")
points(allspiderspcoa1$points[dionycha,1],allspiderspcoa1$points[dionycha,2], cex=0.9,pch=dionychapch,col="#9ad0f3")
legend(-13,-2.8,c("Mygalomorphae","Synspermiata","Araneoidea", "Other Entelygynae","Marronoid clade","Oval Calamistrum","Dionycha"), col=c("#000000","#D55E00","#0072b2","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=16,cex=1,pt.cex=1,text.font=3)
dev.off()

png('Fig 01B Phylo Morphospace.png', width=800, height=800)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,2], cex=0,xlab="Relative size(Axis 1)", ylab="Body/Leg length ratio(Axis 2)")
points(allspiderspcoa1$points[mygalo,1],allspiderspcoa1$points[mygalo,2], cex=0.9,pch=16)
points(allspiderspcoa1$points[synsperm,1],allspiderspcoa1$points[synsperm,2], cex=0.9,pch=16,col="#D55E00")
points(allspiderspcoa1$points[araneio,1],allspiderspcoa1$points[araneio,2], cex=0.9,pch=16,col="#0072b2")
points(allspiderspcoa1$points[entely,1],allspiderspcoa1$points[entely,2], cex=0.9,col="#009e73",pch=16)
points(allspiderspcoa1$points[marro,1],allspiderspcoa1$points[marro,2], cex=0.9,pch=16,col="#e79f00")
points(allspiderspcoa1$points[calamis,1],allspiderspcoa1$points[calamis,2], cex=0.9,pch=16,col="#cc79a7")
points(allspiderspcoa1$points[dionycha,1],allspiderspcoa1$points[dionycha,2], cex=0.9,pch=16,col="#9ad0f3")
legend(-13,-2.8,c("Mygalomorph","Synspermiata","Araneoidea", "Other Entelygynae","Marronoid clade","Oval Calamistrum clade","Dionycha"), col=c("#000000","#D55E00","#0072b2","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=16,cex=1,pt.cex=1,text.font=3)
dev.off()

png('Fig 02 Phylo Morphospace.png', width=800, height=800)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,3], cex=0,xlab="Relative size", ylab="Proximal/distal leg ratio")
points(allspiderspcoa1$points[mygalo,1],allspiderspcoa1$points[mygalo,3], cex=0.9,pch=mygalopch)
points(allspiderspcoa1$points[synsperm,1],allspiderspcoa1$points[synsperm,3], cex=0.9,pch=synspermpch,col="#D55E00")
points(allspiderspcoa1$points[araneio,1],allspiderspcoa1$points[araneio,3], cex=0.9,pch=araneiopch,col="#0072b2")
points(allspiderspcoa1$points[entely,1],allspiderspcoa1$points[entely,3], cex=0.9,col="#009e73",pch=entelypch)
points(allspiderspcoa1$points[marro,1],allspiderspcoa1$points[marro,3], cex=0.9,pch=marropch,col="#e79f00")
points(allspiderspcoa1$points[calamis,1],allspiderspcoa1$points[calamis,3], cex=0.9,pch=calamispch,col="#cc79a7")
points(allspiderspcoa1$points[dionycha,1],allspiderspcoa1$points[dionycha,3], cex=0.9,pch=dionychapch,col="#9ad0f3")
legend(-13,-1.75,c("Mygalomorphae","Synspermiata","Araneiodea", "Other Entelygynae","Marronoid clade","Oval Calamistrum","Dionycha"), col=c("#000000","#D55E00","#0072b2","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=16,cex=0.9,pt.cex=0.9,text.font=3)
dev.off()

png('Fig 02B Phylo Morphospace.png', width=800, height=800)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,3], cex=0,xlab="Relative size(Axis 1)", ylab="Proximal/distal leg ratio(Axis 3)")
points(allspiderspcoa1$points[mygalo,1],allspiderspcoa1$points[mygalo,3], cex=0.9,pch=16)
points(allspiderspcoa1$points[synsperm,1],allspiderspcoa1$points[synsperm,3], cex=0.9,pch=16,col="#D55E00")
points(allspiderspcoa1$points[araneio,1],allspiderspcoa1$points[araneio,3], cex=0.9,pch=16,col="#0072b2")
points(allspiderspcoa1$points[entely,1],allspiderspcoa1$points[entely,3], cex=0.9,col="#009e73",pch=16)
points(allspiderspcoa1$points[marro,1],allspiderspcoa1$points[marro,3], cex=0.9,pch=16,col="#e79f00")
points(allspiderspcoa1$points[calamis,1],allspiderspcoa1$points[calamis,3], cex=0.9,pch=16,col="#cc79a7")
points(allspiderspcoa1$points[dionycha,1],allspiderspcoa1$points[dionycha,3], cex=0.9,pch=16,col="#9ad0f3")
legend(-13,-1.75,c("Mygalomorph","Synspermiata","Araneoidea", "Other Entelygynae","Marronoid clade","Oval Calamistrum","Dionycha"), col=c("#000000","#D55E00","#0072b2","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=16,cex=0.9,pt.cex=0.9,text.font=3)
dev.off()

png('Fig 03 Hunt Morphospace.png', width=800, height=800)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,2], cex=0,xlab="Relative size", ylab="Body/Leg length ratio")
points(allspiderspcoa1$points[ambushers,1],allspiderspcoa1$points[ambushers,2], cex=0.9,pch=ambushpch)
points(allspiderspcoa1$points[hunters,1],allspiderspcoa1$points[hunters,2], cex=0.9,pch=16,col="#D55E00")
points(allspiderspcoa1$points[orbweaver,1],allspiderspcoa1$points[orbweaver,2], cex=0.9,pch=16,col="#0072b2")
points(allspiderspcoa1$points[sedentary,1],allspiderspcoa1$points[sedentary,2], cex=0.9,col="#009e73",pch=sedentpch)
points(allspiderspcoa1$points[irregularweb,1],allspiderspcoa1$points[irregularweb,2], cex=0.9,pch=16,col="#e79f00")
points(allspiderspcoa1$points[sheetwebs,1],allspiderspcoa1$points[sheetwebs,2], cex=0.9,pch=16,col="#cc79a7")
points(allspiderspcoa1$points[mimic,1],allspiderspcoa1$points[mimic,2], cex=0.9,pch=16,col="#9ad0f3")
legend(-13,-2.7,c("Ambush","Ambush silk","Ground hunter","Orb web","Burrow","Retreat","Irregular web","Sheet web","Mimicry"), col=c("#000000","#000000","#D55E00","#0072b2","#009e73","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=c(17,16,16,16,17,16,16,16,16),cex=1,pt.cex=1,text.font=3)
dev.off()

png('Fig 03B Hunt Morphospace.png', width=800, height=800)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,2], cex=0,xlab="Relative size(Axis 1)", ylab="Body/Leg length ratio(Axis 2)")
points(allspiderspcoa1$points[ambushers,1],allspiderspcoa1$points[ambushers,2], cex=0.9,pch=ambushpch)
points(allspiderspcoa1$points[hunters,1],allspiderspcoa1$points[hunters,2], cex=0.9,pch=16,col="#D55E00")
points(allspiderspcoa1$points[orbweaver,1],allspiderspcoa1$points[orbweaver,2], cex=0.9,pch=16,col="#0072b2")
points(allspiderspcoa1$points[sedentary,1],allspiderspcoa1$points[sedentary,2], cex=0.9,col="#009e73",pch=sedentpch)
points(allspiderspcoa1$points[irregularweb,1],allspiderspcoa1$points[irregularweb,2], cex=0.9,pch=16,col="#e79f00")
points(allspiderspcoa1$points[sheetwebs,1],allspiderspcoa1$points[sheetwebs,2], cex=0.9,pch=16,col="#cc79a7")
points(allspiderspcoa1$points[mimic,1],allspiderspcoa1$points[mimic,2], cex=0.9,pch=16,col="#9ad0f3")
legend(-13,-2.7,c("Ambush","Ambush silk","Ground hunter","Orb web","Burrow","Retreat","Irregular web","Sheet web","Mimicry"), col=c("#000000","#000000","#D55E00","#0072b2","#009e73","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=c(17,16,16,16,17,16,16,16,16),cex=1,pt.cex=1,text.font=3)
dev.off()

png('Fig 04 Hunt Morphospace.png', width=800, height=800)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,3], cex=0,xlab="Relative size", ylab="Proximal/Distal leg ratio")
points(allspiderspcoa1$points[ambushers,1],allspiderspcoa1$points[ambushers,3], cex=0.9,pch=ambushpch)
points(allspiderspcoa1$points[hunters,1],allspiderspcoa1$points[hunters,3], cex=0.9,pch=16,col="#D55E00")
points(allspiderspcoa1$points[orbweaver,1],allspiderspcoa1$points[orbweaver,3], cex=0.9,pch=16,col="#0072b2")
points(allspiderspcoa1$points[sedentary,1],allspiderspcoa1$points[sedentary,3], cex=0.9,col="#009e73",pch=sedentpch)
points(allspiderspcoa1$points[irregularweb,1],allspiderspcoa1$points[irregularweb,3], cex=0.9,pch=16,col="#e79f00")
points(allspiderspcoa1$points[sheetwebs,1],allspiderspcoa1$points[sheetwebs,3], cex=0.9,pch=16,col="#cc79a7")
points(allspiderspcoa1$points[mimic,1],allspiderspcoa1$points[mimic,3], cex=0.9,pch=16,col="#9ad0f3")
legend(-13,-1.7,c("Ambush","Ambush silk","Ground hunter","Orb web","Burrow","Retreat","Irregular web","Sheet web","Mimicry"), col=c("#000000","#000000","#D55E00","#0072b2","#009e73","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=c(17,16,16,16,17,16,16,16,16),cex=0.8,pt.cex=1,text.font=3)
dev.off()

png('Fig 04B Hunt Morphospace.png', width=800, height=800)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,3], cex=0,xlab="Relative size(Axis 1)", ylab="Proximal/Distal leg ratio(Axis 3)")
points(allspiderspcoa1$points[ambushers,1],allspiderspcoa1$points[ambushers,3], cex=0.9,pch=ambushpch)
points(allspiderspcoa1$points[hunters,1],allspiderspcoa1$points[hunters,3], cex=0.9,pch=16,col="#D55E00")
points(allspiderspcoa1$points[orbweaver,1],allspiderspcoa1$points[orbweaver,3], cex=0.9,pch=16,col="#0072b2")
points(allspiderspcoa1$points[sedentary,1],allspiderspcoa1$points[sedentary,3], cex=0.9,col="#009e73",pch=sedentpch)
points(allspiderspcoa1$points[irregularweb,1],allspiderspcoa1$points[irregularweb,3], cex=0.9,pch=16,col="#e79f00")
points(allspiderspcoa1$points[sheetwebs,1],allspiderspcoa1$points[sheetwebs,3], cex=0.9,pch=16,col="#cc79a7")
points(allspiderspcoa1$points[mimic,1],allspiderspcoa1$points[mimic,3], cex=0.9,pch=16,col="#9ad0f3")
legend(-13,-1.7,c("Ambush","Ambush silk","Ground hunter","Orb web","Burrow","Retreat","Irregular web","Sheet web","Mimicry"), col=c("#000000","#000000","#D55E00","#0072b2","#009e73","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=c(17,16,16,16,17,16,16,16,16),cex=0.8,pt.cex=1,text.font=3)
dev.off()

png('Fig 05A Phylo Cluster.png', width=800, height=1100)
par( mgp=c(3,2,1), cex=1.3, cex.lab=1.4, cex.axis=1.2,las=1, oma=c(0,2,0,0),mar=c(5,4.5,4,1),mfrow=c(2,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,2], cex=0,xlab="Relative size", ylab="Body/Leg length ratio",xlim=c(-18,18),ylim=c(-7,5))
points(allspiderspcoa1$points[mygalo,1],allspiderspcoa1$points[mygalo,2], cex=0.9,pch=mygalopch)
points(allspiderspcoa1$points[synsperm,1],allspiderspcoa1$points[synsperm,2], cex=0.9,pch=synspermpch,col="#D55E00")
points(allspiderspcoa1$points[araneio,1],allspiderspcoa1$points[araneio,2], cex=0.9,pch=araneiopch,col="#0072b2")
points(allspiderspcoa1$points[entely,1],allspiderspcoa1$points[entely,2], cex=0.9,col="#009e73",pch=entelypch)
points(allspiderspcoa1$points[marro,1],allspiderspcoa1$points[marro,2], cex=0.9,pch=marropch,col="#e79f00")
points(allspiderspcoa1$points[calamis,1],allspiderspcoa1$points[calamis,2], cex=0.9,pch=calamispch,col="#cc79a7")
points(allspiderspcoa1$points[dionycha,1],allspiderspcoa1$points[dionycha,2], cex=0.9,pch=dionychapch,col="#9ad0f3")
text(-15, 4, labels="a",cex=2.2)
legend(-13,-2.8,c("Mygalomorphae","Synspermiata","Araneiodea", "Other Entelygynae","Marronoid clade","Oval Calamistrum","Dionycha"), col=c("#000000","#D55E00","#0072b2","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=16,cex=1,pt.cex=1,text.font=3)
dataEllipse(allspiderspcoa1$points[pamspidfit1$clustering==1,1],allspiderspcoa1$points[pamspidfit1$clustering==1,2], add = T, levels = 0.95,center.cex = 0,col = "#000000", fill.alpha = 0.2,fill = T, plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit1$clustering==2,1],allspiderspcoa1$points[pamspidfit1$clustering==2,2], add = T,levels = 0.95,center.cex = 0,col = "#D55E00",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit1$clustering==3,1],allspiderspcoa1$points[pamspidfit1$clustering==3,2], add = T,levels = 0.95,center.cex = 0,col = "#0072b2",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit1$clustering==4,1],allspiderspcoa1$points[pamspidfit1$clustering==4,2], add = T,levels = 0.95,center.cex = 0,col = "#009e73",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit1$clustering==5,1],allspiderspcoa1$points[pamspidfit1$clustering==5,2], add = T,levels = 0.95,center.cex = 0,col = "#e79f00",fill.alpha = 0.2,fill = T,plot.points = F)
plot(adjspiderpcoa$points[,1],adjspiderpcoa$points[,2], cex=0,xlab="Body/Leg length ratio", ylab="Proximal/distal leg ratio",xlim=c(-15,20),ylim=c(-10,10))
points(adjspiderpcoa$points[mygalo,1],adjspiderpcoa$points[mygalo,2], cex=0.9,pch=mygalopch)
points(adjspiderpcoa$points[synsperm,1],adjspiderpcoa$points[synsperm,2], cex=0.9,pch=synspermpch,col="#D55E00")
points(adjspiderpcoa$points[araneio,1],adjspiderpcoa$points[araneio,2], cex=0.9,pch=araneiopch,col="#0072b2")
points(adjspiderpcoa$points[entely,1],adjspiderpcoa$points[entely,2], cex=0.9,col="#009e73",pch=entelypch)
points(adjspiderpcoa$points[marro,1],adjspiderpcoa$points[marro,2], cex=0.9,pch=marropch,col="#e79f00")
points(adjspiderpcoa$points[calamis,1],adjspiderpcoa$points[calamis,2], cex=0.9,pch=calamispch,col="#cc79a7")
points(adjspiderpcoa$points[dionycha,1],adjspiderpcoa$points[dionycha,2], cex=0.9,pch=dionychapch,col="#9ad0f3")
text(-12, 7, labels="b",cex=2.2)
legend(-13,-4.5,c("Mygalomorphae","Synspermiata","Araneiodea", "Other Entelygynae","Marronoid clade","Oval Calamistrum","Dionycha"), col=c("#000000","#D55E00","#0072b2","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=16,cex=1,pt.cex=1,text.font=3)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==1,1],adjspiderpcoa$points[adjspidfit$clustering==1,2], add = T, levels = 0.95,center.cex = 0,col = "#000000", fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==2,1],adjspiderpcoa$points[adjspidfit$clustering==2,2], add = T,levels = 0.95,center.cex = 0,col = "#D55E00",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==3,1],adjspiderpcoa$points[adjspidfit$clustering==3,2], add = T,levels = 0.95,center.cex = 0,col = "#0072b2",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==4,1],adjspiderpcoa$points[adjspidfit$clustering==4,2], add = T,levels = 0.95,center.cex = 0,col = "#009e73",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==5,1],adjspiderpcoa$points[adjspidfit$clustering==5,2], add = T,levels = 0.95,center.cex = 0,col = "#e79f00",fill.alpha = 0.2,fill = T,plot.points = F)
dev.off()




png('Fig 05B Phylo Cluster.png', width=700, height=800)
par( mgp=c(3,2,1), cex=1.3, cex.lab=1.4, cex.axis=1.2,las=1, oma=c(0,2,0,0),mar=c(5,4.5,4,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,3], cex=0,xlab="Relative size", ylab="Proximal/distal leg ratio", xlim=c(-18,14),ylim=c(-7,5))
points(allspiderspcoa1$points[mygalo,1],allspiderspcoa1$points[mygalo,3], cex=0.9,pch=mygalopch)
points(allspiderspcoa1$points[synsperm,1],allspiderspcoa1$points[synsperm,3], cex=0.9,pch=synspermpch,col="#D55E00")
points(allspiderspcoa1$points[araneio,1],allspiderspcoa1$points[araneio,3], cex=0.9,pch=araneiopch,col="#0072b2")
points(allspiderspcoa1$points[entely,1],allspiderspcoa1$points[entely,3], cex=0.9,col="#009e73",pch=entelypch)
points(allspiderspcoa1$points[marro,1],allspiderspcoa1$points[marro,3], cex=0.9,pch=marropch,col="#e79f00")
points(allspiderspcoa1$points[calamis,1],allspiderspcoa1$points[calamis,3], cex=0.9,pch=calamispch,col="#cc79a7")
points(allspiderspcoa1$points[dionycha,1],allspiderspcoa1$points[dionycha,3], cex=0.9,pch=dionychapch,col="#9ad0f3")
legend(-13,-3.75,c("Mygalomorphae","Synspermiata","Araneoidea", "Other Entelygynae","Marronoid clade","Oval Calamistrum","Dionycha"), col=c("#000000","#D55E00","#0072b2","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=16,cex=0.9,pt.cex=0.9,text.font=3)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==1,1],allspiderspcoa1$points[pamspidfit$clustering==1,3], add = T, levels = 0.95,center.cex = 0,col = "#000000", fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==2,1],allspiderspcoa1$points[pamspidfit$clustering==2,3], add = T,levels = 0.95,center.cex = 0,col = "#D55E00",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==3,1],allspiderspcoa1$points[pamspidfit$clustering==3,3], add = T,levels = 0.95,center.cex = 0,col = "#0072b2",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==4,1],allspiderspcoa1$points[pamspidfit$clustering==4,3], add = T,levels = 0.95,center.cex = 0,col = "#009e73",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==5,1],allspiderspcoa1$points[pamspidfit$clustering==5,3], add = T,levels = 0.95,center.cex = 0,col = "#e79f00",fill.alpha = 0.2,fill = T,plot.points = F)
dev.off()

png('Fig 06A Hunt Cluster.png', width=800, height=1100)
par( mgp=c(3,2,1), cex=1.3, cex.lab=1.4, cex.axis=1.2,las=1, oma=c(0,2,0,0),mar=c(5,4.5,4,1),mfrow=c(2,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,2], cex=0,xlab="Relative size", ylab="Body/Leg length ratio",xlim=c(-18,18),ylim=c(-7,5))
points(allspiderspcoa1$points[ambushers,1],allspiderspcoa1$points[ambushers,2], cex=0.9,pch=ambushpch)
points(allspiderspcoa1$points[hunters,1],allspiderspcoa1$points[hunters,2], cex=0.9,pch=16,col="#D55E00")
points(allspiderspcoa1$points[orbweaver,1],allspiderspcoa1$points[orbweaver,2], cex=0.9,pch=16,col="#0072b2")
points(allspiderspcoa1$points[sedentary,1],allspiderspcoa1$points[sedentary,2], cex=0.9,col="#009e73",pch=sedentpch)
points(allspiderspcoa1$points[irregularweb,1],allspiderspcoa1$points[irregularweb,2], cex=0.9,pch=16,col="#e79f00")
points(allspiderspcoa1$points[sheetwebs,1],allspiderspcoa1$points[sheetwebs,2], cex=0.9,pch=16,col="#cc79a7")
points(allspiderspcoa1$points[mimic,1],allspiderspcoa1$points[mimic,2], cex=0.9,pch=16,col="#9ad0f3")
text(-17, 4, labels="a",cex=2.2)
legend(-13,-2.7,c("Ambush","Ambush silk","Ground hunter","Orb web","Burrow","Retreat","Irregular web","Sheet web","Mimicry"), col=c("#000000","#000000","#D55E00","#0072b2","#009e73","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=c(17,16,16,16,17,16,16,16,16),cex=1,pt.cex=1,text.font=3)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==1,1],allspiderspcoa1$points[pamspidfit$clustering==1,2], add = T, levels = 0.95,center.cex = 0,col = "#000000", fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==2,1],allspiderspcoa1$points[pamspidfit$clustering==2,2], add = T,levels = 0.95,center.cex = 0,col = "#D55E00",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==3,1],allspiderspcoa1$points[pamspidfit$clustering==3,2], add = T,levels = 0.95,center.cex = 0,col = "#0072b2",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==4,1],allspiderspcoa1$points[pamspidfit$clustering==4,2], add = T,levels = 0.95,center.cex = 0,col = "#009e73",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==5,1],allspiderspcoa1$points[pamspidfit$clustering==5,2], add = T,levels = 0.95,center.cex = 0,col = "#e79f00",fill.alpha = 0.2,fill = T,plot.points = F)
plot(adjspiderpcoa$points[,1],adjspiderpcoa$points[,2], cex=0,xlab="Body/Leg length ratio", ylab="Proximal/Distal leg ratio")
points(adjspiderpcoa$points[ambushers,1],adjspiderpcoa$points[ambushers,2], cex=0.9,pch=ambushpch)
points(adjspiderpcoa$points[hunters,1],adjspiderpcoa$points[hunters,2], cex=0.9,pch=16,col="#D55E00")
points(adjspiderpcoa$points[orbweaver,1],adjspiderpcoa$points[orbweaver,2], cex=0.9,pch=16,col="#0072b2")
points(adjspiderpcoa$points[sedentary,1],adjspiderpcoa$points[sedentary,2], cex=0.9,col="#009e73",pch=sedentpch)
points(adjspiderpcoa$points[irregularweb,1],adjspiderpcoa$points[irregularweb,2], cex=0.9,pch=16,col="#e79f00")
points(adjspiderpcoa$points[sheetwebs,1],adjspiderpcoa$points[sheetwebs,2], cex=0.9,pch=16,col="#cc79a7")
points(adjspiderpcoa$points[mimic,1],adjspiderpcoa$points[mimic,2], cex=0.9,pch=16,col="#9ad0f3")
text(-8, 6, labels="b",cex=2.2)
legend(14,-2,c("Ambush","Ambush silk","Ground hunter","Orb web","Burrow","Retreat","Irregular web","Sheet web","Mimicry"), col=c("#000000","#000000","#D55E00","#0072b2","#009e73","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=c(17,16,16,16,17,16,16,16,16),cex=0.8,pt.cex=1,text.font=3)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==1,1],adjspiderpcoa$points[adjspidfit$clustering==1,2], add = T, levels = 0.95,center.cex = 0,col = "#000000", fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==2,1],adjspiderpcoa$points[adjspidfit$clustering==2,2], add = T,levels = 0.95,center.cex = 0,col = "#D55E00",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==3,1],adjspiderpcoa$points[adjspidfit$clustering==3,2], add = T,levels = 0.95,center.cex = 0,col = "#0072b2",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==4,1],adjspiderpcoa$points[adjspidfit$clustering==4,2], add = T,levels = 0.95,center.cex = 0,col = "#009e73",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(adjspiderpcoa$points[adjspidfit$clustering==5,1],adjspiderpcoa$points[adjspidfit$clustering==5,2], add = T,levels = 0.95,center.cex = 0,col = "#e79f00",fill.alpha = 0.2,fill = T,plot.points = F)
dev.off()
png('Fig 06B Hunt Cluster.png', width=800, height=800)
par( mgp=c(3,2,1), cex=1.3, cex.lab=1.4, cex.axis=1.2,las=1, oma=c(0,2,0,0),mar=c(5,4.5,4,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,3], cex=0,xlab="Relative size", ylab="Proximal/distal leg ratio", xlim=c(-18,15),ylim=c(-8,6))
points(allspiderspcoa1$points[ambushers,1],allspiderspcoa1$points[ambushers,3], cex=0.9,pch=ambushpch)
points(allspiderspcoa1$points[hunters,1],allspiderspcoa1$points[hunters,3], cex=0.9,pch=16,col="#D55E00")
points(allspiderspcoa1$points[orbweaver,1],allspiderspcoa1$points[orbweaver,3], cex=0.9,pch=16,col="#0072b2")
points(allspiderspcoa1$points[sedentary,1],allspiderspcoa1$points[sedentary,3], cex=0.9,col="#009e73",pch=sedentpch)
points(allspiderspcoa1$points[irregularweb,1],allspiderspcoa1$points[irregularweb,3], cex=0.9,pch=16,col="#e79f00")
points(allspiderspcoa1$points[sheetwebs,1],allspiderspcoa1$points[sheetwebs,3], cex=0.9,pch=16,col="#cc79a7")
points(allspiderspcoa1$points[mimic,1],allspiderspcoa1$points[mimic,3], cex=0.9,pch=16,col="#9ad0f3")
legend(-15,-3.5,c("Ambush","Ambush silk","Ground hunter","Orb web","Burrow","Retreat","Irregular web","Sheet web","Mimicry"), col=c("#000000","#000000","#D55E00","#0072b2","#009e73","#009e73","#e79f00","#cc79a7","#9ad0f3"),pch=c(17,16,16,16,17,16,16,16,16),cex=1,pt.cex=1,text.font=3)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==1,1],allspiderspcoa1$points[pamspidfit$clustering==1,3], add = T, levels = 0.95,center.cex = 0,col = "#000000", fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==2,1],allspiderspcoa1$points[pamspidfit$clustering==2,3], add = T,levels = 0.95,center.cex = 0,col = "#D55E00",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==3,1],allspiderspcoa1$points[pamspidfit$clustering==3,3], add = T,levels = 0.95,center.cex = 0,col = "#0072b2",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==4,1],allspiderspcoa1$points[pamspidfit$clustering==4,3], add = T,levels = 0.95,center.cex = 0,col = "#009e73",fill.alpha = 0.2,fill = T,plot.points = F)
dataEllipse(allspiderspcoa1$points[pamspidfit$clustering==5,1],allspiderspcoa1$points[pamspidfit$clustering==5,3], add = T,levels = 0.95,center.cex = 0,col = "#e79f00",fill.alpha = 0.2,fill = T,plot.points = F)
dev.off()