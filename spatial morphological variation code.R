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
#environmental gradients
sites=matrix(nrow=39, ncol=3)
sites[,1]=sitenames 
sites[,2]=sitelatclim
sites[,3]=sitelongclim 
library(raster)
library(sp)
totalclim <-getData("worldclim",var="bio",res=10)
sitecoord <- SpatialPoints(cbind(as.numeric(sites[,3]),as.numeric(sites[,2])))
eastclim <- totalclim[[c(1,12)]]
names(eastclim) <- c("Temp","Prec")
sitevalues <- extract(eastclim,sitecoord)
eastdry <- totalclim[[c(15)]]
sitedry<- extract(eastdry,sitecoord)
eastalt <-getData('alt',country= 'AUS')
sitealt <-extract(eastalt,sitecoord)
siteseason= totalclim [[c(4)]]
siteseasonvalues= extract(siteseason,sitecoord)
eastquarter= totalclim[[c(17)]]
sitequarter= extract(eastquarter,sitecoord)
sitevalues=cbind(sitevalues,sitedry,sitequarter,siteseasonvalues,sitealt,sitelats,sitelongs)
row.names(sitevalues)=sitenames
colnames(sitevalues)= c("Temp","Prec","Prec Season","Prec Dry Quarter","Temp Season","Elevation","Lat","Long")
#creation of mean values for each site
library(psych)
sitemean1=array()
sitemeanlist=list()
sitemeancommon=list()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]),1]
  sitemean1[i]<- geometric.mean(sitescores, na.rm=T)
  sitemeanlist[[i]]=sitescores
  sitemeancommon[[i]]=sitespecies[which(sitespecies%in%AllSpiders[,2])]
}

sitemed1=array()
for (i in 1:20) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]),1]
  sitemed1[i]<- median(sitescores, na.rm=T)
}

sitesd1=array()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]),1]
  sitesd1[i]<- sd(sitescores, na.rm=T)
}
sitemed2=array()
for (i in 1:20) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]),2]
  sitemed2[i]<- median(sitescores, na.rm=T)
}

sitemean2=array()
sitemean2list=list()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]),2]
  sitemean2[i]<- geometric.mean(sitescores, na.rm=T)
  sitemean2list[[i]] = sitescores
}
sitemed3=array()
for (i in 1:20) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]),3]
  sitemed3[i]<- median(sitescores, na.rm=T)
}

sitesd2=array()
for (i in 1:20) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]),2]
  sitesd2[i]<- sd(sitescores, na.rm=T)
}
sitesd3=array()
for (i in 1:20) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]),3]
  sitesd3[i]<- sd(sitescores, na.rm=T)
}
sitemean3=array()
sitemean3list=list()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]),3]
  sitemean3[i]<- geometric.mean(sitescores, na.rm=T)
  sitemean3list[[i]] = sitescores
}
#araneidae and salticidae only scores
sitemeanaran1=array()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]&AllSpiders[,1]=="Araneidae"),1]
  sitemeanaran1[i]<- mean(sitescores, na.rm=T)
}

sitemeansalt1=array()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]&AllSpiders[,1]=="Salticidae"),1]
  sitemeansalt1[i]<- mean(sitescores, na.rm=T)
}

sitemeanaran2=array()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]&AllSpiders[,1]=="Araneidae"),2]
  sitemeanaran2[i]<- mean(sitescores, na.rm=T)
}

sitemeanaran3=array()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]&AllSpiders[,1]=="Araneidae"),3]
  sitemeanaran3[i]<- mean(sitescores, na.rm=T)
}

sitemeansalt2=array()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]&AllSpiders[,1]=="Salticidae"),2]
  sitemeansalt2[i]<- geometric.mean(sitescores, na.rm=T)
}
sitemeansalt3=array()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- allspiderspcoa1$points[which(sitespecies%in%AllSpiders[,2]&AllSpiders[,1]=="Salticidae"),3]
  sitemeansalt3[i]<- geometric.mean(sitescores, na.rm=T)
}
sitemeanadj1=array()
sitemeanadj1list=list()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- adjspiderpcoa$points[which(sitespecies%in%AllSpiders[,2]),1]
  sitemeanadj1[i]<- geometric.mean(sitescores, na.rm=T)
  sitemeanadj1list[[i]] = sitescores
}
sitemeanadj2=array()
sitemeanadj2list=list()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- adjspiderpcoa$points[which(sitespecies%in%AllSpiders[,2]),2]
  sitemeanadj2[i]<- geometric.mean(sitescores, na.rm=T)
  sitemeanadj2list[[i]] = sitescores
}
sitemeanadj3=array()
sitemeanadj3list=list()
for (i in 1:39) {
  sitespecies<- unlist(specieslist[i])
  sitescores <- adjspiderpcoa$points[which(sitespecies%in%AllSpiders[,2]),3]
  sitemeanadj3[i]<- geometric.mean(sitescores, na.rm=T)
  sitemeanadj3list[[i]] = sitescores
}
#spatial autoregression
library(raster)
library(spdep)
library(spatialreg)

regcoord= cbind(as.numeric(sites[,3]),as.numeric(sites[,2]))

sitereg = which(! duplicated(regcoord))

nb = knn2nb(knearneigh(regcoord[sitereg,],k=2,longlat=T))
nb2 = knn2nb(knearneigh(regcoord,k=2,longlat=T))
summary(spautolm(sitemean1[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemean2[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemean3[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
#medians
summary(spautolm(sitemed1[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemed2[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemed3[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
#standard deviation
summary(spautolm(sitesd1[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitesd2[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitesd3[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
#araneidae only
summary(spautolm(sitemeanaran1[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemeanaran2[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemeanaran3[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
#salticidae only
summary(spautolm(sitemeansalt1[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemeansalt2[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemeansalt3[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
#size adjusted#
summary(spautolm(sitemeanadj1[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemeanadj2[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemeanadj3[sitereg] ~ sitevalues[sitereg,1:6],listw=nb2listw(nb)),Nagelkerke=T)
#site morphological richness & dispersion
summary(spautolm(sitemorphrich[sitereg] ~ sitevalues[sitereg,1:6]+ sitecounts,listw=nb2listw(nb)),Nagelkerke=T)
summary(spautolm(sitemorphdisp[sitereg] ~ sitevalues[sitereg,1:6]+ sitecounts,listw=nb2listw(nb)),Nagelkerke=T)
#figures
png('Fig 01 Total Morphospace.png', width=600, height=500)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,2], cex=0.8, pch= 16, xlab='Axis 1', ylab='Axis 2')
dev.off()
png('Fig 02 Total Morphospace.png', width=600, height=500)
par( mgp=c(3,1,0), cex=1.3, cex.lab=2, cex.axis=1.8,las=1, oma=c(0,2,0,0),mar=c(4,4.5,0.5,1))
plot(allspiderspcoa1$points[,1],allspiderspcoa1$points[,3], cex=0.8, pch= 16, xlab='Axis 1',ylab='Axis 3')
dev.off()