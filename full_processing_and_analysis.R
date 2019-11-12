# #################################################################### ####
# Title: Comparing the recovery of richness, structure and biomass in  ####
#        naturally regrowing and planted reforestation                 #### 
# Author: Timothy L Staples                                            ####
# Details: Full analysis (data not available due to data agreements)   ####  
# #################################################################### ####
# Libraries ####
rm(list=ls())
setwd("/home/timothy/Dropbox/Tim/PhD/Data/Chapter IV - regrowth")
library(raster) # for environmental rasters
library(rgdal) # for environmental rasters
library(maptools) # manipulating environmental rasters and plotting
library(usdm) # for variance-inflation factor comparison
library(vegan) # for rarefied richness etc
library(nlme) # mixed-effect modelling
library(lme4) # mixed-effect modelling
library(MuMIn) # model selection and averaging
library(gamm4) # additive modelling
library(parallel) # for parallel computation of some apply and model selection functions
library(FD) # to calculate functional diversity indices
library(fields) # for k-means clustering
library(multcomp) # for glht
library(geoR) # for variogram
library(plotrix) # for gradient rectangles in plots
library(rgeos) # for mapping
library(dismo) # for circles in long/lat coords
library(gdistance) # for getting distance between spatial objects
library(geosphere) # area of spatial polygons
library(vioplot) # for violin plots
library(car) # to logit transform proportions
library(gamm4) # for mixed-effect gam
library(funrar) # functional rarity package
library(dendextend) # for making adjustments to dendrograms
library(lsmeans)
library(grImport)
library(merTools)

# Functions ####
gc.dist<-function(lat1, lon1, lat2, lon2){
  
  # custom function to turn degrees to radians
  deg2rad <- function(deg) {(deg * pi) / (180)} 
  R<-6371 # Radius of the earth in km
  dLat<-deg2rad(lat2-lat1) # deg2rad below
  dLon<-deg2rad(lon2-lon1); 
  a<-sin(dLat/2) * sin(dLat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * 
    sin(dLon/2) * sin(dLon/2)
  
  c = 2 * atan2(sqrt(a), sqrt(1-a)); 
  d = R * c # Distance in km
  return(d)
} 

sapply(list.files(path="./Functions", pattern=".R", full.names=TRUE),
       source)

# 0. IMPORT DATA ####

plant<-read.csv("./Data/plantscale.csv", header=TRUE)
site<-read.csv("./Data/plotscale.csv", header=TRUE)
species<-read.csv("./Data/species.list.csv", header=TRUE)

head(site)

#             EXTRACT ENVIRONMENTAL CONDITIONS AT PLOTS ####

shape.file.dir<-paste(ifelse(Sys.info()['sysname']=="Linux", "/home/Storage HDD/", "F:/"),
                      "University files/Shape files", sep="")

data.coords<-as.data.frame(cbind(site$long, site$lat))
colnames(data.coords)<-c("Longitude","Latitude")

point.values<-sapply(list.files(path=shape.file.dir,
                                pattern=".tif"), function(x){
                                  extract(raster(paste0(shape.file.dir,"/",x)), data.coords, method="simple")
                                })

colnames(point.values)<-substr(colnames(point.values), 
                               1, regexpr(".tif", colnames(point.values))-1)

site<-cbind(site, point.values)

pre.RE<-readShapeSpatial("/home/Storage HDD/University files/Shape files/Shape files/Queensland regional ecosystems/Biodiversity status of pre-clearing and remnant regional ecosystems - South East Qld/sth_east/pre.shp")
cur.RE<-readShapeSpatial("/home/Storage HDD/University files/Shape files/Shape files/Queensland regional ecosystems/Vegetation management regional ecosystem and remnant map - version 8.0 non-coastal/Cropped RE.shp")

IBRA<-readShapeSpatial(paste0(shape.file.dir,
                              "/Shape files/",
                              "IBRA7_regions/ibra7_regions.shp"))
IBRAsub<-readShapeSpatial(paste0(shape.file.dir,
                                 "/Shape files/",
                                 "IBRA7_subregions_states/IBRA7_subregions_states.shp"))

data.coords.nona<-data.coords[rowSums(is.na(data.coords))==0,]
data.coords.plp<-site$pl.p[rowSums(is.na(data.coords))==0]
coordinates(data.coords.nona)<-c("Longitude","Latitude")

RE.points<-data.frame(curRE=over(data.coords.nona, cur.RE)$LANDZONE,
                      preRE=over(data.coords.nona, pre.RE)$LANDZONE,
                      pl.p=site$pl.p[rowSums(is.na(data.coords))==0])

IBRA.points<-data.frame(IBRA=over(data.coords.nona,IBRA)$REG_CODE_7,
                        IBRAsub=over(data.coords.nona,IBRAsub)$SUB_CODE_7,
                        pl.p=site$pl.p[rowSums(is.na(data.coords))==0])

temp<-data.coords.nona[IBRA.points$IBRA %in% c("SEQ","BBS"),]


site1<-merge(site, IBRA.points, all.x=TRUE, all.y=FALSE, by.x="pl.p", by.y="pl.p")
site1<-merge(site1, RE.points, all.x=TRUE, all.y=FALSE, by.x="pl.p", by.y="pl.p")

site1<-cbind(pl.p=site1$pl.p,
             site1[,c("IBRA", "IBRAsub")],
             site1[,!colnames(site1) %in% c("pl.p","IBRA", "IBRAsub")])

site<-site1

# 1. FIND STUDY SITES ####

table(site$IBRA, site$plot.type)

# The only two areas with decent regrowth values aer in BBS, SEQ and maybe MDD
# We have fairly good planting and remnant representation in the region too.

# let's cut down to these three groups and look at sub-regions
reg.site<-droplevels(site[site$IBRA %in% c("BBS", "SEQ"),])

reg.site.count<-as.matrix(table(droplevels(reg.site$IBRAsub), reg.site$plot.type))

reg.site.count<-reg.site.count[reg.site.count[,"Natural regrowth"] > 0,]
colSums(reg.site.count)

# SEQ06 is a good option, as is BS17. Maybe MDD05. We might want to look at these
# plantings and build up a map of their distances, to see if they fall into a couple of
# clusters
ausmap<-readShapeSpatial("/home/Storage HDD/University files/Shape files/Shape files/IBRA7_regions/ibra7_regions.shp")

with(reg.site,
     plot(lat ~ long, xlim=c(151,153), ylim=c(-30,-23),
          col=plot.type))
plot(ausmap, add=TRUE)

# Looks like 2 clusters. One in SA, and one in SEQ. Let's get each one and look for
# distances. Cluster 1 is between 150 and 153 longitude, cluster 2 is less than 142.
cluster1<-reg.site[reg.site$long >= 150 & reg.site$lat > -28,]
dim(cluster1)

dim(combn(c(1:86), 2))

cluster1.dist<-apply(combn(c(1:dim(cluster1)[1]), 2), 2, function(x){
  gc.dist(lat1=cluster1$lat[x[1]],
          lon1=cluster1$long[x[1]],
          lat2=cluster1$lat[x[2]],
          lon2=cluster1$long[x[2]])
})

hist(cluster1.dist)
summary(cluster1.dist)

# So the maximum distance we have is ~300km. Let's see what proportion of plots
# we exclude using a circle extending out a certain radius from the coordinate means
# turn degrees into radians
deg2rad <- function(deg) {(deg * pi) / (180)}

dist.to.mean<-gc.dist(lat1=mean(cluster1$lat),
                      lon1=mean(cluster1$long),
                      lat2=cluster1$lat,
                      lon2=cluster1$long)

dist.diff<-sapply(1:300, function(x){
  sum(dist.to.mean <= x) / length(dist.to.mean)
})

dist.diff[123]

# Let's set it to 125km, that gives us almost all of the sites in this region,
# and is a reasonable distance.

sub.site<-cluster1[dist.to.mean <=125,]
table(sub.site$plot.type)

# and we have almost equal sizes of regrowth, planting and remnant data

table(sub.site$plot.type, sub.site$age)
rowSums(table(sub.site$planting, 
              paste0(sub.site$long, 
                     sub.site$lat))) / table(sub.site$planting)

site<-sub.site
plant<-droplevels(plant[plant$pl.p %in% site$pl.p, ])

# 2. DATA PREP ####

#             CALCULATE BASIC STATS (DENSITY, RICHNESS) ####

# biomass
target.plots<-site$pl.p[site$plot.type %in% c("Planting", "Natural regrowth")]

# some regrowth sites have old pasture trees in them - we don't want to count
# them as part of growth rates, they'll skew things

# first off, exclude really small plants
plant<-plant[plant$tot.ag.mass > 0.01, ]

mass.sd<-with(droplevels(plant[plant$pl.p %in% target.plots,]),
              tapply(tot.ag.mass, pl.p, function(x){sd(x, na.rm=TRUE)}))

# find plants that are outside the 99th quantile
mean(mass.sd, na.rm=TRUE)
sd(mass.sd, na.rm=TRUE)
mass.sd.centered<-(mass.sd-mean(mass.sd, na.rm=TRUE)) / sd(mass.sd, na.rm=TRUE)
hist(mass.sd.centered, breaks=20)

varying.plots<-names(mass.sd.centered[mass.sd.centered >=3])
varying.plots<-varying.plots[!is.na(varying.plots)]

error.rows<-lapply(varying.plots, function(x){
  print(x)
  temp<-plant[plant$pl.p==x,]
  
  iqr<-quantile(temp$tot.ag.mass, 0.75)-quantile(temp$tot.ag.mass, 0.25)
  rownames(temp)[which(temp$tot.ag.mass > 3*iqr)]
})

plant<-plant[-as.numeric(error.rows),]

plot.biomass<-with(plant, tapply(tot.ag.mass, pl.p, sum))

site$plot.biomass<-plot.biomass[match(site$pl.p, names(plot.biomass))]

site$biomass.area<-site$plot.biomass/site$plot.area

hist(log(site$biomass.area))

# density
stem.count<-with(plant, tapply(species, pl.p, length))
site$stem.count<-stem.count[match(site$pl.p, names(stem.count))]

site$density<-site$stem.count/site$plot.area

# alpha diversity
plant.alive<-plant[plant$dead==0,]
alpha<-with(plant.alive, tapply(species, pl.p, 
                                function(x){length(unique(x))}))
site$alpha<-alpha[match(site$pl.p, names(alpha))]

# rarefied richness

with(site, tapply(stem.count, plot.type, summary))
hist(site$stem.count[site$plot.type=="Remnant" & site$stem.count<100])
hist(site$stem.count[site$plot.type=="Natural regrowth" & site$stem.count<100])
hist(site$stem.count[site$plot.type=="Planting" & site$stem.count<100])

# how many sites
plot(y=sapply(1:50, 
            function(x){length(site$pl.p[site$stem.count>=x])/
                        length(site$pl.p)}),
     x=1:50)

alive.ssmat<-unclass(table(plant.alive$pl.p, plant.alive$species))

rare.rich<-t(rarefy(alive.ssmat[rowSums(alive.ssmat)>=12,], 
                  sample=12, se=TRUE))

site$rare.rich=NA
site$rare.rich.se=NA
site$rare.rich[match(rownames(rare.rich),
                     site$pl.p)]<-rare.rich[,1]
site$rare.rich.se[match(rownames(rare.rich),
                     site$pl.p)]<-rare.rich[,2]

#             CALCULATE TRAITS ####
#                               CONTINUOUS TRAIT IMPORT #####

# import all lists of trait values


raw.trait<-list(height=read.csv("./Data/max.height.csv", header=T),
                sla=read.csv("./Data/sla.csv", header=T),
                seed.mass=read.csv("./Data/seed.mass.csv", header=T),
                wood.density=read.csv("./Data/wood.density.csv", header=T))

trait.mean.list<-lapply(raw.trait, function(x){
  sapply(split(x[,2], f=x$species), function(y){mean(y, na.rm=TRUE)})
})

trait.means<-data.frame(species=species$species)

trait.mean.data<-sapply(trait.mean.list, function(x){
  temp.trait<-data.frame(species=names(x),
                         trait=x)
  temp<-merge(trait.means, temp.trait, all.x=TRUE)
  return(as.vector(temp[,2]))
})
trait.means<-cbind(trait.means, trait.mean.data)
colnames(trait.means)<-c("species", "MH","SLA","SM","WD")

# now we want to get average values for the genus only IDs
# let's add genus to our trait list
trait.means<-merge(trait.means, species[,colnames(species) %in% 
                                          c("species", "genus", "genus.only")], 
                   all.x=TRUE)

# Now we can apply these values to the plantscale dataframe
# (giving each individual plant it's relevant trait score)
plant<-merge(plant, trait.means[,colnames(trait.means) %in%
                                                      c("species","SLA","WD","SM","MH")], 
                       all.x=TRUE)

# Now I want to set up an indicator for each trait telling me how many species it was calculated from
# 0 means it is a species-level value
trait.n<-as.data.frame(matrix(0, ncol=4, nrow=length(plant[,1]), 
                              dimnames=list(NULL, paste(c("SLA","WD","SM","MH"), ".n", sep=""))))
# I also want the radius of the circle that was used to calculate genus means
trait.dist<-as.data.frame(matrix(NA, ncol=4, nrow=length(plant[,1]), 
                                 dimnames=list(NULL, paste(c("SLA","WD","SM","MH"), ".dist", sep=""))))

plant<-cbind(plant, trait.n, trait.dist)

#                               CALCULATE GENUS MEANS ####

# now here we want to set genus means, but universal means are a bad idea, especially in genera with a lot of
# variation (e.g., Eucalyptus). So what we'll do is a hierarchical process. First, I'm assuming if there's another species
# in the sample plot as a genus ID, the data collector recognised the genus ID was NOT the same species, so I don't want
# to apply the species trait values to the genus ID. Instead what I'll do is apply a mean from a discreet region.

# for starters we want all genus IDs, and species without trait values
genus<-unique(plant$species[plant$genus.only==1])

species.missing.trait<-
  unique(plant$species[
    rowSums(is.na(plant[,colnames(plant) %in% 
                                    c("SLA","WD","SM","MH")]))>0 &
      plant$genus.only==0 &
      plant$unknown==0])


# next thing we need is a distance matrix for our plantings
grid.plp<-expand.grid(site$pl.p, site$pl.p, stringsAsFactors=FALSE)
grid.longlat<-data.frame(long1=site$long[match(grid.plp[,1], site$pl.p)],
                         lat1=site$lat[match(grid.plp[,1], site$pl.p)],
                         long2=site$long[match(grid.plp[,2], site$pl.p)],
                         lat2=site$lat[match(grid.plp[,2], site$pl.p)])
head(grid.longlat)

# custom function to turn degrees to radians
deg2rad <- function(deg) {(deg * pi) / (180)} 

R = 6371e3 # earth circumference in metres
lat1<-deg2rad(grid.longlat$lat1) # φ1 (lat of plot 1 in radians)
lat2<-deg2rad(grid.longlat$lat2) # φ2 (lat of plot 2 in radians)
dlat<-with(grid.longlat, deg2rad(lat2 - lat1)) # Δφ (lat difference in radians)
dlong<-with(grid.longlat, deg2rad(long2 - long1)) # Δλ (long difference in radians)

# funky trigonometry (converted from javascript function found at
# http://www.movable-type.co.uk/scripts/latlong.html)
a<-sin(dlat/2) * sin (dlat/2) +
  cos(lat2) * cos (lat2) *
  sin(dlong/2) * sin(dlong/2)
c<- 2 * atan2(sqrt(a), sqrt(1-a))
d <- (R * c) /1000 # distance in km

pl.p.dist<-matrix(d, ncol=length(site$pl.p), 
                  nrow=length(site$pl.p),
                  dimnames=list(site$pl.p, site$pl.p),
                  byrow=TRUE)

# now for each plot we need to find each species missing a trait value, then attempt to
# find other congeners in a certain radius
species.by.pl.p<-split(plant, f=plant$pl.p)

# identify which genus IDs / species are missing traits
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
setDefaultCluster(cl)
# export data and load libraries in cluster
clusterExport(varlist=c("plant", "pl.p.dist", "species", "species.by.pl.p"))
genus.means.list<-parLapply(cl=cl, 1:length(species.by.pl.p), function(row.num){
  
  ### 1. GET PL.P ###
  pl.p<-species.by.pl.p[[row.num]]
  
  ### 2. GET PL.P CHARACTERISTICS ###
  
  # get species with missing traits
  species.traits<-pl.p[!duplicated(pl.p$species),
                       colnames(pl.p) %in% 
                         c("species", "SLA","WD","SM","MH")]
  missing.traits<-is.na(species.traits[,-1])
  missing.species<-species.traits[rowSums(is.na(species.traits[,-1]))>0,1]
  
  if(length(missing.species)==0){return(NULL)}
  
  # get distance scores for pl.p
  pl.p.distances<-pl.p.dist[row.num,]
  
  # identify genus of species
  missing.genus<-species$genus[species$species %in% missing.species]
  
  # now for each species, attempt to find other members of genus within a given radius
  
  # 3. IDENTIFY TRAIT VALUES FOR EACH GENUS WITHIN GIVEN RADIUS ###
  
  pl.p.mean.traits<-lapply(missing.genus, function(genus){
    
    candidate.species<-species$species[species$genus==genus]
    
    if(length(candidate.species)==1){return(rep(NA,12))}
    
    # for a number of radii, calculate mean traits and sample size
    cand.traits<-t(sapply(c(50,100,200,500,1000,10000), function(radius){
      
      # subset pl.ps within radii
      pl.p.radii<-plant[
        plant$pl.p %in%
          names(pl.p.distances)[pl.p.distances<radius] &
          plant$species %in% candidate.species,]
      if(dim(pl.p.radii)[1]==0){return(c(rep(NA,4),rep(0,4), rep(NA,4)))} # if no species from genus, return NAs
      
      if(dim(pl.p.radii)[1]>0){
        
        pl.p.radii<-pl.p.radii[!duplicated(pl.p.radii$species),]
        
        trait.means<-colMeans(pl.p.radii[,colnames(pl.p.radii) %in%
                                           c("SLA","WD","SM","MH")], 
                              na.rm=TRUE)
        n<-apply(pl.p.radii[,colnames(pl.p.radii) %in% 
                              c("SLA","WD","SM","MH")], 
                 MARGIN=2, function(x){length(x[!is.na(x)])})
        
        return(c(trait.means, n, rep(radius,4)))
      }
    }))
    
    # now accept the first non-NA value
    smallest.distance<-apply(cand.traits[,1:4], MARGIN=2, function(x){match("TRUE", !is.na(x))})
    
    if(sum(!is.na(smallest.distance))>0){ # if there's a value for all traits
      
      genus.means<-cand.traits[cbind(smallest.distance, 1:4)]
      genus.n<-cand.traits[cbind(smallest.distance, 5:8)]
      genus.distance<-cand.traits[cbind(smallest.distance, 9:12)]
      return(c(genus.means, genus.n, genus.distance))
    }
    
    if(sum(!is.na(smallest.distance))==0){ # if there's no value for any traits (like a genus with no species)
      return(rep(NA,12))
    }
    
  })
  
  pl.p.mean.traits.df<-do.call("rbind", pl.p.mean.traits)
  pl.p.mean.traits.df<-as.data.frame(matrix(pl.p.mean.traits.df, ncol=12))
  colnames(pl.p.mean.traits.df)<-c("SLA","WD","SM","MH",
                                   "SLA.n","WD.n","SM.n","MH.n",
                                   "SLA.rad","WD.rad","SM.rad","MH.rad")
  rownames(pl.p.mean.traits.df)<-missing.species 
  return(pl.p.mean.traits.df)
  
})
stopCluster(cl=cl)

# now we need to overwrite each pl.p's genus traits with new values
# except where traits already exist (for some species IDs)
pl.p.missing.traits<-mapply(pl.p=species.by.pl.p, 
                            trait.means=genus.means.list, 
                            function(pl.p, trait.means){
                              
                              if(is.null(trait.means)){return(NULL)}
                              
                              species<-split(pl.p, f=droplevels(pl.p$species))
                              
                              # first match up the rows of the pl.p to the missing traits
                              trait.match<-match(rownames(trait.means), names(species))
                              
                              # now we need to know whether any of the pl.p species 
                              # already have species level traits
                              trait.cols<-match(c("SLA","WD","SM","MH"), colnames(pl.p))
                              
                              # now we need to override trait values that are NA
                              species.missing.trait<-mapply(missing.species=species[trait.match], 
                                                            traits=split(trait.means, rownames(trait.means)), 
                                                            function(missing.species, traits){
                                                              
                                                              traits.df<-traits[rep(seq_len(nrow(traits)), 
                                                                                    each=length(missing.species[,1])),]
                                                              
                                                              #over-write trait values
                                                              missing.traits<-is.na(missing.species[1,
                                                                                                    grepl("SLA|WD|SM|MH", 
                                                                                                          colnames(missing.species))][,1:4])
                                                              
                                                              trait.cols<-colnames(missing.species) %in% colnames(missing.traits)
                                                              n.cols<-grepl("\\.n", colnames(missing.species))
                                                              dist.cols<-grepl("\\.dist", colnames(missing.species))
                                                              
                                                              missing.species[,trait.cols][,missing.traits]=traits.df[,1:4][,missing.traits]
                                                              missing.species[,n.cols][,missing.traits]=traits.df[,5:8][,missing.traits]
                                                              missing.species[,dist.cols][,missing.traits]=traits.df[,9:12][,missing.traits]
                                                              
                                                              return(missing.species)  
                                                            }, SIMPLIFY=FALSE)
                              
                              
                              pl.p.new<-rbind(do.call("rbind", species[-trait.match]), 
                                              do.call("rbind", species.missing.trait))
                              return(pl.p.new)
                            })

# now we need to combine our list (excluding the Nulls, with the pl.ps with no missing traits)
# first let's get the NULLs out
null.pl.p<-sapply(pl.p.missing.traits, 
                  function(x){if(is.null(x)){return(1)} else{return(0)}})

null.pl.p.names<-names(null.pl.p[null.pl.p==1])

plant<-rbind(plant[plant$pl.p %in% 
                                         null.pl.p.names,],
                       do.call("rbind", pl.p.missing.traits[null.pl.p==0]))

#                               CATEGORICAL TRAIT IMPORT ####

# I've got a list of which species are nitrogen fixers, and categorised
# species as either shrubs (max height of <6m) and trees (max height >6m)
# so we can calculate a proportion of each of these for each plot and 
# analyse them using a logistic regression

cat.trait.species<-species[,colnames(species) %in% c("species", "n.fixer", "tree")]

plant<-merge(plant, cat.trait.species)

cat.traits.pl.p<-as.data.frame(t(sapply(split(plant, plant$pl.p), 
                                        
                                        function(x){
                                          c(mean(x$n.fixer, na.rm=TRUE),
                                            mean(x$tree, na.rm=TRUE))
                                        }
)))

colnames(cat.traits.pl.p)<-c("nfixer.prop", "tree.prop")
cat.traits.pl.p$pl.p<-rownames(cat.traits.pl.p)

head(cat.traits.pl.p)

site<-merge(site, cat.traits.pl.p)

#                               TRANSFORM TRAITS ####

par(mfrow=c(2,2))
sapply(c("SLA","WD","SM","MH"), function(x){
  hist(unique(plant[,colnames(plant) %in% x]), 
       xlab="", main="")
  mtext(side=3, text=x)
})

# traits are pretty skewed, especially seed mass. We'll need to log-transform first
# otherwise our standardisation will be off by orders of magnitude.
log.traits<-sapply(c("SLA","WD","SM","MH"), function(x){
  log(plant[,colnames(plant) %in% x])
})

par(mfrow=c(2,2))
sapply(c("SLA","WD","SM","MH"), function(x){
  hist(unique(log.traits[,colnames(log.traits) %in% x]), 
       xlab="", main="")
  mtext(side=3, text=x)
})

# save our raw traits
plant.raw<-plant

# override raw traits with standardised values
plant[,match(colnames(log.traits), colnames(plant))]=
  log.traits

#             SAVE PROCESSED DATA ####
write.csv(site, "./Data/mixed site sub.csv")
write.csv(plant, "./Data/mixed plant sub.csv")
write.csv(plant.raw, "./Data/mixed plant sub (raw traits).csv")

# ####
# DATA ANALYSIS ####
#             IMPORT PROCESSED DATA ####

site<-read.csv("./Data/mixed site sub.csv")

plot.type2<-site$plot.type
levels(plot.type2)<-c(levels(plot.type2), "Young regrowth", "Old regrowth")

plot.type2[site$plot.type=="Natural regrowth" & site$age <=20]="Young regrowth"
plot.type2[site$plot.type=="Natural regrowth" & site$age >20]="Old regrowth"
plot.type2<-relevel(plot.type2, "Remnant")

site$plot.type<-plot.type2
site<-droplevels(site)
summary(site$plot.type)

table(site$plot.type, site$alpha, site$age)
site<-droplevels(site[site$alpha >1, ])

plant<-read.csv("./Data/mixed plant sub.csv")
summary(exp(plant$SLA))

plant<-droplevels(plant[plant$pl.p %in% site$pl.p,])
plant<-merge(plant, site[,c(colnames(site)[!colnames(site) %in% colnames(plant)],
                            "pl.p"),],
             by.x="pl.p", by.y="pl.p", all.x=TRUE, all.y=FALSE)
species<-read.csv("./Data/species.list.csv")
plant<-merge(plant, species[,c("species","fam","genus")], all.x=TRUE, all.y=FALSE)
plant<-droplevels(plant)
head(plant)

# Colours for plotting
colours<-read.csv("./Data/plot.colours.csv")

#             COMPARISONS BETWEEN PLOT TYPES ####
#                               data prep ####

site$log.biomass.area<-log(site$biomass.area)
site$log.rare.rich<-log(site$rare.rich)

#site$pl.p <- as.numeric(site$pl.p)
#site$planting <- as.numeric(site$planting)
#write.csv(site, "./Outputs/plot_attribute_data.csv")

response.vars<-c("log.biomass.area",
                 "log.density",
                 "log.rare.rich")

#                               models ####

plot.type.models<-lapply(response.vars, function(x){
  
  temp.data<-site[!is.na(site[,x]),]

  lmer(as.formula(paste0(x, "~ plot.type + (1|preRE/planting)")), 
       data=temp.data)
})

spat.auto.test <- do.call("rbind", lapply(plot.type.models, function(x){
  
  temp.lat <- site$lat + rnorm(nrow(site), 0, 1e-3)
  temp.long <- site$long + rnorm(nrow(site), 0, 1e-3)
  
  inv.mat <- 1/as.matrix(dist(cbind(temp.long,temp.lat)))
  diag(inv.mat) = 0
  summary(inv.mat)
  
  rbind(Moran.I(x=resid(x),
          weight= inv.mat))

  gm0 <- gearymoran(inv.mat, as.data.frame(resid(x)), 999, alter="greater")
  data.frame(obs = gm0$obs,
          exp = gm0$expvar$Expectation,
          exp.var =gm0$expvar$Variance,
          alter=gm0$alter,
          pvalue=gm0$pvalue)
  
  }))

#                               glht ####

plot.type.models[[1]]
plot.type.glht<-lapply(plot.type.models, function(x){
  glht(x, linfct = mcp(plot.type="Tukey"))
})

plot.glht.summ<-lapply(plot.type.glht, function(x){
  do.call("cbind", summary(x)$test[-(c(1,2,7))])
})
names(plot.glht.summ)<-response.vars

# which models show a difference between two groups?
plot.glht.sub<-plot.glht.summ[which(sapply(plot.glht.summ, function(x){
  ifelse(sum(x[,4] <=0.05) > 0, TRUE, FALSE)
}))]

glht.output<-do.call("rbind", lapply(1:length(plot.glht.sub), function(x){
  temp<-as.data.frame(plot.glht.sub[[x]])
  temp<-round(temp, 3)
  temp$var<-rep(names(plot.glht.sub)[x], dim(temp)[1])
  return(temp)
}))

write.csv(glht.output, "./Outputs/plot comparison glht.csv")

#                               basic coef plot ####

pdf("./Plots/plot comparison models.pdf", height=4.5, width=1.9, useDingbats=FALSE)
par(mfrow=c(3,1), oma=c(2.5,3,1,0.5), mar=c(0,0,0,0), 
    las=1, tck=-0.02, ps=8, mgp=c(3,0,0))

names(plot.type.models)<-response.vars
temp.names<-cbind(response.vars,
                  c("Plot aboveground biomass ln(kg/ha",
                    "Plant density ln(plants/ha",
                    "Rarefied species richness (n = 12)"))

lapply(1:3, function(x){
  
  iless<-update(plot.type.models[[x]], .~. -1)
  raw.data<-plot.type.models[[x]]@frame
  head(raw.data)
  
  temp.coefs<-summary(iless)$coefficients
  levels<-substr(rownames(temp.coefs),
                 nchar("plot.type")+1, nchar(rownames(temp.coefs)))
  
  ylims<-rbind(log(c(1000,600000)),
               log(c(100,9000)),
               log(c(1,9)))[x,]
  
  plot(x=NULL, y=NULL, xlim=c(0.5,4.5), ylim=ylims,
       axes=FALSE, xlab="", ylab="")
  box()
  
# Y-axes
  if(x==1){
    
    axis(side=2, at=log(c(1000,10000,100000,500000)), labels=c(1,10,100,500), mgp=c(3,0.5,0))
    axis(side=2, at=log(c(seq(1000,10000,1000),seq(10000,100000,10000), 
                          seq(100000,1000000,100000))), tck=-0.01,
         labels=NA, mgp=c(3,0.5,0))
  }
  
  if(x==2){
    
    axis(side=2, at=log(c(100,1000,5000)), labels=c(100,1000,5000), mgp=c(3,0.5,0))
    axis(side=2, at=log(c(seq(100,1000,100),seq(1000,10000,1000))), tck=-0.01,
         labels=NA, mgp=c(3,0.5,0))
  }
  
  if(x==3){
    
    axis(side=2, at=log(c(1,2,5)), labels=c(1,2,5), mgp=c(3,0.5,0))
    axis(side=2, at=log(seq(1,10,1)), tck=-0.01,
         labels=NA, mgp=c(3,0.5,0))
  }
  
    mtext(side=2, text= c("Plot aboveground biomass (Mg/ha)",
                          "Plant density (plants/ha)",
                          "Rarefied species richness (n = 12)")[x], line=c(1.75,2,1.5)[x], 
          las=0, cex=0.7)
    
    if(x==3){
      axis(side=1, at=1:4, 
           labels=c("Remnant", NA, NA, NA), cex=0.8)
      
      mtext(side=1, at=2, line=0, text="Planting", cex=0.6)
      
      mtext(side=1, at=3:4, line=0.6, 
            text=c("Young\nregrowth", "Old\nregrowth"), cex=0.6)
    } else {
      axis(side=1, at=1:4, 
           labels=NA, cex=0.8)
    }

    # raw points
  sapply(1:4, function(y){
    
    box.data<-plot.type.models[[x]]@frame[plot.type.models[[x]]@frame[,2]==levels[y],]
    
    #boxplot(box.data[,1], at=y, add=TRUE)
    points(y=box.data[,1], x=jitter(rep(y, dim(box.data)[1]), amount=0.2), 
           col=rgb(colours[y,5],
                   colours[y,6],
                   colours[y,7], 0.4),
           pch=16)
    
  })
  
  points(x=1:4, y=temp.coefs[,1], pch=16,
         col="black")
  
  segments(x0=1:4, x1=1:4,
           y0=temp.coefs[,1] + 1.96*temp.coefs[,2],
           y1=temp.coefs[,1] - 1.96*temp.coefs[,2], lwd=1.5,
           col="black")  
  
  # significance labels
  sig.labs <- rbind(c("A", "C", "B", "AB"),
                    c("A", "B", "A", "AB"),
                    c("A", "B", "B", "B"))[x,]
  
  text.adj <- c(0.35,0.2,0.1)[x]
  
  text(x=1:4, y=temp.coefs[,1] + 1.96*temp.coefs[,2] + text.adj,
       labels=sig.labs, font=2)

  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"), 
       labels=paste0("(",letters[1:3],")")[x], font=2,
         adj=0, cex=1.1)
  
})
dev.off()
#             SIZE CLASS DISTRIBUTIONS ####

#                               data prep ####
  
# cut by arbitrary number of points
plant$mass.cut<-cut(log(plant$tot.ag.mass), breaks=20)

cut.points<-cbind(as.numeric(substr(levels(plant$mass.cut), 
                                    2, 
                                    regexpr("\\,", levels(plant$mass.cut))-1)),
                  as.numeric(substr(levels(plant$mass.cut), 
                                    regexpr("\\,", levels(plant$mass.cut))+1, 
                                    nchar(levels(plant$mass.cut))-1)))
mid.point<-cut.points[,1] + 0.5*(cut.points[,2]-cut.points[,1])
levels(plant$mass.cut)<-mid.point

cut.data<-do.call("rbind", lapply(split(plant, f=plant$pl.p), function(x){
  
  temp<-as.data.frame(table(x$mass.cut))
  temp$Freq<-temp$Freq / x$plot.area[1]
  temp$pl.p<-rep(x$pl.p[1], dim(temp)[1])
  temp$planting<-plant$planting[plant$pl.p==temp$pl.p[1]][1]
  temp$plot.type<-plant$plot.type[plant$pl.p==temp$pl.p[1]][1]
  temp$preRE <- plant$preRE[plant$pl.p==temp$pl.p[1]][1]

  return(temp)
}))
cut.data$pl.p <- as.numeric(cut.data$pl.p)
cut.data$planting <- as.numeric(cut.data$planting)

write.csv(cut.data, "./Outputs/size_class_data.csv")

#                             models ####
cut.data$Var1<-as.numeric(as.character(cut.data$Var1))

cut.data.bin<-cut.data
cut.data.bin$Freq <- ifelse(cut.data.bin$Freq > 0 , 1, 0)
cut.data.nozeros<-cut.data[cut.data$Freq > 0,]

bin.counts<-with(cut.data.bin[cut.data.bin$Freq==1,], as.data.frame(table(Var1, plot.type)))
bin.counts<-bin.counts[bin.counts$Freq==0,]
# model by mass and height as thin plate spline

# model will only converge if we truncate splines in each category so there are no
# 0 plot.type x mass bins
cut.data.bin<-cut.data.bin[!paste0(cut.data.bin$Var1, cut.data.bin$plot.type) %in%
               paste0(bin.counts[,1], bin.counts[,2]),]

cut.model.occur<-gamm(Freq ~ plot.type + s(Var1, bs="cr", k=6, by=plot.type),
                       random=list(preRE = ~1,
                                   planting=~1,
                                   pl.p=~1),
                       family=binomial,
                       data=cut.data.bin)

summary(cut.model.occur$gam)
write.csv(cbind(summary(cut.model.occur$gam)$p.table,
                summary(cut.model.occur$gam)$s.table), "./Outputs/size class occurrence coefs.csv")

cut.model.count<-gamm(log(Freq) ~ plot.type + s(Var1, bs="cr", k=6, by=plot.type),
                random=list(preRE = ~1,
                            planting=~1,
                            pl.p=~1),
                family=gaussian,
                data=cut.data.nozeros)

summary(cut.model.count$lme)
write.csv(cbind(summary(cut.model.count$gam)$p.table,
                summary(cut.model.count$gam)$s.table), "./Outputs/size class density coefs.csv")

write.csv(summary(cut.model.occur$gam)$p.table, "size class occur parametric coefs.csv")
write.csv(summary(cut.model.occur$gam)$s.table, "size class occur smoothing coefs.csv")

#                             Plot ####

count.gam.predict<-lapply(levels(cut.data.nozeros$plot.type), function(x){
  
  temp.data<-cut.data.nozeros[cut.data.nozeros$plot.type==x,]
  
  pred.data=data.frame(Var1=rep(seq(min(temp.data$Var1),
                                    max(temp.data$Var1),
                                    length=200)),
                       plot.type=rep(x, 200))
  
  temp.predict<-do.call("cbind", predict(cut.model.count$gam, 
                        newdata=pred.data,
                        se.fit=TRUE, type="response"))

  return(cbind(pred.data, temp.predict))
  
})

occur.gam.predict<-lapply(levels(cut.data$plot.type), function(x){
  
  temp.data<-cut.data.bin[cut.data.bin$plot.type==x,]
  
  pred.data=data.frame(Var1=rep(seq(min(temp.data$Var1),
                                    max(temp.data$Var1),
                                    length=200)),
                       plot.type=rep(x, 200))
  
  temp.predict<-do.call("cbind", predict(cut.model.occur$gam, 
                                         newdata=pred.data,
                                         se.fit=TRUE, type="link"))
  
  return(cbind(pred.data, temp.predict))
  
})
names(count.gam.predict)<-levels(cut.data$plot.type)
names(occur.gam.predict)<-levels(cut.data$plot.type)

pdf("./Plots/new size bin plot (no height).pdf", height=3.75, width=4, useDingbats=FALSE)
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(2.5,3,1.5,1), las=1, 
    tck=-0.025, mgp=c(3,0.5,0), ps=8)

#                                       Occurrence ####

ylims<-c(0,1)
xlims<-summary(cut.data$Var1)[c(1,6)]+c(-0.2,0.25)

#                                             Plantings ####

plot(y=NULL, x=NULL, type="n",
     ylim=ylims, xlim=xlims, yaxs="i", xlab="", ylab="", axes=FALSE)

axis(side=1, at=log(c(0.01,0.1,1,10,100,1000)), 
     labels=NA, mgp=c(3,0.2,0))
box()

axis(side=1, at=log(c(seq(0.1,1,0.1),
                      seq(1,10,1), seq(10,100,10), seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck=-0.01)

axis(side=2, at=seq(0,1,0.2), labels=seq(0,1,0.2))

mtext(side=2, line=1.5, text=expression("Probability of occurrence"), las=0, cex=0.8)

temp.predict<-occur.gam.predict[["Young regrowth"]]

temp.data<-droplevels(cut.data.bin[cut.data.bin$plot.type=="Young regrowth" &
                                     cut.data.bin$Var2 %in% temp.predict$Var2,])

points(jitter(temp.data$Freq, amount=0.05) ~ jitter(temp.data$Var1, amount=0.1),
       col=rgb(colours[3,2],
               colours[3,3],
               colours[3,4], 0.7),
       pch=c(1, 16)[as.factor(temp.data$Var2)])

# Young regrowth
with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(plogis(fit + 1.96*se.fit),
                 rev(plogis(fit - 1.96*se.fit))),
             col=rgb(colours[3,2],
                     colours[3,3],
                     colours[3,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(plogis(fit) ~ Var1, type="l",
            col=rgb(colours[3,2],
                    colours[3,3],
                    colours[3,4]), lwd=1.5))


# Plantings
temp.predict<-occur.gam.predict[["Planting"]]
temp.data<-droplevels(cut.data.bin[cut.data.bin$plot.type=="Planting" &
                                     cut.data.bin$Var2 %in% temp.predict$Var2,])

points(jitter(temp.data$Freq, amount=0.05) ~ jitter(temp.data$Var1, amount=0.1),
       col=rgb(colours[2,2],
               colours[2,3],
               colours[2,4], 0.7),
       pch=c(1, 16)[as.factor(temp.data$Var2)])

with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(plogis(fit + 1.96*se.fit),
                 rev(plogis(fit - 1.96*se.fit))),
             col=rgb(colours[2,2],
                     colours[2,3],
                     colours[2,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(plogis(fit) ~ Var1, type="l",
            col=rgb(colours[2,2],
                    colours[2,3],
                    colours[2,4]), lwd=1.5))


text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(a)", font=2, cex=1, adj=0)

mtext(side=3, line=0.2, font=2, text="Young stands")

text(x=relative.axis.point(0.5, "x"), y=relative.axis.point(0.5, "y"),
     labels=c("Young regrowth", "Planting"), col=rgb(colours[3:2,2],
                                                     colours[3:2,3],
                                                     colours[3:2,4],1))

#                                             Old regrowth ####

plot(y=NULL, x=NULL, type="n",
     ylim=ylims, xlim=xlims, yaxs="i", xlab="", ylab="", axes=FALSE)

axis(side=1, at=log(c(0.01,0.1,1,10,100,1000)), 
     labels=NA, mgp=c(3,0.2,0))
box()

axis(side=1, at=log(c(seq(0.1,1,0.1),
                      seq(1,10,1), seq(10,100,10), seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck=-0.01)

axis(side=2, at=seq(0,1,0.2), labels=NA)

temp.predict<-occur.gam.predict[["Old regrowth"]]

temp.data<-droplevels(cut.data.bin[cut.data.bin$plot.type=="Young regrowth" &
                                     cut.data.bin$Var2 %in% temp.predict$Var2,])

points(jitter(temp.data$Freq, amount=0.05) ~ jitter(temp.data$Var1, amount=0.1),
       col=rgb(colours[4,2],
               colours[4,3],
               colours[4,4], 0.7),
       pch=c(1, 16)[as.factor(temp.data$Var2)])

# Planting
with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(plogis(fit + 1.96*se.fit),
                 rev(plogis(fit - 1.96*se.fit))),
             col=rgb(colours[4,2],
                     colours[4,3],
                     colours[4,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(plogis(fit) ~ Var1, type="l",
            col=rgb(colours[4,2],
                    colours[4,3],
                    colours[4,4]), lwd=1.5))

#                                             Remnants ####

temp.predict<-occur.gam.predict[["Remnant"]]

# Young regrowth
with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(plogis(fit + 1.96*se.fit),
                 rev(plogis(fit - 1.96*se.fit))),
             col=rgb(colours[1,2],
                     colours[1,3],
                     colours[1,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(plogis(fit) ~ Var1, type="l",
            col=rgb(colours[1,2],
                    colours[1,3],
                    colours[1,4]), lwd=1.5))

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(b)", font=2, cex=1, adj=0)

mtext(side=3, line=0.2, font=2, text="Old stands")

text(x=relative.axis.point(0.5, "x"), y=relative.axis.point(0.5, "y"),
     labels=c("Old regrowth", "Remnant"), col=rgb(colours[c(4,1),2],
                                                  colours[c(4,1),3],
                                                  colours[c(4,1),4],1))

#                                       Abundance ####

ylims<-log(c(4,400))
xlims<-summary(cut.data$Var1)[c(1,6)]+c(-0.2,0.25)

#                                             Plantings ####

plot(y=NULL, x=NULL, type="n",
     ylim=ylims, xlim=xlims, yaxs="i", xlab="", ylab="", axes=FALSE)

axis(side=1, at=log(c(0.01,0.1,1,10,100,1000)), 
     labels=c(0.01,0.1,1,10,100,1000), mgp=c(3,0.2,0))
box()

axis(side=1, at=log(c(seq(0.1,1,0.1),
                      seq(1,10,1), seq(10,100,10), seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck=-0.01)
mtext(side=1, line=1.2, at=par("usr")[2], text="Plant aboveground biomass (kg)", cex=0.8)

axis(side=2, at=log(c(0,1,10,100,800)+1), labels=c(0,1,10,100,800))
axis(side=2, at=log(c(seq(1,10,1), seq(10,100,10), seq(100,800,100))+1), labels=NA, tck=-0.01)

mtext(side=2, line=1.5, text="Plant density (plants/ha)", las=0, cex=0.8)

temp.predict<-count.gam.predict[["Planting"]]

with(temp.predict[1:200,],
polygon(x=c(Var1, rev(Var1)),
        y=c(fit + 1.96*se.fit,
            rev(fit - 1.96*se.fit)),
        col=rgb(colours[2,2],
                colours[2,3],
                colours[2,4], 0.5), border=NA))

with(temp.predict[1:200,],
points(fit ~ Var1, type="l",
       col=rgb(colours[2,2],
               colours[2,3],
               colours[2,4]), lwd=1.5))

with(temp.predict[201:400,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(fit + 1.96*se.fit,
                 rev(fit - 1.96*se.fit)),
             border=rgb(colours[2,2],
                     colours[2,3],
                     colours[2,4], 1), col=rgb(1,1,1,0)), lty="dashed")

with(temp.predict[201:400,],
     points(fit ~ Var1, type="l",
            col=rgb(colours[2,2],
                    colours[2,3],
                    colours[2,4]), lwd=1.5, lty="dashed"))


#                                             Young regrowth ####

temp.predict<-count.gam.predict[["Young regrowth"]]

with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(fit + 1.96*se.fit,
                 rev(fit - 1.96*se.fit)),
             col=rgb(colours[3,2],
                     colours[3,3],
                     colours[3,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(fit ~ Var1, type="l",
            col=rgb(colours[3,2],
                    colours[3,3],
                    colours[3,4]), lwd=1.5))

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.92, "y"),
     labels="(c)", font=2, cex=1, adj=0)

text(x=relative.axis.point(0.5, "x"), y=relative.axis.point(0.5, "y"),
     labels=c("Young regrowth", "Planting"), col=rgb(colours[3:2,2],
                                                     colours[3:2,3],
                                                     colours[3:2,4],1))

#                                             Old regrowth  ####
plot(y=NULL, x=NULL, type="n",
     ylim=ylims, xlim=xlims, yaxs="i", xlab="", ylab="", axes=FALSE)

axis(side=1, at=log(c(0.01,0.1,1,10,100,1000)), 
     labels=c(0.01,0.1,1,10,100,1000), mgp=c(3,0.2,0))
box()

axis(side=1, at=log(c(seq(0.1,1,0.1),
                      seq(1,10,1), seq(10,100,10), seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck=-0.01)

temp.predict<-count.gam.predict[["Old regrowth"]]

with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(fit + 1.96*se.fit,
                 rev(fit - 1.96*se.fit)),
             col=rgb(colours[4,2],
                     colours[4,3],
                     colours[4,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(fit ~ Var1, type="l",
            col=rgb(colours[4,2],
                    colours[4,3],
                    colours[4,4]), lwd=1.5))

#                                             Remnants ####
 
temp.predict<-count.gam.predict[["Remnant"]]

with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(fit + 1.96*se.fit,
                 rev(fit - 1.96*se.fit)),
             col=rgb(colours[1,2],
                     colours[1,3],
                     colours[1,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(fit ~ Var1, type="l",
            col=rgb(colours[1,2],
                    colours[1,3],
                    colours[1,4]), lwd=1.5))

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.92, "y"),
     labels="(d)", font=2, cex=1, adj=0)

text(x=relative.axis.point(0.5, "x"), y=relative.axis.point(0.5, "y"),
     labels=c("Old regrowth", "Remnant"), col=rgb(colours[c(4,1),2],
                                                  colours[c(4,1),3],
                                                  colours[c(4,1),4],1))

dev.off()


#             PLOT AGE-BIOMASS TRENDS OVER TIME ####
#                                           DATA PREP ####

model.data<-site[complete.cases(site[,c("plot.type","biomass.area","age", "log.aridity.index",
                                        "planting","preRE","elevation.relief",
                                        "silt","carbon")]),]

plot.uncent<-model.data[!is.na(model.data$log.density),]

center.data.output<-sapply(plot.uncent[,c("log.age", "log.density", 
                                          "log.aridity.index",
                                          "SLA.cwm", "WD.cwm", "SM.cwm", "MH.cwm",
                                          "FRv", "FEm", "FDm")],
                           function(x){
                             (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
                           })

plot.cent<-cbind(plot.uncent[,c("biomass.area", "preRE", "planting", 
                                "plot.type")],
                 center.data.output)

#                                           MODELS ####

plot.uncent$plot.type<-as.character(plot.uncent$plot.type)
plot.uncent$plot.type[plot.uncent$plot.type %in% 
                        c("Young regrowth", "Old regrowth")] =  "Natural regrowth"
plot.uncent$plot.type<-as.factor(plot.uncent$plot.type)
plot.uncent$plot.type<-relevel(plot.uncent$plot.type, "Remnant")

model.data<-droplevels(plot.uncent[plot.uncent$plot.type !="Remnant", ])

null.model<-lmer(log(biomass.area) ~ plot.type + (1|preRE/planting),
                 data=model.data, REML=FALSE)
summary(null.model)

age.model<-lmer(log(biomass.area) ~ log.age * plot.type + (1|preRE/planting),
                data=model.data,
                REML=FALSE)
summary(age.model)
anova(null.model, age.model)

age.quad.model<-lmer(log(biomass.area) ~ (log.age + I(log.age^2)) * plot.type + (1|preRE/planting),
                     data=model.data,
                     REML=FALSE,
                     control=lmerControl(optimizer="bobyqa"))
summary(age.quad.model)
anova(age.model, age.quad.model)
plot(age.quad.model)

biomass.base.model<-lmer(log(biomass.area) ~ plot.type + (1|planting),
                         data=plot.uncent, REML=FALSE)
summary(biomass.base.model)
plot(biomass.base.model)

#                                           PLOT ####

predict.points<-lapply(c("Planting","Natural regrowth"), function(x){
  
  temp.data<-model.data[model.data$plot.type==x,]
  factor.level<-unique(as.numeric(model.data$plot.type[model.data$plot.type == x]))
  
  new.data<-data.frame(log.age=seq(min(temp.data$log.age),
                                   max(temp.data$log.age), 
                                   length=200),
                       plot.type=factor(rep(x, 200),
                                        levels=levels(temp.data$plot.type)),
                       planting = "a",
                       preRE = "a")
  
  points <- predictInterval(age.quad.model, 
                            newdata = new.data,
                            n.sims = 3999,
                            which="fixed",
                            level=0.95,
                            include.resid.var=FALSE,
                            type="linear.prediction")

  return(cbind(new.data, points))
  
})

pdf(paste0("./Plots/SEQ-BBS log(mass)",
           Sys.Date(), ".pdf"), height=2.2, width=2.5, useDingbats=FALSE)
par(mar=c(1.5,2,1,1), ps=6, tck=-0.01, mgp=c(3,0.2,0), las=1)
plot(x=NULL, y=NULL, 
     xlim=log(c(4, 150)), 
     ylim=c(log(1200), log(500000)), axes=FALSE)
box()
abline(v=log(85), lty="dashed")
axis(side=2, at=log(c(1000,5000,10000,50000,100000,500000)),
     labels=c(1,5,10,50,100,500))
axis(side=2, at=log(c(seq(1000,10000,1000),
                      seq(10000,100000,10000),
                      seq(100000,700000,100000))), tck=-0.005, labels=NA)
mtext(side=2, line=1, text="Plot aboveground biomass (Mg/ha)", las=0)

axis(side=1, at=log(c(1, 2, 5, 10, 20, 50)),
     labels=c(1, 2, 5, 10, 20, 50), mgp=c(3,-0.2,0))
axis(side=1, at=log(c(seq(1,20,1),
                      seq(20,80,10))), tck=-0.005, labels=NA)
axis(side=1, at=log(120), label="Remnant", mgp=c(3,-0.2,0))
mtext(side=1, line=0.5, text=expression("Regrowing forest age (years)"))

# Raw points
with(model.data,
     points(log(biomass.area) ~ jitter(log.age, amount=0.025),
            col=rgb(colours[c(4,2),5],
                    colours[c(4,2),6],
                    colours[c(4,2),7], 0.6)[plot.type], 
            cex=0.6, pch=c(17,15)[plot.type], lwd=0.5))

# separate regrowth polygon into 3 segments, based on where we have data and 
# where we don't. Not running at the moment
if(FALSE){
  reg.1<-max(which(predict.points[[2]]$log.age <= log(5)))
  reg.2<-max(which(predict.points[[2]]$log.age < log(10)))
  
  with(predict.points[[2]][1:reg.1,],
       polygon(x=c(log.age, rev(log.age)),
               y=c(fit + 1.96* se.fit,
                   rev(fit - 1.96* se.fit)),
               border=NA, col=rgb(0,0.5,0,0.5)))
  
  with(predict.points[[2]][reg.1:reg.2,],
       polygon(x=c(log.age, rev(log.age)),
               y=c(fit + 1.96* se.fit,
                   rev(fit - 1.96* se.fit)),
               border=NA, col=rgb(0,0.5,0,0.2)))
  
  with(predict.points[[2]][reg.2:200,],
       polygon(x=c(log.age, rev(log.age)),
               y=c(fit + 1.96* se.fit,
                   rev(fit - 1.96* se.fit)),
               border=NA, col=rgb(0,0.5,0,0.5)))
  
  with(predict.points[[2]][1:reg.1,],
       points(fit ~ log.age, type='l', col="darkgreen", lwd=2))
  
  with(predict.points[[2]][reg.1:reg.2,],
       points(fit ~ log.age, type='l', col="darkgreen", lwd=2, lty="22"))
  
  with(predict.points[[2]][reg.2:200,],
       points(fit ~ log.age, type='l', col="darkgreen", lwd=2))
}

# Remnant points
rem.points<-site[site$plot.type=="Remnant",]
with(rem.points, 
     points(log(biomass.area) ~ jitter(rep(log(120), dim(rem.points)[1]),
                                       amount=0.2),
            col=rgb(colours[1,5], 
                    colours[1,6], 
                    colours[1,7], 0.4), 
            cex=0.6, pch=16, lwd=0.5))

base.coefs<-summary(biomass.base.model)$coefficients

points(y=base.coefs[1,1], x=log(120), pch=16, col=rgb(colours[1,2],
                                                      colours[1,3],
                                                      colours[1,4]), cex=0.8)
segments(y0=base.coefs[1,1] + 1.96* base.coefs[1,2],
         y1=base.coefs[1,1] - 1.96* base.coefs[1,2],
         x0=log(120), x1=log(120),
         col=rgb(colours[1,2],
                 colours[1,3],
                 colours[1,4]), lwd=1.5)

# CI Polygons
polygon(x=c(predict.points[[2]]$log.age, rev(predict.points[[2]]$log.age)),
        y=c(predict.points[[2]]$upr,
            rev(predict.points[[2]]$lwr)),
        border=NA, col=rgb(colours[4,2],
                           colours[4,3],
                           colours[4,4], 0.5))

polygon(x=c(predict.points[[1]]$log.age, rev(predict.points[[1]]$log.age)),
        y=c(predict.points[[1]]$upr,
            rev(predict.points[[1]]$lwr)),
        border=NA, col=rgb(colours[2,2],
                           colours[2,3],
                           colours[2,4], 0.5))

# Slopes
points(predict.points[[1]]$fit ~ predict.points[[1]]$log.age,
       type='l', lwd=1.5, col=rgb(colours[2,2],
                                  colours[2,3],
                                  colours[2,4]))

points(predict.points[[2]]$fit ~ predict.points[[2]]$log.age,
       type='l', lwd=1.5, col=rgb(colours[4,2],
                                  colours[4,3],
                                  colours[4,4]))

rect(xleft=log(25), 
     xright=log(85), 
     ybottom=par("usr")[3], 
     ytop=log(3500))

legend(x=log(85),
       y=relative.axis.point(-0.05, "y"),
       pch=c(15,17,16), legend=c("Planting", "Regrowth", "Remnant"),
       col=rgb(colours[c(2,4,1),2],
               colours[c(2,4,1),3],
               colours[c(2,4,1),4]),
       y.intersp=0.5, x.intersp=0.5, xjust=1, yjust=0,
       pt.cex=0.6, bty="n")

box()

dev.off()

summary(age.quad.model)
write.csv(round(summary(age.quad.model)$coefficients, 3), 
          "./Outputs/age quad model coefs.csv")

temp<-plant[!duplicated(plant$species), c("species", "SM")]
temp$SM<-exp(temp$SM)

summary(temp$SM)



#             FUNCTIONAL GROUP ANALYSIS ####
#                               SPECIES GROUPING (GOWER's DISTANCE) ####

plot.species<-do.call("rbind", lapply(split(plant, f=plant$pl.p), function(x){
  x[!duplicated(x$species),]
}))
plot.species$plot.type<-as.character(plot.species$plot.type)
plot.species$plot.type[plot.species$plot.type=="Young regrowth"]<-"Youngregrowth"
plot.species$plot.type[plot.species$plot.type=="Old regrowth"]<-"Oldregrowth"
plot.species$plot.type<-as.factor(plot.species$plot.type)
plot.species$plot.type<-relevel(plot.species$plot.type, "Remnant")

unique.species<-plot.species[!duplicated(plot.species$species),c("SLA","WD","SM","MH")]
rownames(unique.species)<-plot.species$species[!duplicated(plot.species$species)]

unique.species <- exp(unique.species)
#unique.species <- t(apply(unique.species, 1, scale))
colnames(unique.species)<-c("SLA","WD","SM","MH")

write.csv(unique.species, "./Data/subset species.csv")

cat.traits<-read.csv("./Data/subset species traits.csv", row.names=1)
cat.traits$species <- rownames(cat.traits)

unique.species<-as.data.frame(cbind(species=rownames(unique.species),
                      unique.species))

unique.species<-merge(unique.species, cat.traits[,c("species", "n.fixer",
                                                    "seed.class")],
                      all.x=TRUE, all.y=FALSE, by.x="species", by.y="species")

unique.species<-unique.species[unique.species$species != "Parsonsia straminea",]

unique.species[,2:5] <- sapply(unique.species[,2:5], as.numeric)
unique.species$species <- as.character(unique.species$species)
rownames(unique.species) <- unique.species$species

unique.species$SM <- log(unique.species$SM)

# SEED DENDROGRAM
seed.dist<-gowdis(unique.species[,c("SM","seed.class")])
seed.groups<-hclust(seed.dist, method="ward.D")

# GROWTH/CARBON DENDROGRAM
growth.dist<-gowdis(unique.species[,c("SLA","WD","MH", "n.fixer")])
growth.groups<-hclust(growth.dist, method="ward.D")

pdf(paste0("./Plots/seed and growth dendrograms 1 ", Sys.Date(), ".pdf"), 
           height=5, width=5)
par(mfrow=c(1,2), mar=c(2,2,1,5), ps=5)

plot(dendrapply(as.dendrogram(seed.groups), 
           function(x) { attr(x, "height") <- log(attr(x, "height")+1); x }),
     horiz=TRUE, main="Seed mass and dispersal")

plot(dendrapply(as.dendrogram(growth.groups), 
                function(x) { attr(x, "height") <- log(attr(x, "height")+1); x }),
     horiz=TRUE, main="SLA, wood density, height, N fixer")

dev.off()

# if we split our species into n groups, how many plot*group interactions can
# we not estimate (because there's never a species from that group in that plot type)
seed.split.test<-sapply(2:50,
                 function(x){
  
  x<-cutree(seed.groups, k=x)
  
  plot.species$group<-x[match(plot.species$species, names(x))]
  
  group.count<-as.data.frame(t(table(plot.species$group, plot.species$pl.p)))
  colnames(group.count)<-c("pl.p", "group", "count")
  group.count<-merge(group.count, site[,c("pl.p", "planting", "plot.type")])
  
  mm<-model.matrix(object=as.formula("count ~ group*plot.type"),
                   data=group.count)
  
  # model can't handle categories with all 0s, so remove these by identifying
  # which groups have a mean of 0
  means<-with(group.count, 
              tapply(count, paste0(group,":plot.type", plot.type), mean))
  
  names(means[means==0])
  
})
sapply(seed.split.test, length)
# It doesn't take long for us to lose our ability to model these interactions. 

growth.split.test<-sapply(2:50,
                        function(x){
                          
                          x<-cutree(growth.groups, k=x)
                          
                          plot.species$group<-x[match(plot.species$species, names(x))]
                          
                          group.count<-as.data.frame(t(table(plot.species$group, plot.species$pl.p)))
                          colnames(group.count)<-c("pl.p", "group", "count")
                          group.count<-merge(group.count, site[,c("pl.p", "planting", "plot.type")])
                          
                          mm<-model.matrix(object=as.formula("count ~ group*plot.type"),
                                           data=group.count)
                          
                          # model can't handle categories with all 0s, so remove these by identifying
                          # which groups have a mean of 0
                          means<-with(group.count, tapply(count, paste0(group,":plot.type", plot.type), mean))
                          
                          names(means[means==0])
                          
                        })
sapply(growth.split.test, length)

plot(y=sapply(growth.split.test, length), x=2:50, type="l", col="darkgreen")
points(y=sapply(seed.split.test, length), x=2:50, type="l", col="red")

#                               SEED MODEL ####

seed.split<-cutree(seed.groups, k=8)
seed.split[order(seed.split)]

temp<-cbind(c(seq(0.01,2, length.out=50), 0.5),
            sapply(c(seq(0.01,2, length.out=50), 0.5), function(x){
              length(unique(cutree(seed.groups, h=x)))
            }))

seed.dend<-as.dendrogram(seed.groups)
seed.dend <- set(seed.dend, "labels_cex", 0.5)

dend.nodes<-as.data.frame(cbind(get_nodes_xy(seed.dend, type = c("rectangle"), center = FALSE,
                                             horiz = TRUE),
                                get_nodes_attr(seed.dend, attribute="label")))

dend.nodes$V1<-as.numeric(as.character(dend.nodes$V1))
dend.nodes<-dend.nodes[!is.na(dend.nodes$V3),]
dend.nodes$group<-seed.split[match(dend.nodes$V3, names(seed.split))]
dend.groups<-do.call("rbind", lapply(split(dend.nodes, f=dend.nodes$group), 
                                     function(x){c(min(x$V1),max(x$V1))}))
group.centers<-tapply(dend.nodes$V1, dend.nodes$group, mean)

seed.split.names<-names(seed.split)
seed.split<-(8:1)[match(seed.split, order(group.centers))]
names(seed.split) <- seed.split.names

# group centers to rename groups

plot.species$seed.group<-as.factor(seed.split[match(plot.species$species, names(seed.split))])

plot.species[plot.species$species=="Cassinia laevis",]
seed.count<-as.data.frame(t(table(plot.species$seed.group, plot.species$pl.p)))
colnames(seed.count)<-c("pl.p", "seed.group", "count")
seed.count<-merge(seed.count, site[,c("pl.p", "planting", "plot.type", "plot.area", "preRE")])
seed.count$plot.area<-scale(seed.count$plot.area)

t(table(plot.species$seed.group, plot.species$species))
seed.mm<-model.matrix(object=as.formula("count ~ seed.group*plot.type + plot.area"),
                 data=seed.count)
# model can't handle categories with all 0s, so remove these by identifying
# which groups have a mean of 0
means<-with(seed.count, tapply(count, paste0(seed.group,":plot.type", plot.type), mean))
seed.removal<-names(means[means==0])

# Then cut them from dataset
if(length(seed.removal)>0){
  seed.mm<-seed.mm[,!grepl(paste0(seed.removal, collapse="|"), colnames(seed.mm))]
}

seed.count$pl.p <- as.numeric(seed.count$pl.p)
seed.count$planting <- as.numeric(seed.count$planting)
write.csv(seed.count, "./Outputs/seed_fungroup_count.csv")
saveRDS(seed.mm, "./Outputs/seed_model_matrix.RDS")
saveRDS(seed.dend, "./Outputs/seed_dendrogram.RDS")
saveRDS(seed.removal, "./Outputs/seed_missing_groups.RDS")
saveRDS(seed.split, "./Outputs/seed_group_membership.RDS")

seed.model<-glmmPQL(count ~ 0 + seed.mm, random=list(~1|preRE,
                                                     ~1|planting,
                                                     ~1|pl.p),
            family=negative.binomial(theta = 1),
            data=seed.count)

write.csv(summary(seed.model)$tTable,
          "./Outputs/seed group model coefs.csv")

#                               GROWTH MODEL ####
growth.split<-cutree(growth.groups, k=8)
saveRDS(growth.split, "./Outputs/growth_group_membership.RDS")

temp<-cbind(c(seq(0.01,2, length.out=50), 0.5),
            sapply(c(seq(0.01,2, length.out=50), 0.5), function(x){
length(unique(cutree(growth.groups, h=x)))
}))

growth.split[order(names(growth.split))]

growth.dend<-as.dendrogram(growth.groups)
growth.dend <- set(growth.dend, "labels_cex", 0.5)

dend.nodes<-as.data.frame(cbind(get_nodes_xy(growth.dend, type = c("rectangle"), center = FALSE,
                                             horiz = TRUE),
                                get_nodes_attr(growth.dend, attribute="label")))

dend.nodes$V1<-as.numeric(as.character(dend.nodes$V1))
dend.nodes<-dend.nodes[!is.na(dend.nodes$V3),]
dend.nodes$group<-growth.split[match(dend.nodes$V3, names(growth.split))]
dend.groups<-do.call("rbind", lapply(split(dend.nodes, f=dend.nodes$group), 
                                     function(x){c(min(x$V1),max(x$V1))}))
group.centers<-tapply(dend.nodes$V1, dend.nodes$group, mean)

growth.split.names<-names(growth.split)
growth.split<-(8:1)[match(growth.split, order(group.centers))]
names(growth.split) <- growth.split.names

plot.species$growth.group<-growth.split[match(plot.species$species, names(growth.split))]
table(plot.species$growth.group, plot.species$plot.type)
t(table(plot.species$growth.group, plot.species$species))

with(droplevels(plot.species[plot.species$growth.group==2,]),
     table(species, plot.type))

growth.count<-as.data.frame(t(table(plot.species$growth.group, plot.species$pl.p)))
colnames(growth.count)<-c("pl.p", "growth.group", "count")
growth.count<-merge(growth.count, site[,c("pl.p", "planting", "plot.type", "plot.area", "preRE")])
growth.count$plot.area<-scale(growth.count$plot.area)

with(growth.count[growth.count$plot.type=="Planting",], table(growth.group, count))
table(growth.count$count, growth.count$growth.group, growth.count$plot.type)

growth.mm<-model.matrix(object=as.formula("count ~ growth.group*plot.type + plot.area"),
                 data=growth.count)

# model can't handle categories with all 0s, so remove these by identifying
# which groups have a mean of 0
means<-with(growth.count, tapply(count, paste0(growth.group,":plot.type", plot.type), mean))
growth.removal<-names(means[means==0])
# Then cut them from dataset
if(length(growth.removal)>0){
  growth.mm<-growth.mm[,!grepl(paste0(growth.removal, collapse="|"), colnames(growth.mm))]
}

growth.count$pl.p <- as.numeric(growth.count$pl.p)
growth.count$planting <- as.numeric(growth.count$planting)
write.csv(growth.count, "./Outputs/growth_fungroup_count.csv")
saveRDS(growth.mm, "./Outputs/growth_model_matrix.RDS")
saveRDS(growth.dend, "./Outputs/growth_dendrogram.RDS")
saveRDS(growth.removal, "./Outputs/growth_missing_groups.RDS")

growth.model<-glmmPQL(count ~ -1 + growth.mm, random=list(~1|preRE,
                                                          ~1|planting,
                                                          ~1|pl.p),
                    family=negative.binomial(theta = 1),
                    data=growth.count)

summary(growth.model)
write.csv(summary(growth.model)$tTable,
          "./Outputs/growth group model coefs.csv")

#                               PLOT ####
#                                     seed data prep ####

# rerun seed models to get reference level coefs and SEs
seed.coef.df<-as.data.frame(summary(seed.model)$tTable)

seed.coef.df$plot.type<-substr(rownames(seed.coef.df),
                             regexpr("R|P|O|Y", rownames(seed.coef.df)),
                             nchar(rownames(seed.coef.df)))
seed.coef.df$plot.type[!seed.coef.df$plot.type %in% c("Planting", 
                                                  "Young regrowth", 
                                                  "Old regrowth")] = "Remnant"
seed.coef.df$group<-substr(rownames(seed.coef.df),
                           regexpr("seed\\.group", rownames(seed.coef.df))+10,
                           regexpr("seed\\.group", rownames(seed.coef.df))+10)
seed.coef.df$group[is.na(as.numeric(seed.coef.df$group))]=1

# calculate and write comparison to remnants 
seed.cut.mat<-model.matrix(seed.model)
colnames(seed.cut.mat)<-gsub("seed.mm", "", colnames(seed.cut.mat))
colnames(seed.cut.mat)
seed.cut.mat[,"plot.area"]=0
seed.cut.mat<-unique.matrix(seed.cut.mat, MARGIN=1)
seed.marg.mat<-seed.cut.mat # save this to predict actual coefficients

# create contrast matrix to compare group estimates to remnants
seed.cont.mat<-seed.cut.mat
seed.cont.mat[,1:8]<-0
seed.cont.mat<-seed.cont.mat[rowSums(seed.cont.mat)>0,]

# get group names and plot types for each row of marginal matrix
ind<-which(seed.cont.mat !=0, arr.ind=TRUE)
ind<-ind[order(ind[,1]),]

colnames(seed.cont.mat)
plot.type.ind<-ind[ind[,2] %in% 9:11,]
group.ind<-ind[!ind[,2] %in% 9:11,]

colnames(seed.cont.mat)
plot.types<-gsub("plot.type", "", colnames(seed.cont.mat)[plot.type.ind[,2]])

groups<-rep(1, length(plot.types))
groups[group.ind[,1]]<-substr(colnames(seed.cont.mat)[group.ind[,2]],
                  nchar("seed.group")+1,
                  regexpr(":", colnames(seed.cont.mat)[group.ind[,2]])-1)

seed.marg.groups<-cbind(groups, plot.types)
seed.keeps<-!duplicated(paste0(seed.marg.groups[,1],
                               seed.marg.groups[,2]))

seed.marg.groups<-seed.marg.groups[seed.keeps,]
seed.cont.mat<-seed.cont.mat[seed.keeps,]

seed.contrast<-summary(glht(seed.model, linfct=seed.cont.mat))

seed.contrast.df<-data.frame(group=seed.marg.groups[,1],
                             plot.type=seed.marg.groups[,2],
                             estimate=seed.contrast$test$coefficients,
                             se=seed.contrast$test$sigma,
                             t=seed.contrast$test$tstat,
                             p=seed.contrast$test$pvalues)

seed.contrast.df[,-(1:2)] = round(seed.contrast.df[,-(1:2)], 3)
write.csv(seed.contrast.df, "./Outputs/seed contrast test.csv")

# functional group mean data for boxplots and pictures
temp.species<-unique.species
temp.species$group<-seed.split[match(names(seed.split), rownames(temp.species))]

seed.SM<-tapply(temp.species$SM, temp.species$group, mean)
seed.SM<- 3 * (seed.SM / max(seed.SM))
saveRDS(seed.SM, "./Outputs/seed_mass_group_means.RDS")

seed.class<-table(temp.species$seed.class, temp.species$group)
seed.class<-rownames(which(seed.class>0, arr.ind=TRUE))
seed.class<-c("El","Fl","Un","Wi")[as.factor(seed.class)]
saveRDS(seed.class, "./Outputs/seed_disp_class.RDS")

# get plot coefficients and error

seed.marg.df<-data.frame(seed.mm=I(seed.marg.mat))
gsub("seed\\.mm\\.", "seed\\.mm", colnames(seed.marg.df))

seed.rem.df<-data.frame(group=substr(seed.removal, 1,1),
                        plot.type=substr(seed.removal, 12, nchar(seed.removal)),
                        stringsAsFactors = FALSE)

seed.rem.rows<-apply(seed.rem.df, 1, function(x){

  #group column
  if(x[1]==1){
  gr.col<-which(rowSums(seed.marg.mat[,2:8])==0)
  } else {gr.col<-which(seed.marg.mat[,2:8][,grepl(x[1], 
                                                   colnames(seed.marg.mat[,2:8]))] == 1)}
  
  #plot type column
  pl.col<-which(seed.marg.mat[,9:11][,grepl(x[2], colnames(seed.marg.mat[,9:11]))] == 1)
  
  return(gr.col[gr.col %in% pl.col])
  
})

seed.marg.mat<-seed.marg.mat[-seed.rem.rows,]

seed.preds<-data.frame(estimate=seed.marg.mat%*%fixef(seed.model),
                       se=sqrt(diag(seed.marg.mat %*% tcrossprod(vcov(seed.model),seed.marg.mat))))
seed.preds$plot.type= "Remnant" 

plot.type.ind<-which(seed.marg.mat[,9:11] == 1, arr.ind=TRUE)
seed.preds$plot.type[plot.type.ind[,1]] =
  c("Planting", "Young regrowth", "Old regrowth")[plot.type.ind[,2]]

seed.preds$group<-1
colnames(seed.marg.mat)
group.ind<-which(seed.marg.mat[,2:8] == 1, arr.ind=TRUE)
seed.preds$group[group.ind[,1]] = (2:8)[group.ind[,2]]

#                                     growth data prep ####
growth.coef.df<-as.data.frame(summary(growth.model)$tTable)

growth.coef.df$plot.type<-substr(rownames(growth.coef.df),
                               regexpr("R|P|O|Y", rownames(growth.coef.df)),
                               nchar(rownames(growth.coef.df)))
growth.coef.df$plot.type[!growth.coef.df$plot.type %in% c("Planting", 
                                                      "Young regrowth", 
                                                      "Old regrowth")] = "Remnant"

growth.coef.df$group<-substr(rownames(growth.coef.df),
                           regexpr("growth\\.group", rownames(growth.coef.df))+12,
                           regexpr("growth\\.group", rownames(growth.coef.df))+12)
growth.coef.df$group[is.na(as.numeric(growth.coef.df$group))]=1

# calculate and write comparison to remnants 
growth.cut.mat<-model.matrix(growth.model)
colnames(growth.cut.mat)<-gsub("growth.mm", "", colnames(growth.cut.mat))
growth.cut.mat[,"plot.area"]=0
growth.cut.mat<-unique.matrix(growth.cut.mat, MARGIN=1)
growth.marg.mat<-growth.cut.mat # save this to predict actual coefficients

# create contrast matrix to compare group estimates to remnants
# rows
growth.cont.mat<-growth.cut.mat
growth.cont.mat[,1:8]<-0
growth.cont.mat<-growth.cont.mat[rowSums(growth.cont.mat)>0,]

# get group names and plot types for each row of marginal matrix
ind<-which(growth.cont.mat !=0, arr.ind=TRUE)
ind<-ind[order(ind[,1]),]
plot.type.ind<-ind[ind[,2] %in% 9:11,]
group.ind<-ind[!ind[,2] %in% 9:11,]

plot.types<-gsub("plot.type", "", colnames(growth.cont.mat)[plot.type.ind[,2]])

groups<-rep(1, length(plot.types))
groups[group.ind[,1]]<-substr(colnames(growth.cont.mat)[group.ind[,2]],
                              nchar("growth.group")+1,
                              regexpr(":", colnames(growth.cont.mat)[group.ind[,2]])-1)

growth.marg.groups<-cbind(groups, plot.types)
growth.keeps<-!duplicated(paste0(growth.marg.groups[,1],
                               growth.marg.groups[,2]))

growth.marg.groups<-growth.marg.groups[growth.keeps,]
growth.cont.mat<-growth.cont.mat[growth.keeps,]

growth.contrast<-summary(glht(growth.model, linfct=growth.cont.mat))

growth.contrast.df<-data.frame(group=growth.marg.groups[,1],
                             plot.type=growth.marg.groups[,2],
                             estimate=growth.contrast$test$coefficients,
                             se=growth.contrast$test$sigma,
                             t=growth.contrast$test$tstat,
                             p=growth.contrast$test$pvalues)

growth.contrast.df[,-(1:2)] = round(growth.contrast.df[,-(1:2)], 3)
write.csv(growth.contrast.df, "./Outputs/growth contrast test.csv")

# get plot coefficients and error

growth.rem.df<-data.frame(group=substr(growth.removal, 1,1),
                        plot.type=substr(growth.removal, 12, nchar(growth.removal)),
                        stringsAsFactors = FALSE)

growth.rem.rows<-apply(growth.rem.df, 1, function(x){
  
  #group column
  if(x[1]==1){
    gr.col<-which(rowSums(growth.marg.mat[,2:8])==0)
  } else {gr.col<-which(growth.marg.mat[,2:8][,grepl(x[1], 
                                                   colnames(growth.marg.mat[,2:8]))] == 1)}
  
  #plot type column
  pl.col<-which(growth.marg.mat[,9:11][,grepl(x[2], colnames(growth.marg.mat[,9:11]))] == 1)
  
  return(gr.col[gr.col %in% pl.col])
  
})

growth.marg.mat<-growth.marg.mat[-growth.rem.rows,]

growth.preds<-data.frame(estimate=growth.marg.mat%*%fixef(growth.model),
                       se=sqrt(diag(growth.marg.mat %*% tcrossprod(vcov(growth.model), 
                                                              growth.marg.mat))))

growth.preds$plot.type= "Remnant" 
colnames(growth.marg.mat)
plot.type.ind<-which(growth.marg.mat[,9:11] == 1, arr.ind=TRUE)
growth.preds$plot.type[plot.type.ind[,1]] =
  c("Planting", "Young regrowth", "Old regrowth")[plot.type.ind[,2]]

growth.preds$group<-1
group.ind<-which(growth.marg.mat[,2:8] == 1, arr.ind=TRUE)
growth.preds$group[group.ind[,1]] = (2:8)[group.ind[,2]]

# functional group mean data for plotting groups
temp.species<-unique.species
temp.species$group<-growth.split[match(names(growth.split), rownames(temp.species))]

growth.SLA<-tapply(temp.species$SLA, temp.species$group, mean)
growth.SLA<-3 * (growth.SLA/max(growth.SLA))
saveRDS(growth.SLA, "./Outputs/SLA_group_means.RDS")

growth.WD<-tapply(temp.species$WD, temp.species$group, mean)
growth.WD<- 3 *(growth.WD)/max(growth.WD)
saveRDS(growth.WD, "./Outputs/wood_dens_group_means.RDS")

growth.MH<-tapply(temp.species$MH, temp.species$group, mean)
growth.MH<- 3*(growth.MH/max(growth.MH))
saveRDS(growth.MH, "./Outputs/max_height_group_means.RDS")

#                                     plot setup ####
pdf(paste0("./Plots/trait seed plot ",
           Sys.Date(),
           ".pdf"), height=9, width=3.74, useDingbats=FALSE)
# 
# split.screen(rbind(c(0.05,0.65,0.555,0.98),
#                    c(0.5,0.99,0.555,0.98),
#                    
#                    c(0.05,0.65, 0.05, 0.475),
#                    c(0.5,0.99,0.05,0.475),
#                    
#                    c(0.05,0.99,0.555,0.98),
#                    c(0.05,0.99,0.05,0.475)))

split.screen(rbind(c(0.05,0.65,0.65,0.98),
                   c(0.5,0.99,0.65,0.98),
                   
                   c(0.05,0.65, 0.16, 0.57),
                   c(0.5,0.99,0.16,0.57),
                   
                   c(0.05,0.99,0.65,0.98),
                   c(0.05,0.99,0.16,0.57),
                   
                   c(0.05,0.99,0.01,0.13)))

#                                     seed dendro ####
screen(1)
par(mar=c(0,0,0,5), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))

plot(dendrapply(as.dendrogram(seed.dend), 
                function(x) { attr(x, "height") <- log(attr(x, "height")+1); x }), 
     horiz=TRUE, yaxs="i", leaflab="none", axes=FALSE, ylab="", xlab="",
     ylim=c(-1,68.50), xlim=c(3,0), lwd=0.5, hang=1)

dend.pars<-par("usr")
dend.nodes<-as.data.frame(cbind(get_nodes_xy(seed.dend, type = c("rectangle"), center = FALSE,
                                             horiz = TRUE),
                                get_nodes_attr(seed.dend, attribute="label")))

dend.nodes$V1<-as.numeric(as.character(dend.nodes$V1))
dend.nodes<-dend.nodes[!is.na(dend.nodes$V3),]
dend.nodes$group<-seed.split[match(dend.nodes$V3, names(seed.split))]
dend.groups<-do.call("rbind", lapply(split(dend.nodes, f=dend.nodes$group), 
                                    function(x){c(min(x$V1),max(x$V1))}))
group.centers<-tapply(dend.nodes$V1, dend.nodes$group, mean)

text(x=relative.axis.point(0.08, "x"), 
     y=sort(group.centers), labels=8:1, cex=1.2, font=2,
     col=c("grey80", "grey50"))

seed.group.decode<-cbind(names(group.centers[order(group.centers)]),
                         8:1)
write.csv(seed.group.decode, "./Outputs/seed.groups.decoded.csv")

par(xpd=NA)
rect(xleft=par("usr")[1], xright=relative.axis.point(2.76, "x"),
     ybottom=dend.groups[names(sort(group.centers))[c(2,4,6,8)],1]-0.5, 
     ytop=dend.groups[names(sort(group.centers))[c(2,4,6,8)],2]+0.5,
     col=rgb(0,0,0,0.1), border=NA)

rect(xleft=-0.2, xright=-0.3, ybottom=group.centers-0.5, ytop=group.centers+2.5, lty="11")
rect(xleft=-0.2, xright=-0.3, ybottom=group.centers-0.5, ytop=(group.centers-0.5)+seed.SM, 
     col="black")
segments(x0=-0.15, x1=-0.35, y0=group.centers-0.5, y1=group.centers-0.5)
text(x=-0.25, y=group.centers-1.5, labels="M", adj=0.5, cex=0.8)
text(x=-0.6, y=group.centers-1.5, labels=seed.class, adj=0.5, cex=0.8)

picture.ys<-0.555 + ((par("usr")[3] + group.centers)/par("usr")[4])*(0.98-0.555)

PostScriptTrace("./Data/wind.ps")
wind<-readPicture("wind.ps.xml")
PostScriptTrace("./Data/fleshy.ps")
fleshy<-readPicture("fleshy.ps.xml")
PostScriptTrace("./Data/elaisome.ps")
elaisome<-readPicture("elaisome.ps.xml")
PostScriptTrace("./Data/unassisted.ps")
unassisted<-readPicture("unassisted.ps.xml")

lapply(picture.ys[seed.class=="Wi"], function(x){
  grid.picture(wind, 
               x=0.435, 
               y=x+0.03, 
               width=0.03, height=0.03)
})

lapply(picture.ys[seed.class=="Fl"], function(x){
  grid.picture(fleshy, 
               x=0.435, 
               y=x+0.055, 
               width=0.04, height=0.04)
})

lapply(picture.ys[seed.class=="El"], function(x){
  grid.picture(elaisome, 
               x=0.435, 
               y=x+0.055, 
               width=0.04, height=0.04)
})

lapply(picture.ys[seed.class=="Un"], function(x){
grid.picture(unassisted, 
             x=0.435, 
             y=x+0.08, 
             width=0.03, height=0.03)
})

par(xpd=FALSE)

axis(side=1, mgp=c(3,0,0))
dendro.pars<-par("usr")
mtext(side=1, line=0.75, text="ln(Ward's minimum variance)")
close.screen(1)

#                                     seed points ####
screen(2)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot(x=NULL, y=NULL, xlim=c(-5.6, 1.25), ylim=dend.pars[3:4], 
     yaxs="i", axes=FALSE, ylab="", xlab="")

axis(side=1, at=log(c(0.001,0.01,0.1,1,2,3)), labels=c(0.001,0.01,0.1,1,2,3),
     mgp=c(3,0,0))
axis(side=1, at=log(c(seq(0.001,0.01,0.001), seq(0.01,0.1,0.01), seq(0.1,1,0.1),
                      1,3,1)), labels=NA, tck=-0.01)
mtext(side=1, line=0.75, text="Mean species per plot")

plot.type.pos<-group.centers[seed.preds$group] + 
                      c(-2.25, -0.75, 0.75, 2.25)[as.factor(seed.preds$plot.type)]

segments(x0=seed.preds$estimate + 1.96*seed.preds$se,
         x1=seed.preds$estimate - 1.96*seed.preds$se,
         y0=plot.type.pos, 
         y1=plot.type.pos, lwd=1.5,
         col=rgb(colours[match(seed.preds$plot.type, colours[,1]),2],
                 colours[match(seed.preds$plot.type, colours[,1]),3],
                 colours[match(seed.preds$plot.type, colours[,1]),4]))

points(x=seed.preds$estimate,
       y=plot.type.pos,
       pch=c(25,22,21,24)[as.factor(seed.preds$plot.type)],
       cex=0.8,
       bg=rgb(colours[match(seed.preds$plot.type, colours[,1]),2],
              colours[match(seed.preds$plot.type, colours[,1]),3],
              colours[match(seed.preds$plot.type, colours[,1]),4]))
close.screen(2)

#                                     growth dendro ####
screen(3)

par(mar=c(0,0,0,5), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))

growth.dend<-as.dendrogram(growth.groups)
growth.dend <- set(growth.dend, "labels_cex", 0.5)

plot(dendrapply(as.dendrogram(growth.dend), 
                function(x) { attr(x, "height") <- log(attr(x, "height")+1); x }), 
     horiz=TRUE, axes=FALSE, yaxs="i", leaflab="none", ylab="", xlab="",
     ylim=c(-1,68.50), xlim=c(2.5,0), lwd=0.5)

dend.pars<-par("usr")
dend.nodes<-as.data.frame(cbind(get_nodes_xy(growth.dend, type = c("rectangle"), center = FALSE,
                                             horiz = TRUE),
                                get_nodes_attr(growth.dend, attribute="label")))

dend.nodes$V1<-as.numeric(as.character(dend.nodes$V1))
dend.nodes<-dend.nodes[!is.na(dend.nodes$V3),]
dend.nodes$group<-growth.split[match(dend.nodes$V3, names(growth.split))]
dend.groups<-do.call("rbind", lapply(split(dend.nodes, f=dend.nodes$group), 
                                     function(x){c(min(x$V1),max(x$V1))}))
group.centers<-tapply(dend.nodes$V1, dend.nodes$group, mean)

growth.group.decode<-cbind(names(group.centers[order(group.centers)]),
                         8:1)
write.csv(growth.group.decode, "./Outputs/growth.groups.decoded.csv")

par(xpd=NA)
rect(xleft=par("usr")[1], xright=relative.axis.point(2.76, "x"),
     ybottom=dend.groups[names(sort(group.centers))[c(2,4,6,8)],1]-0.5, 
     ytop=dend.groups[names(sort(group.centers))[c(2,4,6,8)],2]+0.5,
     col=rgb(0,0,0,0.1), border=NA)

text(x=relative.axis.point(0.08, "x"),
     y=sort(group.centers), labels=8:1, cex=1.2, font=2,
     col=c("grey80", "grey50"))

# SLA
rect(xleft=-0.2, xright=-0.3, ybottom=group.centers-0.5, 
     ytop=group.centers+2.5,
     lty="11")
rect(xleft=-0.2, xright=-0.3, ybottom=group.centers-0.5, 
     ytop=(group.centers-0.5)+growth.SLA, 
     col="black")

# WD
rect(xleft=-0.35, xright=-0.45, ybottom=group.centers-0.5, ytop=group.centers+2.5, lty="11")
rect(xleft=-0.35, xright=-0.45, ybottom=group.centers-0.5, 
     ytop=(group.centers-0.5)+growth.WD, 
     col="black")

# MH
rect(xleft=-0.5, xright=-0.6, ybottom=group.centers-0.5, ytop=group.centers+2.5, lty="11")
rect(xleft=-0.5, xright=-0.6, ybottom=group.centers-0.5, 
     ytop=(group.centers-0.5)+growth.MH, 
     col="black")

n.groups<-unique(plot.species$growth.group[plot.species$n.fixer==1])
saveRDS(n.groups, "./Outputs/n_fixer.RDS")

draw.ellipse(x=rep(-0.775, sum(names(group.centers) %in% n.groups)), 
             y=group.centers[names(group.centers) %in% n.groups]+1,
            a=0.1, b=1, lwd=1)
text(x=-0.775, y=group.centers[names(group.centers) %in% n.groups]+1, labels="N", adj=0.5, cex=0.8)

segments(x0=-0.15, x1=-0.65, y0=group.centers-0.5, y1=group.centers-0.5)
text(x=-0.25, y=group.centers-1.5, labels="S", adj=0.5, cex=0.8)
text(x=-0.4, y=group.centers-1.5, labels="W", adj=0.5, cex=0.8)
text(x=-0.55, y=group.centers-1.5, labels="H", adj=0.5, cex=0.8)
               
par(xpd=FALSE)

axis(side=1, mgp=c(3,0,0))
dendro.pars<-par("usr")
mtext(side=1, line=0.75, text="ln(Ward's minimum variance)")
close.screen(3)

#                                     growth points ####
screen(4)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot(x=NULL, y=NULL, xlim=c(-5.6, 1.25), ylim=dend.pars[3:4], 
     yaxs="i", axes=FALSE, xlab="", ylab="")

axis(side=1, at=log(c(0.001,0.01,0.1,1,2,3)), labels=c(0.001,0.01,0.1,1,2,3),
     mgp=c(3,0,0))
axis(side=1, at=log(c(seq(0.001,0.01,0.001), seq(0.01,0.1,0.01), seq(0.1,1,0.1),
                      1,3,1)), labels=NA, tck=-0.01)
mtext(side=1, line=0.75, text="Mean species per plot")

plot.type.pos<-group.centers[growth.preds$group] + 
  c(-2.25, -0.75, 0.75, 2.25)[as.factor(growth.preds$plot.type)]
#plot.type.pos[growth$group==5]<-1.5

segments(x0=growth.preds$estimate + 1.96*growth.preds$se,
         x1=growth.preds$estimate - 1.96*growth.preds$se,
         y0=plot.type.pos, y1=plot.type.pos, lwd=1.5,
         col=rgb(colours[match(growth.preds$plot.type, colours[,1]),2],
                 colours[match(growth.preds$plot.type, colours[,1]),3],
                 colours[match(growth.preds$plot.type, colours[,1]),4]))

points(x=growth.preds$estimate,
       y=plot.type.pos,
       pch=c(25,22,21,24)[as.factor(growth.preds$plot.type)],
       cex=0.8,
       bg=rgb(colours[match(growth.preds$plot.type, colours[,1]),2],
              colours[match(growth.preds$plot.type, colours[,1]),3],
              colours[match(growth.preds$plot.type, colours[,1]),4]))
close.screen(4)

#                                     other screens ####
screen(5)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot.new()
#box()
par(xpd=NA)
text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(1.01, "y"), labels="(a) Seed dispersal", font=2, adj=0)
par(xpd=FALSE)
close.screen(5)

screen(6)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot.new()
#box()
par(xpd=NA)
text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(1.01, "y"), labels="(b) Growth and structure", 
     font=2, adj=0)
par(xpd=FALSE)
close.screen(6)

screen(7)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot(x=NULL, y=NULL, xlim=c(0,12), ylim=c(0,1), axes=FALSE)

gap<-seq(0,0.15, length=2)[2]

text(x=1.25, y=0.825, labels="Categorial traits", font=2)

text(x=0.5, y=0.8 - rep(gap, 5)*c(1:5),
     labels=c("Wind dispersed", "Elaisome", "Fleshy", "Unassisted",
              "Nitrogen fixer"), adj=0)

PostScriptTrace("./Data/wind.ps")
wind<-readPicture("wind.ps.xml")
grid.picture(wind, 
             x=0.1, 
             y=0.09, 
             width=0.03, height=0.03)

PostScriptTrace("./Data/elaisome.ps")
elaisome<-readPicture("elaisome.ps.xml")
grid.picture(elaisome, 
             x=0.1, 
             y=0.07, 
             width=0.03, height=0.03)

PostScriptTrace("./Data/fleshy.ps")
 fleshy<-readPicture("fleshy.ps.xml")
 grid.picture(fleshy, 
              x=0.1, 
              y=0.055, 
              width=0.03, height=0.03)
 
PostScriptTrace("./Data/unassisted.ps")
unassisted<-readPicture("unassisted.ps.xml")
grid.picture(unassisted, 
                x=0.1, 
                y=0.035, 
                width=0.03, height=0.03)

draw.ellipse(x=0.2, 
             y=0.8-gap*5,
             a=0.2, b=0.06, lwd=1)
text(x=0.2, y=0.8-gap*5, labels="N", adj=0.5, cex=0.8)

text(x=5.5, y=0.825, labels="Continuous traits", font=2)

text(x=4.5, y=0.8 - rep(gap, 4)*c(1:4),
     labels=c("Maximum height", "Seed mass", 
              "Specific leaf area", "Wood density"), adj=0)

text(x=4, y=0.8 - rep(gap, 4)*c(1:4),
     labels=c("H", "M", "S", "W"), font=2)

text(x=10, y=0.825, labels="Stand categories", font=2)

text(x=9.5, y=0.8 - rep(gap, 4)*c(1:4),
     labels=c("Remnant", "Planting",
              "Young regrowth", "Old regrowth"), adj=0)

points(x=rep(9, 4), y=0.8 - rep(gap, 4)*c(1:4),
       pch=c(21,22,24,25),
       bg=rgb(colours[1:4,5],
                 colours[1:4,6],
                 colours[1:4,7]),
       cex=1.2)
close.screen(7)

close.screen(all.screens=TRUE)
dev.off()

#                                     Species group table ####

sp.table<-data.frame(species=names(seed.split),
                     seed.group.raw=seed.split)
sp.table$growth.group.raw<-growth.split[match(sp.table$species,
                                              names(growth.split))]

sp.table$family<-species$fam[match(sp.table$species, species$species)]

write.csv(sp.table,
          "./Outputs/species functional group table.csv")

#                               SMALL PLANTS IN PLANTINGS PLOT ####
smallplants<-droplevels(plant[plant$plot.type=="Planting" &
                                plant$tot.ag.mass <= 40, ])

smallplants$growth.group<-plot.species$growth.group[match(smallplants$species,
                                                          plot.species$species)]
  
small.groups<-as.matrix(table(smallplants$species, smallplants$growth.group))
saveRDS(small.groups, "./Outputs/small_planting_groups_mat.RDS")


pdf("./Plots/small plants breakdown.pdf", height=3, width=6)

par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(3,3.5,1,1), ps=8, tck=-0.025, mgp=c(3,0.5,0),
    las=1)

plot(x=NULL, y=NULL, xlim=c(0.5,8.5), ylim=c(0,1000), axes=FALSE,
     xlab="", ylab="")
axis(side=1, at=1:8, mgp=c(3,0.2,0))
axis(side=2)
mtext(side=1, at=par("usr")[2], line=1.25, text="Growth functional group")
mtext(side=2, line=2, text="Count of small plants (< 40 kg)", las=0)
box()

x<-"5"
sapply(colnames(small.groups), function(x){
  
  temp.sort<-sort(small.groups[,x])
  temp.sort<-temp.sort[temp.sort > 0]
  
  temp.mat<-c(0, temp.sort)
  
  for(i in 1:length(temp.mat)){
    
    rect(ybottom=ifelse(i>1, sum(temp.mat[1:(i-1)]), 0),
         ytop=sum(temp.mat[1:i]),
         xleft=as.numeric(x)-0.3,
         xright=as.numeric(x)+0.3,
         col=c("grey80", "grey50")[(i %% 2) + 1])

  }
  
  top.sort<-temp.sort>=20
  
  if(sum(top.sort)>0){
  text(labels=names(temp.sort[top.sort]),
       x=as.numeric(x),
       y=sapply(names(temp.sort[top.sort]), function(y){
         sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort],
       pos=2, font=3, cex=0.8)
    
    segments(x0=as.numeric(x), x1=as.numeric(x)-0.3,
             y0=sapply(names(temp.sort[top.sort]), function(y){
               sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort],
             y1=sapply(names(temp.sort[top.sort]), function(y){
               sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort])
  }
  
  })

segments(x0=5-0.3, x1=5+0.3, y0=0, y1=0)
text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(a) Plantings", font=2, adj=0)

# ANd compared to young regrowth....
smallreg<-droplevels(plant[plant$plot.type=="Young regrowth" &
                                plant$tot.ag.mass <= 40, ])

smallreg$growth.group<-plot.species$growth.group[match(smallreg$species,
                                                       plot.species$species)]

small.reg.groups<-as.matrix(table(smallreg$species, smallreg$growth.group))
small.reg.roups <- saveRDS(small.reg.groups, "./Outputs/small_regrowth_groups_mat.RDS")
colSums(small.reg.groups)

plot(x=NULL, y=NULL, xlim=c(0.5,8.5), ylim=c(0,1000), axes=FALSE)
axis(side=1, at=1:8, mgp=c(3,0.2,0))
axis(side=2, labels=NA)
box()

x<-"5"
sapply(colnames(small.reg.groups), function(x){
  
  print(x)
  
  temp.sort<-sort(small.reg.groups[,x])
  temp.sort<-temp.sort[temp.sort > 0]
  
  temp.mat<-c(0, temp.sort)
  
  for(i in 1:length(temp.mat)){
    
    rect(ybottom=ifelse(i>1, sum(temp.mat[1:(i-1)]), 0),
         ytop=sum(temp.mat[1:i]),
         xleft=as.numeric(x)-0.3,
         xright=as.numeric(x)+0.3,
         col=c("grey80", "grey50")[(i %% 2) + 1])
    
  }
  
  top.sort<-temp.sort>=20
  
  if(sum(top.sort)>0){
    text(labels=names(temp.sort[top.sort]),
         x=as.numeric(x),
         y=sapply(names(temp.sort[top.sort]), function(y){
           sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort],
         pos=2, font=3, cex=0.8)
    
    segments(x0=as.numeric(x), x1=as.numeric(x)-0.3,
             y0=sapply(names(temp.sort[top.sort]), function(y){
               sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort],
             y1=sapply(names(temp.sort[top.sort]), function(y){
               sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort])
  }
  
})

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(b) Young regrowth", font=2, adj=0)

dev.off()

#                               TABLE OF LAND ZONES ####

RE.table <- as.matrix(table(site$preRE, site$plot.type))
write.csv(RE.table, "./Outputs/RE table.csv")

#               Write data files ####
write.csv(cut.data, "./Outputs/size class data.csv")
write.csv(growth.count, "./Outputs/growth functional group data.csv")
write.csv(seed.count, "./Outputs/seed functional group data.csv")
write.csv(site, "./Outputs/plot-level data.csv")
#             REGION MAP ####

ausmap<-readShapeSpatial("/home/timothy/University files - offline/PhD - offline/Shape files/Australia outline/AUS_adm/AUS_adm1.shp")
rainfall<-raster("/home/timothy/University files - offline/PhD - offline/Shape files/prec.annual.tif")

# cur.RE<-readShapeSpatial("/home/timothy/University files - offline/PhD - offline/Shape files/Queensland RE/croppedRE.shp")
# rem.areas <- cur.RE[cur.RE@data$RE != "non-rem",]
# rem.rast <- rasterize(rem.areas, rain.sub)
# writeRaster(rem.rast, "rasterize remnants.tif", format="GTiff", overwrite=TRUE)
rem.rast <- raster("rasterize remnants.tif")

library(rworldmap)
worldmap<-getMap(resolution="low")
smallaus<-worldmap[worldmap@data$ADMIN.1=="Australia" &
                     !is.na(worldmap@data$ADMIN.1),]
plot(smallaus)
summary(site[,c("long","lat")])

bbox<-extent(150, 154, -29, -24)
ausmap.sub<-crop(ausmap, bbox)

ausmap.simp<-disaggregate(ausmap)

sort(sapply(ausmap.simp@polygons,
            function(x){slot(x, "area")}), decreasing=TRUE)[1:10]

ausmap.simp1<-ausmap.simp[sapply(ausmap.simp@polygons,
                                 function(x){slot(x, "area")}) >1, ]

rain.sub<-crop(rainfall, bbox)
rain.smooth<- focal(rain.sub, w=matrix(1, 9, 9), mean)

data(world.cities)
towns<-world.cities[world.cities$lat > -29 &
                      world.cities$lat < -24 &
                      world.cities$long > 150 &
                      world.cities$long < 154, ]

pdf(paste0("./Plots/map trial ", Sys.Date(), ".pdf"), height=4, width=4)

split.screen(rbind(c(0.125,0.95,0.1,0.95),
                   c(0.125,0.95,0.1,0.95),
                   c(0.125,0.95,0.1,0.95),
                   c(0.7,0.95,0.7,0.95)))

screen(1)
par(mar=c(0,0,0,0), ps=6, las=1, tck=-0.01, mgp=c(3,0.3,0))

plot(x=NULL, y=NULL, xlim=c(150.75, 153.5), ylim=c(-27.65, -25.25), 
     axes=FALSE, xlab="", ylab="", asp=1)

rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
     border=NA, col="grey90")

plot(ausmap.sub, axes=FALSE, bg="grey90", col="white",
     xaxs="i", yaxs="i", add=TRUE)

axis(side=1, mgp=c(3,0,0),
     at=seq(151,153, 0.5),
     labels=parse(text=paste(seq(151,153, 0.5), "*degree~E", sep="")))

axis(side=2, at=seq(-27.5,-25.5, 0.5),
     labels=parse(text=paste(seq(-27.5,-25.5, 0.5), "*degree~S", sep="")))

usrs<-par("usr")

close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=6, las=1, tck=-0.01, mgp=c(3,0.3,0))

plot(x=NULL, y=NULL, xlim=c(150.75, 153.5), ylim=c(-27.65, -25.25), 
     axes=FALSE, xlab="", ylab="", asp=1)

plot(rem.rast, asp=1,
     col="grey70", border=NA,
     legend=FALSE, axes=FALSE, add=TRUE, xlim=usrs[1:2], ylim=usrs[3:4],
     maxpixels=1000000)

close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=6, las=1, tck=-0.01, mgp=c(3,0.3,0))

plot(x=NULL, y=NULL, xlim=c(150.75, 153.5), ylim=c(-27.65, -25.25), 
     axes=FALSE, xlab="", ylab="", asp=1)

plot(ausmap.sub, axes=FALSE,
     xaxs="i", yaxs="i", add=TRUE)

contour(crop(rain.smooth, extent(par("usr")[1],
                                 par("usr")[2],
                                 par("usr")[3],
                                 par("usr")[4])),
        levels=c(600,800,1000,1200,1400),
        col="black", method="edge", add=TRUE, lty="31")


with(site[!duplicated(site$planting),], 
     points(lat ~ long,
            pch=c(21,22,24,25)[plot.type],
            bg=rgb(colours[1:4,5],
                   colours[1:4,6],
                   colours[1:4,7])[match(plot.type,
                                         colours$plot.type)]))

with(towns[towns$name %in% c("Kingaroy", "Toowoomba", "Dalby", "Gympie", "Brisbane"),],
     points(lat ~ long, pch=16, cex=0.8))

with(towns[towns$name %in% c("Kingaroy", "Toowoomba", "Dalby", "Gympie"),],
     text(long, lat, name, pos=c(2,2,2,4), offset=0.25))

with(towns[towns$name == "Brisbane",],
     text(long-0.075, lat+0.07, labels=name))

rect(xleft=par("usr")[1],
     xright=151.6,
     ybottom=par("usr")[3],
     ytop=-27.35,
     col=rgb(1,1,1,0.9))

legend(x=150.6,
       y=-27.28,
       legend=c("Remnant", "Planting", "Young regrowth", "Old regrowth", "Uncleared vegetation"),
       pch=c(21,22,24,25,15),
       pt.bg=rgb(colours[1:4,5],
                 colours[1:4,6],
                 colours[1:4,7]),
       col = c(rep("black", 4), "grey70"),
       pt.cex=1.2,
       y.intersp = 0.65,
       x.intersp = 0.8,
       bty="n")

rect(xleft=relative.axis.point(0.515, "x"),
     xright=relative.axis.point(0.765, "x"),
     ybottom=par("usr")[3],
     ytop=relative.axis.point(0.07, "y"),
     border=NA, col=rgb(1,1,1,0.9))

scalebar(d=50, xy=c(relative.axis.point(0.55, "x"),
                    relative.axis.point(0.025, "y")),
         type="bar", label=c("0 km", NA, "50 km"))


box()
close.screen(3)

screen(4)
par(mar=c(0,0,0,0), ps=8, las=1, tck=-0.01, mgp=c(3,0.3,0))
plot(ausmap.simp1, col="white", bg="grey90", lwd=0.5)
rect(xleft=usrs[1], xright=usrs[2]+1, ybottom=usrs[3], ytop=usrs[4],
     border=rgb(1,1,1,0.9), col="black")
plot(ausmap.simp1, lwd=0.5, add=TRUE)
box()
close.screen(4)

close.screen(all.screens=TRUE)
dev.off()
