library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(raster)
library(fasterize)
#Outside of Jasper National Park
beetles<-st_read("0_data/mpb/largerBeetleLeftStanding.shp")
beetlesTMERC<-beetles%>%
  st_transform(crs="+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

points<-st_read("0_data/mpb/All_MPB_Locations.shp")
pointsTMERC<-points%>%
  st_transform(crs="+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")

st_bbox(pointsTMERC)
studyArea<-st_as_sfc(st_bbox(pointsTMERC))
st_write(studyArea, "0_data/mpb/studyArea.shp", append=FALSE)

studyArea<-st_read("0_data/mpb/studyArea.shp")
studyArea<-studyArea%>%
  st_transform(crs=crs(beetlesTMERC))

beetles_studyArea<-st_crop(beetlesTMERC, studyArea)
plot(beetlesTMERC["betl_yr"])
plot(beetles_studyArea["betl_yr"])


#corrupt geometries detected the first time loops were run for this file
nrow(beetles_studyArea[st_is_valid(beetles_studyArea)==FALSE,])#13
#beetles_studyArea<-beetles_studyArea[st_is_valid(beetles_studyArea)==TRUE,]

#convert to raster. Assign most recent beetle year to each pixel

#f <- system.file("ex/lux.shp", package="terra")
byear<-beetles_studyArea["betl_yr"]
v <- vect(byear)#instead of f
#extent      : 493946.1, 813585.9, 6056702, 6401906
#x range=319639.8 y range=345204
r <- rast(v, res=1000)
zBeetle <- rasterize(v, r, field="betl_yr")
plot(zBeetle)#what to do about zero values?

fyear<-beetles_studyArea["YEAR_1"]
v <- vect(fyear)
zFire<- rasterize(v, r, field="YEAR_1")
plot(zFire)
zFire[zFire<1800]<-1800#set default maximum age to 200?

hyear<-beetles_studyArea["HrvstYr"]
v <- vect(hyear)
zHarvest <- rasterize(v, r, field="HrvstYr")
plot(zHarvest)
zHarvest[zHarvest<1800]<-1800#set default maximum age to 200?

#create buffers around points
b150<-st_buffer(pointsTMERC,150)#150 m circle buffer
b300<-st_buffer(pointsTMERC,300)#300 m circle buffer
#b150.300<-st_difference(b300, b150)#150-300 m ring

#extract mean harvest year value in each buffer
pointList<-list()
for (i in 1:nrow(pointsTMERC)){
  pointData<-pointsTMERC[i,]
  #Extracted and summarized within 50 m of point
  try(meanBeetleYear_150m <- terra::extract(zBeetle, b150[i,], fun="max", na.rm=TRUE))
  try(pointData$BY150m<-as.integer(meanBeetleYear_150m$last))
  try(meanHarvestYear_150m <- terra::extract(zHarvest, b150[i,], fun="max", na.rm=TRUE))
  try(pointData$HY150m<-as.integer(meanHarvestYear_150m$last))
  try(meanFireYear_150m <- terra::extract(zFire, b150[i,], fun="max", na.rm=TRUE))
  try(pointData$FY150m<-as.integer(meanFireYear_150m$last))
  
  try(meanBeetleYear_300m <- terra::extract(zBeetle, b300[i,], fun="max", na.rm=TRUE))
  try(pointData$BY300m<-as.integer(meanBeetleYear_300m$last))
  try(meanHarvestYear_300m <- terra::extract(zHarvest, b300[i,], fun="max", na.rm=TRUE))
  try(pointData$HY300m<-as.integer(meanHarvestYear_300m$last))
  try(meanFireYear_300m <- terra::extract(zFire, b300[i,], fun="max", na.rm=TRUE))
  try(pointData$FY300m<-as.integer(meanFireYear_300m$last))
  
  try(pointList[[i]]<-pointData)
  print(paste0("Line ",i," done."))
}
pointBeetleYear<-do.call(rbind,pointList)
summaries<-data.frame(pointBeetleYear)
write.csv(summaries, file="2_outputs/meanBeetleHarvestFireYear_NotJasper.csv")

#Inside Jasper National Park
st_layers("0_data/jnp-vri.gdb-20250116T164242Z-001/jnp-vri.gdb")
#VRI
VRIpoly<-st_read("0_data/jnp-vri.gdb-20250116T164242Z-001/jnp-vri.gdb", layer="JNP_VRI_DBO_VRI_Polygons")
#shapefile / simple features data frame
str(VRIpoly)

VRIdisturb<-st_read("0_data/jnp-vri.gdb-20250116T164242Z-001/jnp-vri.gdb", layer="JNP_VRI_DBO_Disturbance")
#regular data frame
str(VRIdisturb)
levels(as.factor(VRIdisturb$DISTURBANCE_START_YEAR))

VRIpoly2<-left_join(VRIpoly, VRIdisturb, by=c("VRI_ID"))
str(VRIpoly2)

levels(as.factor(VRIpoly2$DISTURBANCE_TYPE_CODE))#IBM = beetles
#filter to beetle-infested polygons

Beetlepoly<-VRIpoly2%>%
  filter(DISTURBANCE_TYPE_CODE=="IBM")
plot(Beetlepoly["DISTURBANCE_TYPE_CODE"])

#Now fire
JasperFire<-st_read("0_data/jnp-fire-perimeters-shapefile-20250116T164417Z-001/jnp-fire-perimeters-shapefile/JNP_Fire.DBO.shp")
str(JasperFire)
JasperFire<-JasperFire%>%
  st_transform(crs=crs(VRIpoly))
JasperFire$YEAR<-JasperFire$FIRE_YEAR
JasperFire$BurnSever<-as.factor(JasperFire$Burn_Sever)

#Now read in harvest
st_layers("0_data/WRR Areas_Emily Swerdfager_2024-01-19.gdb")

JasperHarvest<-st_read("0_data/WRR Areas_Emily Swerdfager_2024-01-19.gdb", layer="WRR_Complete")
str(JasperHarvest)
JasperHarvest<-JasperHarvest%>%
  st_transform(crs=crs(VRIpoly))

#Read in point counts
points<-st_read("0_data/MPB_points-20250116T164310Z-001/MPB_points/jnp-sites/jnp_sites.shp")
str(points)

#change CRS of points to that of VRI
pointsVRI<-points%>%
  st_transform(crs=crs(VRIpoly)) 

pointsVRI_150<-st_buffer(pointsVRI, 150)
pointsVRI_300<-st_buffer(pointsVRI, 300)

#Get area of beetle disturbance,
#mean and maximum disturbance year,
#and typical severity level

#functions
beetleAreaAndCategory<-function(pointfile, polyfile, maxdist){
  ssyr.sf.summary<-list()
  #ssyr.sf.demo<-pointfile[1:100,]
  for (i in 1:nrow(pointfile)){
    ssyr.sf.i<-pointfile[i,]#ssyr.sf.demo[i,]
    SS.i<-ssyr.sf.i$location
    YEAR.i<-ssyr.sf.i$year
    #LF.i<-ssyr.sf.i$LF
    #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
    polyfile.y<-polyfile%>%
      filter(YEAR<YEAR.i)#should this filtering occur on the polygons after isolating those within 150 and 565 m of each point count?
    #1b. also after time-filtering, create 150-m and 565-m buffers
    b150<-st_buffer(ssyr.sf.i,150)
    p.150<-st_intersection(polyfile.y, b150)
    p.150$Recalc_Area<-st_area(p.150)#recalculate areas of clipped polygons
    
    dom.150<-ifelse(nrow(p.150)==0, 0, as.numeric(p.150 %>% st_drop_geometry() %>% 
                                                    group_by(ATTACK_CODE) %>% 
                                                    summarise(area=sum(Recalc_Area)%>% as_tibble) %>%
                                                    arrange(desc(area)) %>% dplyr::select(ATTACK_CODE)%>%
                                                    filter(row_number()==1))-1)

    sev.150<-ifelse(nrow(p.150)==0, 0, as.numeric(p.150 %>% st_drop_geometry() %>% 
                                                    group_by(PEST_SEVERITY_CODE) %>% 
                                                    summarise(area=sum(Recalc_Area)%>% as_tibble) %>%
                                                    arrange(desc(area)) %>% dplyr::select(PEST_SEVERITY_CODE)%>%
                                                    filter(row_number()==1))-1)
    
    #calculate area amount after dissolving polygons
    #(in case of polygons overlapping before dissolution, causing overestimates of harvest)
    area.150 <- ifelse((nrow(p.150)==0), 0, (p.150 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    
    b300<-st_buffer(ssyr.sf.i,300)
    p.300<-st_intersection(polyfile.y, b300)
    p.300$Recalc_Area<-st_area(p.300)#recalculate areas of clipped polygons
    dom.300<-ifelse(nrow(p.300)==0, 0, as.numeric(p.300 %>% st_drop_geometry() %>% 
                                                    group_by(ATTACK_CODE) %>% 
                                                    summarise(area=sum(Recalc_Area)%>%  as_tibble) %>%
                                                    arrange(desc(area)) %>% dplyr::select(ATTACK_CODE)%>%
                                                    filter(row_number()==1))-1)
    
    sev.300<-ifelse(nrow(p.300)==0, 0, as.numeric(p.300 %>% st_drop_geometry() %>% 
                                                    group_by(PEST_SEVERITY_CODE) %>% 
                                                    summarise(area=sum(Recalc_Area)%>% as_tibble) %>%
                                                    arrange(desc(area)) %>% dplyr::select(PEST_SEVERITY_CODE)%>%
                                                    filter(row_number()==1))-1)
    
    
    #calculate area amount after dissolving polygons
    #(in case of polygons overlapping before dissolution, causing overestimates of harvest)
    area.300 <- ifelse((nrow(p.300)==0), 0, (p.300 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    b500<-st_buffer(ssyr.sf.i,500)#a bit under 1 square kilometre
    p.500<-st_intersection(polyfile.y, b500)
    p.500$Recalc_Area<-st_area(p.500)#recalculate areas of clipped polygons
    dom.500<-ifelse(nrow(p.500)==0, 0, as.numeric(p.500 %>% st_drop_geometry() %>% 
                                                    group_by(ATTACK_CODE) %>% 
                                                    summarise(area=sum(Recalc_Area)%>%  as_tibble) %>%
                                                    arrange(desc(area)) %>% dplyr::select(ATTACK_CODE)%>%
                                                    filter(row_number()==1))-1)
    
    sev.500<-ifelse(nrow(p.500)==0, 0, as.numeric(p.500 %>% st_drop_geometry() %>% 
                                                    group_by(PEST_SEVERITY_CODE) %>% 
                                                    summarise(area=sum(Recalc_Area)%>% as_tibble) %>%
                                                    arrange(desc(area)) %>% dplyr::select(PEST_SEVERITY_CODE)%>%
                                                    filter(row_number()==1))-1)
    
    #calculate area amount after dissolving polygons
    #(in case of polygons overlapping before dissolution, causing overestimates of harvest)
    area.500 <- ifelse((nrow(p.500)==0), 0, (p.500 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    
    
    ssyr.sf.summary[[i]]<-data.frame(SS=SS.i,#PKEY=PKEY.i,
                                     YEAR=YEAR.i,
                                     BEETLEAREA150=area.150,
                                     BEETLEAREA300=area.300,
                                     BEETLEAREA500=area.500,
                                     SEVERITY150=sev.150,
                                     DOMCLASS150=dom.150,
                                     SEVERITY300=sev.300,
                                     DOMCLASS300=dom.300,
                                     SEVERITY500=sev.500,
                                     DOMCLASS500=dom.500)
    print(paste0("point ",i," done"))
  }
  summaries<-do.call(bind_rows,ssyr.sf.summary)
  return(summaries)
}

Beetlepoly$YEAR<-2015
Beetlepoly$ATTACK_CODE<-as.factor(Beetlepoly$ATTACK_CODE)
Beetlepoly$PEST_SEVERITY_CODE<-as.factor(Beetlepoly$PEST_SEVERITY_CODE)
summary1<-beetleAreaAndCategory(pointsVRI,Beetlepoly,500)
#write.csv(summaries, "2_outputs/MPB_Damage_Jasper.csv")


FireAreaAndCategory<-function(pointfile, polyfile, maxdist){
  ssyr.sf.summary<-list()
  #ssyr.sf.demo<-pointfile[1:100,]
  for (i in 1:nrow(pointfile)){
    ssyr.sf.i<-pointfile[i,]#ssyr.sf.demo[i,]
    SS.i<-ssyr.sf.i$location
    YEAR.i<-ssyr.sf.i$year
    #LF.i<-ssyr.sf.i$LF
    #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
    polyfile.y<-polyfile%>%
      filter(YEAR<YEAR.i)#should this filtering occur on the polygons after isolating those within 150 and 565 m of each point count?
    #1b. also after time-filtering, create 150-m and 565-m buffers
    b150<-st_buffer(ssyr.sf.i,150)
    p.150<-st_intersection(polyfile.y, b150)
    p.150$Recalc_Area<-st_area(p.150)#recalculate areas of clipped polygons
    
    dom.150<-ifelse(nrow(p.150)==0, 0, as.numeric(p.150 %>% st_drop_geometry() %>% 
                                                    group_by(BurnSever) %>% 
                                                    summarise(area=sum(Recalc_Area)%>% as_tibble) %>%
                                                    arrange(desc(area)) %>% dplyr::select(BurnSever)%>%
                                                    filter(row_number()==1))-1)
    #summarise(prop.table(table(Recalc_Area)) %>% as.list %>% as_tibble) %>% 
    #mutate(
    #  across(-1, coalesce, 0),
    #  Categ_dominant = across(-1) %>% {names(.)[max.col(.)]}
    #)  
    #calculate area amount after dissolving polygons
    #(in case of polygons overlapping before dissolution, causing overestimates of harvest)
    area.150 <- ifelse((nrow(p.150)==0), 0, (p.150 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    
    b300<-st_buffer(ssyr.sf.i,300)
    p.300<-st_intersection(polyfile.y, b300)
    p.300$Recalc_Area<-st_area(p.300)#recalculate areas of clipped polygons
    dom.300<-ifelse(nrow(p.300)==0, 0, as.numeric(p.300 %>% st_drop_geometry() %>% 
                                                    group_by(BurnSever) %>% 
                                                    summarise(area=sum(Recalc_Area)%>%  as_tibble) %>%
                                                    arrange(desc(area)) %>% dplyr::select(BurnSever)%>%
                                                    filter(row_number()==1))-1)
    
    #calculate area amount after dissolving polygons
    #(in case of polygons overlapping before dissolution, causing overestimates of harvest)
    area.300 <- ifelse((nrow(p.300)==0), 0, (p.300 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    b500<-st_buffer(ssyr.sf.i,500)#a bit under 1 square kilometre
    p.500<-st_intersection(polyfile.y, b500)
    p.500$Recalc_Area<-st_area(p.500)#recalculate areas of clipped polygons
    dom.500<-ifelse(nrow(p.500)==0, 0, as.numeric(p.500 %>% st_drop_geometry() %>% 
                                                    group_by(BurnSever) %>% 
                                                    summarise(area=sum(Recalc_Area)%>%  as_tibble) %>%
                                                    arrange(desc(area)) %>% dplyr::select(BurnSever)%>%
                                                    filter(row_number()==1))-1)
    
    #calculate area amount after dissolving polygons
    #(in case of polygons overlapping before dissolution, causing overestimates of harvest)
    area.500 <- ifelse((nrow(p.500)==0), 0, (p.500 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    
    
    ssyr.sf.summary[[i]]<-data.frame(SS=SS.i,#PKEY=PKEY.i,
                                     YEAR=YEAR.i,
                                     BURNAREA150=area.150,
                                     BURNAREA300=area.300,
                                     BURNAREA500=area.500,
                                     DOMBURNCLASS150=dom.150,
                                     DOMBURNCLASS300=dom.300,
                                     DOMBURNCLASS500=dom.500)
    print(paste0("point ",i," done"))
  }
  summaries<-do.call(bind_rows,ssyr.sf.summary)
  return(summaries)
}
summaries2<-FireAreaAndCategory(pointsVRI,JasperFire,500)
#write.csv(summaries, "2_outputs/Fire_Damage_Jasper.csv")

#get mean and maximum fire year
#create fire-year raster. 
byear<-JasperFire["FIRE_YEAR"]
v <- vect(byear)#instead of f
#extent      : 493946.1, 813585.9, 6056702, 6401906
#x range=319639.8 y range=345204
r <- rast(v, res=30)
zFire <- rasterize(v, r, field="FIRE_YEAR")
plot(zFire)#what to do about zero values?

b150<-st_buffer(pointsVRI, 150)
b300<-st_buffer(pointsVRI, 300)
b500<-st_buffer(pointsVRI, 500)

#extract mean fire year value in each buffer
pointList<-list()
for (i in 1:nrow(pointsVRI)){
  pointData<-pointsVRI[i,]
  #Extracted and summarized within 50 m of point
  pointData$SS<-pointData$location
  try(maxFireYear_150m <- terra::extract(zFire, b150[i,], fun="max", na.rm=TRUE))
  try(pointData$maxFY150m<-as.integer(maxFireYear_150m$FIRE_YEAR))
  try(meanFireYear_150m <- terra::extract(zFire, b150[i,], fun="mean", na.rm=TRUE))
  try(pointData$meanFY150m<-as.integer(meanFireYear_150m$FIRE_YEAR))
  
  try(maxFireYear_300m <- terra::extract(zFire, b300[i,], fun="max", na.rm=TRUE))
  try(pointData$maxFY300m<-as.integer(maxFireYear_300m$FIRE_YEAR))
  try(meanFireYear_300m <- terra::extract(zFire, b300[i,], fun="mean", na.rm=TRUE))
  try(pointData$meanFY300m<-as.integer(meanFireYear_300m$FIRE_YEAR))

  try(maxFireYear_500m <- terra::extract(zFire, b500[i,], fun="max", na.rm=TRUE))
  try(pointData$maxFY500m<-as.integer(maxFireYear_500m$FIRE_YEAR))
  try(meanFireYear_500m <- terra::extract(zFire, b500[i,], fun="mean", na.rm=TRUE))
  try(pointData$meanFY500m<-as.integer(meanFireYear_500m$FIRE_YEAR))
  
  try(pointList[[i]]<-pointData)
  print(paste0("Line ",i," done."))
}
pointFireYear<-do.call(rbind,pointList)
summaries3<-data.frame(pointFireYear)
summaries3$geometry<-NULL
#write.csv(summaries, file="2_outputs/meanandmaxFireYear_Jasper.csv")

#Harvest Area
HarvestArea<-function(pointfile, polyfile, maxdist){
  ssyr.sf.summary<-list()
  #ssyr.sf.demo<-pointfile[1:100,]
  for (i in 1:nrow(pointfile)){
    ssyr.sf.i<-pointfile[i,]#ssyr.sf.demo[i,]
    SS.i<-ssyr.sf.i$location
    YEAR.i<-ssyr.sf.i$year
    #LF.i<-ssyr.sf.i$LF
    #1. use the year of the survey then filter to only the footprint polygons that are older than that point count survey.
    polyfile.y<-polyfile%>%
      filter(Treatment_Year<YEAR.i)#should this filtering occur on the polygons after isolating those within 150 and 565 m of each point count?
    #1b. also after time-filtering, create 150-m and 565-m buffers
    b150<-st_buffer(ssyr.sf.i,150)
    p.150<-st_intersection(polyfile.y, b150)
    p.150$Recalc_Area<-st_area(p.150)#recalculate areas of clipped polygons
    area.150 <- ifelse((nrow(p.150)==0), 0, (p.150 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    
    b300<-st_buffer(ssyr.sf.i,300)
    p.300<-st_intersection(polyfile.y, b300)
    p.300$Recalc_Area<-st_area(p.300)#recalculate areas of clipped polygons
    area.300 <- ifelse((nrow(p.300)==0), 0, (p.300 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    b500<-st_buffer(ssyr.sf.i,500)#a bit under 1 square kilometre
    p.500<-st_intersection(polyfile.y, b500)
    p.500$Recalc_Area<-st_area(p.500)#recalculate areas of clipped polygons
    area.500 <- ifelse((nrow(p.500)==0), 0, (p.500 %>% st_union() %>% st_area()))# select just the polygons in the intersection
    
    
    
    ssyr.sf.summary[[i]]<-data.frame(SS=SS.i,#PKEY=PKEY.i,
                                     YEAR=YEAR.i,
                                     HARVESTAREA150=area.150,
                                     HARVESTAREA300=area.300,
                                     HARVESTAREA500=area.500)
    print(paste0("point ",i," done"))
  }
  summaries<-do.call(bind_rows,ssyr.sf.summary)
  return(summaries)
}
summaries4<-HarvestArea(pointsVRI,JasperHarvest,500)


#Harvest Year
#get mean and maximum harvest year
#create harvest-year raster. 
JasperHarvest<-JasperHarvest[st_is_valid(JasperHarvest)==TRUE,]
hyear<-JasperHarvest["Treatment_Year"]
v <- vect(hyear)#instead of f
#extent      : 493946.1, 813585.9, 6056702, 6401906
#x range=319639.8 y range=345204
r <- rast(v, res=30)
zHarvest <- rasterize(v, r, field="Treatment_Year")
plot(zHarvest)#what to do about zero values?

#extract mean and max harvest year value in each buffer
pointList<-list()
for (i in 1:nrow(pointsVRI)){
  pointData<-pointsVRI[i,]
  #Extracted and summarized within 50 m of point
  try(maxHarvestYear_150m <- terra::extract(zHarvest, b150[i,], fun="max", na.rm=TRUE))
  try(pointData$maxHY150m<-as.integer(maxHarvestYear_150m$Treatment_Year))
  try(meanHarvestYear_150m <- terra::extract(zHarvest, b150[i,], fun="mean", na.rm=TRUE))
  try(pointData$meanHY150m<-as.integer(meanHarvestYear_150m$Treatment_Year))
  
  try(maxHarvestYear_300m <- terra::extract(zHarvest, b300[i,], fun="max", na.rm=TRUE))
  try(pointData$maxHY300m<-as.integer(maxHarvestYear_300m$Treatment_Year))
  try(meanHarvestYear_300m <- terra::extract(zHarvest, b300[i,], fun="mean", na.rm=TRUE))
  try(pointData$meanHY300m<-as.integer(meanHarvestYear_300m$Treatment_Year))
  
  try(maxHarvestYear_500m <- terra::extract(zHarvest, b500[i,], fun="max", na.rm=TRUE))
  try(pointData$maxHY500m<-as.integer(maxHarvestYear_500m$Treatment_Year))
  try(meanHarvestYear_500m <- terra::extract(zHarvest, b500[i,], fun="mean", na.rm=TRUE))
  try(pointData$meanHY500m<-as.integer(meanHarvestYear_500m$Treatment_Year))
  
  try(pointList[[i]]<-pointData)
  print(paste0("Line ",i," done."))
}
pointHarvestYear<-do.call(rbind,pointList)
summaries5<-data.frame(pointHarvestYear)
summaries5$geometry<-NULL

m1<-merge(summary1,summaries2,by=c("SS","YEAR"))
str(m1)
summaries3$SS<-summaries3$location
summaries3$YEAR<-summaries3$year
m2<-merge(m1,summaries3,by=c("SS","YEAR"))
m3<-merge(m2,summaries4,by=c("SS","YEAR"))
summaries5$SS<-summaries5$location
summaries5$YEAR<-summaries5$year
m4<-merge(m3,summaries5,by=c("SS","YEAR"))
write.csv(m4, file="2_outputs/allFireHarvestBeetleSummaries.csv")
