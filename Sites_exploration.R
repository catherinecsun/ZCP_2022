### exploring proposed 2022 CT sites ###

#library
library(ggplot2)
library(ggmap)
library(gridExtra)
library(raster)
library(rgdal)
library(maps)
library(ggsn)
library(FNN)
library(sp)
library(sf)
library(dplyr)
library(lubridate)


#data
sites<-read.csv("../Data/Sites.csv")
data_2016<-read.table("../MilanVinks/captfile.txt")
obs_all<-read.csv("../MilanVinks/tblSightings_7_17_19.csv")
cats_suggestions<-read.csv("../Data/ZCP_2022_suggestedChanges_Cat.csv")

#change from degrees to UTM Zone 35S
sites_latLong <- SpatialPoints(cbind(sites$long, sites$lat), proj4string = CRS("+proj=longlat"))
sites_UTM <- spTransform(sites_latLong, CRS("+init=epsg:20935"))
sites<-cbind(sites,sites_UTM@coords)
sites$Status[sites$Status=="New"]<-'2022'
sites$Status[sites$Status=="Previous"]<-'2016'
colnames(sites)<-c("Site_ID","CT Year","lat","long","UTM_35S_Easting","UTM_35S_Northing")

par(mfrow = c(1, 2))
plot(sites_latLong, axes = TRUE, main = "Lat-Long", cex.axis = 0.95)
plot(sites_UTM, axes = TRUE, main = "UTM 35S", col = "red", cex.axis = 0.95)
dev.off()

ggplot(as.data.frame(sites_UTM@coords),aes(x=coords.x1,y=coords.x2))+
  geom_point()

###Plot ####
# courtesy R Lovelace
ggmap_rast <- function(map){
  map_bbox <- attr(map, 'bb') 
  .extent <- extent(as.numeric(map_bbox[c(2,4,1,3)]))
  my_map <- raster(.extent, nrow= nrow(map), ncol = ncol(map))
  rgb_cols <- setNames(as.data.frame(t(col2rgb(map))), c('red','green','blue'))
  red <- my_map
  values(red) <- rgb_cols[['red']]
  green <- my_map
  values(green) <- rgb_cols[['green']]
  blue <- my_map
  values(blue) <- rgb_cols[['blue']]
  stack(red,green,blue)
}

ZCP_CTregion_Sat<-get_map(location=c(lon=26.15,lat=-14.39080),
                            #c(left = 25.95, bottom = -14.55, 
                            #  right = 26.3, top = -14.3),
                            zoom=11,#"auto",
                            #maptype = "satellite",
                            source = "google")
ZCP_CTregion_SatRast<-ggmap_rast(map = ZCP_CTregion_Sat)
ZCP_CTregion_SatRast<-crop(ZCP_CTregion_SatRast,as(extent(25.9, 26.3, -14.55, -14.3), 'SpatialPolygons'))
ZCP_CTregion_SatRast <- data.frame(rasterToPoints(ZCP_CTregion_SatRast))

ZCP_CTregion_Hybrid<-get_map(location=c(lon=26.15,lat=-14.39080),
                          #c(left = 25.95, bottom = -14.55, 
                          #  right = 26.3, top = -14.3),
                          zoom=11,#"auto",
                          maptype = "hybrid",
                          source = "google")
ZCP_CTregion_HybRast<-ggmap_rast(map = ZCP_CTregion_Hybrid)
ZCP_CTregion_HybRast<-crop(ZCP_CTregion_HybRast,as(extent(25.9, 26.3, -14.55, -14.3), 'SpatialPolygons'))
ZCP_CTregion_HybRast <- data.frame(rasterToPoints(ZCP_CTregion_HybRast))

#this doesnt work very well
# theMap<-ggplot() + 
#   geom_point(alpha=0.7,data=ZCP_CTregion_SatRast,aes(x=x, y=y, col=rgb(layer.1/255, layer.2/255, layer.3/255))) + 
#   geom_point(alpha=0.7,data=ZCP_CTregion_HybRast,aes(x=x, y=y, col=rgb(layer.1/255, layer.2/255, layer.3/255))) + 
#   scale_color_identity()
# 
# theMap+
#   geom_point(data=sites,aes(x = long, y = lat),
#              size = 1,show.legend = TRUE) +
#   ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
#                  st.color="white", st.dist=0.05,
#                  height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
#   xlab("Longitude") +
#   ylab("Latitude") +
#   theme_bw()
                          
map_proposedFor2022<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`,
                            fill=`CT Year`),
             size = 2,pch=24,show.legend = TRUE) +
  scale_color_manual(values=c( "black", "red"))+
  scale_fill_manual(values=c( "black", "red"))+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
            st.color="white", st.dist=0.1,
            height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()+
  theme(legend.position="bottom")
   
map_proposedFor2022

####Distances ####
distances_nearest<-knn.dist(sites[,(ncol(sites)-1):ncol(sites)],k=1)
distances_nearest_old<-knn.dist(sites[sites$`CT Year`=="2016",(ncol(sites)-1):ncol(sites)],k=1)
summary(distances_nearest)
summary(distances_nearest_old)

distances<-data.frame(Period=c(rep(2022,nrow(distances_nearest)),rep("2016",nrow(distances_nearest_old))),
                      Meters=c(distances_nearest,distances_nearest_old))

plot_Dists<-ggplot(data=distances[distances$Period==2022,],aes(x=Meters))+
  geom_histogram(alpha=0.5,position="identity")+
  xlim(1000,2700)+
  geom_vline(aes(xintercept=mean(Meters)),
              linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot_Dists_v2<-ggplot(data=distances,aes(x=Meters,fill=Period))+
  geom_histogram(alpha=0.5,position="identity")+
  scale_fill_manual(values=c("#56B4E9","#E69F00"))+
  geom_vline(aes(xintercept=mean(Meters[Period=='2022'])),
             color="darkblue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=mean(Meters[Period=='Previous'])),
             color="darkorange", linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# plot_multi_histogram <- function(df, feature, label_column) {
#   plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
#     geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
#     geom_density(alpha=0.7) +
#     geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
#     labs(x=feature, y = "Density")
#   plt + guides(fill=guide_legend(title=label_column))
# }
# plot_multi_histogram(distances, 'Meters', 'Period')


### previous detections ####
data_2016
colnames(data_2016)<-c("Year","Ind","Occ","Site_ID")
table(data_2016$Site_ID)
sites<-left_join(sites,as.data.frame(table(data_2016$Site_ID)),
          by=c("Site_ID"="Var1"))
colnames(sites)[ncol(sites)]<-"CT Dets in 2016"
sites$`CT Dets in 2016`[is.na(sites$`CT Dets in 2016`)]<-0


ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "#56B4E9", "black"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +

  scale_size(range = c(3, 6))+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()+
  theme(legend.position="bottom")

### previous observations ####
obs_all$SightingDate<-as.Date(obs_all$SightingDate,"%m/%d/%Y")
obs_all_drySeason<-obs_all[month(obs_all$SightingDate)<12&month(obs_all$SightingDate)>5,]
obs_all_drySeason<-obs_all_drySeason[obs_all_drySeason$Lat!=0,]
obs_all_drySeason_Leopard<-obs_all_drySeason[obs_all_drySeason$LeopardSighting==TRUE,]
obs_all_drySeason_Leopard<-obs_all_drySeason_Leopard[!is.na(obs_all_drySeason_Leopard$ID),]
obs_all_drySeason_Hyena<-obs_all_drySeason[obs_all_drySeason$HyenaSighting==TRUE,]
obs_all_drySeason_Hyena<-obs_all_drySeason_Hyena[!is.na(obs_all_drySeason_Hyena$ID),]

obs_all_drySeason_Leopard%>%group_by(year(obs_all_drySeason_Leopard$SightingDate))%>%summarize(n())
obs_all_drySeason_Hyena%>%group_by(year(obs_all_drySeason_Hyena$SightingDate))%>%summarize(n())

### Leopard Plots ####

map_sitesDetsObs_Leop_2012<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Leopard[year(obs_all_drySeason_Leopard$SightingDate)==2012,],
             aes(x=Lon,y=Lat,fill=LeopardSighting),pch=21)+
   scale_fill_manual(values="yellow")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Leopards: Jun-Nov 2012",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))

map_sitesDetsObs_Leop_2013<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Leopard[year(obs_all_drySeason_Leopard$SightingDate)==2013,],
             aes(x=Lon,y=Lat,fill=LeopardSighting),pch=21)+
  scale_fill_manual(values="yellow")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Leopards: Jun-Nov 2013",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))

map_sitesDetsObs_Leop_2014<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Leopard[year(obs_all_drySeason_Leopard$SightingDate)==2014,],
             aes(x=Lon,y=Lat,fill=LeopardSighting),pch=21)+
  scale_fill_manual(values="yellow")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Leopards: Jun-Nov 2014",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))

map_sitesDetsObs_Leop_2015<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Leopard[year(obs_all_drySeason_Leopard$SightingDate)==2015,],
             aes(x=Lon,y=Lat,fill=LeopardSighting),pch=21)+
  scale_fill_manual(values="yellow")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Leopards: Jun-Nov 2015",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))


map_sitesDetsObs_Leop_2016<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Leopard[year(obs_all_drySeason_Leopard$SightingDate)==2016,],
             aes(x=Lon,y=Lat,fill=LeopardSighting),pch=21)+
  scale_fill_manual(values="yellow")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Leopards: Jun-Nov 2016",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))

map_sitesDetsObs_Leop_2017<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Leopard[year(obs_all_drySeason_Leopard$SightingDate)==2017,],
             aes(x=Lon,y=Lat,fill=LeopardSighting),pch=21)+
  scale_fill_manual(values="yellow")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Leopards: Jun-Nov 2017",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))


map_sitesDetsObs_Leop_2018<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Leopard[year(obs_all_drySeason_Leopard$SightingDate)==2018,],
             aes(x=Lon,y=Lat,fill=LeopardSighting),pch=21)+
  scale_fill_manual(values="yellow")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Leopards: Jun-Nov 2018",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))


map_proposedFor2022
map_sitesDetsObs_Leop_2012
map_sitesDetsObs_Leop_2013
map_sitesDetsObs_Leop_2014
map_sitesDetsObs_Leop_2015
map_sitesDetsObs_Leop_2016
map_sitesDetsObs_Leop_2017
map_sitesDetsObs_Leop_2018

#### Hyena Plots ####
map_sitesDetsObs_Hyena_2012<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Hyena[year(obs_all_drySeason_Hyena$SightingDate)==2012,],
             aes(x=Lon,y=Lat,fill=HyenaSighting),pch=21)+
  scale_fill_manual(values="white")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Hyena: Jun-Nov 2012",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))

map_sitesDetsObs_Hyena_2013<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Hyena[year(obs_all_drySeason_Hyena$SightingDate)==2013,],
             aes(x=Lon,y=Lat,fill=HyenaSighting),pch=21)+
  scale_fill_manual(values="white")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Hyena: Jun-Nov 2013",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))

map_sitesDetsObs_Hyena_2014<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Hyena[year(obs_all_drySeason_Hyena$SightingDate)==2014,],
             aes(x=Lon,y=Lat,fill=HyenaSighting),pch=21)+
  scale_fill_manual(values="white")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Hyena: Jun-Nov 2014",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))

map_sitesDetsObs_Hyena_2015<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Hyena[year(obs_all_drySeason_Hyena$SightingDate)==2015,],
             aes(x=Lon,y=Lat,fill=HyenaSighting),pch=21)+
  scale_fill_manual(values="white")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Hyena: Jun-Nov 2015",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))


map_sitesDetsObs_Hyena_2016<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Hyena[year(obs_all_drySeason_Hyena$SightingDate)==2016,],
             aes(x=Lon,y=Lat,fill=HyenaSighting),pch=21)+
  scale_fill_manual(values="white")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Hyena: Jun-Nov 2016",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))

map_sitesDetsObs_Hyena_2017<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Hyena[year(obs_all_drySeason_Hyena$SightingDate)==2017,],
             aes(x=Lon,y=Lat,fill=HyenaSighting),pch=21)+
  scale_fill_manual(values="white")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Hyena: Jun-Nov 2017",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))


map_sitesDetsObs_Hyena_2018<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),pch=21,fill=NA,size=2) +
  scale_color_manual(values=c( "black", "red"))+
  geom_point(data=sites[sites$`CT Dets in 2016`>0,],
             aes(x = long, y = lat,size=`CT Dets in 2016`),pch=21,
             fill="black",color="black") +
  scale_size(range = c(3, 6))+
  geom_point(data=obs_all_drySeason_Hyena[year(obs_all_drySeason_Hyena$SightingDate)==2018,],
             aes(x=Lon,y=Lat,fill=HyenaSighting),pch=21)+
  scale_fill_manual(values="white")+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(title="Hyena: Jun-Nov 2018",x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))


map_proposedFor2022
map_sitesDetsObs_Hyena_2012
map_sitesDetsObs_Hyena_2013
map_sitesDetsObs_Hyena_2014
map_sitesDetsObs_Hyena_2015
map_sitesDetsObs_Hyena_2016
map_sitesDetsObs_Hyena_2017
map_sitesDetsObs_Hyena_2018

grid.arrange(map_sitesDetsObs_Leop_2012, map_sitesDetsObs_Hyena_2012, ncol=2)
grid.arrange(map_sitesDetsObs_Leop_2013, map_sitesDetsObs_Hyena_2013, ncol=2)
grid.arrange(map_sitesDetsObs_Leop_2014, map_sitesDetsObs_Hyena_2014, ncol=2)
grid.arrange(map_sitesDetsObs_Leop_2015, map_sitesDetsObs_Hyena_2015, ncol=2)
grid.arrange(map_sitesDetsObs_Leop_2016, map_sitesDetsObs_Hyena_2016, ncol=2)
grid.arrange(map_sitesDetsObs_Leop_2017, map_sitesDetsObs_Hyena_2017, ncol=2)
grid.arrange(map_sitesDetsObs_Leop_2018, map_sitesDetsObs_Hyena_2018, ncol=2)

#concerns
# the trap spacing may still be too large
# how were sites chosen? 
#    trails, trees, known detections, or, etc?
#    why the focus on the northeast?
#    
#

####new proposition
cats_suggestions_UTM <- SpatialPoints(cbind(cats_suggestions$long, cats_suggestions$lat), proj4string = CRS("+proj=longlat"))
cats_suggestions_UTM <- spTransform(cats_suggestions_UTM, CRS("+init=epsg:20935"))

cats_suggestions$UTM_35S_Easting<-cats_suggestions_UTM@coords[,1]
cats_suggestions$UTM_35S_Northing<-cats_suggestions_UTM@coords[,2]
cats_suggestions$nCams<-2
cats_suggestions$nCams[c(12,16,17)]<-1

sites_closer<-sites[c(3,4,11,14,15,16),] #2022 site locations im suggesting we keep
sites_closer$nCams<-2
sites_closer$nCams[c(3,5)]<-1 # but reduce to 1 cam for 2 sites
sites_closer<-plyr::rbind.fill(sites_closer,sites[sites$`CT Year`=="2016",])
sites_closer$nCams[is.na(sites_closer$nCams)]<-2
sites_closer<-plyr::rbind.fill(sites_closer,cats_suggestions)
sites_closer$`CT Year`[is.na(sites_closer$`CT Year`)]<-2022

plot_proposedLocs<-ggmap(ZCP_CTregion_Hybrid)+
  geom_point(data=sites_closer,
             aes(x = long, y = lat,color=`CT Year`, fill=`CT Year`,
                 pch=as.factor(nCams),show.legend=FALSE),
             size=2) +
  #scale_color_continuous(guide = "none") +
  scale_color_manual(values=c( "black", "red"))+
  ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
                 st.color="white", st.dist=0.1,
                 height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  labs(x="Longitude",y="Latitude") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  guides(size=guide_legend(nrow=2, byrow=TRUE),
         color=guide_legend(nrow=2, byrow=TRUE))

#distances
distances_nearest_proposed<-knn.dist(sites_closer[,5:6],k=1)
distances_nearest_proposed<-data.frame(Period=rep("Newly Proposed",nrow(distances_nearest_proposed)),
                                       Meters=distances_nearest_proposed)
dist_comparison<-rbind(distances_nearest_proposed,distances[distances$Period=="2022",])
dist_comparison$Period[dist_comparison$Period=="2022"]<-"Originally Proposed"

plot_proposedDists<-ggplot(data=dist_comparison[dist_comparison$Period=="Newly Proposed",],aes(x=Meters))+ #,fill=Period
  geom_histogram(alpha=0.5,fill="red",position="identity")+
  #scale_fill_manual(values=c("red"))+ #,"gray"
  xlim(1000,2700)+
  geom_vline(aes(xintercept=mean(Meters[Period=='Newly Proposed'])),
             color="darkred", linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


### points to polygons for convex hulls
chull_wBuff_proposedPrev2022<-st_buffer(dist=5000,st_as_sf(sites, coords=c("UTM_35S_Easting","UTM_35S_Northing")))
chull_wBuff_proposedPrev2022<-st_union(chull_wBuff_proposedPrev2022)
st_crs(chull_wBuff_proposedPrev2022) <- "+init=epsg:20935" #set the CRS
chull_wBuff_proposedPrev2022<-st_transform(chull_wBuff_proposedPrev2022, "+proj=longlat") # change to lat logn
st_area(chull_wBuff_proposedPrev2022) #622.5 km2

map_proposedFor2022+ # cat's proposed
  geom_sf(data = chull_wBuff_proposedPrev2022, colour = "black", fill = NA,inherit.aes = FALSE)+
  geom_sf(data = chull_wBuff_proposed, color =alpha("red",0.3), fill = NA,inherit.aes = FALSE)

  

chull_wBuff_proposed<-st_buffer(dist=5000,st_as_sf(sites_closer, coords=c("UTM_35S_Easting","UTM_35S_Northing")))
chull_wBuff_proposed<-st_union(chull_wBuff_proposed)
st_crs(chull_wBuff_proposed) <- "+init=epsg:20935" #set the CRS
chull_wBuff_proposed<-st_transform(chull_wBuff_proposed, "+proj=longlat") # change to lat logn
st_area(chull_wBuff_proposed)#565.9 km2

plot_proposedLocs+
  geom_sf(data = chull_wBuff_proposed, colour = "red", fill = NA,inherit.aes = FALSE)

####comparison of newly proposed by Cat vs originally proposed
plot_Dists
plot_proposedDists

plot_proposedLocs
map_proposedFor2022

grid.arrange(map_proposedFor2022+
               geom_sf(data = chull_wBuff_proposedPrev2022, colour = "black", fill = NA,inherit.aes = FALSE)+
              theme(legend.position="none"),
             plot_proposedLocs+
               geom_sf(data = chull_wBuff_proposedPrev2022, colour = "black", fill = NA,inherit.aes = FALSE)+
               geom_sf(data = chull_wBuff_proposed, color =alpha("red",0.3), fill = NA,inherit.aes = FALSE)+
               theme(legend.position="none"),
             plot_Dists,plot_proposedDists+ theme(legend.position="none"),
             ncol = 2, nrow = 2)
