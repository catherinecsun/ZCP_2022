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
  geom_point(data=sites,aes(x = long, y = lat,color=`CT Year`),
             size = 1,show.legend = TRUE) +
   ggsn::scalebar(sites[,3:4], location = "bottomleft", dist = 5,
            st.color="white", st.dist=0.1,
            height = 0.05, transform = TRUE, model = "WGS84", dist_unit = "km")+
  xlim(c(25.95,26.27))+
  ylim(c(-14.55,-14.3))+
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()+
  theme(legend.position="bottom")
   
 

#### distances ####
distances_nearest<-knn.dist(sites[,(ncol(sites)-1):ncol(sites)],k=1)
distances_nearest_old<-knn.dist(sites[sites$`CT Year`=="2016",(ncol(sites)-1):ncol(sites)],k=1)
summary(distances_nearest)
summary(distances_nearest_old)

distances<-data.frame(Period=c(rep(2022,nrow(distances_nearest)),rep("2016",nrow(distances_nearest_old))),
                      Meters=c(distances_nearest,distances_nearest_old))

ggplot(data=distances,aes(x=Meters,fill=Period))+
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
