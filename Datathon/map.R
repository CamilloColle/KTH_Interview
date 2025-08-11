library(mapview)
library(sp)
library(sf)
library(dplyr)
library(arrow)
library(ggplot2)
library(tmap)

aed_locations<-read_parquet('Data/aed_locations.parquet.gzip')
belgium.sf<-st_read("4326/postaldistricts.shp")

aedCountByPostalCode<-aed_locations %>% count(postal_code)
aedCountByPostalCode$postal_code<-as.character(aedCountByPostalCode$postal_code)
names(belgium.sf)[names(belgium.sf)=='nouveau_PO']='postal_code'


belgium.sf<-left_join(belgium.sf,aedCountByPostalCode,by='postal_code')

belgium.sf<-left_join(belgium.sf,new_aed_1,by='postal_code')
df$postal_code<-as.character(df$postal_code)
belgium.sf<-left_join(belgium.sf,df,by='postal_code')

Breaks <- c(0,1,5,10,)
#Labels <- c("0 - 2000", "2000 - 4000", "4000 - 6000", ">6000")
MyPalette <- c("#9C3848","#FFAD69", "#F5E663", "#5DA271", "#204E4A")

tm_shape(belgium.sf) + tm_polygons(style='cont',
                                   palette= 'Greens',
                                   col="add_vector_2", 
                                   border.col="black"
                                   )
tm_shape(belgium.sf) + tm_polygons(style='cont',
                                   palette= 'Reds',
                                   col="DOA_cardiac_count", 
                                   border.col="black"
)
trans = "reverse"
scale_fill_continuous(guide = guide_colourbar())
library(ggthemes)
library(colorspace)
p1<-ggplot(belgium.sf) + 
  geom_sf(aes(fill = DOA_cardiac_count),color="white")+
  scale_fill_gradient2(mid='white',high='#7c0202ff', na.value = '#292929',midpoint=-70,
                        guide = guide_colorbar(title = "Dead on Arrival Ambulance Calls (Cardiac)",
                                               title.position = "top",
                                               title.theme = element_text(size = 10,face = "bold", colour = "white",angle = 0),
                                               label.theme= element_text(size = 10,face = "bold",colour = "white",angle = 0))
                        )+
  theme_set(theme_map()+theme(legend.background = element_rect(fill = 'transparent')))
p1
p1<-ggplot(belgium.sf) + 
  geom_sf(aes(fill = DOA_cardiac_count.x),color="grey30")+
  scale_fill_gradient2(mid ='#eebfbfff',high='#7c0202ff', na.value = 'white',midpoint=0,
                       guide = guide_colorbar(title = "Dead on Arrival Ambulance Calls (Cardiac)",
                                              title.position = "top",
                                              title.theme = element_text(size = 10,face = "bold", colour = "white",angle = 0),
                                              label.theme= element_text(size = 10,face = "bold",colour = "white",angle = 0))
  )+
  theme_set(theme_map()+theme(legend.background = element_rect(fill = 'transparent')))
p1
ggsave('cardiac_DOA_map.png', p1, bg='transparent',width = 7, height = 7,dpi = 500)
df$DOA_cardiac_count_norm<-df$DOA_cardiac_count/df$population
library(scales)
p2<-ggplot(belgium.sf) + 
  geom_sf(aes(fill = aed_count.x),color="grey30")+
  scale_fill_gradient2(mid ='#eebfbfff',high='#7c0202ff', na.value = 'white', limits=c(0, 150), oob=squish,midpoint=0,labels = c(0, 50, 100,">150"),
                       guide = guide_colorbar(title = "AED Counts by Postal Code (original)",
                                              title.position = "top",
                                              title.theme = element_text(size = 10,face = "bold", colour = "white",angle = 0),
                                              label.theme= element_text(size = 10,face = "bold",colour = "white",angle = 0))
  )+
  theme_set(theme_map()+theme(legend.background = element_rect(fill = 'transparent')))
p2
ggsave('cardiac_DOA_map_w2.png', p2, bg='transparent',width = 7, height = 7,dpi = 500)