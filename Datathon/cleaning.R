library(ggplot2)
library(mapview)
library(arrow)
aed_locations<-read_parquet('Data/aed_locations.parquet.gzip')
interventions_bxl<-read_parquet('Data/interventions_bxl.parquet.gzip')
interventions_bxl2<-read_parquet('Data/interventions_bxl2.parquet.gzip')
interventions1<-read_parquet('Data/interventions1.parquet')
interventions2<-read_parquet('Data/interventions2.parquet')
interventions3<-read_parquet('Data/interventions3.parquet')


#Cleaning up Intervention bxl
## Removing any ambulance ride where they abandoned the trip
## Except for ones where they're dead on arrival (probably useful for AED research)
interventions_bxl<-as.data.frame(subset(interventions_bxl,(is.na(interventions_bxl$abandon_reason)|interventions_bxl$abandon_reason=='Overleden')))
interventions_bxl$abandon_reason<-as.factor(interventions_bxl$abandon_reason)
interventions_bxl$eventtype_trip<-as.factor(interventions_bxl$eventtype_trip)
## Removing rows where the problem doesn't seem relevent (imo)
aed_reasons<-c("P003 - Cardiac arrest",
               "P001 - Traffic accident",
               "P010 - Respiratory problems",
               "P011 - Chest pain",
               "P014 - Electrocution - electrification",
               "P019 - Unconscious - syncope",
               "P024 - CO intoxication",
               "P025 - Non-traumatic headache",
               "P026 - Unclear problem",
               "P039 - Cardiac problem (other than thoracic pain)",
               "P038 - Person does not answer the call")
interventions_bxl<-interventions_bxl[interventions_bxl$eventtype_trip %in% aed_reasons,]
## Cleaning the times up and getting the time it takes to intervention
interventions_bxl$t0<-as.POSIXct(strptime(substr(interventions_bxl$t0,1,nchar(interventions_bxl$t0)-15),format ="%Y-%m-%d %H:%M:%S"),tz="EET")
interventions_bxl$t3<-as.POSIXct(strptime(substr(interventions_bxl$t3,1,nchar(interventions_bxl$t3)-15),format ="%Y-%m-%d %H:%M:%S"),tz="UTC")
timeToIntervention_bxl<-interventions_bxl$t3-interventions_bxl$t0
length(na.omit(timeToIntervention_bxl))/nrow(interventions_bxl)

#Cleaning up intervention bxl2
## Removing any ambulance ride where they abandoned the trip
## Except for ones where they're dead on arrival (probably useful for AED research)
interventions_bxl2<-as.data.frame(subset(interventions_bxl2,(is.na(interventions_bxl2$`Abandon reason NL`)|interventions_bxl2$`Abandon reason NL` =='Dood Ter Plaatse')))
interventions_bxl2$`Abandon reason NL`<-as.factor(interventions_bxl2$`Abandon reason NL`)
## Removing rows where the problem doesn't seem relevent (imo)
aed_reasons_2<-c("P004 N01 - CVA-TIA",
                 "P004 N02 - CVA-TIA",
                 "P004 N03 - CVA-TIA",
                 "P004 N05 - CVA-TIA",
                 "P003  N05 - HARTSTILSTAND - DOOD - OVERLEDEN",
                 "P003  N05 - HARTSTILSTAND - DOOD - OVERLEDEN",
                 "P011 N01 - PIJN OP DE BORST",
                 "P011 N02 - PIJN OP DE BORST",
                 "P011 N03 - PIJN OP DE BORST",
                 "P011 N04 - PIJN OP DE BORST",
                 "P011 N05 - PIJN OP DE BORST",
                 "P019 N01 - Bewusteloos - coma - syncope",
                 "P019 N03 - Bewusteloos - coma - syncope",
                 "P019 N04 - Bewusteloos - coma - syncope",
                 "P019 N05 - Bewusteloos - coma - syncope",
                 "P019 N06 - Bewusteloos - coma - syncope",
                 "P024 N05 - CO-intoxicatie",
                 "P026 N01 - ONWEL ZONDER DUIDELIJKE REDEN",
                 "P026 N02 - ONWEL ZONDER DUIDELIJKE REDEN",
                 "P026 N04 - ONWEL ZONDER DUIDELIJKE REDEN",
                 "P026 N05 - ONWEL ZONDER DUIDELIJKE REDEN",
                 "P026 N06 - ONWEL ZONDER DUIDELIJKE REDEN",
                 "P038 N01 - PATIËNT BEANTWOORDT DE OPROEP NIET",
                 "P038 N03 - PATIËNT BEANTWOORDT DE OPROEP NIET",
                 "P038 N05 - PATIËNT BEANTWOORDT DE OPROEP NIET",
                 "P039 N01 - CARDIAAL PROBLEEM (ANDERE DAN PIJN AAN DE BORST)",
                 "P039 N03 - CARDIAAL PROBLEEM (ANDERE DAN PIJN AAN DE BORST)",
                 "P039 N04 - CARDIAAL PROBLEEM (ANDERE DAN PIJN AAN DE BORST)",
                 "P039 N05 - CARDIAAL PROBLEEM (ANDERE DAN PIJN AAN DE BORST)"
                 )
interventions_bxl2<-interventions_bxl2[interventions_bxl2$`EventType and EventLevel` %in% aed_reasons_2,]
## Cleaning the times up and getting the time it takes to intervention
interventions_bxl2$T0<-as.POSIXct(strptime(interventions_bxl2$T0,format ="%d%b%y:%H:%M:%S"),tz="UTC")
interventions_bxl2$T3<-as.POSIXct(strptime(interventions_bxl2$T3,format ="%d%b%y:%H:%M:%S"),tz="UTC")
timeToIntervention_bxl2<-interventions_bxl2$T3-interventions_bxl2$T0
length(na.omit(timeToIntervention_bxl2))/nrow(interventions_bxl2)

#Cleaning up Intervention 1
## Removing any ambulance ride where they abandoned the trip
## Except for ones where they're dead on arrival (probably useful for AED research)
interventions1<-as.data.frame(subset(interventions1,(is.na(interventions1$`Abandon reason`)|interventions1$`Abandon reason`=='Overleden')))
interventions1$`Abandon reason`<-as.factor(interventions1$`Abandon reason`)
interventions1$`EventType Trip`<-as.factor(interventions1$`EventType Trip`)
## Removing rows where the problem doesn't seem relevent (imo)
interventions1<-interventions1[interventions1$`EventType Trip` %in% aed_reasons,]
## Cleaning the times up and getting the time it takes to intervention
interventions1$T0<-as.POSIXct(strptime(interventions1$T0,format ="%d%b%y:%H:%M:%S"),tz="UTC")
interventions1$T3<-as.POSIXct(strptime(substr(interventions1$T3,1,nchar(interventions1$T3)-4),format ="%Y-%m-%d %H:%M:%S"),tz="UTC")
timeToIntervention1<-interventions1$T3-interventions1$T0
length(na.omit(timeToIntervention1))/nrow(interventions1)

#Cleaning up Intervention 2
## Removing any ambulance ride where they abandoned the trip
## Except for ones where they're dead on arrival (probably useful for AED research)
interventions2<-as.data.frame(subset(interventions2,(is.na(interventions2$`Abandon reason`)|interventions2$`Abandon reason`=='Overleden')))
interventions2$`Abandon reason`<-as.factor(interventions2$`Abandon reason`)
interventions2$`EventType Trip`<-as.factor(interventions2$`EventType Trip`)
## Removing rows where the problem doesn't seem relevent (imo)
interventions2<-interventions2[interventions2$`EventType Trip` %in% aed_reasons,]
## Cleaning the times up and getting the time it takes to intervention
interventions2$T0<-as.POSIXct(strptime(interventions2$T0,format ="%d%b%y:%H:%M:%S"),tz="UTC")
interventions2$T3<-as.POSIXct(strptime(substr(interventions2$T3,1,nchar(interventions2$T3)-4),format ="%Y-%m-%d %H:%M:%S"),tz="UTC")
timeToIntervention2<-interventions2$T3-interventions2$T0
length(na.omit(timeToIntervention2))/nrow(interventions2)

#Cleaning up Intervention 3
## Removing any ambulance ride where they abandoned the trip
## Except for ones where they're dead on arrival (probably useful for AED research)
interventions3<-as.data.frame(subset(interventions3,(is.na(interventions3$`Abandon reason`)|interventions3$`Abandon reason`=='Overleden')))
interventions3$`Abandon reason`<-as.factor(interventions3$`Abandon reason`)
interventions3$`EventType Trip`<-as.factor(interventions3$`EventType Trip`)
## Removing rows where the problem doesn't seem relevent (imo)
interventions3<-interventions3[interventions3$`EventType Trip` %in% aed_reasons,]
## Cleaning the times up and getting the time it takes to intervention
interventions3$T0<-as.POSIXct(strptime(interventions3$T0,format ="%d%b%y:%H:%M:%S"),tz="UTC")
interventions3$T3<-as.POSIXct(strptime(substr(interventions3$T3,1,nchar(interventions3$T3)-4),format ="%Y-%m-%d %H:%M:%S"),tz="UTC")
timeToIntervention3<-interventions3$T3-interventions3$T0
length(na.omit(timeToIntervention3))/nrow(interventions3)

