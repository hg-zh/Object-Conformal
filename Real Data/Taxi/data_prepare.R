library(arrow)
library(lubridate)
library(dplyr)
library(readxl)
library(timeDate)
Get_ad_mat<-function(dataset){
  out<-matrix(0,13,13)
  for (i in 1:length(dataset$putime_d)) {
    out[mapping_vec[which(taxi_zone_mh_id ==as.numeric(dataset[i,1]))],mapping_vec[which(taxi_zone_mh_id ==as.numeric(dataset[i,2]))]]<-
      out[mapping_vec[which(taxi_zone_mh_id ==as.numeric(dataset[i,1]))],mapping_vec[which(taxi_zone_mh_id ==as.numeric(dataset[i,2]))]]+1
  }
  out<-out/max(out)
  return(out)
}
taxi_zone_mh_id <-read_excel("taxi_zone_mh_id.xlsx")$ID
taxi_zone_mh_id <- taxi_zone_mh_id[!taxi_zone_mh_id %in% c(103,104,105)]
num_in_bin<-1000
#### ####
mapping_vec<-rep(0,length(taxi_zone_mh_id))
for (i in 1:length(taxi_zone_mh_id)) {
  if(taxi_zone_mh_id[i] %in% c(12,13,87,88,209,231,261)){mapping_vec[i]<-1}
  if(taxi_zone_mh_id[i]  %in% c(113,114,125,144,158,211,249)){mapping_vec[i]<-2}
  if(taxi_zone_mh_id[i]  %in% c(4,45,79,148,232)){mapping_vec[i]<-3}
  if(taxi_zone_mh_id[i]  %in% c(48,50,68,90,246)){mapping_vec[i]<-4}
  if(taxi_zone_mh_id[i]  %in% c(100,161,163,164,186,230,234)){mapping_vec[i]<-5}
  if(taxi_zone_mh_id[i]  %in% c(107,137,162,170,224,229,233)){mapping_vec[i]<-6}
  if(taxi_zone_mh_id[i]  %in% c(24,142,143,151,238,239)){mapping_vec[i]<-7}
  if(taxi_zone_mh_id[i]  %in% c(140,141,202,236,237,262,263)){mapping_vec[i]<-8}
  if(taxi_zone_mh_id[i]  %in% c(116,152,166)){mapping_vec[i]<-9}
  if(taxi_zone_mh_id[i]  %in% c(41,42)){mapping_vec[i]<-10}
  if(taxi_zone_mh_id[i]  %in% c(74,75,194)){mapping_vec[i]<-11}
  if(taxi_zone_mh_id[i]  %in% c(120,127,128,153,243,244)){mapping_vec[i]<-12}
  if(taxi_zone_mh_id[i]  %in% c(43)){mapping_vec[i]<-13}
}
Data_wd<-list();Data_ho<-list()
#### ####
Mon_of_Int<-1:12
for(m in 1:length(Mon_of_Int)){
  mon<-Mon_of_Int[m]
  if(mon<10){
    yetrips <- read_parquet(paste0('yellow_2023/yellow_tripdata_2023-0',mon,'.parquet'))
  }else{
    yetrips <- read_parquet(paste0('yellow_2023/yellow_tripdata_2023-',mon,'.parquet'))
  }
  yetrips <- yetrips %>%
    mutate(tpep_pickup_datetime = with_tz(tpep_pickup_datetime, tzone = "America/Los_Angeles") %>%
             with_tz(tzone = "America/New_York"))
  yetrips <- yetrips %>%
    mutate(tpep_dropoff_datetime = with_tz(tpep_dropoff_datetime, tzone = "America/Los_Angeles") %>%
             with_tz(tzone = "America/New_York"))
  yetrips <- yetrips %>%
    mutate(pudate = as.Date(tpep_pickup_datetime,tz='EST'),
           putime = format(tpep_pickup_datetime, "%H:%M:%S"))
  yetrips <- yetrips %>%
    mutate(dodate = as.Date(tpep_dropoff_datetime,tz='EST'),
           dotime = format(tpep_dropoff_datetime, "%H:%M:%S"))
  yetrips <- yetrips %>% mutate(putime_d = hour(tpep_pickup_datetime) + minute(tpep_pickup_datetime) / 60
                                + second(tpep_pickup_datetime) / 3600)
  yetrips <- yetrips %>% mutate(dotime_d = hour(tpep_dropoff_datetime) + minute(tpep_dropoff_datetime) / 60
                                + second(tpep_dropoff_datetime) / 3600)
  days_of_mon<-unique(yetrips$pudate )
  days_of_mon <- days_of_mon[year(days_of_mon) == 2023]
  days_of_mon <- days_of_mon[month(days_of_mon) == mon]
  for (d in 1:length(days_of_mon)) {
    day_of_Int<-days_of_mon[d]
    D_day<-select(yetrips %>%  filter(pudate == day_of_Int ), PULocationID, DOLocationID,putime_d)
    D_day<-D_day %>% filter(PULocationID %in% taxi_zone_mh_id)
    D_day<-D_day %>% filter(DOLocationID %in% taxi_zone_mh_id)
    num_of_bin<-floor(length(D_day$putime_d)/num_in_bin)
    attach(D_day) # -> no need to add "dataVol$" when calling variables in dataVol
    # define bins s.t. each bin has same/similar numbers of measurements
    timebreaks <- quantile(putime_d , seq(0,1,length.out = num_of_bin + 1))
    timebreaks[1] <- min(putime_d)
    timebreaks[length(timebreaks)] <- max(putime_d)
    timeBinMids <- timebreaks[-1] - diff(timebreaks)/2 # mid points of each bin
    names(timebreaks) <- names(timeBinMids) <- NULL
    timeBinIdx<- findInterval(x =putime_d, vec = timebreaks, rightmost.closed = TRUE) # indices of bins where each observation lies
    detach(D_day)
    VperBin_D1<-list()
    for (i in 1:length(timebreaks)) {
      VperBin_D1[[i]]<-D_day[which(timeBinIdx==i) ,]
      if( isHoliday(timeDate(day_of_Int), holidayNYSE()) ){
        Data_ho[[length(Data_ho)+1]]<-list(x=timebreaks[i],y=Get_ad_mat(VperBin_D1[[i]]),date=day_of_Int)
      }else{
        Data_wd[[length(Data_wd)+1]]<-list(x=timebreaks[i],y=Get_ad_mat(VperBin_D1[[i]]),date=day_of_Int)
      }
    }
  }
  print(m)
}

save(Data_ho,file = "Holiday_2023.RData")
save(Data_wd,file = "Weekday_2023.RData")



