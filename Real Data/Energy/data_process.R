###
library(readxl)
library(tidyverse)
orginal_data<-read_excel("existcapacity_annual.xlsx",skip = 1)
orginal_data = orginal_data[,!(names(orginal_data) %in% 'Producer Type')]
source_type<-unique(orginal_data$`Fuel Source`)
type1<-c("Coal","Petroleum","Wood and Wood Derived Fuels")
type2<-c("Natural Gas","Other Gases"  )
type3<-c("Hydroelectric","Wind","Nuclear","Geothermal" ,"Solar Thermal and Photovoltaic")
range(orginal_data$Year)
x_grid<-1990:2021
Data_year<-list()
for (i in 1:length(x_grid)) {
  year<-x_grid[i]
  Data_i<-orginal_data[orginal_data$Year==year,]
  states<-unique(Data_i$`State Code`)
  Data_state<-list()
  for (j in 1:length(states)) {
    state<-states[j]
    Data_i_state<-Data_i[Data_i$`State Code`==state,]
    Data_i_state<-Data_i_state[!duplicated(Data_i_state), ]
    G_typ1<-Data_i_state %>% filter(`Fuel Source` %in% type1)
    G_typ2<-Data_i_state %>% filter(`Fuel Source` %in% type2)
    G_typ3<-Data_i_state %>% filter(`Fuel Source` %in% type3)
    prop<-c(sum(G_typ1$`Nameplate Capacity (Megawatts)`),sum(G_typ2$`Nameplate Capacity (Megawatts)`),
            sum(G_typ3$`Nameplate Capacity (Megawatts)`))
    Data_state[[state]]<-sqrt(prop/sum(prop))
  }
  Data_year[[i]]<-Data_state
}
save(Data_year,file = "pre-processed.RData")
