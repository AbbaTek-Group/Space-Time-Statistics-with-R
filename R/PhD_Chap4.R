library("readxl")
library("tidyverse")
library("ggplot2")
library("gstat")
library("maps")
library("STRbook")
library("grid")
library("gridExtra")
library("lubridate")
## ------------------------------------------------------------------------
setwd('/home/ivo/Downloads/PhD-Data-code/')

## load the soil moisture data 

Siloam_orig <- readRDS('SPdata')

## load the atmospheric data and wrangle date and time column

Atmos_orig <- read_excel('2012-2017Atmosph.xlsx')

##----------------------------------------------------------------------

## clean up data and remove missing rows

Siloam_clean <- Siloam_orig[complete.cases(Siloam_orig),]

Atmos_clean <- Atmos_orig[complete.cases(Atmos_orig),]
##------------------------------------------------------------------------

## add columns for year, month, and day to dataset

Siloam <- Siloam_clean %>% mutate (Day = day(DateTime), Month = month(DateTime), Year = year(DateTime))  
Atmos <- Atmos_clean %>% mutate (Day = day(DateTime), Month = month(DateTime), Year = year(DateTime))

Atmos <- Atmos_clean %>% mutate (Day = day(DateTime), Month = month(DateTime), Year = year(DateTime)) %>% gather (Proc, z, -Day,-Month, -Year, -DateTime)

##-------------------------------------------------------------------------------

## summarise mean rainfall and ET 

Atmos %>% group_by(Year, Proc) %>% summarise(mean_proc = mean(z))

## summarise mean soil moisture

Siloam_summary_time <- Siloam %>% filter(!(`soil_moisture(v/v)`<= 0)) %>% group_by(Probe, Year, Month, Day) %>% summarise(mean_soil_moisture = mean(`soil_moisture(v/v)`)) 

Siloam_summary_space <- Siloam %>% filter(!(`soil_moisture(v/v)`<= 0)) %>% group_by(Probe, Year, Month, Day, `Depth(cm)`) %>% summarise(mean_soil_moisture = mean(`soil_moisture(v/v)`)) 




## -----------------------------------------------------------

## Temporal Stability Analysis

# Reload the data and rename some columns
Siloam <- readRDS('SPdata')
Siloam <- Siloam %>% rename (Depth_cm = `Depth(cm)`, soil_moisture =`soil_moisture(v/v)`,
Elevation =`Elevation(m)`) %>% 
mutate ( Year = year(DateTime), Month = month(DateTime), Day = day(DateTime) )

Siloam <- Siloam %>% mutate (Probe = as.factor(Probe))

saveRDS(Siloam, 'Siloam1')

Siloam <- readRDS('Siloam1')     # Actually just load tidy data set here to start Analysis

## compute the mean relative difference at various time scales and depths for all locations

Siloam <-  Siloam[complete.cases(Siloam),] 

## spatio-temporal xtics for the root zone and beyond

Siloam %>% select(Year, Probe, Depth_cm, soil_moisture) %>% filter(Depth_cm >= 90) %>% 
  mutate(mean_VSM = mean(soil_moisture), sd_VSM = sd(soil_moisture))


## relative mean difference at various emporal and spatial intsances

Siloam %>% filter(!(soil_moisture >= 100) & !(soil_moisture <= 0)) %>%        # change here for topsoil, root and below root zone
  group_by(Probe, Depth_cm, Month, Year) %>% 
  summarise(relDiff_VSM =(soil_moisture - mean(soil_moisture, na.rm=T))/mean(soil_moisture, na.rm = T), sd_VSM = sd(relDiff_VSM, na.rm=T), rank_relDiff = rank(relDiff_VSM)) %>% 
  #filter( !(Month == 4) & !(Month == 5) & !(Month == 6) & !(Month == 7) & !(Month == 8) & !(Month == 9)) %>%  # for the wet summer months
  filter(Depth_cm == 180 & !(Month == 10) & !(Month == 11) & !(Month == 12) & !(Month == 1) & !(Month == 2) & !(Month == 3)) %>%  # for the dry winter months
  ggplot() + geom_boxplot(aes(y = relDiff_VSM, x = rank_relDiff, fill = Probe))


## fixed code for relDiff for aggregated data (daily) and correlogram

Siloam$date<- as.Date(Siloam$DateTime, format = '%Y-%m-%d')
Siloam %>% group_by(Day, Probe, Depth_cm,lat, long)%>% summarise(soil_moisture = mean(soil_moisture)) -> Siloam_daily
Siloam_daily <- Siloam_daily[complete.cases(Siloam_daily),]    # edit from here
Siloam_daily %>% ungroup() %>% select(Day, Depth_cm, daily_mean_VSM )-> M

Siloam %>% filter(!(soil_moisture >= 100) & !(soil_moisture <= 0)) %>% 
group_by(Depth_cm,Probe,Day,Month, Year)%>% summarise(mean_time = mean(soil_moisture, nam.rm=T), delta_mean = soil_moisture - mean_time) %>% 
  ungroup() %>% group_by(Day,Month,Year,Probe)%>% summarise(field_mean = mean(mean_time), field_sd = sd(mean_time), 
                                                        relDiff_VSM = delta_mean/field_mean, rank_relDiff = rank(relDiff_VSM)) %>%
#filter( !(Month == 4) & !(Month == 5) & !(Month == 6) & !(Month == 7) & !(Month == 8) & !(Month == 9)) %>%  # for the wet summer months
filter(!(Month == 10) & !(Month == 11) & !(Month == 12) & !(Month == 1) & !(Month == 2) & !(Month == 3)) %>%  # for the dry winter months
  ggplot() + geom_boxplot(aes(y = relDiff_VSM, x = rank_relDiff, fill = Probe))


 # plot of temporal variability of soil moisture 

Siloam %>% filter(!(soil_moisture >= 100) & !(soil_moisture <= 0)) %>%        # change here for topsoil, root and below root zone
  group_by(Probe, Year) %>% 
  summarise(mean_VSM = mean(soil_moisture), sd_VSM = sd(soil_moisture)) %>% 
  ggplot() + geom_point(aes(y=mean_VSM, x=Year)) + geom_line(aes(x=Year, y=mean_VSM, group=Probe, colour=Probe)) 

### spatial plots with mean and standard deviation

Siloam %>% filter(!(soil_moisture <= 0) & !(soil_moisture >= 100)) %>% filter (Depth_cm <= 10 & Year == 2012) %>%  # change here for topsoil, root and below root zone
  group_by(Probe, Depth_cm, lat, long) %>% summarise(mean_VSM=mean(soil_moisture), sd_VSM=sd(soil_moisture)) %>% ungroup() %>%
  ggplot(aes(long, lat)) + geom_point(aes(long, lat, colour = mean_VSM, size = sd_VSM)) + geom_text(aes(label = Probe),nudge_x = 0.00005) +
  xlab('Longitude (deg)') + labs(title = '10 cm top soil variability-mean soil moisture')

### compute other statistical parameters

# CoV in space
Siloam %>% group_by(Day, Month, Year, Probe) %>% filter(!(soil_moisture <= 0) & !(soil_moisture >= 100))%>% group_by(Probe, Depth_cm, Day,Month, Year) %>% 
  summarise(mean_VSM = mean(soil_moisture, na.rm=T), sd_VSM = sd(soil_moisture, na.rm=T), CoV = (sd_VSM/mean_VSM))

## change to wide format to compute CoV in time

Siloam %>% group_by(Day, Month, Year, Probe) %>% filter(!(soil_moisture <= 0) & !(soil_moisture >= 100))%>% group_by(Day,Month, Year, Probe, Depth_cm) %>%   
    summarise(mean_VSM = mean(soil_moisture)) %>% spread(key = Probe, value = mean_VSM) -> tmp


tmp$temp_means <- rowMeans(tmp[,c(-1:-4)], na.rm = T)
gather(tmp, av_VSM,Probe, -Day,-Month,-Year,-Depth_cm)


# CoV in time
tmp <- tmp %>% mutate(sd_VSM = sd(c(`12699`,`17433`, `17434`, `17435`, `17437`, `17440`, `20916`, `20918`, `20922`), na.rm=T), 
              mean_VSM = mean(c(`12699`,`17433`, `17434`, `17435`, `17437`, `17440`, `20916`, `20918`, `20922`), na.rm=T), CoV = sd_VSM/mean_VSM)

# CoV in space

Siloam %>% group_by(Day, Month, Year, Probe) %>% filter(!(soil_moisture <= 0) & !(soil_moisture >= 100))%>% group_by(Probe, Depth_cm, Day,Month, Year) %>% 
  summarise(mean_VSM = mean(soil_moisture)) %>% spread(key = Day, value = mean_VSM) %>% 
  mutate(sd_VSM = sd(c(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`10`,`11`,`12`,`13`,`14`,`15`,`16`,`17`,`18`,`19`,`20`,`21`,`22`,`23`,`24`,`25`,`26`,`27`), na.rm = T), 
         mean_VSM = mean(c(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`10`,`11`,`12`,`13`,`14`,`15`,`16`,`17`,`18`,`19`,`20`,`21`,`22`,`23`,`24`,`25`,`26`,`27`), na.rm = T), 
         CoV =sd_VSM/mean_VSM) %>% select(CoV)%>% arrange(CoV)%>% tail()


## Spatial Anomaly plots

Siloam %>% group_by(Year, Day, Probe, lat, long, Depth_cm, soil_moisture) %>% summarise(watershed_VSM = mean(soil_moisture), Anomaly_VSM = soil_moisture - watershed_VSM) -> tmp

lim_lat <- range(tmp$lat)
lim_t <- range(tmp$Day)            # time range
lat_axis <- seq(lim_lat[1],       # latitude axis
                lim_lat[2],
                length=25)
t_axis <- seq(lim_t[1],           # time axis
              lim_t[2],
              length=100)
lat_t_grid <- expand.grid(lat = lat_axis,
                          t = t_axis)

tmp_grid <- tmp
dists <- abs(outer(tmp$lat, lat_axis, "-"))
tmp_grid$lat <- lat_axis[apply(dists, 1, which.min)]

tmp_lat_Hov <- group_by(tmp_grid, lat, Day) %>%
  summarise(Anomaly_VSM = mean(Anomaly_VSM))

Hovmoller_lat <- ggplot(tmp_lat_Hov) +            # take data
  geom_tile(aes(x = lat, y = Day, fill = Anomaly_VSM)) + # plot
  fill_scale(name = "VSM") +     # add color scale
  scale_y_reverse() +             # rev y scale
  ylab("Day number (days)") +     # add y label
  xlab("Latitude (degrees)") +    # add x label
  theme_bw()                      # change theme





### empirical Spatial and Temporal Means calculations and plots

Siloam1 <- filter(Siloam,     # subset the data
               Month %in% 5:9 &   # May to September
                 Year == 2012)      # year of 1993

Siloam1 <- Siloam          # Full data set
Siloam1$t <- hour(Siloam1$DateTime)
Siloam1$date <- as.Date(Siloam1$DateTime, format = '%Y-%m-%d')                              # added a date column
Siloam1_daily <- group_by(Siloam1, Day)%>% summarise(soil_moisture = mean(soil_moisture)) 

#Siloam1 <- Siloam1 %>% group_by(DateTime, Probe) %>% spread(DateTime, soil_moisture) 
#Siloam1 %>% group_by(Day, Probe, lat, long)%>% summarise(soil_moisture = mean(soil_moisture)) %>% ungroup()-> Siloam1


spat_av <- Siloam1 %>% filter(!(soil_moisture >= 100) & !(soil_moisture <= 0)) %>%      # group soil moisture by lat and long to compute the empirical spatial means at each location
  group_by(soil_moisture, lat, long) %>% summarise(mu_VSM = mean(soil_moisture))

lat_means <- ggplot(spat_av) +
  geom_point(aes(lat, mu_VSM)) +                                                       # plot of avg VSM at the lat and long 
  xlab("Latitude (deg)") +
  ylab("soil moisture (VSM)") + theme_bw()


lon_means <- ggplot(spat_av) +
  geom_point(aes(long, mu_VSM)) +
  xlab("Longitude (deg)") +
  ylab("soil moisture (VSM)") + theme_bw()



## computation and plotting of the empirical temporal means

Siloam1 <-  filter(Siloam1, !(soil_moisture >= 100) & !(soil_moisture <= 0)) 
  
  Tmax_av <- group_by(Siloam1,date) %>%
  summarise(meanVSM = mean(soil_moisture))

gTmaxav <-
  ggplot() +
  geom_line(data = Siloam1,aes(x = date, y = soil_moisture, group = Probe),
            colour = "blue", alpha = 0.04) +
  geom_line(data = Tmax_av, aes(x = date, y = meanVSM)) + facet_wrap(~Probe)
  xlab("Year") + ylab("soil moisture (VSM)") +
  theme_bw()


## computationand plotting of Empirical covariances (uncomment to subset)

#Siloam1 <- filter(Siloam,   !(soil_moisture >= 100) & !(soil_moisture <= 0) &  # subset the data
#                  Month %in% 5:9 &   # May to September
#                    Year == 2012)      # year of 1993
#Siloam1$t <- hour(Siloam1$DateTime)

Siloam1 <- Siloam1 %>% filter(Depth_cm == 60 & Year == 2012)      # subset to specific depth and Year

lm1 <- lm(soil_moisture ~ lat + t + I(t^2), data = Siloam1) # fit a linear model with daily data for Siloam
Siloam1$residuals <- residuals(lm1)      # store the residuals


spat_df <- filter(Siloam1, t == 0) %>%      #  lon/lat coords of stations
  select(long, lat) %>%                         #select lon/lat only
  arrange(long, lat)           #sort ascending by lon/lat 
m <- nrow(spat_av)            #number of stations

#  converting data into space wide format

X <- select(Siloam1, long, lat, residuals, t, DateTime) %>%       # select columns
  group_by(t, DateTime) %>% spread(DateTime, residuals) %>% ungroup() %>%                       # make time-wide
  select(-long, -lat, -DateTime) %>%                     # drop coord info
  t()                                         # make space-wide

# create lags

Lag0_cov <- cov(X, use = 'complete.obs')
Lag1_cov <- cov(X[-1, ], X[-nrow(X),], use = 'complete.obs')


spat_df$n <- 1:nrow(spat_df)       # assign an index to each station
lim_lon <- range(spat_df$long)        # range of lon coordinates
lon_strips <- seq(lim_lon[1],                  # create 4 long. strip boundaries
                  lim_lon[2],
                  length = 5)                        # bin the lon into
spat_df$lon_strip <- cut(spat_df$long, lon_strips,                           # their respective bins
                         labels = FALSE, # don't assign labels
                         include.lowest = TRUE) # include edges

plot_cov_strips(Lag0_cov, spat_df)   # plot the lag-0 matrices
plot_cov_strips(Lag1_cov, spat_df)   # plot the lag-1 matrices
