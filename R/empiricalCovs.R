library("readxl")
library("tidyverse")
library("ggplot2")
library("gstat")
library("maps")
library("STRbook")
library("grid")
library("gridExtra")
library("lubridate")
library("sp")
library("spacetime")

setwd('/home/ivo/Downloads/PhD-Data-code/')
## ------------------------------------------------------------------------


Siloam <- readRDS('SPdata')
Siloam <- Siloam %>% rename (Depth_cm = `Depth(cm)`, soil_moisture =`soil_moisture(v/v)`,
                             Elevation =`Elevation(m)`) %>% 
  mutate ( Year = year(DateTime), Month = month(DateTime), Day = day(DateTime) )

Siloam <- Siloam %>% mutate (Probe = as.factor(Probe))

saveRDS(Siloam, 'Siloam1')
Siloam <- readRDS('Siloam1')

## compute the mean relative difference at various time scales and depths for all locations

Siloam <- filter(Siloam,!(soil_moisture >= 100) & !(soil_moisture <= 0)) #& Depth_cm == 30 & Year == 2012)

Siloam$date <- as.Date(Siloam$DateTime, format = '%Y-%m-%d')                              # added a date column

Siloam <-  Siloam[complete.cases(Siloam),] 

Siloam <- group_by(Siloam,Day, Month, Year, Probe, lat, long )%>% summarise(soil_moisture = mean(soil_moisture, na.rm=T))

spat_av <- Siloam %>%                                                  # group soil moisture by lat and long to compute the empirical spatial means at each location
  group_by(soil_moisture, lat, long) %>% summarise(mu_VSM = mean(soil_moisture))


lm1 <- lm(soil_moisture ~ lat + Day + I(Day^2), data = Siloam)     # fit a linear model with daily data for Siloam

#Siloam <- Siloam[c(1:2790),]

Siloam$residuals <- residuals(lm1)      # store the residuals

spat_df <- filter(Siloam, Day ==2) %>%  ungroup() %>%    #  lon/lat coords of stations
  select(long, lat) %>%                         #select lon/lat only
  arrange(long, lat)           #sort ascending by lon/lat 

#spat_df <- read.csv('spat_df.csv', header = T) %>% select(long, lat)    # manualy created spat_df table
m <- nrow(spat_av)            #number of stations


X <- Siloam %>% filter(Depth_cm == 10 &  Day== 2 & Month %in% 6:12 & Year == 2012) %>% select(long, lat, residuals, Day) %>%       # select columns
  spread(Day, residuals) %>%                  # make time-wide
  select(-long, -lat) %>%                     # drop coord info
  t()                                         # make space-wide


# create lags

Lag0_cov <- cov(X, use = 'complete.obs')
Lag1_cov <- cov(X[-1, ], X[-nrow(X),], use = 'complete.obs')

spat_df$n <- 1:nrow(spat_df)   # assign an index to each station
lim_lon <- range(spat_df$long)  # range of lon coordinates
lon_strips <- seq(lim_lon[1],                    # create 4 long. strip boundaries
                  lim_lon[2],
                  length = 6)                  # bin the lon into

spat_df$lon_strip <- cut(spat_df$long,
                         lon_strips,                           # their respective bins
                         labels = FALSE,                 # don't assign labels
                         include.lowest = TRUE)        # include edges


plot_cov_strips(Lag0_cov, spat_df)  # plot the lag-0 matrices
plot_cov_strips(Lag1_cov, spat_df)  # plot the lag-1 matrices

##############################################################################################################################################

Siloam <- readRDS('SPdata')
Siloam <- Siloam %>% rename (Depth_cm = `Depth(cm)`, soil_moisture =`soil_moisture(v/v)`,
                             Elevation =`Elevation(m)`) %>% 
  mutate ( Year = year(DateTime), Month = month(DateTime), Day = day(DateTime) )

Siloam <- Siloam %>% mutate (Probe = as.factor(Probe))
saveRDS(Siloam, 'Siloam1')
Siloam <- readRDS('Siloam1')

Siloam <- filter(Siloam,!(soil_moisture >= 100) & !(soil_moisture <= 0)) 
Siloam$date <- as.Date(Siloam$DateTime, format = '%Y-%m-%d')                              # added a date column
Siloam <-  Siloam[complete.cases(Siloam),] 

# Integrate measurements over average profile depth

Siloam_daily <- Siloam %>% group_by(Probe, Depth_cm, date, lat, long) %>% summarise(daily_VSM = mean(soil_moisture)) %>% ungroup()         # aggregate soil moisture to daily for each site

Siloam_monthly <- Siloam %>% group_by(Month, Day, date, Probe, Depth_cm) %>% summarise(monthly_VSM = mean(soil_moisture)) %>% ungroup()

 
Siloam_daily %>%  mutate(Day = day(date), Month = month(date), Year = year(date)) %>%                                                                                  
  group_by(Probe, Depth_cm, date, lat, long) %>%                                               # aggegate the daily soil moisture for each monitoring site
  spread(Depth_cm,daily_VSM) %>% select(-lat,-long,-Probe,-date, -Month, -Year) %>% rowwise() %>% 
  summarise(m = sum(c(`10`*10,`20`*10,`30`*10,`40`*10,`60`*20,`80`*10,`90`*10,`120`*30,`150`*20,`180`*30)/120, na.rm=T)) -> profile_daily_VSM

profile_daily_VSM <- profile_daily_VSM %>% ungroup() %>% mutate(Month = month(date), Year = year(date))


 ### calculation of the Random Combination Method for NRL for long probes 

mu_bar_12699 <- Siloam %>% filter(Probe == 12699) %>% group_by (Day, Month, Year, Depth_cm) %>% summarise(mu = mean(soil_moisture)) %>% ungroup() %>% spread(Depth_cm, mu)%>%    # daily average VSM at location 12699
   select(-Day,-Month,-Year) %>% rowwise() %>% mutate(mu_12699 = 1/6*sum(`30`,`60`,`90`,`120`,`150`,`180`, na.rm = T)) 
tmp <- Siloam %>% filter(Probe == 12699) %>% group_by (Day, Month, Year, Depth_cm) %>% summarise(mu = mean(soil_moisture)) %>% ungroup() %>% spread(Depth_cm, mu) %>% select(Day,Month,Year)
mu_bar_12699 <- cbind(mu_bar_12699,tmp) %>% select(Day,Month,Year, mu_12699) 

mu_bar_20922 <- Siloam %>% filter(Probe == 20922) %>% group_by (Day, Month, Year, Depth_cm) %>% summarise(mu = mean(soil_moisture)) %>% ungroup() %>% spread(Depth_cm, mu) %>%    # daily average VSM at location 12699
  select(-Day,-Month,-Year) %>% rowwise() %>% mutate(mu_20922 = 1/6*sum(`30`,`60`,`90`,`120`,`150`,`180`, na.rm = T))
tmp <- Siloam %>% filter(Probe == 20922) %>% group_by (Day, Month, Year, Depth_cm) %>% summarise(mu = mean(soil_moisture)) %>% ungroup() %>% spread(Depth_cm, mu) %>% select(Day,Month,Year)
mu_bar_20922 <- cbind(tmp, mu_bar_20922) %>% select(Day,Month,Year, mu_20922)

mu_bar_20916 <- Siloam %>% filter(Probe == 20916) %>% group_by (Day, Month, Year, Depth_cm) %>% summarise(mu = mean(soil_moisture)) %>% ungroup() %>% spread(Depth_cm, mu) %>%    # daily average VSM at location 12699
  select(-Day,-Month,-Year) %>% rowwise() %>% mutate(mu_20916 = 1/6*sum(`30`,`60`,`90`,`120`,`150`,`180`, na.rm = T))
tmp <- Siloam %>% filter(Probe == 20916) %>% group_by (Day, Month, Year, Depth_cm) %>% summarise(mu = mean(soil_moisture)) %>% ungroup() %>% spread(Depth_cm, mu) %>% select(Day,Month,Year)
mu_bar_20916 <- cbind(tmp, mu_bar_20916)%>% select(Day,Month,Year, mu_20916)

mu_bar_20918 <- Siloam %>% filter(Probe == 20918) %>% group_by (Day, Month, Year, Depth_cm) %>% summarise(mu = mean(soil_moisture)) %>% ungroup() %>% spread(Depth_cm, mu) %>%    # daily average VSM at location 12699
  select(-Day,-Month,-Year) %>% rowwise() %>% mutate(mu_20918 = 1/6*sum(`30`,`60`,`90`,`120`,`150`,`180`, na.rm = T))
tmp <- Siloam %>% filter(Probe == 20918) %>% group_by (Day, Month, Year, Depth_cm) %>% summarise(mu = mean(soil_moisture)) %>% ungroup() %>% spread(Depth_cm, mu) %>% select(Day,Month,Year)
mu_bar_20918 <- cbind(tmp, mu_bar_20918) %>% select(Day,Month,Year, mu_20918)

### computing the RMSD

read.csv('longprobes_dailymean_profile.csv', header=T, stringsAsFactors = F) %>% group_by(Day,Month,Year) %>% summarise(field_mu = mean(soil_moisture, na.rm=T))%>% ungroup() -> fieldmean_VSM   # read the field mean data

P_12699 <- merge(mu_bar_12699, fieldmean_VSM)
P_20922 <- merge(mu_bar_20922, fieldmean_VSM)
P_20916 <- merge(mu_bar_20916, fieldmean_VSM)
P_20918 <- merge(mu_bar_20918, fieldmean_VSM)

cbind(P_12699[c(1:477),],P_20916[c(1:477),c(-1:-3)],P_20918[,c(-1:-3)],P_20922[c(1:477),c(-1:-3)]) ->tmp    # combine new data frame with all locations
tmp <- tmp[,c(-5,-7,-9)]
tmp <- tmp %>% select(-Day,-Month,-Year)%>% rowwise() %>% mutate(rel_diff = 1/4*sum((mu_12699-field_mu)^2, (mu_20918-field_mu)^2, (mu_20916-field_mu)^2, (mu_20922-field_mu)^2)) %>% mutate(rmsd = sqrt(sum(rel_diff)/477))   # RMSD for long probes
tmp <- gather(tmp, id, z, -field_mu, -rmsd) 
arrange(tmp)  # to get the location with highest value of RMSD









# Visualize the temporal variation of VSM

profile_daily_VSM %>% select (m, Month, Year) %>%
           mutate(Month2 = as.Date(paste0('2013-' , profile_daily_VSM$Month, '-01'), '%Y-%m-%d')) %>% 
  
ggplot(aes(x = Month2, y = m)) + 
  geom_bar(stat = 'identity', fill = 'darkorchid4') + 
  facet_wrap(~Year, ncol = 3) +
  labs(title= 'Daily Soil Profile', 
       subtitle = 'plotted by year', 
       y = 'volumetric soil moisture',
       x = 'Month') + theme_bw(base_size = 15) + scale_x_date(date_labels = '%b')



ggplot(data = profile_daily_VSM, aes(x = date, y = m, group = Probe, colour = Probe)) + geom_line() + xlab('Year') + ylab("Daily Profile Soil Moisture (VSM)") + theme_bw()


# Empirical Temporal Means with Daily Data

field_daily_av <- group_by(profile_daily_VSM, date) %>% summarise(mean_profile_VSM = mean(m))

# plot of field average and profile daily average
  ggplot() +
  geom_line(data = profile_daily_VSM,aes(x = date, y = m, group = Probe),
            colour = "blue", alpha = 0.04) +
  geom_line(data = field_daily_av, aes(x = date, y = mean_profile_VSM)) +
  xlab("Year") + ylab("Field daily Average Soil Moisture (VSM)") +
  theme_bw()


# Empirical Spatial Means with Daily Data  

  spat_av <- group_by(profile_daily_VSM, lat, long) %>%
    summarise(mu_emp = mean(m))

  lat_means <- ggplot(spat_av) +
    geom_point(aes(lat, mu_emp)) +
    xlab("Latitude (deg)") +
    ylab("Daily Soil Moisture (VSM)") + theme_bw()
  lon_means <- ggplot(spat_av) +
    geom_point(aes(long, mu_emp)) +
    xlab("Longitude (deg)") +
    ylab("Daily Soil Moisture (VSM)") + theme_bw()
  
  
  
# Empirical covariances
  
  profile_daily_VSM <- mutate(profile_daily_VSM, t = day(date))
  
  lm1 <- lm(m ~ lat + t + I(t^2), data = profile_daily_VSM) # fit a linear model
  profile_daily_VSM$residuals <- residuals(lm1)
  
  spat_df <- filter(profile_daily_VSM, t == 1) %>%
    select(long, lat) %>%
    arrange(long, lat)
  s <- nrow(spat_av)

  profile_daily_VSM <- read.csv('profile_daily.csv', header = T, stringsAsFactors=F)     # manually fixed file
  profile_daily_VSM <- profile_daily_VSM[,c(-1,-9)]
  
  X <- select(profile_daily_VSM, long, lat, residuals, t) %>% # select columns
    spread(t, residuals) %>%      # make time-wide
    select(-long, -lat) %>%          # drop coord info
    t()
  
  Lag0_cov <- cov(X, use = 'complete.obs')
  Lag1_cov <- cov(X[-1, ], X[-nrow(X),], use = 'complete.obs')

  
  spat_df$n <- 1:nrow(spat_df)        # assign an index to each station
  lim_lon <- range(spat_df$long)         # range of lon coordinates
  lon_strips <- seq(lim_lon[1],                         # create 4 long. strip boundaries
                    lim_lon[2],
                    length = 5)                   # bin the lon into
  spat_df$lon_strip <- cut(spat_df$long,
                           lon_strips,                 # their respective bins
                           labels = FALSE,               # don't assign labels
                           include.lowest = TRUE)
  
  head(spat_df)
  
  plot_cov_strips(Lag0_cov, spat_df)
  plot_cov_strips(Lag1_cov, spat_df)
  

      
#  Semivariogram Analysis

Siloam %>% group_by(Probe, Depth_cm, date, lat, long) %>% summarise(soil_moisture_mu = mean(soil_moisture)) %>%         # aggregate to daily 
  mutate(Day = day(date), Month = month(date), Year = year(date))  %>% write.csv('mu_all_probes_daily.csv')                                                 

tmp <- read.csv('mu_all_probes_daily.csv', header = T, stringsAsFactors = F)
tmp$date <- as.Date(tmp$date, "%Y-%m-%d")
tmp$Probe <- as.factor(tmp$Probe)
tmp <- select(tmp, -X)

###### test covariance
tmp %>% filter(Depth_cm ==  30 & Month %in% 5:9) %>% write_csv('compare.csv')
tmp <- read.csv('compare.csv', header=T, stringsAsFactors=F)

tmp <- rename(tmp, t='count', z='soil_moisture_mu',)
spat_av <- group_by(tmp, lat, long) %>%
  summarise(mu_emp = mean(z))


lm1 <- lm(z~ lat + t + I(t^2), data = tmp)
tmp$residuals <- residuals(lm1)

spat_df <- filter(tmp, t==1) %>%
  select(long, lat) %>%
  arrange(long, lat)
m <- nrow(spat_av)

X <- select(tmp, long, lat, residuals, t) %>%    #   select columns
  spread(t, residuals) %>%                       # make time-wide
  select(-long, -lat) %>%
  t()


Lag0_cov <- cov(X, use = 'complete.obs')
Lag1_cov <- cov(X[-1, ], X[-nrow(X),], use = 'complete.obs')


spat_df$n <- 1:nrow(spat_df)   # assign an index to each station
lim_lon <- range(spat_df$lon)  # range of lon coordinates
lon_strips <- seq(lim_lon[1],
                  lim_lon[2],
                  length = 5)   # create 4 long. strip boundaries

spat_df$lon_strip <- cut(spat_df$lon,                         # bin the lon into
                         lon_strips,                          # their respective bins
                         labels = FALSE,                      # don't assign labels
                         include.lowest = TRUE)

head(spat_df)  

plot_cov_strips(Lag0_cov, spat_df)
plot_cov_strips(Lag1_cov, spat_df)



#####################################################################################

Siloam_daily   # use for aggregated daily VSM 

# creation of the STFDF

#spatial part
spat_part <- read.csv('spatial_location_probes.csv', header = T, stringsAsFactors = F)
spat_part <- spat_part[c(-5, -6,-8, -13),] %>%select(lat,long) 
spat_part <- SpatialPoints(coords = spat_part[, c("long", "lat")])
length(spat_part)

#temporal part
tmp %>% filter(Depth_cm == 120 & Month %in% 5:9) -> temp_part                #filter(Probe == 20918 & Depth_cm == 120)-> temp_part
temp_part <- temp_part$date
temp_part <- as.Date(temp_part, '%Y-%m-%d')
temp_part <- temp_part %>% sort()
length(temp_part)

# data
data <- tmp %>% select(-Depth_cm,-date,-lat,-long)
nrow(data)



STObj3 <- STFDF(sp = spat_part,
                time = temp_part, data = data[c(1:6201),])

project4string(STObj3) <- CRS('+proj=longlat + ellps=WGS84')

vv <- variogram(object = z ~ 1 + lat,
          data = STObj3,
          width = 10,
          cutoff = 1000,
          tlags = 0.01:6.01)

vv <- vv[complete.cases(vv),]
plot(vv)


#### EOFs and CCA

# read the data 

Siloam %>% group_by(Probe, Depth_cm, date, lat, long) %>% summarise(soil_moisture_mu = mean(soil_moisture)) %>%         # aggregate to daily 
  mutate(Day = day(date), Month = month(date), Year = year(date)) -> tmp

## Put data into space-wide form

Siloam_10 <- tmp %>% filter(Depth_cm == 10) %>% ungroup()%>% select(date,soil_moisture_mu, Probe)%>% spread(Probe, soil_moisture_mu) %>% select(-date)    #data
Siloam_latlong <- tmp %>% filter(Depth_cm == 10) %>% ungroup()%>% select(date,soil_moisture_mu, Probe, lat, long)%>% spread(Probe, soil_moisture_mu) %>% select(lat, long)
Z <- Siloam_10
dim(Z)

## First find the matrix we need to subtract:
spat_mean <- tmp %>% filter(Depth_cm == 10) %>% ungroup() %>% select(date, Probe, soil_moisture_mu) %>% spread(Probe, soil_moisture_mu) %>% select(-date) %>% apply(1,mean)
nT <- tmp %>%filter(Depth_cm == 10) %>% ungroup() %>% select(date, Probe, soil_moisture_mu) %>% spread(Probe, soil_moisture_mu) %>% select(-date) %>% ncol()

## Then subtract and standardize:
Zspat_detrend <- Z - outer(rep(1, nT), spat_mean)
Zt <- 1/sqrt(nT - 1)*Zspat_detrend
E <- Zt[complete.cases(Zt),] %>% svd(Zt)


V <- E$v
colnames(E$v) <- paste0("EOF", 1:ncol(Siloam_10)) # label columns
EOFs <- cbind(Siloam_latlong[c(1:2261),], E$v)

head(EOFs[, 1:6])
TS <- data.frame(E$u) %>%
  mutate(t = 1:nrow(E$u)) %>%
  gather(EOF, PC, -t) 

TS$nPC <- TS$PC * sqrt(nT-1)
ggplot(EOFs) + geom_tile(aes(x = long, y = lat, fill = EOF1)) +
  fill_scale(name = "VSM") + theme_bw() +
  xlab("Longitude (deg)") + ylab("Latitude (deg)")
