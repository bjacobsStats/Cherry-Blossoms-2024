#### team members: benjamin jacobs, spencer wadsworth, kaushiki singh, aditya ranade,
## phd students at iowa state university
#### title: cherry blossom random coefficients model
#### last update: feb 29, 2024
#### description: in this document, I read in data from the github, make some graphics
# to explore it, and propose a random coefficients model. The model is motivated by the graphics,
# and my (weak) understanding of the science involved, 
# which I think would help the story-telling element of our submission.
#
# the model assumes a linear relationship between altitude and bloom date, and between 
# latitude and bloom date; these relationships change randomly year to year, with an ar(1) process.
# regions (DC vs other) also have their own effect, and there is a random longitude effect which 
# is correlated to account for correlation between data with similar longitudes.
#
# the model is not totally adequate, the residuals are too heavily tailed for a gaussian distribution
#
# this model is still not as accurate as i'd like. i have some notions of how to improve it:
# i.  add in weather-related covariates. i've played around with this, and it hasn't seemed to 
#     help. still, it should be considered.
# ii. fix the longitude correlation, which i have currently as ar(1) because it was easy to code.
#     we'll want something that takes the circular shape of the earth into account...
# 
# another idea is to re-implement this in stan to see if it gives the same results, and get 
# prediction intervals which i'm having trouble getting from inla. we can just use the credible intervals
# for the fitted values as prediction intervals, they seem pretty wide as is...

#### dependencies ####
library(ggplot2)
library(dplyr)
library(INLA)


##### read in data #####
japan <- read.csv("https://raw.githubusercontent.com/GMU-CherryBlossomCompetition/peak-bloom-prediction/main/data/japan.csv",header=T)
japan <- japan[japan$bloom_doy>=70,] 
kyoto <- read.csv("https://raw.githubusercontent.com/GMU-CherryBlossomCompetition/peak-bloom-prediction/main/data/kyoto.csv",header=T)
liestal <- read.csv("https://raw.githubusercontent.com/GMU-CherryBlossomCompetition/peak-bloom-prediction/main/data/liestal.csv",header=T)
meteoswiss <- read.csv("https://raw.githubusercontent.com/GMU-CherryBlossomCompetition/peak-bloom-prediction/main/data/meteoswiss.csv",header=T)
south.korea <- read.csv("https://raw.githubusercontent.com/GMU-CherryBlossomCompetition/peak-bloom-prediction/main/data/south_korea.csv",header=T)
vancouver <- read.csv("https://raw.githubusercontent.com/GMU-CherryBlossomCompetition/peak-bloom-prediction/main/data/vancouver.csv",header=T)
dc <- read.csv("https://raw.githubusercontent.com/GMU-CherryBlossomCompetition/peak-bloom-prediction/main/data/washingtondc.csv",header=T)

# useful info about the prediction locations:

# Location	Latitude (°)	Longitude (°)	Altitude (m)	Years available	Bloom definition	Species
# Kyoto (Japan)	35.0120	135.6761	44	801–2023	80%	Prunus jamasakura
# Liestal-Weideli (Switzerland)	47.4814	7.730519	350	1895–2023	25%	Prunus avium
# Washington, D.C. (USA)	38.8853	–77.0386	0	1921–2023	70%	Prunus × yedoensis ‘Somei-yoshino’
# Vancouver, BC (Canada)	49.2237	–123.1636	24	2022–20231	70%	Prunus × yedoensis ‘Akebono’
# New York City, NY (USA)	40.73040	–73.99809	8.5	2019–20232	70%	Prunus × yedoensis


# 'region' is an indicator for washington dc; i was having trouble predicting dc but this seems to have 
# helped by giving it its own effect. i'm expecting we'll want to code new york as TRUE
japan$region <- F 
kyoto$region <- F 
liestal$region <- F 
meteoswiss$region <- F 
south.korea$region <- F  
dc$region <- T # New york will also be this category
vancouver$region <- F

vancouver$long <- 236.8364 # adjust vancouver's longitude so it will regress on japan, not washington dc

all.data <- rbind(japan,kyoto,liestal,meteoswiss,south.korea,dc,vancouver)
all.data$year2 <- all.data$year # needed for INLA's syntax
all.data <- filter(all.data,year>=1950) # throwing out old data on the theory that trends have changed

##### visualizing the data in various and sundry ways #####

hist(japan$bloom_doy) 
hist(kyoto$bloom_doy)
hist(liestal$bloom_doy)
hist(meteoswiss$bloom_doy)
hist(south.korea$bloom_doy)
hist(dc$bloom_doy)
hist(vancouver$bloom_doy) # obvisously not a lot of data, but both fall in the bottom half of the other data sets
plot(all.data$alt,all.data$bloom_doy) # clear linear trend 
plot(all.data$lat,all.data$bloom_doy) # clear linear trend
plot(all.data$long,all.data$bloom_doy) # not so much of an effect
plot(all.data$year,all.data$bloom_doy) # honestly not seeing a pattern, if anything it just gets more varied, but that probably has to do with more reporting in recent years
plot(all.data$lat,all.data$alt) # clearly not independent, northernmost points also have highest elevation. this may lead to confounding?
plot(all.data$long,all.data$alt) # also clearly not independent
ggplot(all.data)+geom_boxplot(aes(x=bloom_doy,y=region))

# here's a sort of topographic map of our observations
ggplot(all.data)+geom_point(aes(x=long,y=lat,color=alt))





###### my current best model #####

inlaB <- inla(bloom_doy~ region+
                f(year,alt,model="ar",order=1)+ # higher orders are really slow to fit, and don't help
                f(year2,lat,model="ar",order=1)+ # random altitude and latitude effects for each year, with ar1 correlation; i looked at adding an ineraction but it didn't seem to help
                f(long,model="ar",order=1), # this is the wrong correlation i think, but it works for now; higher orders of three or more mess everything up
              family="gaussian",data=all.data,verbose=T) # gaussian seems to be working well enough

qqnorm(inlaB$summary.fitted.values[,4]-all.data$bloom_doy)
abline(0,sd(inlaB$summary.fitted.values[,4]-all.data$bloom_doy)) # fat tails, how can we deal with this? do we need to?

mean(abs(inlaB$summary.fitted.values[,4]-all.data$bloom_doy)) # mean abs error: 3.971098
median(abs(inlaB$summary.fitted.values[,4]-all.data$bloom_doy)) # median abs error: 2.89277

ggplot()+
  geom_point(aes(x=inlaB$summary.fitted.values[,4],y=all.data$bloom_doy-inlaB$summary.fitted.values[,4]))+
  labs(x="fitted values",y="residuals")# residual plot looks roughly adequate


###### 2020 test ######
# there's no vancouver data from this year

test <- filter(all.data,year!=2020 & year!=2021 & year!=2023 & year!=2022)
test <- rbind(test,
              data.frame(
                location=c("kyoto","liestal","dc","vancouver"),
                bloom_date=rep(NA,4),
                bloom_doy=rep(NA,4),
                lat=c(35.0120,47.4814,38.8853,49.2237),
                long=c(135.6761,7.730519,-77.0386,-123.1636),
                alt= c(44,350,0,24),
                year= c(2020,2020,2020,2020),
                year2= c(2020,2020,2020,2020),
                region=c(F,F,T,F)
              ))
inla20<- inla(bloom_doy~ region+
                f(year,alt,model="ar",order=1)+ # higher orders are really slow to fit, and don't help
                f(year2,lat,model="ar",order=1)+ # random altitude and latitude effects for each year, with ar1 correlation
                f(long,model="ar",order=1), # this is the wrong correlation i think, but it works for now; higher orders of three or more mess everything up
              family="gaussian",data=test,verbose=T) 

# errors:
floor(tail(inla20$summary.fitted.values,4)[,4])[1:3]-c(kyoto$bloom_doy[kyoto$year==2020],
                                                       liestal$bloom_doy[liestal$year==2020],
                                                       dc$bloom_doy[dc$year==2020])
# average error of 6.333 days...



###### 2021 test ######
# there's no vancouver data from this year

test <- filter(all.data,year!=2021 &year!=2023 & year!=2022)
test <- rbind(test,
              data.frame(
                location=c("kyoto","liestal","dc","vancouver"),
                bloom_date=rep(NA,4),
                bloom_doy=rep(NA,4),
                lat=c(35.0120,47.4814,38.8853,49.2237),
                long=c(135.6761,7.730519,-77.0386,-123.1636),
                alt= c(44,350,0,24),
                year= c(2021,2021,2021,2021),
                year2= c(2021,2021,2021,2021),
                region=c(F,F,T,F)
              ))

inla21<- inla(bloom_doy~ region+
                 f(year,alt,model="ar",order=1)+ # higher orders are really slow to fit, and don't help
                 f(year2,lat,model="ar",order=1)+ # random altitude and latitude effects for each year, with ar1 correlation
                 f(long,model="ar",order=1), # this is the wrong correlation i think, but it works for now; higher orders of three or more mess everything up
               family="gaussian",data=test,verbose=T) 

# errors:
floor(tail(inla21$summary.fitted.values,4)[,4])[1:3]-c(kyoto$bloom_doy[kyoto$year==2021],
                                                  liestal$bloom_doy[liestal$year==2021],
                                                  dc$bloom_doy[dc$year==2021])
# average error of 3 days...

###### 2022 test ######
set.seed(432432)
test <- filter(all.data,year!=2023 & year!=2022)
test <- rbind(test,
              data.frame(
                location=c("kyoto","liestal","dc","vancouver"),
                bloom_date=rep(NA,4),
                bloom_doy=rep(NA,4),
                lat=c(35.0120,47.4814,38.8853,49.2237),
                long=c(135.6761,7.730519,-77.0386,-123.1636),
                alt= c(44,350,0,24),
                year= c(2022,2022,2022,2022),
                year2= c(2022,2022,2022,2022),
                region=c(F,F,T,F)
              ))
inla22 <- inla(bloom_doy~ region+
                f(year,alt,model="ar",order=1)+ # higher orders are really slow to fit, and don't help
                f(year2,lat,model="ar",order=1)+ # random altitude and latitude effects for each year, with ar1 correlation
                f(long,model="ar",order=1), # this is the wrong correlation i think, but it works for now; higher orders of three or more mess everything up
              family="gaussian",data=test,verbose=T) 

# errors:
floor(tail(inla22$summary.fitted.values,4)[,4])-c(kyoto$bloom_doy[kyoto$year==2022],
                                                 liestal$bloom_doy[liestal$year==2022],
                                                 dc$bloom_doy[dc$year==2022],
                                                 vancouver$bloom_doy[vancouver$year==2022])

# an average absolute error of 3 days



####### 2023 test ######

test <- filter(all.data,year!=2023)
test <- rbind(test,
              data.frame(
                location=c("kyoto","liestal","dc","vancouver"),
                bloom_date=rep(NA,4),
                bloom_doy=rep(NA,4),
                lat=c(35.0120,47.4814,38.8853,49.2237),
                long=c(135.6761,7.730519,-77.0386,-123.1636),
                alt= c(44,350,0,24),
                year= c(2023,2023,2023,2023),
                year2= c(2023,2023,2023,2023),
                region=c(F,F,T,F)
              ))
inla23 <- inla(bloom_doy~region+
                f(year,alt,model="ar",order=1)+ # higher orders are really slow to fit, and don't help
                f(year2,lat,model="ar",order=1)+ # random altitude and latitude effects for each year, with ar1 correlation
                f(long,model="ar",order=1), # this is the wrong correlation i think, but it works for now; higher orders of three or more mess everything up
              family="gaussian",data=test,verbose=F) 

# errors:
floor(tail(inla23$summary.fitted.values,4)[,4])-c(kyoto$bloom_doy[kyoto$year==2023],
                                                 liestal$bloom_doy[liestal$year==2023],
                                                 dc$bloom_doy[dc$year==2023],
                                                 vancouver$bloom_doy[vancouver$year==2023])
# vancouver is still bad with an error of 9, but the average mean error 
# of the other 3 is 2....



##### 2024 Prediction ######

# New York City, NY (USA)	40.73040	–73.99809	8.5	2019–20232	70%	Prunus × yedoensis
pred <- rbind(all.data,
              data.frame(
                location=c("kyoto","liestal","dc","vancouver","new york"),
                bloom_date=rep(NA,5),
                bloom_doy=rep(NA,5),
                lat=c(35.0120,47.4814,38.8853,49.2237,40.73040),
                long=c(135.6761,7.730519,-77.0386,-123.1636,-73.99809),
                alt= c(44,350,0,24,8.5),
                year= c(2024,2024,2024,2024,2024),
                year2= c(2024,2024,2024,2024,2024),
                region=c(F,F,T,F,T)
              ))
inla24 <- inla(bloom_doy~ region+
                 f(year,alt,model="ar",order=1)+ # higher orders are really slow to fit, and don't help
                 f(year2,lat,model="ar",order=1)+ # random altitude and latitude effects for each year, with ar1 correlation
                 f(long,model="ar",order=1), # this is the wrong correlation i think, but it works for now; higher orders of three or more mess everything up
               family="gaussian",data=pred,verbose=T) 
tail(inla24$summary.fitted.values,5) # predicting later than last year across the board. new york is ten days behind dc, does that sound right?

# comparing the predictions to historical distributions
preds24 <- floor(tail(inla24$summary.fitted.values,5)[,4])
preds24-c(kyoto$bloom_doy[kyoto$year==2023],
          liestal$bloom_doy[liestal$year==2023],
          dc$bloom_doy[dc$year==2023],
          vancouver$bloom_doy[vancouver$year==2023],NA)

hist(kyoto$bloom_doy)
abline(v=preds24[1])

hist(liestal$bloom_doy)
abline(v=preds24[2])

hist(dc$bloom_doy)
abline(v=preds24[3])

hist(vancouver$bloom_doy)
abline(v=preds24[4])

# and there's no data to compare the new york prediction to
