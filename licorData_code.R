library(tidyverse) # for readr, dplyr, ggplot...
library(lubridate) # for dealing with dates and times (e.g., hms)
library(MCMCglmm) # for MCMCglmm stats
library(gridExtra) # for plotting multiple plots in a single frame
library(Rmisc) # for summarySE

setwd("~/Dropbox/andrea_plankton/clean_for_upload/light_temp/")

# Load the file
d0 <- read_csv("licor_data.csv")
head(d0)
summary(d0)

# make depth a factor, not a continuous variable
d0$depth <- as.factor(d0$depth)

# get rid of spaces in site names
d0$siteName <- gsub(" ", "", d0$siteName)
summary(d0)

# plot licor input over time
d0 %>% ggplot(aes(x = localTime, y = licorData, group = siteType, color = siteType)) +
  geom_line()+
  facet_grid(siteName ~ .)+
  theme_bw()

# Get rid of the first 4 measurements of each depth for each site (time during descent)----
d1 <- d0 %>% 
  group_by(siteName,depth) %>% 
  slice(5:n())
head(d1)

# Gid rid of measurements at Bastimentos South at 6 ft after 12:35 
# (meter was taken out but not turned off)
endBS6 <- hms("12:35:00")

d2 <- d1 %>%
  filter(!(siteName=="BASTIMSOUTH" & depth==6 & hms(localTime)>endBS6))
head(d2)

# Get rid of any measurements before 11 to 2
dstart <- hms("11:00:00")
dend <- hms("14:00:00")

d3 <- d2 %>%
  filter(cloudy=="N") %>% # try analysis with and without this line. Removing "cloudy" points removes 1198 observations
  filter(!(hms(localTime)<dstart), 
         !(hms(localTime)>dend) )


table(as.factor(d2$cloudy))
table(as.factor(d3$cloudy))

# Times look good?
# before filtering ...
summary(hms(d2$localTime))
# min is 955 am , max is 238 pm

# after filtering...
summary(hms(d3$localTime))
# min time is 11H 0M 4S
# max time is 13H 59M 50S
# looks good!

# Mine out the top 20 values for each site and depth
head(d3)

top20 <- d3 %>%
  ungroup() %>%
  group_by(siteName,depth) %>%
  top_n(n = 20, wt = licorData)

head(top20)
summary(top20)

# let's look at the top20 values 
licor_data <- top20

# plot each site type
licor_data %>% ggplot(aes(x = siteType, y = licorData)) +
  geom_boxplot(aes(x = siteType, y = licorData, fill = depth)) +
  geom_jitter(width = 0.4, size = 0.1) +
  theme_bw()

licor_data %>% ggplot(aes(x = depth, y = licorData)) +
  geom_boxplot(aes(x = depth, y = licorData, fill = siteType)) +
  # geom_jitter(width = 0.4, size = 0.1) +
  scale_fill_manual(values = c("salmon", "royalblue4")) +
  theme_bw()

licor_data %>% ggplot(aes(x = siteType, y = licorData)) +
  geom_boxplot(aes(x = siteType, y = licorData, fill = siteType)) +
  scale_fill_manual(values = c("salmon", "royalblue4")) +
  theme_bw()

# stats for inshore vs offshore -----
set.seed(1)
model1 <- MCMCglmm(licorData ~ siteType,
                   random = ~depth,
                   data = licor_data)
summary(model1)

set.seed(1)
model2 <- MCMCglmm(licorData ~ siteType,
                   random = ~depth+siteName,
                   data = licor_data)
summary(model2)


set.seed(1)
model3 <- MCMCglmm(licorData ~ siteType+depth,
                   random = ~siteName,
                   data = licor_data)
summary(model3)

set.seed(1)
model4 <- MCMCglmm(licorData ~ siteType+as.numeric(depth),
                   random = ~siteName,
                   data = licor_data)
summary(model4)

# plot mean ± standard deviation for each site and depth
sumstats <- licor_data %>%
  ungroup() %>%
  filter(cloudy=="N") %>%
  summarySE("licorData", c("siteName","depth"))
sumstats
str(sumstats)

ggplot(data = sumstats, 
       aes(x = depth,
           y = licorData, 
           ymin = licorData-sd, 
           ymax = licorData+sd,       
           fill = siteName)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(position = position_dodge(1), colour="black", width = 0.3, size = 0.5) +
  theme_bw()

ggplot(data = sumstats, 
       aes(x = siteName,
           y = licorData, 
           ymin = licorData-sd, 
           ymax = licorData+sd,       
           fill = depth)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(position = position_dodge(1), colour="black", width = 0.3, size = 0.5) +
  theme_bw()
