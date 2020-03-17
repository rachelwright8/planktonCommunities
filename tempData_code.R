setwd("~/Dropbox/andrea_plankton/clean_for_upload/light_temp/")

# Load libraries
library(tidyverse)
library(stringr)

# Load data
temp_data <- read.csv("minmeanmax_tempdata.csv")

# Reformate date
temp_data <- temp_data %>% separate(date, c("month", "day", "year"), "/")

head(temp_data)
summary(temp_data)

# Reformate site_location
temp_data$site_location <- as.factor(ifelse(temp_data$site_location=="IR", "inshore", "offshore"))

# Add site names
temp_data <- temp_data %>% 
  mutate(site_name = ifelse(site_number=="IR1", "PuntaDonato", 
                            ifelse(site_number == "IR2", "STRIPoint", 
                                   ifelse(site_number== "IR3", "Cristobal", 
                                          ifelse(site_number=="OR3", "PopaIsland", "DragoMar")))))

temp_data$site_name <- as.factor(temp_data$site_name)

summary(temp_data)
  
# Subset data for just June 2015
sub <- temp_data %>% filter(month=="6" & year=="15")
head(sub)

# Making colors and shapes
summary(sub)
ggplot_shapescolors <- sub %>% 
  mutate(colors = ifelse(site_location=="inshore", "inshore", "offshore"), 
         shapes = ifelse(site_name=="Cristobal", 15, 
                         ifelse(site_name=="DragoMar", 1, 
                                ifelse(site_name=="PopaIsland", 17, 
                                       ifelse(site_name=="PuntaDonato", 19, 5))))) %>% 
  select(colors, shapes)

str(ggplot_shapescolors$shapes)

# Plot temp over time - different color for each site
ggplot(sub, aes(x = day, y = Max, group = site_name, color = site_name)) +
  geom_line(size=1) +
  geom_point(aes(shape=factor(ggplot_shapescolors$shapes)), size = 3) +
  labs(x = "Day in June", y = "Max Temp °C")+
  geom_line()+
  geom_point()+
  theme_bw()

# Plot temp over time - color by inshore/offshore, symbol by site
ggplot(sub, aes(x = day, y = Mean, group = site_name, color = sub$site_location, shape = sub$site_name)) +
  geom_line(size=1) +
  labs(x = "Day in June 2015", y = "Max Temp °C") +
  geom_ribbon(aes(ymin=Max, ymax=Min, fill= sub$site_location), alpha = 0.5, colour=NA)+
  scale_color_manual(values=c("salmon", "royalblue4")) +
  scale_fill_manual(values=c("salmon", "royalblue4")) +
  geom_point(size = 5) +
  scale_shape_manual(values=c(15,1,5,19,17))+
  theme_bw()
