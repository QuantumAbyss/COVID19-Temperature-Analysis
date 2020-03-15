#'
#' Analysis of COVID-19 transmission rate and its correlation
#' with regional temperature
#'
#' COVID Data: https://coronavirus.jhu.edu/map.html
#' Weather Data: http://shinyapps.admenergy.com/app/getNOAA
#'

### INITIALIZATION --------------------------------------------------------------------------------
.x <- c("data.table", "dplyr", "lubridate", "ggplot2", "broom")
lapply(.x, library, character.only=T)

setwd("~/1) Projects/COVID19-Temperature-Analysis/")

### DATA CLEAN ------------------------------------------------------------------------------------

### WEATHER -----------------
weather.us <- lapply(list.files("./WeatherData/", full.names=T, pattern="2020.csv"), 
                     function(x) {fread(x) %>% mutate(Region=gsub("-.*|.*\\/", "", x))}) %>%
  rbindlist() %>%
  mutate(timestamp = ymd_hms(date_local)) %>%
  select(Region, timestamp, temp, rel_hum) %>%
  mutate(Date=date(timestamp)) %>%
  group_by(Region, Date) %>% 
  summarise(min.temp=min(temp),
            mean.temp=mean(temp),
            max.temp=max(temp),
            min.relhum=100*min(rel_hum),
            mean.relhum=100*mean(rel_hum),
            max.relhum=100*max(rel_hum))

weather.wu <- fread("./WeatherData/weatherundergroundDaily.csv") %>%
  mutate(Date=mdy(Date))

weather <- bind_rows(weather.us, weather.wu)

### TRACKER -----------------
covid <- lapply(list.files("./csse_covid_19_daily_reports/", full.names=T), fread) %>% 
  rbindlist(fill=T)

table(covid$`Province/State`)

covid <- covid %>%
  mutate(Region = ifelse(grepl(", ", `Province/State`), 
                         gsub(".*, ", "", `Province/State`), 
                         `Province/State`)) %>%
  mutate(Region = ifelse(`Country/Region` != "US", `Country/Region`, Region)) %>% 
  merge(data.table(abb=state.abb,
                   StateName=state.name), by.x='Region', by.y='abb',all.x=T) %>%
  mutate(Region = ifelse(!is.na(StateName), StateName, Region)) %>%
  mutate(timestamp = ymd_hms(`Last Update`)) %>%
  mutate(timestamp = ymd_hms(ifelse(is.na(timestamp), as.character(mdy_hm(`Last Update`)), as.character(timestamp)))) %>%
  filter(!(StateName %in% c("Diamond Princess", "Unassigned Location (From Diamond Princess"))) %>%
  mutate(Region = ifelse(Region=="Iran (Islamic Republic of)", "Iran", 
                         ifelse(Region=="Korea, South", "South Korea", 
                                ifelse(Region=="Republic of Korea", "South Korea", Region)))) %>%
  arrange(Region, timestamp) %>%
  filter(Region %in% unique(weather$Region),
         !is.na(Confirmed)) %>%
  unique() %>%
  select(Region, `Province/State`, timestamp, Confirmed, Deaths)

# covid %>% 
#   group_by(Region) %>%
#   summarise(CurrentCases = max(Confirmed)) %>%
#   arrange(desc(CurrentCases))

covid <- covid %>%
  split(.$Region) %>%
  lapply(function(x) {
    timestamp.set <- sort(unique(x$timestamp))
    province.set <- unique(x$`Province/State`)
    if (length(province.set)==1) {return(x)}
    complete <- lapply(1:length(timestamp.set), function(y) {
      curr <- filter(x, timestamp==timestamp.set[y])
      data.table(Region= curr$Region[1],
                 `Province/State` = province.set[which(!(province.set %in% curr$`Province/State`))],
                 timestamp=timestamp.set[y],
                 Confirmed = lapply(province.set[which(!(province.set %in% curr$`Province/State`))], 
                                    function(z) {
                                      prev <- filter(x, timestamp < timestamp.set[y], `Province/State`==z)
                                      if (nrow(prev)==0) {
                                        return(0)
                                      } else {
                                        return(prev$Confirmed[which.max(prev$timestamp)])
                                      }
                                    }) %>% unlist(),
                 Deaths = lapply(province.set[which(!(province.set %in% curr$`Province/State`))], 
                                 function(z) {
                                   prev <- filter(x, timestamp < timestamp.set[y], `Province/State`==z)
                                   if (nrow(prev)==0) {
                                     return(0)
                                   } else {
                                     return(prev$Deaths[which.max(prev$timestamp)])
                                   }
                                 }) %>% unlist()) %>%
        bind_rows(curr) %>%
        return()
    }) %>% rbindlist()
  }) %>% rbindlist()

covid <- covid %>%
  mutate(Date=date(timestamp)) %>%
  group_by(timestamp, Date, Region) %>%
  summarise(Confirmed=sum(Confirmed,na.rm=T),
            Deaths=sum(Deaths,na.rm=T)) %>%
  arrange(Region, timestamp) %>%
  filter(!is.na(Confirmed))

full <- merge(covid, weather, by=c("Region", "Date")) %>%
  mutate(Day = yday(Date)) %>%
  arrange(Region, timestamp) %>%
  mutate(timenum=as.numeric(timestamp))

popdens <- fread("./PopulationDensities.csv")

### ANALYSIS --------------------------------------------------------------------------------------

result <- full %>%
  split(.$Region) %>%
  lapply(function(x) {
    
    model <- nls(Confirmed ~ (1+r)^(Day-s), x, start=list(r=.2,s=40))
    x$pred <- predict(model, x)
    cor(x$Confirmed, x$pred)^2
    
    # ggplot(x) + geom_point(mapping=aes(x=Day,y=Confirmed,color="Reported"),size=2) + geom_line(mapping=aes(x=Day,y=Confirmed,color="Reported"),size=1.5) +
    #   geom_point(mapping=aes(x=Day,y=pred,color="Modeled"),size=2) + geom_line(mapping=aes(x=Day,y=pred,color="Modeled"),size=1.5) +
    #   xlab("Day of the Year") + ylab("Confirmed Cases") + ggtitle("California COVID-19 Cases Over Time") + 
    #   theme_light() + theme(legend.position="top") 

    return(data.table(Region=x$Region[1],
                      GrowthRate=tidy(model)$estimate[1],
                      t(colMeans(x[,6:11]))))
  }) %>% rbindlist() %>%
  merge(popdens, by='Region', all.x=T)

ggplot(result, mapping=aes(x=PeoplePerSqkm, y=GrowthRate)) + geom_label(aes(label=Region)) +
   xlab("Population Density (Person/km^2)") + ylab("Viral Growth Rate") +
  scale_x_continuous(breaks=seq(0,26000,2000)) + 
  scale_y_continuous(breaks=seq(0,.9,.2)) + theme_light() + 
  ggtitle("COVID-19 Viral Growth Rate versus Population Density Across Geographic Region")

result <- result %>%
  mutate(Group=as.factor(ifelse(Region %in% c("New York", "France", "Greece", "Iran", "Japan"), "G2", "G1")))

finalreg <- lm(GrowthRate ~ mean.temp + mean.relhum + PeoplePerSqkm + Group, data=result)
tidy(finalreg)
glance(finalreg)

ggplot(result, mapping=aes(x=mean.relhum,y=GrowthRate)) + geom_label(aes(label=Region)) +
  xlab("Average Relative Humidity") + ylab("Viral Growth Rate") +
  scale_x_continuous(breaks=seq(0,100,10)) + 
  scale_y_continuous(breaks=seq(0,.9,.2)) + theme_light() + 
  ggtitle("COVID-19 Viral Growth Rate versus Average Relative Humidity Across Region")

ggplot(result, mapping=aes(x=mean.temp,y=GrowthRate)) + geom_label(aes(label=Region)) +
  xlab("Average Temperature (Fahrenheit)") + ylab("Viral Growth Rate") +
  scale_x_continuous(breaks=seq(0,100,10)) + 
  scale_y_continuous(breaks=seq(0,.9,.2)) + theme_light() + 
  ggtitle("COVID-19 Viral Growth Rate versus Average Temperature Across Region")































