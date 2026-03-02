
# must install FRED package
#library(devtools)
#install_github("cykbennie/fbi")

library(fbi)
library(tidyverse)
library(TTR)

start_date = "2003-01-01"
end_date = "2025-01-01"

data = fredmd("ORIGINAL_MEDEIROS/data/2025-01.csv")
data_raw = fredmd("ORIGINAL_MEDEIROS/data/2025-01.csv",
                  transform = FALSE)
varlist = fredmd_description
vars = intersect(colnames(data),varlist$fred)

data = data %>% as_tibble()%>%
  select(all_of(c("date",vars)))
varlist = varlist%>%filter(fred%in%vars)
prices_varlist = varlist%>%filter(group=="Prices",tcode==6)

data = data%>% as_tibble()%>%
  select( -all_of(prices_varlist$fred) )

prices = data_raw %>% as_tibble() %>%
  select(all_of(prices_varlist$fred))%>%
  mutate_all(.funs = function(x)100*c(NA,x%>%log()%>%diff()))

data = cbind(data%>%as.data.frame(),prices%>%as.data.frame())

data = data %>%
  filter(date>=start_date)%>%
  select_if(~ !any(is.na(.)))

data = subset(data, date < end_date)

save(data,file = "ORIGINAL_MEDEIROS/data/data.rda")

