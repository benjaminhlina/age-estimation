library(dplyr)
library(here)
library(janitor)
library(lubridate)
library(readr)


# bring in data -----

lt <- read_csv(here("Data", 
                    "IPE_SPE.csv")) %>% 
  clean_names()


meta_data <- read_csv(here("Data", 
                           "IPE_STA.csv")) %>% 
  clean_names()


# look at the structure ----

glimpse(lt)


glimpse(meta_data) 




meta_data <- meta_data %>% 
  mutate(date_pose = mdy(date_pose), 
         date_levee = mdy(date_levee))



lt <- lt %>% 
  left_join(meta_data, by = c("ue", "no_inv", "no_peche", 
                              "no_pose", "commentair", "ind_action"))


glimpse(lt)




openxlsx::write.xlsx(here("Results",
                          "qc_collection_data_lineup.xlsx"), x = lt,  
                     )
