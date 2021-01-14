####
# Summarize by plant to make Table 1
####

require(tidyverse)
require(here)

# read in data
  av <- read_csv(file = here("Data", "aviaries_data.csv"))

  nectar <- read_csv(file = here("Data/nectar_experiments.csv"))

# Relabel one of the treatments
  nectar$Treatment[which(nectar$Pollen_added_after_visit == 1)] <- "NEHP"

# Reduce to necessary columns
  cols <- c("Species", "Plant", "Treatment")

  nectar2 <- nectar[,which(names(nectar)%in%cols)]
  av2 <- av[,which(names(av)%in%cols)]


# summarize by unique values
  
 sumry_av <- group_by(av2, Species, Treatment) %>%
   summarise(plants=n_distinct(Plant))

 sumry_nectar <- group_by(nectar2, Species, Treatment) %>%
   summarise(plants=n_distinct(Plant))











