#Chi Square
#Package used
library(tidyverse)
#Do species composition differ between indoors and outdoors?
tab_loc <- df %>%
  group_by(Location, Species) %>%
  summarise(n = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) %>%
  column_to_rownames("Location") %>%
  as.matrix()

chisq.test(tab_loc)
#Actually used Fisher's test because some expected counts were small
fisher.test(tab_loc)

#Do species composition differ across towns?
tab_town <- df %>%
  group_by(Town, Species) %>%
  summarise(n = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) %>%
  column_to_rownames("Town") %>%
  as.matrix()

chisq.test(tab_town)
#Same here
fisher.test(tab_town)