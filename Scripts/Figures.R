# Set working directory manually first if needed
#Seasonal variation in vector abundance
#Package used
library(tidyverse)

#Enter monthly data
monthly <- tribble(
  ~Month,     ~Count,
  "January",   4,
  "February",  3,
  "March",     8,
  "April",    26,
  "May",      57,
  "June",     99
) %>%
  mutate(
    Month = factor(Month, levels = c("January","February","March","April","May","June")),
    Season = if_else(Month %in% c("January","February","March"), "Dry season", "Rainy season")
  )

#Compute key percentages
season_summary <- monthly %>%
  group_by(Season) %>%
  summarise(n = sum(Count), .groups = "drop") %>%
  mutate(percent = round(100 * n / sum(n), 1))

season_summary

#Time-series style
fig3_line <- ggplot(monthly, aes(x = Month, y = Count, group = 1)) +
  geom_line() +
  geom_point() +
  labs(
    x = NULL,
    y = "Number of Anopheles gambiae s.l.",
    title = "Temporal dynamics of Anopheles gambiae s.l. abundance (January–June 2024)"
  )

fig3_line

#To save it
ggsave("Figure4_monthly_abundance.tiff", fig3_line, width = 18, height = 12, units = "cm", dpi = 300)
#Or Manually with:
#Export on the Plot section
#Save as image
#Choose your width and height


#Knockdown curve
#Package used
library(tidyverse)
kd_long <- bind_rows(
  perm_kd %>% pivot_longer(starts_with("rep"), names_to = "replicate", values_to = "kd") %>%
    mutate(insecticide = "Permethrin 0.75%"),
  delt_kd %>% pivot_longer(starts_with("rep"), names_to = "replicate", values_to = "kd") %>%
    mutate(insecticide = "Deltamethrin 0.05%")
) %>%
  mutate(kd_pct = 100 * kd / 25)

kd_summary <- kd_long %>%
  group_by(insecticide, time) %>%
  summarise(mean_kd_pct = mean(kd_pct), .groups = "drop")

fig_kd <- ggplot(kd_summary, aes(x = time, y = mean_kd_pct, color = insecticide)) +
  geom_line(size = 1.2) +
  labs(x = "Exposure time (min)", y = "Mean knockdown (%)", color = "Insecticide") + 
  theme_minimal()

fig_kd 

ggsave("Figure_knockdown_curves.tiff", fig_kd, width = 16, height = 12, units = "cm", dpi = 300)
#Or Manually with:
#Export on the Plot section
#Save as image
#Choose your width and height

#Zooming + Styling of Figure 1
#Study site map (Nigeria + zoomed inlet)
# --- Packages ---
# Install if needed:
# install.packages(c("sf", "ggplot2", "dplyr", "rnaturalearth", "rnaturalearthdata",
#                    "ggrepel", "patchwork"))
# Nigeria states (ADM1)
#Packages used
library(sf)
library(ggplot2)
library(dplyr)
library(geodata)
library(ggrepel)
library(patchwork)

# Nigeria state boundaries
nga_adm1 <- geodata::gadm(
  country = "NGA",
  level = 1,
  path = tempdir()
)

nga_adm1 <- st_as_sf(nga_adm1)

# Extract Ogun State
ogun <- nga_adm1 %>%
  filter(NAME_1 == "Ogun")

#Nigeria outline
nga <- ne_countries(
  country = "Nigeria",
  returnclass = "sf"
)

#Study site points
sites <- tibble::tribble(
  ~Town,        ~lat,      ~lon,
  "Iju Ota",     6.686434,  3.139768,
  "Ado Odo Ota", 6.681760,  3.192634,
  "Sango Ota",   6.681602,  3.192814
)

sites <- sites %>%
  mutate(lon_plot = case_when(Town == "Ado Odo Ota" ~ lon + 0.008, Town == "Sango Ota" ~ lon - 0.008, TRUE ~ lon), lat_plot = case_when(Town == "Ado Odo Ota" ~ lat + 0.005, Town == "Sango Ota" ~ lat - 0.005, TRUE ~ lat))

sites_sf <- st_as_sf(
  sites,
  coords = c("lon_plot", "lat_plot"),
  crs = 4326,
  remove = FALSE
)

#Connector lines
ogun_centroid <- st_centroid(ogun)

lines_sf <- do.call(
  rbind,
  lapply(1:nrow(sites_sf), function(i) {
    st_sf(
      geometry = st_sfc(
        st_linestring(rbind(
          st_coordinates(ogun_centroid),
          st_coordinates(sites_sf[i, ])
        )),
        crs = 4326
      )
    )
  })
)

#Panel A - Nigeria with Ogun highlighted
p1 <- ggplot() +
  geom_sf(data = nga, fill = "grey90", color = "white") +
  geom_sf(data = ogun, fill = "red") +
  labs(title = "Nigeria showing Ogun State") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#Create zoom bounding box around sites
bb <- st_bbox(sites_sf)
bb_expanded <- bb + c(-0.02, -0.02, 0.02, 0.02)
ogun_zoom <- st_crop(ogun, bb_expanded)

#Panel B - Ogun + study sites
p2 <- ggplot() +
  geom_sf(data = ogun, fill = "grey95", color = "black") +
  geom_sf(data = lines_sf, linetype = "dashed", color = "darkred") +
  geom_sf(data = sites_sf, aes(color = Town, shape = Town), size = 5, stroke = 1.5) +
  ggrepel::geom_text_repel(
    data = sites,
    aes(x = lon_plot, y = lat_plot, label = Town),
    size = 4, 
    box.padding = 0.6,
    segment.color = "black"
  ) +
  scale_color_manual(values = c("red", "blue", "darkgreen")) +
  labs(title = "Mosquito sampling sites in Ado-Odo Ota LGA") +
  theme_minimal() +theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

#Combine panels
final_map <- p1 + p2 + plot_layout(widths = c(1, 1.2))
final_map
#Saved Manually with:
#Export on the Plot section
#Save as image
#Choose your width and height
