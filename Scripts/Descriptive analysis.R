# Data cleaning and creation of tables
# Set working directory manually first if needed
#Package used
library(tidyverse)

df <- tribble(
  ~Town,        ~Location,  ~Species,                  ~Count,
  "Iju",        "Outdoors", "Anopheles gambiae s.l.",  23,
  "Iju",        "Outdoors", "Culex spp.",              87,
  "Iju",        "Outdoors", "Aedes spp.",               0,
  "Iju",        "Indoors",  "Anopheles gambiae s.l.",  74,
  "Iju",        "Indoors",  "Culex spp.",              58,
  "Iju",        "Indoors",  "Aedes spp.",               0,
  
  "Ado Odo Ota","Outdoors", "Anopheles gambiae s.l.",  34,
  "Ado Odo Ota","Outdoors", "Culex spp.",              11,
  "Ado Odo Ota","Outdoors", "Aedes spp.",               3,
  "Ado Odo Ota","Indoors",  "Anopheles gambiae s.l.",  42,
  "Ado Odo Ota","Indoors",  "Culex spp.",              43,
  "Ado Odo Ota","Indoors",  "Aedes spp.",               1,
  
  "Sango Ota",  "Outdoors", "Anopheles gambiae s.l.",   9,
  "Sango Ota",  "Outdoors", "Culex spp.",              96,
  "Sango Ota",  "Outdoors", "Aedes spp.",               2,
  "Sango Ota",  "Indoors",  "Anopheles gambiae s.l.",  15,
  "Sango Ota",  "Indoors",  "Culex spp.",              81,
  "Sango Ota",  "Indoors",  "Aedes spp.",               0
)

#Table creation
table1 <- df %>%
  arrange(Town, Location, Species) %>%
  pivot_wider(names_from = Species, values_from = Count, values_fill = 0) %>%
  mutate(Total = `Anopheles gambiae s.l.` + `Culex spp.` + `Aedes spp.`)

table1

write.csv(table1, "Table1_spatial_distribution.csv", row.names = FALSE)

#Quick numbers you can cite
indoor_outdoor <- df %>%
  group_by(Location, Species) %>%
  summarise(n = sum(Count), .groups = "drop") %>%
  group_by(Location) %>%
  mutate(percent = round(100 * n / sum(n), 1))

indoor_outdoor

#Stacked by town
by_town <- df %>%
  group_by(Town, Location, Species) %>%
  summarise(n = sum(Count), .groups = "drop")


#Overall proportion table
overall <- df %>%
  group_by(Species) %>%
  summarise(n = sum(Count), .groups = "drop") %>%
  mutate(
    total = sum(n),
    percent = round(100 * n / total, 1)
  ) %>%
  arrange(desc(n))

overall


#Insecticide susceptibility test
library(tidyverse)

times <- c(10,15,20,30,40,50,60)  # knockdown minutes

# PERMETHRIN (0.75%) - knockdown counts at each time (per replicate)
perm_kd <- tribble(
  ~time, ~rep1, ~rep2, ~rep3, ~rep4,
  10,    0,     0,     0,     0,
  15,    0,     0,     0,     0,
  20,    0,     0,     0,     0,
  30,    0,     0,     0,     1,
  40,    0,     0,     0,     1,
  50,    0,     0,     3,     1,
  60,    0,     0,     3,     2
)

perm_24h_dead <- c(0, 1, 15, 6)   # deaths at 24h per replicate (out of 25)
perm_control_dead <- c(0, 0)      # control deaths (out of 25) if you used 2 controls

# DELTAMETHRIN (0.05%)
delt_kd <- tribble(
  ~time, ~rep1, ~rep2, ~rep3, ~rep4,
  10,    0,     0,     0,     0,
  15,    0,     0,     0,     0,
  20,    3,     3,     1,     3,
  30,    8,     5,     5,     5,
  40,    14,    11,    12,    15,
  50,    17,    18,    21,    19,
  60,    23,    25,    23,    24
)

delt_24h_dead <- c(25, 25, 24, 24)
delt_control_dead <- c(0, 0)

#Mortality% + WHO classification table
mortality_summary <- tibble(
  insecticide = c("Permethrin 0.75%", "Deltamethrin 0.05%"),
  tested = c(100, 100),
  dead_24h = c(sum(perm_24h_dead), sum(delt_24h_dead)),
  control_dead = c(sum(perm_control_dead), sum(delt_control_dead)),
  control_tested = c(length(perm_control_dead)*25, length(delt_control_dead)*25)
) %>%
  mutate(
    mortality_pct = round(100 * dead_24h / tested, 1),
    control_mortality_pct = round(100 * control_dead / control_tested, 1),
    WHO_status = case_when(
      mortality_pct >= 98 ~ "Susceptible",
      mortality_pct >= 90 ~ "Possible resistance",
      TRUE ~ "Resistant"
    )
  )

mortality_summary


#To calculate allele frequencies + genotype frequencies
kdr <- tribble(
  ~mutation, ~n, ~rr, ~Rr, ~RR,
  "L1014F (kdr-west)", 15, 0, 6, 9,
  "L1014S (kdr-east)", 15, 0, 7, 8
)

kdr_freq <- kdr %>%
  mutate(
    p_R = (2*RR + Rr) / (2*n),
    q_r = 1 - p_R,
    rr_freq = rr/n,
    Rr_freq = Rr/n,
    RR_freq = RR/n
  )

kdr_freq

#Hardy-Weinberg equilibrum test
hw_exact <- function(obs_hom1, obs_hets, obs_hom2) {
  # Exact HWE test (Wigginton et al., 2005 style)
  # hom1 = RR, hets = Rr, hom2 = rr (or vice versa; must be consistent)
  if (any(c(obs_hom1, obs_hets, obs_hom2) < 0)) stop("Counts must be >= 0")
  n <- obs_hom1 + obs_hets + obs_hom2
  if (n == 0) return(NA_real_)
  
  obs_homc <- min(obs_hom1, obs_hom2)
  obs_homr <- max(obs_hom1, obs_hom2)
  rare_copies <- 2*obs_homc + obs_hets
  
  # distribution over possible heterozygote counts conditional on rare_copies
  probs <- rep(0, rare_copies + 1)
  
  mid <- floor(rare_copies * (2*n - rare_copies) / (2*n))
  if ((mid %% 2) != (rare_copies %% 2)) mid <- mid + 1
  
  probs[mid + 1] <- 1
  sum_probs <- 1
  
  # downward recursion
  for (h in seq(mid, 2, by = -2)) {
    probs[h - 1 + 1] <- probs[h + 1] * h * (h - 1) / (4 * ( (rare_copies - h)/2 + 1 ) * ( (2*n - rare_copies - h)/2 + 1 ))
    sum_probs <- sum_probs + probs[h - 1 + 1]
  }
  
  # upward recursion
  for (h in seq(mid, rare_copies - 2, by = 2)) {
    probs[h + 2 + 1] <- probs[h + 1] * 4 * ( (rare_copies - h)/2 ) * ( (2*n - rare_copies - h)/2 ) / ( (h + 2) * (h + 1) )
    sum_probs <- sum_probs + probs[h + 2 + 1]
  }
  
  probs <- probs / sum_probs
  
  # observed het count probability under null
  p_obs <- probs[obs_hets + 1]
  
  # two-sided p-value: sum probabilities <= p_obs
  pval <- sum(probs[probs <= p_obs])
  return(min(1, pval))
}

kdr_hwe <- kdr %>%
  rowwise() %>%
  mutate(
    HWE_p_exact = hw_exact(RR, Rr, rr)
  ) %>%
  ungroup()

kdr_hwe

#Final kdr table
final_kdr_table <- kdr_hwe %>%
  mutate(
    p_R = (2*RR + Rr) / (2*n),
    q_r = 1 - p_R,
    p_R = round(p_R, 3),
    q_r = round(q_r, 3),
    HWE_p_exact = signif(HWE_p_exact, 3)
  ) %>%
  select(mutation, n, rr, Rr, RR, p_R, q_r, HWE_p_exact)
#0R
final_kdr_table <- kdr_hwe %>%
  mutate(
    p_R = (2*RR + Rr) / (2*n),
    q_r = 1 - p_R,
    p_R = round(p_R, 3),
    q_r = round(q_r, 3),
    HWE_p_exact = signif(HWE_p_exact, 3)
  ) %>%
  dplyr::select(mutation, n, rr, Rr, RR, p_R, q_r, HWE_p_exact)

final_kdr_table
