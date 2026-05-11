#Knockdown Time for KDT 50 and KDT95 for Permethrin and Deltamethrin
#Package used
library(tidyverse)
library(tidyverse)
library(MASS)

#KDT kept crashing since it was based on extrapolation for permethrin
fit_kdt_safe <- function(dat){
  
  out <- tryCatch({
    m <- glm(cbind(kd, alive) ~ time,
             family = binomial(link = "probit"),
             data = dat)
    
    kdt <- dose.p(m, p = c(0.5, 0.95))
    se  <- attr(kdt, "SE")
    
    tibble(
      KDT = c("KDT50", "KDT95"),
      estimate_min = as.numeric(kdt),
      se = as.numeric(se),
      lower_95 = estimate_min - 1.96 * se,
      upper_95 = estimate_min + 1.96 * se,
      note = NA_character_
    )
  }, error = function(e){
    tibble(
      KDT = c("KDT50", "KDT95"),
      estimate_min = NA_real_,
      se = NA_real_,
      lower_95 = NA_real_,
      upper_95 = NA_real_,
      note = paste0("Not estimable: ", e$message)
    )
  })
  
  out
}

kdt_table <- kd_long %>%
  group_by(insecticide) %>%
  group_modify(~ fit_kdt_safe(.x)) %>%
  ungroup() %>%
  mutate(across(c(estimate_min, lower_95, upper_95), ~ round(.x, 1)))

kdt_table