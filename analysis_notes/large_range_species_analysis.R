suppressMessages(suppressWarnings({library(qs2); library(dplyr)}))
EX <- "exposureApp/inst/extdata"
rc <- qs_read(file.path(EX,"range_cells_mb.qs"))
rsz <- tibble(spName=names(rc), rangeSize=sapply(rc, length))
yrs <- c(1975, 1980, 1985, 1990)

# per species-year pixel footprint on the map (distinct cells across all vars)
peryr <- lapply(yrs, function(y){
  qs_read(file.path(EX, sprintf("allExpForShiny_by_year/allExpForShiny_%d.qs", y))) %>%
    filter(!is.na(cell)) %>% group_by(spName) %>%
    summarize(pixels = n_distinct(cell), propExp = first(propExposed), .groups="drop") %>%
    mutate(year=y)
}) %>% bind_rows() %>% left_join(rsz, by="spName")

# For each species across the explored years: mean pixels, mean fraction, n years on map
agg <- peryr %>% group_by(spName) %>%
  summarize(years_on_map = n_distinct(year),
            mean_pixels  = mean(pixels),
            mean_frac    = mean(propExp),
            rangeSize    = first(rangeSize), .groups="drop")

# Deliverable table: top 25 pixel contributors in the early series
tab <- agg %>% arrange(desc(mean_pixels)) %>% head(25) %>%
  transmute(species = gsub("_"," ",spName),
            range_size_cells = rangeSize,
            mean_pixels_on_map = round(mean_pixels),
            mean_frac_range_exposed = round(mean_frac,3),
            years_present_of_4 = years_on_map)
write.csv(tab, "/tmp/large_range_species_earlyyears.csv", row.names=FALSE)
cat("=== Top 25 pixel contributors, early series (1975,1980,1985,1990) ===\n")
print(as.data.frame(tab), row.names=FALSE)

# Summary stats for the narrative
cat("\n--- diagnostic summary ---\n")
cat("median range size, all species:", median(rsz$rangeSize), "cells\n")
cat("median range size, top-25 contributors:", median(tab$range_size_cells), "cells\n")
cat("all top-25 mean fraction exposed <= 0.03?:", all(tab$mean_frac_range_exposed<=0.03), "\n")
cat("mean fraction across top-25:", round(mean(tab$mean_frac_range_exposed),3), "(i.e. right at the 1-2% floor)\n")

# correlation across ALL species-years in these years: pixels vs range size
cc <- cor(peryr$pixels, peryr$rangeSize, use="complete.obs")
cat("correlation(pixels on map, range size) across all species-years:", round(cc,3), "\n")
# and pixels vs fraction exposed (should be weak/negative -> pixels driven by range not fraction)
cc2 <- cor(peryr$pixels, peryr$propExp, use="complete.obs")
cat("correlation(pixels on map, fraction exposed):", round(cc2,3), "\n")

saveRDS(peryr, "/tmp/peryr_final.rds")
