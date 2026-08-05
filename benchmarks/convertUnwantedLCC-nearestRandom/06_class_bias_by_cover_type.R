## What the removed lowest-class tie-break cost, in cover-type terms.
##
## The deterministic rule briefly present in 1.2.0.9004 resolved every distance tie to the
## lowest class code. The Canada LCC class codes are ordered roughly non-vegetated ->
## non-forest vegetation -> forest, and within forest coniferous < broadleaf < mixedwood, so
## "lowest code wins" is not a neutral rule: it systematically favours sparse/non-forest
## cover over forest, and coniferous over broadleaf and mixedwood. That is why it was
## removed.
##
## This joins the per-landscape assigned-class compositions from 03_method_comparison.R to
## the cover-type labels and reports, per cover type, how far that removed rule and the two
## surviving abundance-weighted methods land from the old algorithm.
##
## Run (after 03_method_comparison.R):
##   LANDR_SRC=. Rscript benchmarks/convertUnwantedLCC-nearestRandom/06_class_bias_by_cover_type.R

suppressPackageStartupMessages(library(data.table))
OUT <- file.path(Sys.getenv("LANDR_SRC", "."), "benchmarks", "convertUnwantedLCC-nearestRandom")

## Canada LCC class codes, per LandR::prepInputs_NTEMS_LCC_FAO() (R/prepInputs_NTEMS.R).
## NB: the SCANFI "CanadaLCCclassCodes" product used for these landscapes carries a single
## code 30 where the NTEMS scheme splits 31 snow_ice / 32 rock_rubble / 33 exposed_barren;
## 30 is excluded from the available set here and never appears as an assignment, so it is
## left unlabelled rather than guessed at.
lccLabels <- data.table(
  class = c(0L, 20L, 31L, 32L, 33L, 40L, 50L, 80L, 81L, 100L, 210L, 220L, 230L, 240L),
  label = c(
    "no change", "water", "snow/ice", "rock/rubble", "exposed barren land",
    "bryoids", "shrubs", "wetland", "wetland-treed", "herbs",
    "coniferous", "broadleaf", "mixedwood", "disturbed"
  ),
  group = c(
    "non-vegetated", "non-vegetated", "non-vegetated", "non-vegetated", "non-vegetated",
    "non-forest veg.", "non-forest veg.", "wetland", "forest", "non-forest veg.",
    "forest", "forest", "forest", "disturbed"
  )
)

comp <- fread(file.path(OUT, "assigned_class_composition.csv"))
tot <- fread(file.path(OUT, "real_landscapes_methods.csv"))[, .(landscape, unwanted)]

## weight each landscape by its number of unwanted pixels, so the pooled figure is
## "share of all unwanted pixels across the four landscapes"
comp <- merge(comp, tot, by = "landscape")
pooled <- comp[, lapply(.SD, function(p) sum(p * unwanted) / sum(unwanted)),
  by = "class", .SDcols = c("old", "lowestCode", "nearestWeighted", "nearestRandom")
]
pooled <- merge(pooled, lccLabels, by = "class", all.x = TRUE)
setorderv(pooled, "class")

pooled[, `:=`(
  lowestCode_vs_old = round(lowestCode / old, 2),
  nearestWeighted_vs_old = round(nearestWeighted / old, 2),
  nearestRandom_vs_old = round(nearestRandom / old, 2),
  pp_shift_lowestCode = round(100 * (lowestCode - old), 2)
)]
pooled[, `:=`(old = round(100 * old, 2), lowestCode = round(100 * lowestCode, 2),
              nearestWeighted = round(100 * nearestWeighted, 2),
              nearestRandom = round(100 * nearestRandom, 2))]

setcolorder(pooled, c("class", "label", "group", "old", "lowestCode", "nearestWeighted", "nearestRandom",
                      "lowestCode_vs_old", "nearestWeighted_vs_old", "nearestRandom_vs_old", "pp_shift_lowestCode"))
print(as.data.frame(pooled), row.names = FALSE)
fwrite(pooled, file.path(OUT, "class_bias_by_cover_type.csv"))

## forest vs non-forest roll-up: the summary that matters for a succession model
roll <- pooled[!is.na(group), .(
  old = sum(old), lowestCode = sum(lowestCode),
  nearestWeighted = sum(nearestWeighted), nearestRandom = sum(nearestRandom)
), by = "group"]
roll[, `:=`(lowestCode_vs_old = round(lowestCode / old, 2),
            nearestWeighted_vs_old = round(nearestWeighted / old, 2),
            nearestRandom_vs_old = round(nearestRandom / old, 2))]
cat("\n-- rolled up by cover group (% of unwanted pixels) --\n")
print(as.data.frame(roll), row.names = FALSE)
fwrite(roll, file.path(OUT, "class_bias_rollup.csv"))

png(file.path(OUT, "fig3_class_bias.png"), width = 1150, height = 520, res = 115)
par(mar = c(7, 4.5, 3, 1))
m <- t(as.matrix(pooled[, .(old, lowestCode, nearestWeighted, nearestRandom)]))
colnames(m) <- sprintf("%d\n%s", pooled$class, pooled$label)
barplot(m,
  beside = TRUE, col = c("grey35", "firebrick", "steelblue", "darkgreen"), las = 2,
  ylab = "% of unwanted pixels assigned", cex.names = 0.75,
  main = "Assigned cover type: the removed lowest-code rule vs abundance weighting"
)
legend("topright",
  c("spiral (old)", "lowest-code (REMOVED)", "nearestWeighted", "nearestRandom"),
  fill = c("grey35", "firebrick", "steelblue", "darkgreen"), bty = "n"
)
dev.off()
cat("wrote fig3_class_bias.png\n")
