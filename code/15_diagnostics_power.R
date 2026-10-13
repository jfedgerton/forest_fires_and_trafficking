################################################################################
# 15_diagnostics_power.R
#
# Two power diagnostics, written after the lag and intensity results came back
# null and the obvious question was whether those are precise nulls or wide ones.
#
#   Q1  LAG MODELS. Reports the 95% interval of the contemporaneous coefficient
#       when lags are entered alongside it, whether the contemporaneous-only
#       point estimate falls inside that interval, the standard-error inflation,
#       and the minimum detectable effect at 80% power (2.802 x SE for a
#       two-sided 5% test). If the MDE exceeds the contemporaneous estimate, the
#       lag specification could not have detected an effect of the size we
#       estimate even if it were exactly true - which is a statement about power,
#       not about the effect.
#
#   Q2  INTENSITY MARGIN. Reports how far the sample falls when conditioning on
#       y > 0, how many municipalities are dropped as fixed-effect singletons,
#       whether treatment variation survives the conditioning, and the minimum
#       detectable effect. The concern is that "fire does not change how many"
#       may be an artefact of a 300-fold sample reduction rather than a finding.
#
# Reads the estimated grid; estimates nothing new.
# Seed 123.
################################################################################

source("00_config.R")
g <- fread(file.path(DIR_RES,"spec_grid_2sls.csv"))[status=="ok" & window=="2002-2021"]
d <- load_panel()[year <= 2021L]
sdf <- sd(d[[FIRE]],na.rm=TRUE)

cat("\n========== Q4: ARE THE LAG NULLS PRECISE OR WIDE? ==========\n")
z <- g[zset=="Z2_max" & kset=="K4_ag10_gdp" & transformation %in% c("binary","count","log1p")]
z[, `:=`(lo=iv_b-1.96*iv_se, hi=iv_b+1.96*iv_se)]
# minimum detectable effect at 80% power, 5% two-sided: 2.8 * SE
z[, mde := 2.802*iv_se]
L0 <- z[lagid=="L0", .(source,transformation,b0=iv_b,se0=iv_se)]
cmp <- merge(z[lagid=="LALL"], L0, by=c("source","transformation"))
cmp[, b0_in_CI := b0 >= lo & b0 <= hi]
cmp[, se_ratio := iv_se/se0]
print(cmp[,.(source,transformation,
             contemp_b=round(b0,5), lag_b=round(iv_b,5),
             lag_se=round(iv_se,5), se_inflation=round(se_ratio,2),
             CI=sprintf("[%.5f, %.5f]",lo,hi),
             contemp_inside_lagCI=b0_in_CI,
             MDE=round(mde,5), MDE_vs_contemp=round(mde/b0,2))])

cat("\n-- same for the lag-only specs (L1, L2) --\n")
cmp2 <- merge(z[lagid %in% c("L1","L2")], L0, by=c("source","transformation"))
cmp2[, b0_in_CI := b0 >= lo & b0 <= hi]
print(cmp2[,.(source,transformation,lagid,lag_b=round(iv_b,5),lag_se=round(iv_se,5),
              CI=sprintf("[%.5f, %.5f]",lo,hi), contemp_inside_lagCI=b0_in_CI)])

cat("\n========== Q3: THE INTENSITY MEASURE ==========\n")
for (dvn in names(DVS)) {
  yv <- DVS[[dvn]]
  pos <- d[get(yv)>0]
  cat(sprintf("\n%s: %d of %d muni-quarters have y>0 (%.2f%%); %d municipalities\n",
              dvn, nrow(pos), nrow(d), 100*nrow(pos)/nrow(d), uniqueN(pos$code_muni_7)))
  cat(sprintf("   y|y>0: mean %.2f  median %.0f  max %.0f  sd %.2f\n",
              mean(pos[[yv]]), median(pos[[yv]]), max(pos[[yv]]), sd(pos[[yv]])))
  cat(sprintf("   log1p(y)|y>0: mean %.3f  sd %.3f\n",
              mean(log1p(pos[[yv]])), sd(log1p(pos[[yv]]))))
  cat(sprintf("   fire in that subsample: mean %.3f  sd %.3f  (full panel sd %.3f)\n",
              mean(pos[[FIRE]]), sd(pos[[FIRE]]), sdf))
  # how many municipalities contribute >1 positive quarter -> identify within-muni FE
  ct <- pos[, .N, by=code_muni_7]
  cat(sprintf("   municipalities with >1 positive quarter: %d of %d (%.0f%%)\n",
              sum(ct$N>1), nrow(ct), 100*mean(ct$N>1)))
  cat(sprintf("   singleton municipalities dropped by muni FE: %d\n", sum(ct$N==1)))
}
i <- g[transformation=="intensity" & lagid=="L0" & zset=="Z2_max" & kset=="K4_ag10_gdp"]
cat("\nIntensity IV estimates (contemporaneous, vpd_max, 6 controls):\n")
print(i[,.(source,n_obs,iv_b=round(iv_b,5),iv_se=round(iv_se,5),iv_p=round(iv_p,3),
           fs_F=round(fs_F,1), MDE=round(2.802*iv_se,5))])
cat("\nFor comparison, the binary MDE:\n")
print(g[transformation=="binary" & lagid=="L0" & zset=="Z2_max" & kset=="K4_ag10_gdp",
        .(source,n_obs,iv_b=round(iv_b,5),MDE=round(2.802*iv_se,5))])
