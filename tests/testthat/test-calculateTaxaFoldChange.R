test_that("calculateTaxaFoldDifference works", {
  data("FuentesIliGutData")
  ps1 <- subset_samples(FuentesIliGutData, ILI != "L2")
  ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
  tb <- calculateTaxaFoldDifference(ps1, group = "ILI")

  asvids <- tb$FeatureID[1:5]
  exptids <- c("ASV302", "ASV636", "ASV500", "ASV7", "ASV2617")
  expect_equal(asvids, exptids)

  colsid <- colnames(tb)
  expcols <- c(
    rank_names(ps1), "FoldDifference", "FeatureID",
    "Prevalence.C", "Prevalence.L1", "Enriched"
  )

  expect_equal(colsid, expcols)

  pd_vals <- round(tb$Prevalence.C[1:5] * 100)
  exp_vals <- c(85, 22, 30, 99, 23)
  expect_equal(pd_vals, exp_vals)

  pd_vals2 <- round(tb$Prevalence.L1[1:5] * 100)
  exp_vals2 <- c(83, 18, 42, 96, 18)
  expect_equal(pd_vals2, exp_vals2)

  fold1 <- as.character(tb$FoldDifference[1:3])
  expfold1 <- c("-0.167551358688032", "-0.100108612147008", "0.157777020059819")
  expect_equal(fold1, expfold1)
})
