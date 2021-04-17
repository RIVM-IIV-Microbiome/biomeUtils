
test_that("calculatePD works", {
  data("FuentesIliGutData")
  ps1 <- subset_samples(FuentesIliGutData, ILI == "C")
  ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
  meta_tib <- calculatePD(ps1, justDF = TRUE)
  pd_vals <- round(meta_tib$PD[1:3])
  exp_vals <- round(c(967.3003, 1035.1305, 1189.0248))
  expect_equal(pd_vals, exp_vals)

  sam_order <- sample_names(ps1)[1:6]
  new_order <- rownames(meta_tib)[1:6]
  expect_equal(sam_order, new_order)
})
