test_that("uniqueness works", {
  data("FuentesIliGutData")
  ps <- filterSampleData(FuentesIliGutData, ILI != "L2")
  ps <- getProportions(ps)
  dist.mat <- phyloseq::distance(ps, "bray")
  # without ref samples
  uniqueness.df <- uniqueness(ps,
                              dist_mat = dist.mat,
                              reference_samples = NULL)
  col.exp <- 63
  expect_equal(ncol(uniqueness.df), col.exp)

  exp.uniq <- c(0.4275774, 0.4499553, 0.5018299, 0.4744228, 0.5967708)
  expect_equal(uniqueness.df$uniqueness[1:5], exp.uniq, tolerance=1.0e-3)

  # With ref
  ps2 <- filterSampleData(FuentesIliGutData, ILI != "L2")
  ref_samples <- rownames(meta(subset_samples(ps2, ILI == "C")))
  dist.mat.ref <- phyloseq::distance(ps2, "bray")
  uniqueness.df.ref <- uniqueness(ps2,
                                  dist_mat = dist.mat.ref,
                                  reference_samples = ref_samples)
  col.exp <- 63
  expect_equal(ncol(uniqueness.df.ref), col.exp)

  exp.uniq.ref <- c(0.4519785, 0.4678278, 0.5050511, 0.5122802, 0.6030346)
  expect_equal(uniqueness.df.ref$uniqueness[1:5], exp.uniq.ref, tolerance=1.0e-3)

})
