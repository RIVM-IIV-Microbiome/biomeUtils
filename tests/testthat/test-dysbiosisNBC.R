test_that("dysbiosisScore works", {
  data("FuentesIliGutData")
  ps <- filterSampleData(FuentesIliGutData, ILI != "L2")
  ps <- getProportions(ps)
  ref_samples <- rownames(meta(subset_samples(ps, ILI == "C")))
  # median BC
  dys_bc <- dysbiosisScore(ps,
                           reference_samples = ref_samples,
                           method="median-BC")

  col.exp <- 63
  expect_equal(ncol(dys_bc), col.exp)
  expect_equal(nrow(dys_bc), 397)

  exp.uniq <- c(0.6041963, 0.6649651, 0.8163710, 0.7969920, 0.8851318)
  expect_equal(dys_bc$dysbiosis.score[1:5], exp.uniq, tolerance=1.0e-3)

  # BC-ED
  dys_bc_ed <- dysbiosisScore(ps,
                              reference_samples = ref_samples,
                              method="BC-ED")

  col.exp <- 64
  expect_equal(ncol(dys_bc_ed), col.exp)
  expect_equal(nrow(dys_bc_ed), 397)

  exp.uniq <- c(0.12505899, 0.04288149, 0.24758961, 0.10788604, 0.15602434)
  expect_equal(dys_bc_ed$dysbiosis.score[1:5], exp.uniq, tolerance=1.0e-3)

})
