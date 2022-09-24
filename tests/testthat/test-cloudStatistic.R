test_that("cloudStatistic test", {

  library(biomeUtils)
  data("FuentesIliGutData")
  # Keep only two groups. controls and L1.
  ps <- filterSampleData(FuentesIliGutData, ILI != "L2")
  ps <- getProportions(ps)
  # Define controls as reference samples
  ref_samples <- rownames(meta(subset_samples(ps, ILI == "C")))

  cloud_results <- cloudStatistic(ps,
                                  dist_mat = phyloseq::distance(ps, "bray"),
                                  reference_samples = ref_samples,
                                  ndim=-1)
  # check stats match
  stats.ret <- cloud_results$stats[1:5]
  exp.stats <- c(0.9021859, 0.9168544, 1.0071445, 0.9591956, 1.0628701)
  expect_equal(stats.ret, exp.stats, tolerance=1.0e-3)

  # check pvals match
  pval.ret <- cloud_results$pvals[1:5]
  exp.pval <- c(0.9890710, 0.9508197, 0.3989071, 0.6721311, 0.1530055)
  expect_equal(pval.ret, exp.pval, tolerance=1.0e-3)

  # check all samples are returned
  expect_equal(nrow(cloud_results), nsamples(ps))

})
