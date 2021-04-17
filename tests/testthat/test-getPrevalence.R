
test_that("calculatePD works", {
  data("FuentesIliGutData")

  prev_sort <- getPrevalence(FuentesIliGutData,
    return_rank = c("Family", "Genus"),
    return_taxa = c("ASV4", "ASV17", "ASV85", "ASV83"),
    sort = TRUE
  )

  pd_vals <- round(prev_sort$prevalence * 100)
  exp_vals <- c(99, 98, 78, 70)
  expect_equal(pd_vals, exp_vals)

  prev_sort2 <- getPrevalence(FuentesIliGutData,
    return_rank = c("Family", "Genus"),
    return_taxa = c("ASV4", "ASV17", "ASV85", "ASV83"),
    sort = F
  )

  pd_vals2 <- round(prev_sort2$prevalence * 100)
  exp_vals2 <- c(70, 78, 99, 98)
  expect_equal(pd_vals2, exp_vals2)

  cols_nms <- c("Family", "Genus")
  cols_nsm_tib <- colnames(prev_sort2)[3:4]
  expect_equal(cols_nms, cols_nms)
})
