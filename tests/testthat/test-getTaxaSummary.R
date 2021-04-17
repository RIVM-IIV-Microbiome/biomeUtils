test_that("getTaxaSummary works", {

  expect_equal(2 * 2, 4)
  data("SprockettTHData")
  sum_tax <- getTaxaSummary(SprockettTHData, rank="Phylum")
  cl_phy <- sum_tax$Phylum[1:5]
  ex_rn <- c("Firmicutes","Bacteroidetes","Actinobacteria","Proteobacteria","Verrucomicrobia")
  expect_equal(cl_phy,ex_rn)

  cl_count <- sum_tax$Counts[1:5]
  ex_cnt <- c(18062700, 11025765, 9931237, 5613362, 381980)
  expect_equal(cl_count,ex_cnt)

  cl_perc <- sum_tax$Percent[1:5]
  ex_perc <- c(39.3254049,24.0048649,21.6219013,12.2211925,0.8316319)
  expect_equal(cl_perc,ex_perc, 0.0000001)

})
