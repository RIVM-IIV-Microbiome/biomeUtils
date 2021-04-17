test_that("trendAbundanceAssignment", {

  data("SprockettTHData")
  plottrend <- trendAbundanceAssignment(SprockettTHData,
                           quantiles = seq(0, 95, by = 5),
                           plot=T)
  df <- plottrend$data
  levels.ranks <- df$abQuant[c(1:7,50:56)]
  exp_vl <-c("2","2","2","2","2","2","2","327", "327", "327", "327", "327", "327", "327")
  expect_equal(levels.ranks,exp_vl)

  tax_rnk <- as.character(df$Var1[c(1:7,50,56)])
  exp_rnk <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Kingdom","Species")
  expect_equal(tax_rnk,exp_rnk)

  numtax_rnk <- as.character(df$numb_taxa[c(1:7,50,56)])
  exp_numtax_rnk <- c("2319","2319","2319","2319","2319","2319","2319","1507","1507")
  expect_equal(numtax_rnk,exp_numtax_rnk)
})
