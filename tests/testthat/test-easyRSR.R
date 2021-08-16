test_that("multiplication works", {
  # First example
  test <- data.frame(
    yqjclv = c(99.54, 96.52, 99.36, 92.83, 91.71, 95.35, 96.09, 99.27, 94.76, 84.80),
    yfswl = c(60.27, 59.67, 43.91, 58.99, 35.40, 44.71, 49.81, 31.69, 22.91, 81.49),
    wceswl = c(16.15, 20.10, 15.60, 17.04, 15.01, 13.93, 17.43, 13.89, 19.87, 23.63))
  rownames(test) <- c(LETTERS[1:10])
  test <- t(test)
  a <- easyRSR(data=test,weights = c(1,1,1),
               impacts = c("+","-","-"),rank = "non-integer",class = "3")
  a$RSR.result
  a$RSR.plot
  a <- easyRSR(data=test,weights = c(1,1,1),
               impacts = c("+","-","-"),rank = "integer",class = "3")
  a$RSR.result
  a$RSR.plot
  # Second example
  b <- c('CD79B', 'CD79A', 'CD19', 'CD180', 'CD200', 'CD3D', 'CD2', 'CD3E',
  'CD7','CD8A','CD14','CD1C','CD68','CD9','CD247')
  a <- easyRSR(data = pbmc_small, geneset = b,
               group.by="RNA_snn_res.1", impacts = rep("+",15))
  a$RSR.result

})
