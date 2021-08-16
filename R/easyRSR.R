#' Easy to calculate Rank Sum Ratio
#'
#' It is easy to calculate Rank Sum Ratio (RSR) for each evaluation
#' objects base on different evaluation indicators. And at the same time, you
#' are able to rank and group these evaluation objects quickly.
#'
#' @param data A matrix, data frame or seurat object. Whe you input the matrxi
#' or data frame, the rows should include evaluation indicators (eg. genes,
#' clinical indicators), and the columns should include evaluation objects
#' (eg. clinical samples, cell subsets). The data should not include 'NA'value.
#' Values of each row (evaluation indicators) in data should not be all zero.
#' When you input seurat object, you need to think about parameter 'assay',
#' 'slot' and 'group.by'.
#'
#' @param assay Pull out matrix from this assay of the Seurat object. if NULL,
#' it use \code{DefaultAssay(data)}).
#'
#' @param slot Pull out matrix from this slot of the Seurat object. And the
#' 'data' is default.
#'
#' @param group.by You can input the colname of meta.data from seurat object
#' (e.g, ident, celltype, clusters). And the 'ident' is default.
#'
#' @param geneset The evaluation indicators (eg. genes, clinical indicators)
#' which are included in the rownames of the data. if NULL, it use all
#' rownames of data.
#'
#' @param weights A numeric vector for weights of each evaluation indicator.
#' The length of weights must be equal to evaluation indicator. And the each
#' value of weights should be greater than zero. You can enter weights or not.
#'
#' @param impacts A character vector for influences of each evaluation
#' indicator. The length of impacts must be equal to evaluation indicator.
#' The impacts must be '+' or '-' sign. And the '+' sign means the indicator
#' has positive effect on evaluation object. You need to enter the impacts.
#'
#' @param rank A character vector for calculation of the rank of evaluation
#' indicator. The rank must be 'non-integer' or 'integer'. The indicator is
#' arranged by the their values and we will get the integer rank. The indicator
#' is arranged by linear imputation and we will get the non-integer rank which
#' it is able to reflects the real difference of the values of indicator
#' if the rank is 'non-integer'. Sometimes we can not group the evaluation
#' objects when the rank is 'non-integer'. And the 'non-integer' is default.
#'
#' @param class A character vector. The class must be '3', '4', '5', '6', '7'
#' '8', '9' or '10'.The higher the number of class is, the more groups are
#' identified. You can simply recognize 'class' as 'resolution' although
#'  it is not.
#'
#' @return If easyRSR can not group the evaluation objects, it will return a
#' data frame. The data frame includes RSR (the rank sum ratio) and RSR.rank.
#' We can compare the evaluation objects by their RSR and acquire their rank.
#' If easyRSR can group the evaluation objects, it will return a list. List
#' includes two data frames. The first data frame includes RSR, RSR.rank and
#' group. The group already satisfies the following conditions : 1. At least 2
#' observations in each group ; 2. The homogeneity of variance of each group
#' is satisfied ; 3. There is statistically significant difference between
#' the RSR of different groups by ANOVA. The sceond data frame includes the
#' p value of multiple-comparisons are calculated by SNK q Test.
#'
#'
#' @export
#'
#' @references # First example is cited from : Zhenqiu Sun, Yongyong Xu. Medical
#' # Statistics (4nd ed): Beijing, China: People’s Medical Publishing House,
#' # 2014:419 (in Chinese). Sceond example is cited from : Tirosh et al,
#' Science (2016).
#'
#' @examples # First example
#' # creata data frame
#' test <- data.frame(prenatal.visit.rate = c(99.54,96.52,99.36,92.83,91.71,95.35,
#' 96.09,99.27,94.76,84.80), maternal.mortality = c(60.27,59.67,43.91,58.99,35.40,
#' 44.71,49.81,31.69,22.91,81.49), perinatal.mortality = c(16.15,20.10,15.60,17.04,
#' 15.01,13.93,17.43,13.89,19.87,23.63))
#' rownames(test) <- c(LETTERS[1:10])
#' test <- t(test)
#' # Calculate
#' weights = c(1,1,1)
#' impacts = c("+","-","-")
#' a <- easyRSR(data = test, weights = weights, impacts = impacts,
#'              rank = "integer", class = "3")
#' a$RSR.result
#' a$RSR.plot
#' a <- easyRSR(data = test, weights = weights, impacts = impacts,
#'              rank = "non-integer", class = "3")
#' a$RSR.result
#' a$RSR.plot
#'
#'
#' # Second example
#' b <- c('CD79B', 'CD79A', 'CD19', 'CD180', 'CD200', 'CD3D', 'CD2', 'CD3E',
#' 'CD7','CD8A','CD14','CD1C','CD68','CD9','CD247')
#' a <- easyRSR(data = pbmc_small, geneset = b,
#'              group.by="RNA_snn_res.1", impacts = rep("+",15))
#' a$RSR.result
#'
#'

easyRSR = function(data = NULL, assay = NULL, slot = "data",
                   group.by = "ident", geneset = NULL,
                   weights = NULL, impacts = NULL,
                   rank = "non-integer", class = NULL){

  if(missing(data))
    stop("'data' does not exist")
  if (is.null(geneset)){
    geneset <- rownames(data)
  }
  if (! is.null(weights) ) {
    if(length(weights) != length(geneset))
      stop("length of 'weights' is not equal to number of geneset")
    if(! all(weights > 0))
      stop("weights must be positive numbers")
  }
  if(missing(impacts))
    stop("'impacts' must be a character vector")
  if(length(impacts) != length(geneset))
    stop("length of 'impacts' is not equal to number of geneset")
  if(! is.character(impacts))
    stop("'impacts' must be a character vector.")
  if(! all(impacts == "+" | impacts == "-"))
    stop("'impacts' must be only '+' or '-' sign")

  if ( all(methods::is(data)=="Seurat") ) {
    if (is.null(assay)) {assay <- Seurat::DefaultAssay(data)}
    data_seurat <- Seurat::AverageExpression(object = data, slot= slot,
                                             assay = assay,  features = geneset,
                                             group.by = group.by)
    if (is.null(data_seurat[[assay]])) {stop("None of the features of geneset are found in the Seurat object.")}
    if (length(intersect(geneset,rownames(data)))==1) {
      data <- data_seurat[[assay]]
      rownames(data) <- intersect(geneset, rownames(data))
    }else{
      data <- data_seurat[[assay]]}
  }else{
    if(! is.matrix(data) && ! is.data.frame(data))
      stop("'data' must be a matirx or data frame")
    gene_intersect <- intersect(geneset, rownames(data))
    if (length( gene_intersect )==0) {stop("None of the features of geneset are found in the data.")}
    if (all(geneset %in% gene_intersect)) {
      data <- data[geneset, ]
    }else{
      data <- data[gene_intersect, ]
      print(paste0("The following features were not found in the data: ",setdiff(gene_intersect, geneset)))
    }
    if( any(is.na(data)) )
      stop("'data' have NA")
    if(! is.numeric(as.matrix(data)) )
      stop("'data' must be numeric")
    if( any(colSums(data)==0) )
      stop("Some genes of geneset in 'data' are all zero")
  }

  data <- t(data)
  weights <- weights[colnames(data) %in% geneset]
  impacts <- impacts[colnames(data) %in% geneset]
  if(! rank %in% c("non-integer","integer"))
    stop("'rank' must be only 'non-integer' or 'integer' sign")


  # creatre empty matrix
  N <- matrix(nrow = nrow(data), ncol = ncol(data))

  # calculate rank
  if (rank == "non-integer") {
    for (i in seq_along(impacts)) {
      if (impacts[i]=="+") {
        for (j in 1:nrow(data)) {
          N[j,i] <- 1 + (nrow(data) - 1) * (data[j,i] - min(data[,i])) / (max(data[,i]) - min(data[,i]))
        }
      }else{
        for (j in 1:nrow(data)) {
          N[j,i] <- 1 + (nrow(data) - 1) * (max(data[,i]) - data[j,i]) / (max(data[,i]) - min(data[,i]))
        }
      }
    }
  }else{
    for (i in seq_along(impacts)) {
      if (impacts[i]=="+") {
        N[,i] <- rank(data[,i], ties.method = "average")
      }else{
        N[,i] <- rank(-data[,i], ties.method = "average")
      }
    }
  }

  # calculate RSR
  if (is.null(weights)) {
    N <- as.data.frame(N)
    N$RSR <- apply(N, 1, sum) / (nrow(N) * ncol(N))
  }else{

    weights <- weights / sum(weights)
    for (i in 1:ncol(N)) {
      N[,i] <- N[,i]*weights[i]
    }
    N <- as.data.frame(N)
    N$RSR <- apply(N, 1, sum) / nrow(N)
  }
  colnames(N)[1:(ncol(N)-1)] <- colnames(data)
  rownames(N) <- rownames(data)
  # calculate the rank of RSR
  N$RSR <- round(N$RSR, 4)
  N$RSR.rank <- data.table::frankv(N$RSR, order = -1, ties.method = "dense")
  # calculate the distribution of RSR
  RSR.distribution <- N[,c("RSR","RSR.rank")]
  RSR.distribution$RSR.rank.mean <- data.table::frankv(RSR.distribution$RSR, order = 1, ties.method = "average")
  RSR.distribution <- stats::aggregate(RSR.distribution, list(RSR.distribution$RSR), mean)
  RSR.distribution <- RSR.distribution[-1]
  RSR.distribution$Probit <- apply(RSR.distribution["RSR.rank.mean"], 1, function(x){(stats::qnorm(x/nrow(data))) + 5})
  RSR.distribution$Probit[nrow(RSR.distribution)] <- stats::qnorm(1-1/(4*nrow(data))) + 5

  # create plot function
  plotRSR = function(plot.data, plot.fill){
    # create the data of plot
    plot.data <- as.data.frame(plot.data)
    plot.data <- tibble::rownames_to_column(plot.data, var = "cluster")
    if (plot.fill == "RSR") {
      temp <- RColorBrewer::brewer.pal(11, "Spectral")
      temp[6] <- "gold"
      p <- ggpubr::ggbarplot(plot.data, x = "cluster", y = "RSR", fill = plot.fill, color = "white",
                sort.val = "desc",
                sort.by.groups = FALSE,
                x.text.angle = 45,
                xlab = FALSE)+ ggpubr::gradient_fill(temp)
      return(p)
    }
    if (plot.fill == "group" & (length(plot.data$group) <= 30) ) {
      # barplot
      p <- ggpubr::ggbarplot(plot.data, x = "cluster", y = "RSR", fill = plot.fill ,
                color = "white", palette = "simpsons",
                sort.val = "desc", sort.by.groups=FALSE,
                x.text.angle = 45, xlab = FALSE)
      return(p)

    }
    if (plot.fill == "group" & (length(plot.data$group) > 30) ) {
      # boxplot
      p <- ggpubr::ggboxplot(plot.data, x = "group", y = "RSR", color = plot.fill ,
                        palette = "npg", width = 0.3, add.params = list(size = 1),
                        bxp.errorbar = T, bxp.errorbar.width = 0.2,
                        add="jitter", xlab = FALSE)
      return(p)
    }

  }

  # calculate the regression equation
  rsr <- stats::lm(RSR ~ Probit, data = RSR.distribution)
  # save original RSR
  RSR.distribution$RSR.original <- RSR.distribution$RSR

  # test the normality of residuals, perform box-cox convert
  if (length(stats::residuals(rsr)) > 5000) {
    if (stats::ks.test(stats::residuals(rsr), stats::pnorm)$p.value <= 0.05) {
      lambada <- car::powerTransform(RSR ~ Probit, data = RSR.distribution)
      RSR.distribution$RSR <- car::bcPower(RSR.distribution$RSR, lambada$lambda)
      rsr <- stats::lm(RSR ~ Probit, data = RSR.distribution)
    }
  }else{
    if (stats::shapiro.test(stats::residuals(rsr))$p.value <= 0.05) {
      lambada <- car::powerTransform(RSR ~ Probit, data = RSR.distribution)
      RSR.distribution$RSR <- car::bcPower(RSR.distribution$RSR, lambada$lambda)
      rsr <- stats::lm(RSR ~ Probit, data = RSR.distribution)
    }

  }

  # test oral model and regression coefficient by ANOVA
  pvalue <- stats::anova(rsr)[1,5]
  if (pvalue > 0.05) {
    writeLines("The P value of the regression equation is greater than 0.05.\nWe can't calculate the regression equation and classify the samples.")
    results <- list(RSR.result = N[,c("RSR","RSR.rank")])
    results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "RSR")
    return(results)
  }



  # test the assumption of regression equation
  gvmodel <- gvlma::gvlma(rsr)
  if (as.numeric(gvmodel$GlobalTest$GlobalStat4$pvalue) <= 0.05) {
    # test the normality of residuals
    if (length(stats::residuals(rsr)) > 5000) {
      if (stats::ks.test(stats::residuals(rsr), stats::pnorm)$p.value <= 0.05) {
        writeLines("Residuals do not satisfy a normal distribution.\nWe can't calculate the regression equation and classify the samples.")
        results <- list(RSR.result = N[,c("RSR","RSR.rank")])
        results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "RSR")
        return(results)
      }
    }else{
      if (stats::shapiro.test(stats::residuals(rsr))$p.value <= 0.05) {
        writeLines("Residuals do not satisfy a normal distribution.\nWe can't calculate the regression equation and classify the samples.")
        results <- list(RSR.result = N[,c("RSR","RSR.rank")])
        results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "RSR")
        return(results)
      }
    }
    # test the homogeneity of residuals
    if (car::ncvTest(rsr)$p <= 0.05) {
      writeLines("The residuals do not satisfy the homogeneity of variances.\nIn this class, we can't calculate the regression equation and classify the samples.")
      results <- list(RSR.result = N[,c("RSR","RSR.rank")])
      results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "RSR")
      return(results)
    }
  }


  # make a boundary table and classify the samples by adjust RSR
  table <- data.frame(classification = c(rep("3",3),rep("4",4),rep("5",5),rep("6",6),rep("7",7),rep("8",8),rep("9",9),rep("10",10)),
                      level = c(1:3,1:4,1:5,1:6,1:7,1:8,1:9,1:10),
                      Percentile = c("[0% , 15.866%)" , "[15.866% , 84.134%)" , "[84.134% , 100%]" ,
                                     "[0% , 6.681%)" , "[6.681% , 50%)" , "[50% , 93.319%)" , "[93.319% , 100%]" ,
                                     "[0% , 3.593%)" , "[3.593% , 27.425%)" , "[27.425% , 72.575%)" , "[72.575% , 96.407%)" , "[96.407% , 100%]" ,
                                     "[0% , 2.275%)" , "[2.275% , 15.866%)" , "[15.866% , 50%)" , "[50% , 84.134%)" , "[84.134% , 97.725%)" , "[97.725% , 100%]" ,
                                     "[0% , 1.618%)" , "[1.618% , 10.027%)" , "[10.027% , 33.360%)" , "[33.360% , 67.003%)" , "[67.003% , 89.973%)" , "[89.973% , 98.382%)" , "[98.382% , 100%]" ,
                                     "[0% , 1.222%)" , "[1.222& , 8.681%)" , "[8.681% , 22.663%)" , "[22.663% , 50%)" , "[50% , 77.337%)" , "[77.337% , 93.319%)" , "[93.319% , 98.678%)", "[98.678% , 100%)" ,
                                     "[0% , 0.99)" ,"[0.99% , 4.746%)" ,"[4.746% , 15.866%)" ,"[15.866% , 37.070%)" ,"[37.070% , 62.930%)" ,"[62.930% , 84.134%)" ,"[84.134% , 95.254%)" ,"[95.254% , 99.010%)" ,"[99.010, 100%]" ,
                                     "[0% , 0.82%)" ,"[0.82% , 3.593%)" ,"[3.593% , 11.507%)" ,"[11.507% , 27.425%)" ,"[27.425% , 50%)" ,"[50% , 72.575%)" ,"[72.575% , 88.493%)" ,"[88.493% , 96.407%)" ,"[96.407% , 99.180%)" ,"[99.180%, 100%]"
                                     ) ,
                      Probit_down =c(-Inf,4,6,
                                     -Inf,3.5,5,6.5,
                                     -Inf,3.2,4.4,5.6,6.8,
                                     -Inf,3,4,5,6,7,
                                     -Inf,2.86,3.72,4.57,5.44,6.28,7.14,
                                     -Inf,2.78,3.50,4.25,5,5.75,6.50,7.22,
                                     -Inf,2.67,3.33,4,4.67,5.33,6,6.67,7.33,
                                     -Inf,2.6,3.2,3.8,4.4,5,5.6,6.2,6.8,7.4),
                      Probit_up = c(4,6,Inf,
                                    3.5,5,6.5,Inf,
                                    3.2,4.4,5.6,6.8,Inf,
                                    3,4,5,6,7,Inf,
                                    2.86,3.72,4.57,5.44,6.28,7.14,Inf,
                                    2.78,3.50,4.25,5,5.75,6.50,7.22,Inf,
                                    2.67,3.33,4,4.67,5.33,6,6.67,7.33,Inf,
                                    2.6,3.2,3.8,4.4,5,5.6,6.2,6.8,7.4,Inf))

  if (is.null(class)) {

    results <- list(RSR.result = N[,c("RSR","RSR.rank")])
    results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "RSR")
    return(results)

  }
  if(! is.character(class))
    stop("'class' must be a character vector")
  if( length(class) > 1 )
    stop("'class' must be a single number")
  if(! class %in% c("3","4","5","6","7","8","9","10"))
    stop("'class' must be '3', '4', '5', '6' , '7' , '8', '9' or '10' " )

  # filter by class
  classification <- NULL
  table <- table %>% dplyr::filter(classification == class)

  table <- dplyr::filter(table, table$classification == class)

  # caulate threshold value
  BRSR <- stats::coefficients(rsr)[1] + (stats::coefficients(rsr)[2] * unique(c(table$Probit_down,table$Probit_up)) )
  RSR.distribution$levels <- cut(RSR.distribution$RSR, breaks = BRSR, labels = table$level , right = F)
  RSR.distribution$levels <- as.factor(as.character(RSR.distribution$levels))

  # merge
  if (any((RSR.distribution[,c("RSR")]==RSR.distribution[,c("RSR.original")])==F)) {
    RSR.merge <- dplyr::left_join(N[,c("RSR","RSR.rank")],RSR.distribution[,c("RSR","levels","RSR.original")],by = c("RSR" = "RSR.original"))
    colnames(RSR.merge)[which(colnames(RSR.merge)=="RSR")] <- "RSR.original"
    colnames(RSR.merge)[which(colnames(RSR.merge)=="RSR.y")] <- "RSR"

  }else{
    RSR.merge <- dplyr::left_join(N[,c("RSR","RSR.rank")],RSR.distribution[,c("RSR","levels","RSR.original")],'RSR')
  }

  rownames(RSR.merge) <- rownames(N)
  RSR.merge$group <- as.numeric(factor(RSR.merge$levels))
  RSR.merge$group <- paste0("G",RSR.merge$group)
  RSR.merge$group <- as.factor(RSR.merge$group)


  # the observation in each group
  if ( all(table(RSR.merge$group) < 2) ) {
    writeLines("Each group only has 1 observation that we can't perform the test of best classification.\nPerhaps the number of oral samples are too few or the class is not fit.")
    RSR.merge$RSR <- RSR.merge$RSR.original
    results <- list(RSR.result = RSR.merge[,c("RSR","RSR.rank","group")])
    results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "group")
    return(results)
  }

  # the observation in each group
  if ((length(which((table(RSR.merge$group) > 1) == T)) < 2)) {
    writeLines("Only 1 group's observations are over than 1 that we can't perform the test of best classification.\nPerhaps the number of oral samples are too few or the class is not fit.")
    RSR.merge$RSR <- RSR.merge$RSR.original
    results <- list(RSR.result = RSR.merge[,c("RSR","RSR.rank","group")])
    results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "group")
    return(results)
  }

  # filiter some groups which observations are less than 2
  RSR.merge.new <- dplyr::filter(RSR.merge, RSR.merge$group %in% names(which((table(RSR.merge$group) > 1) == T)))

  # test the normality of RSR
  if (length(RSR.merge.new[,"RSR"]) > 5000) {
    if (stats::ks.test(RSR.merge.new[,"RSR"], stats::pnorm)$p.value <= 0.05) {
      writeLines("RSR do not satisfy a normal distribution.\nWe can't perform ANOVA analysis and classify the samples.")
      results <- list(RSR.result = N[,c("RSR","RSR.rank")])
      results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "RSR")
      return(results)
    }
  }else{
    if (stats::shapiro.test(RSR.merge.new[,"RSR"])$p.value <= 0.05) {
      writeLines("RSR do not satisfy a normal distribution.\nWe can't perform ANOVA analysis and classify the samples.")
      results <- list(RSR.result = N[,c("RSR","RSR.rank")])
      results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "RSR")
      return(results)
    }
  }

  # test the homogeneity of variance among groups
  if (car::leveneTest( RSR ~ group, data = RSR.merge.new)$`Pr(>F)`[1] <= 0.05) {
    writeLines("The variances of each group are not equal.\nIn this class, we can't classify the samples exactly.\nMaybe the class is not fit.")
    results <- list(RSR.result = N[,c("RSR","RSR.rank")])
    results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "RSR")
    return(results)
  }

  # perform ANOVA
  if (summary(stats::aov( RSR ~ group, data = RSR.merge.new))[[1]][1,5] > 0.05) {
    writeLines("There is not statistically significant in ANOVA analysis.\nIn this class, we can't classify the samples exactly.\nMaybe the class is not fit.")
    results <- list(RSR.result = N[,c("RSR","RSR.rank")])
    results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "RSR")
    return(results)
  }

  # perform SNK-q test
  RSR.snk <- agricolae::SNK.test(stats::aov( RSR ~ group, data = RSR.merge.new),"group", group=FALSE)
  RSR.snk <- as.data.frame(RSR.snk$comparison)
  RSR.snk <- dplyr::select(RSR.snk, "pvalue")
  RSR.merge$RSR <- RSR.merge$RSR.original
  results <- list(RSR.result = RSR.merge[,c("RSR","RSR.rank","group")],
                  RSR.snk.test = RSR.snk)
  results[["RSR.plot"]] <- plotRSR(plot.data = results[["RSR.result"]], plot.fill = "group")
  return(results)

}
