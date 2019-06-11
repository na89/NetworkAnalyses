#Select matching expression pattern
ecdf_fun <- function(x, perc){
  return(ecdf(x)(perc))
}


PreprocessExpression <- function(ExpressionData){
  all_quant <- c()
  for (i in names(ExpressionData)) {
    all_quant <-
      c(all_quant, seq(5, 95, 10)[data.table::between(
        x = ecdf_fun(ExpressionData, ExpressionData[i]),
        lower = seq(0, 0.9, 0.1),
        upper = seq(0.1, 1, 0.1)
      )])
  }
  names(all_quant) <- names(ExpressionData)
  return(all_quant)
}

SelectMatchingExpression <-
  function(TargetGenes, ExpressionData, Nrand, QuantizedGenes) {
    if (missing(QuantizedGenes)){
      all_quant <- PreprocessExpression(ExpressionData)
    } else{
      all_quant <- QuantizedGenes
    }
    expr_quant <- all_quant[names(all_quant) %in% TargetGenes]
    expr_quant_dist <- table(expr_quant)
    rnd_match_genes <- vector("list", Nrand)
    for (i in 1:Nrand){
      for (j in 1:length(expr_quant_dist)){
        rnd_match_genes[[i]] <- c(rnd_match_genes[[i]],
          sample(names(all_quant[all_quant == names(expr_quant_dist)[j]]),
                  size = expr_quant_dist[j]))
      }
    }
    return(rnd_match_genes)
  }

SelectMatchingExpressionParallel <-
  function(TargetGenes, ExpressionData, Nrand, QuantizedGenes) {
    if (missing(QuantizedGenes)){
      all_quant <- PreprocessExpression(ExpressionData)
    } else{
      all_quant <- QuantizedGenes
    }
    expr_quant <- all_quant[names(all_quant) %in% TargetGenes]
    expr_quant_dist <- table(expr_quant)
    rnd_match_genes <- vector("list", Nrand)
    rnd_match_genes <- mclapply(c(1:Nrand), 
                                FUN = function(x){
                                  tmp.c <- c()
                                  for (j in 1:length(expr_quant_dist)){
                                    tmp.c <- 
                                      c(tmp.c,
                                        sample(names(all_quant[all_quant ==
                                                                 names(expr_quant_dist)[j]]),
                                               size = expr_quant_dist[j]))
                                  }
                                  return(tmp.c)
                                }, mc.cleanup = TRUE, mc.cores = 32)
      
    return(rnd_match_genes)
  }



ParallelRandomGeneSetGenerate <- 
  function(tested_datasets, ExpressionDataInput,
           Nrandom = 100, PreprocessedExprData){
    rnd.data <- vector("list", length(tested_datasets))
    for (i in 1:length(tested_datasets)){
      rnd.data[[i]] <- vector("list", length(ExpressionDataInput))
      message(Sys.time(), 
              " processing ", 
              i," out of ", 
              length(tested_datasets))
      rnd.data[[i]] <- 
        mcmapply(FUN = SelectMatchingExpression, 
                 ExpressionData = ExpressionDataInput,
                 TargetGenes = tested_datasets[[i]],
                 Nrand = Nrandom,
                 QuantizedGenes = PreprocessedExprData
        )
    }
    return(rnd.data)
  }

ParallelEstimateConnectivityP <-
  function(ModuleEdges, Network, Module, RndData, Degrees){
    rnd.conn.est <- 
      (length(E(subNetwork(unique(c(V(Module)$name, 
                                 RndData)), Network))) -
      ModuleEdges -
       length(E(subNetwork(RndData, Network)))) /
      sum(Degrees[names(Degrees) %in% RndData], na.rm = TRUE)
    return(rnd.conn.est)
  }