if(!require(graphkernels)) install.packages("graphkernels")
if(!require(parallel)) install.packages("parallel")
library(parallel)

no_cores <- detectCores()
workers <- makeCluster(6)

config <- 10
posfix <- "config"
load(paste("jointddf_", config, "_", posfix, ".RData", sep=""))
name.list <- c("Graphlet","KStepRandomWalk", "GeometricRandomWalk", "ExponentialRandomWalk", "ShortestPath","WL")
func.list <- paste0("Calculate", name.list, "Kernel")
par.list <- c(4, 2.0, 5e-2, .1, NA, 5)

kerdir <- "../../Datawork/NetworkDataset/Joint-benchmark/exp_20180628/kernel/"

N <- length(name.list)
x <- seq_len(N)

# Export variables and library to workers
clusterExport(workers, "config")
clusterExport(workers, "posfix")
clusterExport(workers, "name.list")
clusterExport(workers, "func.list")
clusterExport(workers, "par.list")
clusterExport(workers, "kerdir")
clusterExport(workers, "jointnet")
clusterEvalQ(workers, library(graphkernels))

### run kernel
runkernel <- function(i){
  outfile <- paste(kerdir, "ph_20180628_Euclid_norm_0_dim_1_kernel_", name.list[i], "_", posfix, "_", config, "_par_", par.list[i],".txt", sep="")
  #cat(func.list[i], "\n")
  #cat(par.list[i], "\n")
  
  if (is.na(par.list[i])) {
    runtime <- system.time(K <- eval(call(func.list[i], jointnet)))
  } else {
    if (name.list[i] == "KStepRandomWalk") {
      runtime <- system.time(K <- eval(call(func.list[i], jointnet, c(1.0, 1.0, 1.0))))
    } else {
      runtime <- system.time(K <- eval(call(func.list[i], jointnet, par.list[i])))
    }
  }

  write.table(K, file=outfile, sep="\t", row.names = FALSE, col.names = FALSE)
  #cat(outfile, "\n")
  #cat(runtime, "\n")
}

cat(system.time(parLapply(workers, x, runkernel)), "\n")
stopCluster(workers)