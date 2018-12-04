if(!require(graphkernels)) install.packages("graphkernels")
if(!require(parallel)) install.packages("parallel")
library(parallel)

# run following command before executing this function 
# datname <- "MUTAG"

name.list <- c("Graphlet","KStepRandomWalk", "GeometricRandomWalk", "ExponentialRandomWalk", "ShortestPath","WL")
func.list <- paste0("Calculate", name.list, "Kernel")
par.list <- c(4, 2.0, 5e-2, .1, NA, 5)
N <- length(name.list)

### calculate kernel
x <- seq_len(N)


### run kernel
rkernel <- function(i){
  load(paste("jointddf_", datname, ".RData", sep=""))
  kerdir <- paste("../../Datawork/NetworkDataset/", datname, "/exp_20180628/kernel/", sep="")
  
  outfile <- paste(kerdir, "graph_kernel_", name.list[i], "_par_", par.list[i],".txt", sep="")
  cat(func.list[i], "\t", par.list[i], "\n")
  
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
  cat(outfile, "\n")
  cat(runtime, "\n")
}

### calculate kernel for multi-cores
calkernel <- function(){
  no_cores <- detectCores()
  workers <- makeCluster(N)
  # Export variables and library to workers
  clusterExport(workers, "name.list")
  clusterExport(workers, "func.list")
  clusterExport(workers, "par.list")
  clusterExport(workers, "datname")
  clusterEvalQ(workers, library(graphkernels))
  
  cat(system.time(parLapply(workers, x, rkernel)), "\n")
  stopCluster(workers)
}

