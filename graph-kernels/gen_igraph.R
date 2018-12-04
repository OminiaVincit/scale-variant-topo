#####################################################################################################
##  File: gen_igraph.R 
##    generate igraph database files from .graphml files
##    example: genigraph("MUTAG")
##      to convert all *.graphml files in data/MUTAG/graphml to list of graphs: jointddf_MUTAG.RData
##
#####################################################################################################

if(!require(graphkernels)) install.packages("graphkernels")

genigraph <- function(datname){
  datdir <- paste("../real-networks/", datname, "/", sep="")
  cat(datdir, "\n")
  filelist <- paste(datdir, "graphlist.txt", sep="")
  
  # Load txt file
  lfiles <- scan(filelist, what="", sep="\n")
  jointnet <- list()
  
  for (i in 1:length(lfiles)){
    fname <- paste(datdir, lfiles[[i]], sep="")
    if (!file.exists(fname)){
      cat(fname)
      break
    }
    g <- read.graph(fname, format="graphml")
    jointnet[[i]] <- g
    if (i %% 50 == 0) {
      # for debug
      cat(i,"\n")
    }
  }
  
  outfile <- paste("jointddf_", datname, ".RData", sep="")
  save(jointnet, file=outfile)
}
