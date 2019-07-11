display.pointset <- function(file, mapfile){
  X <- read.csv(file, header=FALSE)
  map <- read.csv(mapfile, header=FALSE)
  dev.new()
  index <- which( X[,1]=="Fixed" | X[,1]=="Moving" | X[,1]=="FixedTransformed" )
  plot(X[index, 2:3], col=X[index,1])
  
  fixed.index <- which( X[,1]=="Fixed")
  moving.index <- which( X[,1]=="Moving")
  segments( X[ fixed.index[ map[,1] +1 ] , 2], X[ fixed.index[ map[,1]+1] , 3],
            X[ moving.index[ map[,2]+1] , 2], X[ moving.index[ map[,2]+1] , 3],
            col="#00000033" )


  title(file)
}

display.pointset("ot-balanced.csv", "ot-balanced-coupling.csv")
display.pointset("ot-unbalanced.csv", "ot-unbalanced-coupling.csv")

