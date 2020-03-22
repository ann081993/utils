#### fp_util_src

library(rcdk)
library(fingerprint)

jacdis <- function(a1, a2) {
        g1 <- grep(1, a1)
        g2 <- grep(1, a2)
        
        lint <- length(intersect(g1, g2))
        luni <- length(union(g1, g2))
        
        return(lint/luni)
}

# function smi2fp
smi2fp <- function(smi, method = "pubchem") {
        result <- list()
        for (i in smi) {
                result <- append(result, get.fingerprint(parse.smiles(i)[[1]], method))
        }
        
        result <- fp.to.matrix(result)
        return(result)
}

# function fp2jacdis
fp2jacdis <- function(matrix, names = NULL) {
        
        t1 <- Sys.time()
        cat("\n")
        print(t1)
        cat("\n")
        
        result <- matrix(nrow = nrow(matrix), ncol = nrow(matrix))
        for (r1 in 1:nrow(matrix)) { # 2019-10-16 Updated to (O^2)/2 operation, not O^2 
                for (r2 in r1:nrow(matrix)) {
                        a1 <- matrix[r1, ]
                        a2 <- matrix[r2, ]
                        b <- jacdis(a1, a2)
                        result[r1, r2] <- b
                        result[r2, r1] <- b
                }
        }
        
        t2 <- Sys.time()
        cat("\n")
        print(t2)
        cat("\n")
        print(t2-t1)
        cat("\n")
        
        return(result)
}

cat("Loaded:\n",
    " library rcdk, fingerprint\n",
    " function smi2fp(), fp2jacdis()\n")
