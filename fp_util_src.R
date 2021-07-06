#### fp_util_src
# Author: Seungchan An
library(rcdk)
library(fingerprint)

# function jacdis
# returns Jaccard (Tanimoto) similarity between two fingerprints
jacdis <- function(a1, a2) {
        return(as.numeric(1 - dist(rbind(a1, a2), method = "binary")))
}

# function smi2fp
# converts SMILES into fingerprint
# default method "pubchem"
smi2fp <- function(smi, method = "pubchem", as_matrix = TRUE) {
        result <- list()
        for (i in smi) {
                cat(i)
                parsed_smi <- suppressWarnings(parse.smiles(i)[[1]])
                fp <- try(get.fingerprint(parsed_smi, method), silent = TRUE)
                if(!is(fp, 'try-error') | !is.null(parsed_smi)) {
                        result <- append(result, fp)
                        cat("\tO", "\n")
                } else cat("\n")
        }
        if(length(result) > 0 & as_matrix) result <- fp.to.matrix(result)
        return(result)
}

# function fp2jacdis
# returns similarity matrix from the matrix of fingerprints
fp2jacdis <- function(matrix, names = NULL, part = NULL) {
        
        t1 <- Sys.time()
        print(t1)
        
        if(is.null(part)) {
                result <- matrix(nrow = nrow(matrix), ncol = nrow(matrix))
                for (r1 in 1:nrow(matrix)) { 
                        for (r2 in r1:nrow(matrix)) {
                                a1 <- matrix[r1, ]
                                a2 <- matrix[r2, ]
                                b <- jacdis(a1, a2)
                                result[r1, r2] <- b
                                result[r2, r1] <- b
                        }
                }
        } else {
                result <- matrix(nrow = part, ncol = nrow(matrix) - part)
                for (r1 in 1:part) { 
                        for (r2 in 1:(nrow(matrix) - part)) {
                                a1 <- matrix[r1, ]
                                a2 <- matrix[r2 + part, ]
                                b <- jacdis(a1, a2)
                                result[r1, r2] <- b
                        }
                }
        }
        
        t2 <- Sys.time()
        print(t2)
        print(t2-t1)
        cat("\n")
        
        return(result)
}

# function plot_smi
# plots with 2D structures from SMILES
plot_smi <- function(smiles, names = NULL, nrow = 3, ncol = 4, part = 1) {
        from = (nrow * ncol) * (part - 1)
        to = from + ifelse(length(smiles) - from > nrow * ncol,
                           nrow * ncol,
                           length(smiles) - from)
        from = from + 1
        
        cat("... Plotting", length(from:to), "of", length(smiles),"chemical structures\n")
        cat("... Part", part, "of", ceiling(length(smiles) / (nrow * ncol)))
        smiles <- smiles[from:to]
        names <- names[from:to]

        par(mar = c(1,1,1,1), mfrow = c(nrow, ncol))
        for(k in 1:length(smiles)) {
                m <- parse.smiles(smiles[k])[[1]]
                plot(0:10, xaxt = 'n', yaxt = 'n', ann = FALSE, type = 'n', bty = 'n', asp = 1) # asp fixes ratio https://stackoverflow.com/questions/8693558/how-to-define-fixed-aspect-ratio-for-scatter-plot
                img <- view.image.2d(m, get.depictor(width = 400, height = 400, abbr = 'reagents', style = 'cow', zoom = 1.2))
                rasterImage(img, 0, 0, 10, 10)
                text(x = 5, y = 1, names[k], cex = 1)
        }
}


cat("Loaded:\n",
    " library rcdk, fingerprint\n",
    " function smi2fp(), fp2jacdis(), plot_smi()\n")
