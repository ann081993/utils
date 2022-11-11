#### fp_util_src
# Author: Seungchan An
library(rcdk)
library(fingerprint)
library(ChemmineOB)
library(ChemmineR)

# function jacdis
# returns Jaccard (Tanimoto) similarity between two fingerprints
jacdis <- function(a1, a2) {
        return(as.numeric(1 - dist(rbind(a1, a2), method = "binary")))
}

# function smi2fp
# converts SMILES into fingerprint
# default method "pubchem"
fp_types <- c("standard", "extended", "graph", 
              "pubchem", "maccs", "kr", "estate", 
              'ECFP0', 'ECFP2', 'ECFP4', 'ECFP6', 'FCFP0', 'FCFP2', 'FCFP4', 'FCFP6')

smi2fp <- function(smi, type = "pubchem", circular.type = "ECFP6",
                   as_matrix = TRUE, verbose = TRUE) {
        if(!(type %in% fp_types)) stop(paste0("Invalid fingerprint category specified\n\tSelect type arg from { ", paste(fp_types, collapse = " "), " }"))
        result <- list()
        ctypes <- c('ECFP0', 'ECFP2', 'ECFP4', 'ECFP6', 'FCFP0', 'FCFP2', 'FCFP4', 'FCFP6')
        if(type %in% ctypes) {
                circular.type = type
                type = "circular"
        }
        
        for (i in smi) {
                if(verbose) cat(i)
                parsed_smi <- suppressWarnings(parse.smiles(i)[[1]])
                fp <- try(get.fingerprint(parsed_smi, type = type, circular.type = circular.type), silent = TRUE)
                if(!is(fp, 'try-error') | !is.null(parsed_smi)) {
                        result <- append(result, fp)
                        if(verbose) cat("\tO", "\n")
                } else if(verbose) cat("\n")
        }
        if(length(result) > 0 & as_matrix) result <- fp.to.matrix(result)
        rJava::.jgc()
        result
}

# function fp2jacdis
# returns similarity matrix from the matrix of fingerprints
fp2jacdis <- function(matrix, names = NULL, part = NULL, verbose = TRUE) {
        t1 <- Sys.time()
        if(verbose) print(t1)
        
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
        if(verbose) print(t2)
        if(verbose) print(t2-t1)
        if(verbose) cat("\n")
        
        result
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

# function smi_ncomp
# returns number of components in SMILES
smi_ncomp <- function(smi) {
        if(length(smi) > 1) {
                sapply(smi, smi_ncomp, USE.NAMES = F)
        } else {
                sum(grepl("[.]", strsplit(smi, split = "")[[1]])) + 1
        }
}

# function smi_largest
# returns the largest component in SMILES
smi_largest <- function(smi) {
        if(length(smi) > 1) {
                sapply(smi, smi_largest, USE.NAMES = F)
        } else {
                smi <- parse.smiles(smi)[[1]]
                smi <- get.smiles(get.largest.component(smi), flavor = smiles.flavors(c("Isomeric")))
                smi
        }
}

# function smi_inorg
# returns whether the structure is inorganic or not
smi_inorg <- function(smi) {
        smi <- gsub("Ca|Cd|Ce|Cl|Co|Cr|Cs|Cu|Ac|Mc|Sc|Tc", "", smi)
        !grepl("C|c", smi)
}

# function parse.smiles.3d
# converts SMILES into molecules with 3D coordinates
parse.smiles.3d <- function(smi, iterator = FALSE) {
        writeLines(smi, "tmp.smiles")
        convertFormatFile("SMI","SDF","tmp.smiles","tmp.sdf")
        sdf <- read.SDFset("tmp.sdf")
        message("... SDF file loaded")
        sdf <- suppressWarnings(generate3DCoords(sdf))
        message("... 3D coordinates generated")
        cid(sdf) <- paste0("CID", 1:length(sdf))
        write.SDF(sdf, file="tmp.sdf")
        rm(sdf)
        if(iterator) {
                mol <- iload.molecules("tmp.sdf", type = "sdf")
                message("... Iterator prepared for molecules with 3D coordinates")
        } else {
                mol <- load.molecules("tmp.sdf")
                message("... Molecules parsed with 3D coordinates")
        }
        rJava::.jgc()
        mol
}

# function smi2desc
# generates molecular descriptors from input SMILES
desc_types <- c("basic", "hybrid", "constitutional", "topological",  "electronic", "geometrical", "all")

smi2desc <- function(smi, type = "basic", as_matrix = TRUE, verbose = TRUE) {
        if(!(type %in% desc_types)) stop(paste0("Invalid descriptor category specified\n\tSelect type arg from { ", paste(desc_types, collapse = " "), " }"))
        if(type == "basic") {
                mol <- parse.smiles(smi)
                type <- get.desc.names()[c(9,5,14,28,27,12)]
        } else {
                mol <- parse.smiles.3d(smi)
                type <- get.desc.names(type = type)
        }
        
        desc <- suppressWarnings(eval.desc(mol, which.desc = type, verbose = T))
        to_remove <- c("Wgamma", "WG.unity", "TAE") # remove NA descriptors
        to_remove <- grepl(paste(to_remove, collapse = "|"), colnames(desc))
        desc <- desc[, !to_remove]
        
        if(as_matrix) desc <- as.matrix(desc)
        rJava::.jgc()
        desc
}

cat("Loaded:\n",
    " library rcdk, fingerprint, ChemmineR, ChemmineOB \n",
    " function smi2fp(), fp2jacdis(), plot_smi(), smi2desc() ... \n")
