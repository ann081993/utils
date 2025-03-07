#### fp_util_src
# Author: Seungchan An
require(rcdk)
require(fingerprint)
require(ChemmineOB)
require(ChemmineR)
require(rdist)

# function smi2fp
# converts SMILES into fingerprint
# default method "pubchem"
fp_types <- c("standard", "extended", "graph", 
              "pubchem", "maccs", "kr", "estate", 
              'ECFP0', 'ECFP2', 'ECFP4', 'ECFP6', 'FCFP0', 'FCFP2', 'FCFP4', 'FCFP6')

smi2fp <- function(smi, type = "pubchem", as_matrix = TRUE) {
        gc()
        
        parts <- ceiling(length(smi) / 200)
        if(parts > 1) {
                result <- NULL
                for(p in 1:parts) {
                        from = (200 * (p - 1) + 1)
                        to = (200 * (p - 1) + ifelse((parts == p) & length(smi) %% 200 > 0, length(smi) %% 200, 200))
                        
                        result <- rbind(result, smi2fp(smi[from:to], type = type, as_matrix = T))
                        gc()
                        cat("... Fingerprinting ", from, "-", to, ":", p, "of", parts, "\n")
                }
        } else {
                if(!(type %in% fp_types)) stop(paste0("Invalid fingerprint category specified\n\tSelect type arg from { ", paste(fp_types, collapse = " "), " }"))
                result <- list()
                ctypes <- c('ECFP0', 'ECFP2', 'ECFP4', 'ECFP6', 'FCFP0', 'FCFP2', 'FCFP4', 'FCFP6')
                if(type %in% ctypes) {
                        circular.type = type
                        type = "circular"
                }
                
                for(i in smi) {
                        parsed_smi <- suppressWarnings(parse.smiles(i)[[1]])
                        fp <- try(get.fingerprint(parsed_smi, type = type, circular.type = circular.type), silent = TRUE)
                        if(!is(fp, 'try-error') | !is.null(parsed_smi)) {
                                result <- append(result, fp)
                        }
                }
                if(length(result) > 0 & as_matrix) result <- fp.to.matrix(result)
        }
        result
}

# function fp2sim
# returns similarity matrix from the matrix of fingerprints
fp2sim <- function(fp1, fp2 = NULL, verbose = TRUE) {
        if(!is.null(fp2)) {
                sim <- 1 - cdist(fp1, fp2, metric = "jaccard")
        } else {
                sim <- 1 - pdist(fp1, metric = "jaccard")
        }
        sim
}

# function remove_salt
remove_salt <- function(smi) {
        result <- rep(TRUE, length(smi))
        split <- sapply(strsplit(smi, split = "[.]"), "[[", 1)
        for(n in grep("[.]", smi)) {
                if(sum(split == split[n]) > 1) {
                        result[n] <- FALSE
                }
        }
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

# function smi_mcs
# returns maximum common structure (MCS) of given list of SMILES
smi_mcs <- function(smi) {
  if(length(smi) < 2) stop("Please provide 2 or more SMILES strings.")
  #print(smi)
  if(max(nchar(smi)) > 150) return(NA)
  mcs <- get.mcs(parse.smiles(smi[1])[[1]], parse.smiles(smi[2])[[1]])
  if(is.null(mcs)) return(NA)
  smi_mcs <- get.smiles(mcs,
                        flavor = smiles.flavors(c("Isomeric")))
  smi <- smi[-1]
  smi[1] <- smi_mcs
  if(length(smi) == 1) {
    gc()
    return(smi)
  } else {
    smi_mcs(smi)
  }
}

# function smi_largest
# returns the largest component in SMILES
smi_largest <- function(smi) {
        if(length(smi) > 1) {
                smi <- sapply(smi, smi_largest, USE.NAMES = F)
                #rJava::.jgc()
                smi
        } else {
                smi <- parse.smiles(smi)[[1]]
                smi <- get.smiles(get.largest.component(smi), flavor = smiles.flavors(c("Isomeric")))
                #rJava::.jgc()
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
        #rJava::.jgc()
        mol
}

# function smi2desc
# generates molecular descriptors from input SMILES
desc_types <- c("basic", "hybrid", "constitutional", "topological",  "electronic", "geometrical", "all")

smi2desc <- function(smi, type = "basic", as_matrix = TRUE, verbose = TRUE) {
        mol <- NULL

        if(!(type %in% desc_types)) stop(paste0("Invalid descriptor category specified\n\tSelect type arg from { ", paste(desc_types, collapse = " "), " }"))
        if(class(smi) == "list") mol <- smi
        if(type == "basic") {
                if(is.null(mol)) mol <- parse.smiles(smi)
                type <- get.desc.names()[c(9,5,14,28,27,12)]
        } else {
                if(is.null(mol)) mol <- parse.smiles.3d(smi)
                type <- get.desc.names(type = type)
        }
        
        desc <- suppressWarnings(eval.desc(mol, which.desc = type, verbose = T))
        to_remove <- c("Wgamma", "WG.unity", "TAE") # remove NA descriptors
        to_remove <- grepl(paste(to_remove, collapse = "|"), colnames(desc))
        desc <- desc[, !to_remove]
        
        if(as_matrix) desc <- as.matrix(desc)
        #rJava::.jgc()
        desc
}

cat("Loaded:\n",
    " library rcdk, fingerprint, ChemmineR, ChemmineOB \n",
    " function smi2fp(), fp2sim(), plot_smi(), smi2desc() ... \n")
