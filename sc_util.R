#### sc_util
# Author: Seungchan An
library(ggplot2)
library(ggpubr)
library(ggsci)
library(viridis)
library(cowplot)
library(patchwork)
library(dplyr)
library(reshape2)

# function OptiClust
OptiClust <- function(object, idents = NULL, rescale = TRUE, feature.plot = NULL,
                      nvarfeats = c(500, 1000, 2000), ndims = c(5, 10, 20)) {
        if(!is.null(idents)) { object <- subset(object, idents = idents) }
        
        nr <- length(nvarfeats)
        nc <- length(ndims)
        
        plist <- list()
        olist <- list()
        i = 1
        for(nf in nvarfeats) {
                for(nd in ndims) {
                        cat(i, "/", nr * nc, " ", paste0(nf, " features, ", nd, " dims\n"))
                        p <- Subcluster(object, rescale = rescale,
                                        nvarfeat = nf, ndim = nd,
                                        only.plot = T, feature.plot = feature.plot, verbose = F)
                        plist[[i]] <- p
                        #olist[[i]] <- DietSeurat(object, features = VariableFeatures(object),
                        #                         dimreducs = c("pca", "tsne"))
                        i = i + 1
                        gc(reset = TRUE)
                }
        }
        p <- wrap_plots(plist, ncol = nc, nrow = nr)
        print(p)
        return(p)
        #result <- list()
        #result$objects <- olist
        #result$dimplot <- p
        #result
}

# function Subcluster
Subcluster <- function(object, idents = NULL,
                       normalize = TRUE, nvarfeat = 1000, vf.filter = FALSE, 
                       scale = TRUE, scale.all = TRUE, vars.to.regress = NULL,
                       component.analysis = "pca", ndim = 20,
                       reduction =  "tsne", iter = 1000,
                       res = 0.03, fn.reudction = "tsne",
                       seed = 1, only.plot = FALSE, feature.plot = NULL, verbose = TRUE) {
        # scale.all == T | F do not affect reduced dimensions and clusters
        if(!is.null(idents)) { object <- subset(object, idents = idents) }
        if(normalize) { object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = verbose) } # same as default
        object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nvarfeat, verbose = verbose)
        vf <- VariableFeatures(object)
        if(vf.filter) {
                ae <- AverageExpression(object, features = vf)
                ae <- colSums(t(ae$RNA))
                ind <- ae > (mean(ae) / 2)
                print(paste("Using", sum(ind), "/", length(ind), "variable features"))
                VariableFeatures(object) <- vf[ind]
        }
        if(scale.all) { vf <- rownames(object) }
        if(scale) { object <- ScaleData(object, features = vf, vars.to.regress = vars.to.regress, verbose = verbose) } # scale using all genes
        
        if(component.analysis == "ica") {		
                object <- RunICA(object, verbose = verbose)
                object <- RunTSNE(object, dims = 1:ndim, seed.use = seed, num_threads = 39, reduction = "ica", max_iter = iter, verbose = verbose)
        } else {
                object <- RunPCA(object, verbose = verbose)
                object <- RunTSNE(object, dims = 1:ndim, seed.use = seed, num_threads = 39, max_iter = iter, verbose = verbose)
        }
        if(fn.reudction == "tsne") {
                object <- FindNeighbors(object, reduction = "tsne", dims = 1:2, verbose = verbose)
        } else {
                object <- FindNeighbors(object, reduction = component.analysis, dims = 1:ndim, verbose = verbose)
        }
        object <- FindClusters(object, resolution = res, verbose = verbose)
        
        nclust <- as.numeric(tail(levels(object$seurat_clusters), 1))
        
        if(is.null(feature.plot)) {
                p <- suppressMessages(DimPlot(object, reduction = "tsne", label = TRUE) + 
                                              ggtitle(paste0(nvarfeat, " variable features for PCA\n", ndim, " dimensions for t-SNE\n", nclust + 1, " clusters (res=", res, ")")) +
                                              theme(plot.title = element_text(hjust = 0.5, size = 10, face = "plain")) &
                                              NoAxes() & NoLegend())
        } else {
                p <- suppressMessages(FeaturePlot(object, reduction = "tsne", label = TRUE, features = feature.plot, order = TRUE) + 
                                              ggtitle(paste0(nvarfeat, " variable features for PCA\n", ndim, " dimensions for t-SNE\n", nclust + 1, " clusters (res=", res, ")")) +
                                              theme(plot.title = element_text(hjust = 0.5, size = 10, face = "plain")) &
                                              NoAxes() & NoLegend() & FeatureCol())
        }
        if(only.plot) { object <- p }
        gc(reset = TRUE)
        return(object)
}

# function CompositionAnalysis
CompositionAnalysis <- function(object, x, y) {
        meta_names <- names(object@meta.data)
        if(!all(c(x,y) %in% meta_names)) stop(paste0(c("Please check whether input meta name is in:", meta_names), collapse = " "))
        
        batch <- object[[x]][, 1]
        clust <- object[[y]][, 1]
        
        fr <- NULL
        gr <- NULL
        cl <- NULL
        for(b in unique(batch)) {
                batch_fr <- table(clust[batch == b]) / sum(table(clust[batch == b])) * 100
                fr <- c(fr, batch_fr)
                gr <- c(gr, rep(b, length(batch_fr)))
                cl <- c(cl, names(batch_fr))
        }
        composition <- data.frame(fraction = round(fr, 1),
                                  group = factor(gr, levels = levels(batch)),
                                  cluster = factor(cl, levels = levels(clust)))
        return(composition)
}

# function CompositionPlot
CompositionPlot <- function(composition, cols = NULL, ncol = NULL) {
        p <- ggbarplot(data = composition, xlab = "", ylab = "Fraction (%)", 
                       x = "group", y = "fraction",
                       fill = "cluster", facet.by = "cluster", palette = cols, 
                       label = TRUE, lab.pos = "out") +
                scale_y_continuous(expand = expansion(mult = c(0,0.2))) + rotate_x_text(45) +
                facet_wrap(~ cluster, ncol = ncol)
        return(p)
}

# function ClusterDEG
ClusterDEG <- function(object, cont, exp, test.use = "wilcox") {
        barcode_df <- FetchData(object, vars = c("orig.ident", "celltypes"))
        barcode_df$barcode <- rownames(barcode_df)
        
        de_results <- NULL
        for(cl in levels(barcode_df$celltypes)) {
                cells1 <- barcode_df %>%
                        filter(orig.ident == cont) %>%
                        filter(celltypes == cl) %>%
                        pull(barcode)
                cells2 <- barcode_df %>%
                        filter(orig.ident == exp) %>%
                        filter(celltypes == cl) %>%
                        pull(barcode)
                
                if(length(cells1) > 3 & length(cells2) > 3) {
                        de <- FindMarkers(object, ident.1 = cells2, ident.2 = cells1,
                                          logfc.threshold = 0.1, min.pct = 0.05, test.use = test.use) # logfc.threshold = 0.25
                        de$celltype <- cl
                        de$symbol <- rownames(de)
                        de_results <- rbind(de_results, de)
                }
        }
        de_results[de_results$p_val_adj < 0.05, ]
}

# function DotPlot2
DotPlot2 <- function(object, features = g) {
        p <- DotPlot(object, features = features) +
                theme(panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill = NA, size = 1), # element_blank(),
                      #axis.text.x.bottom = element_blank(), # axis.ticks.x = element_blank()
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(face = "italic"),
                      axis.line.x = element_blank(),
                      axis.line.y = element_blank()) & FeatureCol2(ngray = 2) & rotate_x_text(45)
        return(p)
}

# function BarPlot
BarPlot <- function(object, features = g, ncol = NULL, cols = NULL, error = "mean_se",
                    group.by = NULL, split.by = NULL, slot = "data", size = NULL) {
        found <- features %in% rownames(object)
        if(any(!found)) { cat(paste0("The following requested features were not found: ",
                                     paste0(features[!found], collapse = ", "), "\n")) }
        features <- features[found]
        g_ex <- GetAssayData(object = object, slot = slot)[features, , drop = FALSE]
        od <- order(object@reductions$pca@cell.embeddings[, "PC_1"])
        
        ncell <- ncol(g_ex)
        nfeat <- nrow(g_ex)
        
        if(is.null(ncol)) { ncol <- ceiling(sqrt(nfeat)) }
        if(!is.null(group.by)) { Idents(object) <- object[[group.by]] }
        
        g_ex <- g_ex[, od, drop = FALSE]
        df <- melt(t(as.matrix(g_ex)), varnames = c("cell", "gene"))
        
        df$ident <- rep(Idents(object)[od], nrow(df) / ncell)
        if(!is.null(split.by)) { df$split <- rep(unlist(object[[split.by]])[od],  nrow(df) / ncell) }
        
        df <- df[with(df, order(gene, ident)), ]
        df <- df %>% group_by(ident, gene) %>% mutate(med = quantile(value)[4])
        df <- df %>% group_by(ident, gene) %>% mutate(med = median(value, na.rm = TRUE))
        df$cell <- rep(1:length(Idents(object)), nrow(df) / length(Idents(object)))
        
        plist <- list()
        n = 1
        for(g in features) {
                df_subset <- df[df$gene == g, ]
                fill <- ifelse(is.null(split.by), "ident", "split")
                h2 <- ggbarplot(df_subset, x = "ident", y = "value", fill = fill, size = size,
                                palette = cols, color = "black", position = position_dodge(0.9), 
                                add = error, error.plot = "upper_errorbar", add.params = list(width = 0.3, size = size)) +
                        ggtitle(g) + FeatureTitle() +
                        theme(panel.background = element_blank(),
                              panel.border = element_rect(colour = "black", fill = NA, size = 1), # element_blank(),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              axis.line.x = element_blank(),
                              axis.line.y = element_blank()) +
                        scale_y_continuous(expand = expansion(mult = c(0,0.1))) + NoLegend() + rotate_x_text(45)
                #p <- wrap_plots(list(h1, h2), widths = c(3,1))
                plist[[n]] <- h2; n = n + 1
        }
        return(wrap_plots(plist, ncol = ncol))
}

# function FeatureScatter2
FeatureScatter2 <- function(object, f1, f2) {
        df_exp <- t(as.matrix(GetAssayData(object, slot = "data")[c(f1, f2), ]))
        df_exp <- data.frame(df_exp)
        summary(df_exp)
        #df_exp <- df_exp[df_exp[, 1] * df_exp[, 2] > 0, ]
        dim(df_exp)
        p <- ggscatter(df_exp, x = colnames(df_exp)[1], y = colnames(df_exp)[2], size = 0.2,
                       add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, # add = "reg.line",
                       cor.coef = TRUE, cor.method = "pearson") +#  + stat_density_2d(aes(fill = ..level..), geom = "polygon") + gradient_fill(c("white", "steelblue"))
                coord_cartesian(expand = c(0), xlim = c(0, max(df_exp[, 1])), ylim = c(0, max(df_exp[, 2])))
        return(p)
}

# function FeatureCol
FeatureCol <- function() {
        return(scale_color_gradientn(colors = c("gray", "cyan4", "yellow", "red", "darkred")))
}

# function FeatureCol2
FeatureCol2 <- function(mode = "red", ngray = 4) {
        grays <- c("gray90", "gray90", "gray85", "gray80", "gray70", "gray60" )
        grays <- grays[1:ngray]
        if(mode == "red") {
                col <- scale_color_gradientn(colors = c(grays, "red", "darkred"))
        } else {
                col <- scale_color_gradientn(colors = c(grays, "blue", "darkblue"))
        }
        return(col)
}

# function FeatureTtile
FeatureTitle <- function(size = 16, face = "bold.italic", color = "black") {
        return(theme(plot.title = element_text(hjust = 0.5, size = size, face = face, color = color)))
}

# function NoTitle
NoTitle <- function(no_text_x = TRUE) {
        return(theme(plot.title = element_blank()))
}

# function NoAxesTitle
NoAxesTitle <- function(no_text_x = TRUE) {
        if (no_text_x) {
                return(theme(axis.title = element_blank(), axis.text.x = element_blank()))
        } else {
                return(theme(axis.title = element_blank()))
        }
}

cat("Loaded:\n",
    " library ggplot2, ggpubr, patchwork, cowplot, reshape2 ...\n",
    " function BarPlot(), CompositionAnalysis(), Subcluster() ...\n")
