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

# function ClusterCor
ClusterCor <- function(object) {
        avex <- AverageExpression(object)$RNA
        avex_cor <- cor(as.matrix(avex), method = "spearman")
        heatmap(avex_cor, col = colorRampPalette(c("green", "white", "red"))(20), symm = T)
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
        if(!is.null(group.by)) { Idents(object) <- unlist(object[[group.by]]) }
        
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

# functional enrichment analysis

string_analysis <- function(gene_list, suffix = NULL, thr = 300) {
  require(STRINGdb)
  require(igraph)
  require(export)
  
  string_db <- STRINGdb$new(species = 10090, version = "12.0", score_threshold = thr) # human 9606
  mapped <- string_db$map(data.frame(gene = gene_list), "gene"); mapped <- mapped[complete.cases(mapped), ]
  
  annotations <- string_db$get_proteins()
  mapped <- merge(mapped, annotations, by.x = "STRING_id", by.y = "protein_external_id",
                  all.x = TRUE)
  
  interactions <- suppressWarnings(string_db$get_interactions(mapped$STRING_id))
  interactions <- interactions[!duplicated(interactions), ]
  edges <- interactions[, c("from", "to")]
  edge_weights <- interactions$combined_score / 1000 
  
  net <- suppressWarnings(graph_from_data_frame(edges, directed = F))
  V(net)$name <- mapped$preferred_name[match(V(net)$name, mapped$STRING_id)]
  E(net)$weight <- edge_weights
  
  isol <- which(igraph::degree(net) == 0)
  net <- delete_vertices(net, isol)
  
  de = igraph::degree(net)
  st = igraph::strength(net)
  
  set.seed(1)
  l <- igraph::layout_nicely(net)
  
  wc <- suppressWarnings(cluster_edge_betweenness(net, weights = E(net)$weight, directed = FALSE, bridges = TRUE))
  mb <- membership(wc)

  top_cluster <- which(table(mb) / length(mb) > 0.05)
  vcolors <- rep("darkgray", length(V(net)))
  for(tc in top_cluster) vcolors[mb == tc] <- lighten(categorical_pal(7), 0.1)[tc]
  V(net)$color <- vcolors
  
  print(plot(net, layout = l,
       edge.color = adjustcolor("black", 0.3), 
       edge.width = (E(net)$weight^2)*3,
       vertex.frame.color = "white", vertex.frame.width = 0.4,
       vertex.size = st^0.5 * 2, 
       vertex.label.dist = ifelse(st^0.5 > 2.5, 0, sample(c(-0.6, 0.6), length(V(net)), replace = T)),
       mark.col = adjustcolor("gray", 0.3), mark.border = NA,
       vertex.label.cex = 6 / 12, vertex.label.font = 3, vertex.label.color = "black", 
       vertex.label.family = "Arial"))
  
  psize = max(as.integer(log2(length(V(net))) * 3) * 0.394, 2)
  graph2ppt(file = paste("Fig STRING network", suffix),
            width = psize,
            height = psize)
  print("STRING network was constructed.")
  
  require(gprofiler2)
  res <- NULL
  for(i in top_cluster) {
    g = names(mb)[mb == i]
    if(length(g) < 5) next
    gostres <- gost(query = g, ordered_query = F, significant = T, evcodes = T,
                    exclude_iea = T, sources = c("GO:BP"),
                    organism = "mmusculus")
    gostres$result <- gostres$result[!gostres$result$term_id %in% unlist(gostres$result$parents) & gostres$result$term_size < 500, ]
    if(!is.null(gostres$result)) {
      if(nrow(gostres$result) > 0) {
        gostres$result$cluster <- i
        res <- rbind(res, gostres$result)
      }
    }
  }
  if(is.null(res)) return()
  
  print("Enrichment analysis was performed.")
  res2 <- res
  res2$mlp <- -log10(res2$p_value)
  res2$cluster <- factor(res2$cluster)
  bl <- grep("protozoan|virus|bacter|fungu|built from|antimicrobial|vertebrate eye| poly|yeast", res2$term_name) # blacklist
  if(length(bl) > 0) res2 <- res2[-bl, ]
  res2$term_name[duplicated(res2$term_name)] <- paste0(" ", res2$term_name[duplicated(res2$term_name)])
  p <- ggbarplot(res2, x = "term_name", xlab = "", ylab = "-log P",
                 legend = "none",
                 y = "mlp", fill = "cluster", color = NA,
                 palette = colorspace::lighten(categorical_pal(7)[top_cluster], amount = 0.5)) +
    geom_hline(yintercept = -log10(0.05),
               linewidth = 3 / 8 / 1.07, linetype = 5,
               color = "darkgray") +
    rotate() + theme_pub() + theme(axis.text.y = element_text(size = 7))
  
  require(grid)
  grob <- ggplotGrob(p)
  panel_index <- which(grob$layout$name == "panel")
  grob$widths[grob$layout[panel_index, ]$l] <- unit(1.5, "cm")
  grob$heights[grob$layout[panel_index, ]$t] <- unit(nrow(p$data) * 0.22, "cm")
  grid.newpage()
  grid.draw(grob)
  
  graph2ppt(file = paste("Fig STRING enrichment", suffix),
            width = (max(nchar(as.character(p$data$term_name))) / 10 + 3) * 0.394,
            height = (nrow(p$data) * 0.3 + 1) * 0.394)
  table2ppt(x = res2[, -c(1,2,12,14,15)], file = paste("Fig STRING enrichment", suffix), append = T, font = "Arial", pointsize = 6)
}

fgsea_analysis <- function(loading, suffix = "",
                           category = c("C5", "H"), pathway_grep = "GOBP_|HALLMARK_|HP_",
                           max_size = 300, nes_cutoff = 2) {
  require(fgsea)
  require(msigdbr)
  require(data.table)
  require(export)
  
  msigdbr_collections() %>% mutate(included = ifelse(gs_cat %in% category, "O", "")) %>% as.data.frame() %>% print()  
  
  sigs <- rbindlist(lapply(category, function(x) msigdbr(species = "mouse", category = x)))
  sigs <- split(sigs$gene_symbol, sigs$gs_name)
  sl <- sapply(sigs, length)
  sigs <- sigs[sl > 5 & sl < max_size]
  
  # perform fgsea
  set.seed(1)
  fgsea_result <- fgsea(pathways = sigs, stats = loading, minSize = 3)
  print("fgsea was successfully performed:")
  
  fgsea_result_sorted <- fgsea_result %>% mutate(leadingEdge = sapply(leadingEdge, toString)) %>% arrange(desc(NES))
  write.csv(fgsea_result_sorted,
            file = paste0("fgsea result ", suffix, ".csv"))
  
  print(paste("    total", nrow(fgsea_result), "pathways"))
  fgsea_result <- fgsea_result[fgsea_result$padj < 0.05, ]
  print(paste("   ", nrow(fgsea_result), "enriched pathways"))
  
  fgsea_result <- fgsea_result[grepl(pathway_grep, fgsea_result$pathway), ]
  fgsea_result <- fgsea_result[abs(fgsea_result$NES) > nes_cutoff, ]
  print(paste("   ", nrow(fgsea_result), "filtered pathways"))
  
  if(nrow(fgsea_result) > 30) {
    print("Too many pathways to plot! Please adjust grep pattern (pathway_grep) or cutoff (nes_cutoff).")
    print(head(fgsea_result_sorted, 20))
    return()
  }
  require(ggrastr)
  p2 <- data.frame(rank = 1:length(loading),
                   metric = sort(loading, decreasing = T)) %>%
    ggplot(aes(x = rank, y = metric, fill = metric)) +
    rasterise(geom_bar(stat = "identity", size = 0), dpi = 600) +
    scale_fill_gradient2(low = "blue", high = "red") +
    geom_hline(yintercept = 0, linewidth = 3 / 8 / 1.07, color = "darkgray") +
    xlab("Rank") + ylab("Metric") + y_zero(0) +
    theme_pub() & NoLegend()
  
  for(pw in fgsea_result$pathway) {
    pname <- pw; print(pw)
    nes <- fgsea_result[fgsea_result$pathway == pname, "NES"]
    pv <- fgsea_result[fgsea_result$pathway == pname, "padj"]; pv <- scales::scientific(pv, digits = 2)
    gs <- paste0(unique(names(head(sort(loading[sigs[[pname]]], decreasing = T), 6))), collapse = ",")
    p1 <- plotEnrichment(sigs[[pname]],
                         loading, ticksSize = 2 / 8 / 1.07)
    p1$layers <- p1$layers[-5]
    
    if(as.numeric(nes) > 0) {
      p1$layers <- p1$layers[-4]
    } else {
      p1$layers <- p1$layers[-3]
    }
    
    p1 <- p1 + ylab("Enrichment score") + 
      ggtitle(paste0(c(pname,
                       paste0("NES = ", round(nes, 3), " / adj. P = ", pv),
                       gs), collapse = "\n")) +
      geom_hline(yintercept = 0, linewidth = 3 / 8 / 1.07, color = "darkgray") +
      theme_pub() + y_zero(0.05) +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank())

    print(plot_grid(plotlist = list(p1, p2), align = "v", ncol = 1, rel_heights = c(2,1)))
    
    graph2ppt(file = paste("Fig fgsea", suffix), width = 5.4 * 0.394, height = 4.5 * 0.394, append = T)
  }
}

gprofiler_analysis <- function(gene_list, suffix = "", terms_of_interest = NULL, term_sources = "GO:BP") {
  require(gprofiler2)
  require(export)
  
  gostres <- gost(query = gene_list, ordered_query = F, significant = T, evcodes = T,
                  exclude_iea = T, sources = term_sources,
                  organism = "mmusculus")
  
  gostres$result <- gostres$result[!gostres$result$term_id %in% unlist(gostres$result$parents) & gostres$result$term_size < 1000, ]
  
  if(is.null(gostres$result)) return()
  if(nrow(gostres$result) == 0) return()

  res2 <- gostres$result
  res2$mlp <- -log10(res2$p_value)
  
  bl <- grep("protozoan|virus|bacter|fungu|built from|antimicrobial|vertebrate eye| poly|yeast", res2$term_name) # blacklist
  if(length(bl) > 0) res2 <- res2[-bl, ]
  res2$term_name[duplicated(res2$term_name)] <- paste0(" ", res2$term_name[duplicated(res2$term_name)])
  
  if(!is.null(terms_of_interest)) {
    if(sum(res2$term_name %in% terms_of_interest) > 0) {
      res2 <- res2[res2$term_name %in% terms_of_interest, ]
    } else {
      res2 <- res2[grep(terms_of_interest, res2$term_name), ]
    }
  }
  
  res2 <- res2[nchar(res2$term_name) > 10, ]
  
  if(is.null(res2)) return()
  p <- ggbarplot(res2, x = "term_name", xlab = "", ylab = "-log P",
                 legend = "none",
                 y = "mlp", fill = "darkgray", color = NA, amount = 0.5) +
    geom_hline(yintercept = -log10(0.05),
               linewidth = 3 / 8 / 1.07, linetype = 5,
               color = "darkgray") +
    rotate() + theme_pub() + theme(axis.text.y = element_text(size = 7))
  
  require(grid)
  grob <- ggplotGrob(p)
  panel_index <- which(grob$layout$name == "panel")
  grob$widths[grob$layout[panel_index, ]$l] <- unit(1.5, "cm")
  grob$heights[grob$layout[panel_index, ]$t] <- unit(nrow(p$data) * 0.22, "cm")
  grid.newpage()
  grid.draw(grob)
  
  graph2ppt(file = paste("Fig gprofiler enrichment", suffix),
            width = (max(nchar(as.character(p$data$term_name))) / 10 + 3) * 0.394,
            height = (nrow(p$data) * 0.3 + 1) * 0.394, append = T)
  table2ppt(x = res2[, -c(1,2,12,14,15)], file = paste("Fig gprofiler enrichment", suffix),
            append = T, font = "Arial", pointsize = 6) # remove parents ... etc
  table2ppt(x = data.frame(V1 = c(paste("n =", length(gene_list)), gene_list), stringsAsFactors = FALSE),
            file = paste("Fig gprofiler enrichment", suffix), append = T, font = "Arial", pointsize = 6, width = 5)
}


cat("Loaded:\n",
    " library ggplot2, ggpubr, patchwork, cowplot, reshape2 ...\n",
    " function BarPlot(), CompositionAnalysis(), Subcluster() ...\n",
    " functions for functional enrichment analyses gprofiler_analysis(), fgsea_analysis(), string_analysis() ...\n")
