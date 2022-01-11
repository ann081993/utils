
library(reshape2)
BarPlot <- function(object, features = g, ncol = 3, cols = NULL, error = "mean_sd") {
        g_ex <- GetAssayData(object = object)[features, ]
        od <- order(object@reductions$pca@cell.embeddings[, "PC_1"])
        
        if(length(features) > 1) {
                g_ex <- g_ex[, od]
                df <- melt(t(as.matrix(g_ex)), varnames = c("cell", "gene"))
        } else {
                g_ex <- g_ex[od]
                df <- melt(as.matrix(g_ex), varnames = c("cell", "unused"))
                df$gene <- features
                df$unused <- NULL
        }
        df$ident <- rep(Idents(object)[od], nrow(df) / length(Idents(object)))
        
        df <- df[with(df, order(gene, ident)), ]
        df <- df %>% group_by(ident, gene) %>% mutate(med = quantile(value)[4])
        df <- df %>% group_by(ident, gene) %>% mutate(med = median(value, na.rm = TRUE))
        df$cell <- rep(1:length(Idents(object)), nrow(df) / length(Idents(object)))
        
        plist <- list()
        n = 1
        for(g in features) {
                df_subset <- df[df$gene == g, ]
                h2 <- ggbarplot(df_subset, x = "ident", y = "value", fill = "ident",
                                palette = cols, color = "black", add = error, error.plot = "upper_errorbar", add.params = list(width = 0.3)) +
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
