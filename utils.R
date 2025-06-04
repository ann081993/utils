#### utils
# Author: Seungchan An

# function tocb()
# copies vectors into line-separated list
tocb <- function(x) {
        x <- paste(x, collapse = "\n")
        cat(x)
        writeClipboard(x)
}

# function fromcb()
# returns a vector from line-separated list
fromcb <- function() {
        x <- readClipboard()
        return(unlist(strsplit(x, split = "\n")))
}

# function ul()
# returns length of unique elements of a vector
ul <- function(x) {
        return(length(unique(x)))
}

# function htable()
# returns a head or tail of ordered frequency table
htable <- function(x, n = 6) {
        tab <- table(x)
        tab <- tab[order(tab, decreasing = TRUE)]
        if(n > 0) {
                tab <- head(tab, n)
        } else {
                tab <- tail(tab, -n)
        }
        return(tab)
}

# function pbar()
pbar <- function() {
        # library(progress)
        l <- 10 # length of for loop
        pb <- progress_bar$new(format = "Progress: [:bar] :percent, Elapsed time :elapsedfull",
                               total = l, width = 80, clear = F, force = T)
        # add tick in for loop
        for(n in 1:l) {
                pb$tick()
        }
}

# function list2df()
# returns data.frame from the named list input
list2df <- function(l) {
        # Create the data frame
        df <- data.frame(
          value = unlist(l),  # Flatten all values
          group = rep(names(l), times = lengths(l))  # Repeat names
        )
        df
}

# plot utils
library(ggplot2)
library(ggpubr)
library(ggsci)
library(cowplot)
library(grid)

if(.Platform$OS.type == "windows") {
  if(!"Arial" %in% names(windowsFonts())) {
    windowsFonts(Arial = windowsFont("Arial"))
    message("'Arial' was registered in windowsFonts()")
  }
}

font_color = "black"
font_size = 8
line_width = 0.75 / 2.14
plot_factor = 1
print(paste("plot factor:", plot_factor))

theme_pub <- function() {
  font_size = font_size * plot_factor
  line_width = line_width * plot_factor
          
  theme(
    plot.background = element_blank(),
    plot.margin = margin(rep(0.03, 4), unit = "cm"),
    strip.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = line_width, fill = NA, color = font_color),
    axis.title.x = element_text(size = font_size, family = "Arial", color = font_color),
    axis.title.y = element_text(size = font_size, family = "Arial", color = font_color, angle = 90),
    axis.text.x = element_text(size = font_size, family = "Arial", color = font_color),
    axis.text.y = element_text(size = font_size, family = "Arial", color = font_color),
    strip.text = element_text(size = font_size, family = "Arial", color = font_color),
    plot.title = element_text(size = font_size, family = "Arial", face = "plain",
                              color = font_color, hjust = 0.5, vjust = 0.5),
    axis.line = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.ticks = element_line(linewidth = line_width, color = font_color),
    axis.ticks.length = unit(0.05 * plot_factor, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.margin = margin(rep(0.01, 4), unit = "cm"),
    legend.text = element_text(size = font_size - 1, family = "Arial", color = font_color),
    legend.title = element_text(size = font_size - 1, family = "Arial", color = font_color),
    legend.frame = element_rect(linewidth = line_width, color = font_color),
    legend.ticks = element_line(linewidth = line_width, color = font_color),
    legend.ticks.length = unit(0.1 * plot_factor, "cm"),
    legend.key.width = unit(0.3 * plot_factor, "cm"), 
    legend.key.height = unit(0.4 * plot_factor, "cm"),
    legend.position = "right")
}

theme_no_axes <- function() {
  return(
    theme_pub() + theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank())
  )
}

theme_no_legend <- function() {
  return(
    theme(legend.position = "none")
  )
}

theme_pub_white <- function() {
  font_color = "white"
  theme_pub()
}

plot_panel_resize <- function(gp, w, h) {
  grob <- ggplotGrob(gp)
  panel_index <- which(grob$layout$name == "panel")
  grob$widths[grob$layout[panel_index, ]$l] <- unit(w, "cm")
  grob$heights[grob$layout[panel_index, ]$t] <- unit(h, "cm")
  grid.newpage()
  grid.draw(grob)
}

y_zero <- function(y_expension = 0.1) {
  scale_y_continuous(expand = expansion(mult = c(0, y_expension)), labels = function(l) {
    # Check if there are any non-integer values
    if (all(l == floor(l), na.rm = TRUE)) {
      return(as.character(l))  # Keep integers as they are
    } else {
      max.decimals <- max(nchar(stringr::str_extract(as.character(l), "\\.[0-9]+")), na.rm = TRUE) - 1
      lnew <- formatC(l, replace.zero = TRUE, zero.print = "0",
                      digits = max.decimals, format = "f", preserve.width = TRUE)
      return(lnew)
    }
  })
}

theme_pub_dens <- function() {
  theme_pub() + theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())
}

theme_dot <- function() {
  theme_pub() + theme(panel.grid.major = element_line(color = "darkgray",
                                                      linewidth = line_width, linetype = 3),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank())
}
