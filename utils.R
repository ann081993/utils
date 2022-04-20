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
        
        
