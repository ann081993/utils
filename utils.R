#### utils
# Author: Seungchan An

# function tocb
# copies vectors into line-separated list
tocb <- function(x) {
        x <- paste(x, collapse = "\n")
        cat(x)
        writeClipboard(x)
}

# function fromcb
# returns a vector from line-separated list
fromcb <- function() {
        x <- readClipboard()
        return(unlist(strsplit(x, split = "\n")))
}

# function ul
# returns length of unique elements of a vector
ul <- function(x) {
        returns(unique(length(x)))
}
