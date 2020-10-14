#### pubchem_util
# Author: Seungchan An
# Ref: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest

# function name2smi()
# returns isomeric SMILES from chemical name using PubChem PUG REST API
name2smi <- function(name) {
        result <- NA
        data <- try(read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                                name, "/property/IsomericSMILES/csv"),
                         stringsAsFactors = F), silent = T)
        if(class(data) != "try-error") { result <- data$IsomericSMILES[1] }
        return(result)
}

# function cid2smi()
# returns isomeric SMILES from chemical name using PubChem PUG REST API
cid2smi <- function(cid) {
        result <- NA
        data <- try(read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                                    cid, "/property/IsomericSMILES/csv"),
                             stringsAsFactors = F), silent = T)
        if(class(data) != "try-error") { result <- data$IsomericSMILES[1] }
        return(result)
} 

# function smi2name()
# returns name (preferred name) from SMILES using PubChem PUG REST API
smi2name <- function(smi) {
        result <- NA
        data <- try(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
                                     smi, "/synonyms/txt"), n = 1), silent = T)
        if(class(data) != "try-error") { result <- data }
        return(result)
}

# function cid2name()
# returns name (preferred name) from cid using PubChem PUG REST API
cid2name <- function(cid) {
        result <- NA
        data <- try(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                                     cid, "/synonyms/txt"), n = 1), silent = T)
        if(class(data) != "try-error") { result <- data }
        return(result)
}

# function name2cid()
# returns cid from name using PubChem PUG REST API
# derived from kegg_flavonoids.R written 2019-07-16
name2cid <- function(name) {
        cid <- tryCatch({ readLines(paste0("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                                           name, "/cids/txt")) },
                        error = function(e) { return(NA) }
        )
        if(length(cid) > 1) {
                cid <- cid[which.min(cid)] # this line was added to prevent get 2 or more cids
        }
        return(cid)
}

# function get_bioassay()
# returns "biological test results" from cid
# derived from kegg_flavonoids.R written 2019-07-16
get_bioassay <- function(cid) {
        bioassay_results <- read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=jsonp&query={%22download%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22cid%22:%22",
                                            cid, "%22}]},%22order%22:[%22acvalue,asc%22],%22start%22:1,%22limit%22:1000000,%22nullatbottom%22:1}"),
                                     stringsAsFactors = FALSE, encoding = "UTF-8")
        if (names(bioassay_results)[1] == "X.U.FEFF.Warning..no.hits.found.satisfying.your.input.query.criteria.") {
                return(NA) } else { return(bioassay_results) }
}

cat("Loaded:\n",
    " function name2smi(), cid2smi(), smi2name(), cid2name(), name2cid(), get_bioassay() \n")
