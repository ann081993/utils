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

cat("Loaded:\n",
    " function name2smi(), cid2smi()\n")