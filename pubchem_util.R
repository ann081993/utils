#### pubchem_util
# Author: Seungchan An
# Ref: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest

# function name2smi()
# returns isomeric SMILES from a chemical name using PubChem PUG REST API
name2smi <- function(name) {
        result <- NA
        data <- try(read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                                name, "/property/IsomericSMILES/csv"),
                         stringsAsFactors = F), silent = T)
        if(class(data) != "try-error") { result <- data$IsomericSMILES[1] }
        return(result)
}

# function cid2smi()
# returns isomeric SMILES from CID using PubChem PUG REST API
cid2smi <- function(cid) {
        result <- NA
        data <- try(read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                                    cid, "/property/IsomericSMILES/csv"),
                             stringsAsFactors = F), silent = T)
        if(class(data) != "try-error") { result <- data$IsomericSMILES[1] }
        return(result)
} 

# function cid2smi_bulk()
# returns isomeric SMILES from list of CIDs using PubChem PUG REST API
# handling 300 CIDs per request
# written 2021-02-26
cid2smi_bulk <- function(cid) {
        result <- NULL
        parts <- ceiling(length(cid) / 300)
        for(p in 1:parts) {
                from = (300 * (p - 1) + 1)
                to = (300 * (p - 1) + ifelse((parts == p) & length(cid) %% 300 > 0, length(cid) %% 300, 300))
                part_cid <- cid[from:to]
                part_cid <- paste(unique(part_cid), collapse = ",")
                data <- try(read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                                            part_cid, "/property/IsomericSMILES/csv"),
                                     stringsAsFactors = F), silent = T)
                result <- rbind(result, data)
                cat("... Fetching SMILES", from, "-", to, ":", p, "of", parts, "\n")
        }
        return(result)
}

# function smi2name()
# returns a name (preferred name) from SMILES using PubChem PUG REST API
smi2name <- function(smi) {
        result <- NA
        data <- try(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
                                     smi, "/synonyms/txt"), n = 1), silent = T)
        if(class(data) != "try-error") { result <- data }
        return(result)
}

# function smi2cid()
# returns PubChem CID from SMILES using PubChem PUG REST API
smi2cid <- function(smi) {
        result <- NA
        data <- try(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
                               smi, "/cids/txt"), n = 1), silent = T)
        if(class(data) != "try-error") { result <- data }
        if(length(result) > 1) result <- result[which.min(result)]
        if(is.na(result)) result <- -999
        return(result)
}

# function cid2name()
# returns a name (preferred name) from CID using PubChem PUG REST API
cid2name <- function(cid) {
        result <- NA
        data <- try(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                                     cid, "/synonyms/txt"), n = 1), silent = T)
        if(class(data) != "try-error") { result <- data }
        return(result)
}

# function cid2cas()
# returns a CAS no. from CID using PubChem PUG REST API
cid2cas <- function(cid) {
        result <- NA
        data <- try(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                                     cid, "/xrefs/RN/txt"), n = 1), silent = T)
        if(class(data) != "try-error") { result <- data }
        return(result)
}

# function name2cid()
# returns CID from a name using PubChem PUG REST API
# derived from kegg_flavonoids.R written 2019-07-16
name2cid <- function(name) {
        cid <- tryCatch({ readLines(URLencode(paste0("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                                                     name, "/cids/txt"))) },
                        error = function(e) { return(NA) }
        )
        if(length(cid) > 1) {
                cid <- cid[which.min(cid)] # this line was added to prevent get 2 or more cids
        }
        return(cid)
}

# function get_bioassay()
# returns "biological test results" from CID
# derived from kegg_flavonoids.R written 2019-07-16
get_bioassay <- function(cid) {
        bioassay_results <- read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=jsonp&query={%22download%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22cid%22:%22",
                                            cid, "%22}]},%22order%22:[%22acvalue,asc%22],%22start%22:1,%22limit%22:1000000,%22nullatbottom%22:1}"),
                                     stringsAsFactors = FALSE, encoding = "UTF-8", quote = "")
        if (names(bioassay_results)[1] == "X.U.FEFF.Warning..no.hits.found.satisfying.your.input.query.criteria.") {
                return(NA) } else { return(bioassay_results) }
}

# function get_assaysummary()
# returns assaysummary from a vector list of CIDs using PubChem PUG REST API
# written 2020-10-14
get_assaysummary <- function(cid) {
        result <- NA
        assaysummary <- try(read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                                            paste(cid, collapse = ","), "/assaysummary/CSV"),
                                     stringsAsFactors = FALSE, encoding = "UTF-8", quote = ""), silent = T)
        if(class(assaysummary) != "try-error") { result <- assaysummary }
        return(result)
}

# function cid2smi_bulk()
# returns assaysummary from a vector list of CIDs using PubChem PUG REST API
# handling 300 CIDs per request
# written 2021-02-26
get_assaysummary_bulk <- function(cid) {
        result <- NULL
        parts <- ceiling(length(cid) / 50)
        for(p in 1:parts) {
                from = (50 * (p - 1) + 1)
                to = (50 * (p - 1) + ifelse((parts == p) & length(cid) %% 50 > 0, length(cid) %% 50, 50))
                part_cid <- cid[from:to]
                part_cid <- paste(unique(part_cid), collapse = ",")
                data <- try(read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                                            part_cid, "/assaysummary/CSV"),
                                     stringsAsFactors = F, encoding = "UTF-8"), silent = T)
                if(class(data) != "try-error") result <- rbind(result, data)
                cat("... Fetching assaysummary", from, "-", to, ":", p, "of", parts, "\n")
        }
        part_cid <- cid[!(cid %in% unique(result$CID))]
        if(length(part_cid) > 0) {
                
                cat("... Retrying to fetch assaysummary for", length(part_cid), "missings\n")
                part_cid <- paste(unique(part_cid), collapse = ",")
                data <- try(read.csv(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                                            part_cid, "/assaysummary/CSV"),
                                     stringsAsFactors = F, encoding = "UTF-8"), silent = T)
                if(class(data) != "try-error") result <- rbind(result, data)
        }
        return(result)
}
# function get_similar()
# returns a vector of CIDs of similar structure from a CID using PubChem PUG REST API
# default Threshold = 99, MaxRecords = 50
# written 2020-10-14
get_similar <- function(cid, Threshold = 99, MaxRecords = 50) {
        result <- NA
        data <- try(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/",
                                     cid, "/cids/txt?Threshold=", Threshold, "&MaxRecords=", MaxRecords)), silent = T)
        if(class(data) != "try-error") { result <- data }
        return(result)
}

# function get_superstructure()
# returns a vector of CIDs with superstructure from a CID using PubChem PUG REST API
# default MaxRecords = 50
# written 2020-10-14
get_superstructure <- function(cid, MaxRecords = 50) {
        result <- NA
        data <- try(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsuperstructure/cid/",
                                     cid, "/cids/txt?StripHydrogen=true&MaxRecords=", MaxRecords)), silent = T) # a parameter StripHydrogen=true is a default in PubChem Web
        if(class(data) != "try-error") { result <- data }
        return(result)
}

# function get_substructure()
# returns a vector of CIDs with substructure from a CID using PubChem PUG REST API
# default MaxRecords = 50
# written 2020-10-14
get_substructure <- function(cid, MaxRecords = 50) {
        result <- NA
        data <- try(readLines(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/cid/",
                                     cid, "/cids/txt?StripHydrogen=true&MaxRecords=", MaxRecords)), silent = T) # a parameter StripHydrogen=true is a default in PubChem Web
        if(class(data) != "try-error") { result <- data }
        return(result)
}

# function filter_smi()
# returns filtered SMILES
filter_smi <- function(smiles) {
        ind_filter <- grepl("2H|3H|13CH|14CH|[.]", smiles)
        return(smiles[!ind_filter])
}

cat("Loaded:\n",
    " functions name2smi(), cid2smi(), cid2smi_bulk(), smi2name(), cid2name(), name2cid(), get_bioassay() \n",
    " functions get_similar(), get_superstructure(), get_substructure(), get_assaysummary() \n")
