.readBam <- function(file, region, fragment.size, unique = TRUE,
                     paired.end = FALSE, swap.strand = FALSE, param = NULL){

        if (!file.exists(file)){stop("[ERROR] 'file' doesn't exist.")}

        if (!missing(region)){
                if (!is(region, "GRanges")){stop("[ERROR] 'region' needs to be of class GRanges.")}
        }

        if (!missing(fragment.size)){
                if ( is.numeric(fragment.size) ){
                        if (fragment.size < 0) {
                                stop("[ERROR] 'fragment.size' needs to be a positive number.")
                        }
                } else {
                        stop("[ERROR] 'fragment.size' needs to be numeric.")
                }
        }

        if (!is.logical(unique)){stop("[ERROR] 'unique' needs to be a boolean.")}

        if (!is.logical(paired.end)){stop("[ERROR] 'paired.end' needs to to be a boolean.")}

        if (!is.logical(swap.strand)){stop("[ERROR] 'swap.strand' needs to be a boolean.")}

        if (!is.null(param)){
                if (!is(param, 'ScanBamParam')){stop("[ERROR] 'param' needs to be an object of class ScanBamParam.")}
        }

        if (is.null(param)){
                if (paired.end){
                        flag <- Rsamtools::scanBamFlag(isPaired = TRUE,
                                                       isProperPair = TRUE,
                                                       hasUnmappedMate=FALSE,
                                                       isDuplicate = FALSE,
                                                       isNotPassingQualityControls = FALSE,
                                                       isUnmappedQuery = FALSE)
                } else {
                        flag <- Rsamtools::scanBamFlag(isPaired = FALSE,
                                                       isDuplicate = FALSE,
                                                       isNotPassingQualityControls = FALSE,
                                                       isUnmappedQuery = FALSE,
                                                       isSecondaryAlignment = FALSE)
                }

                # To avoid 'duplicated record selection', reduce repeated or overlapping regions.
                if (missing(region)){
                        param <- Rsamtools::ScanBamParam(flag = flag)
                } else {
                        param <- Rsamtools::ScanBamParam(which = reduce(region, ignore.strand = FALSE), flag = flag)
                }
        }

        # Retrieve reads
        reads <- granges(GenomicAlignments::readGAlignments(file, param = param))
        # Select unique
        if (unique){reads <- unique(reads)}
        # Swap strand
        if (swap.strand){strand(reads) <- ifelse(strand(reads) == "+", "-", "+")}
        # resize reads
        if (!missing(fragment.size)){
                reads <- resize(reads, width = fragment.size, fix = 'start')
        }

        return(reads)
}
