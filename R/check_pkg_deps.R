#' check_pkg_deps()
#' This function import the libraries that are required
#'
#' @param
#'
#' @return If required packages are not installed, return a message
#'
#' @examples
#'   check_pkg_deps()
#'
#' @export

check_pkg_deps <- function(){
    # import methylKit
    if (!require(methylKit)){
        message('the methylKit package needs to be installed')
    }
    # import genomation
    if (!require(genomation)){
        message('the genomation package needs to be installed')
    }
    # import genomicranges
    if (!require(GenomicRanges)){
        message('the GenomicRanges package needs to be installed')
    }
    # import dplyr
    if (!require(dplyr)){
        message('the dplyr package needs to be installed')
    }
}
