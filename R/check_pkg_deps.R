#' check_pkg_deps()
#' This function import the libraries that are required
#'
#' @param
#'
#' @return If required packages are not installed, return a message. If required packages are installed, return nothing
#'
#' @examples
#'   check_pkg_deps()
#'
#' @export

# dependency checking
check_pkg_deps <- function(){
# if methylKit is not installed, install methylKit
    if (!require(methylKit))
        message('the methylKit package needs to be installed')
# genomation
    if (!require(genomation))
        message('the genomation package needs to be installed')
# genomicranges
    if (!require(GenomicRanges))
        message('the GenomicRanges package needs to be installed')
# dplyr
    if (!require(dplyr))
        message('the dplyr package needs to be installed')
}
