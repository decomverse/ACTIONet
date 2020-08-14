#' Pre-processed genesets from the gProfiler database (Human)
#'
#' @format:
#'   SYMBOL: A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are gene symbols)
#'   ENSEMBL: A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are ENSEMBL gene IDs)
#'
#' @source \url{https://biit.cs.ut.ee/gprofiler}
"gProfilerDB_human"


#' Pre-processed genesets from the gProfiler database (Mouse)
#'
#' @format:
#'   SYMBOL: A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are gene symbols)
#'   ENSEMBL: A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are ENSEMBL gene IDs)
#'
#' @source \url{https://biit.cs.ut.ee/gprofiler}
"gProfilerDB_mouse"


#' Pre-processed genesets from the nanoString database (Human)
#'
#' @format:
#'   A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are gene symbols)
#'
#' @source \url{https://www.nanostring.com/}
"nanoStringDB_human"


#' Pre-processed genesets from the nanoString database (Mouse)
#'
#' @format:
#'   A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are gene symbols)
#'
#' @source \url{https://www.nanostring.com/}
"nanoStringDB_mouse"


#' Pre-processed cell type markers from the CellMarker database (Human)
#'
#' @format:
#'   A list of marker genes for different cell types
#'
#' @source \url{http://biocc.hrbmu.edu.cn/CellMarker/}
"CellMarkerDB_human_markers"


#' Pre-processed cell type markers from the CellMarker database (Mouse)
#'
#' @format:
#'   A list of marker genes for different cell types
#'
#' @source \url{http://biocc.hrbmu.edu.cn/CellMarker/}
"CellMarkerDB_mouse_markers"


#' Pre-processed cell type markers from the PanglaoDB database (Human)
#'
#' @format:
#'   A list of marker genes for different cell types
#'
#' @source \url{https://panglaodb.se/}
"PanglaoDB_human_markers"



#' Pre-processed cell type markers from the PanglaoDB database (Mouse)
#'
#' @format:
#'   A list of marker genes for different cell types
#'
#' @source \url{https://panglaodb.se/}
"PanglaoDB_mouse_markers"


#' Curated cell type markers for hand-selected cell types (Human)
#'
#' @format:
#'   A list of gene symbols for different tissues/cell-types
#'
#' @source \url{http://compbio.mit.edu/ACTIONet/}
"curatedMarkers_human"



#' Collection of Disease-associated genesets
#'
#' @format:
#'   A list of genesets
#'
#' @source \url{https://amp.pharm.mssm.edu/Enrichr/}
"DiseaseDB"


#' Collection of TF-TG assocations
#'
#' @format:
#'   A list of TF-TG associations
#'
#' @source \url{https://amp.pharm.mssm.edu/chea3/}
"ChEA3plusDB"


#' Human-Mouse orthologs
#'
#' @format:
#'   data.frame with two columns corresponding to human and mouse genes, respectively
#'
#' @source \url{http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt}
"Hs2Mm"
