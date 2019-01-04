#' ---
#' title: "Histone Project Propsal"
#' author: "Nicholas D. Socci"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'

#+ Load modules
require(tidyverse)
require(readxl)
require(magrittr)
require(stringr)

require(gplots)
require(RColorBrewer)

#' ##Summary figures
#' 1. Heatmap #1:
#'   - Datasets:  All
#'   - X-axis: ‘Cancer Type Display’
#'   - Y-axis: Hierarchy:
#'      - Histone type (H2A, H2B, H3, H4 [each data set I will send will be labeled either H2A, H2B, H3, or H4] → Column ‘Mutated residue’
#'   - Range: 1-580

mutationSheets <- list()
mutationExcelFile <- "data/18_7_13_oncohistone_Analysis_for_NS_v4.xlsx"
dataSet="Per Patient TMB 10"
dataSet="Per Patient TMB 2"
sheetNames <- excel_sheets(mutationExcelFile) %>%
    as.tibble %>%
    filter(grepl(paste0(dataSet,"_"),value)) %>%
    pull

for(i in seq(sheetNames)) {

    si=sheetNames[i]
    cat("Sheet =",si,"\n")

    mutationSheets[[si]]=read_xlsx(mutationExcelFile,sheet=si)

    dd=mutationSheets[[si]]
    coreHistone=gsub(".*_H","H",si)

    xx=dd %>%
        count(Protein_Change_Histone_Convention) %>%
        arrange(desc(n)) %>%
        slice(1:20) %>%
        arrange(n)
    xb=xx$n
    names(xb)=xx$Protein_Change_Histone_Convention
    xMax=2*ceiling(max(xb)/2+.25)
    if(!interactive()) pdf(file=cc("barPlot_03_",make.names(dataSet),coreHistone,".pdf"),height=16.5,width=12.75)

    bb=barplot(xb,horiz=T,yaxt="n",main=paste0(coreHistone,"\n",dataSet),xlim=c(0,xMax))
    text(x=0,y=bb,names(xb),xpd=T,pos=2)

    if(!interactive()) dev.off()

}
