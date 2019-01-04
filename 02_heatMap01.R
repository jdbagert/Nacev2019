#' ---
#' title: "Histone Project Propsal"
#' author: "Nicholas D. Socci"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'

#+ Load modules
require(readxl)
require(dplyr)
require(magrittr)
require(tidyr)
require(tibble)
require(stringr)

require(gplots)
require(RColorBrewer)

#+ Utility functions

extractHistoneCountsByTypeAndMutation <- function(dd) {
    dd %>%
        mutate(Protein_Change_Position=paste0(
            Reference_Amino_Acid,Protein_Start_Position_Histone_Convention)) %>%
        count(Protein_Change_Position,Curated_Main_Cancer_Type) %>%
        spread(Curated_Main_Cancer_Type,n) %>%
        mutate_all(funs(replace(., is.na(.), 0))) %>%
        as.data.frame %>%
        column_to_rownames("Protein_Change_Position") %>%
        as.matrix
}

plotHeatMap <- function(xx) {
    breaks=c(0,seq(min(xx[xx>0]),max(xx),length=len(cls2)))
    hh=heatmap.2(xx,trace="none",col=cls2,
        main=coreHistone,dendrogram="none",margins =c(12,5),
        keysize=0.75,lwid=c(1,3.5),lhei=c(1,8),breaks=breaks,
        density.info="none",cexRow=.7)
    hh
}

#' ##Summary figures
#' 1. Heatmap #1:
#'   - Datasets:  All
#'   - X-axis: ‘Cancer Type Display’
#'   - Y-axis: Histone Mutation:

mutationSheets <- list()
mutationExcelFile <- "data/18_7_13_oncohistone_Analysis_for_NS_v4.xlsx"
sheetNames <- excel_sheets(mutationExcelFile) %>%
    as.tibble %>%
    filter(grepl("Per Patient_",value)) %>%
    pull

for(i in seq(sheetNames)) {

    si=sheetNames[i]
    cat("Sheet =",si,"\n")

    mutationSheets[[si]]=read_xlsx(mutationExcelFile,sheet=si)

    dd=mutationSheets[[si]]
    coreHistone=gsub(".*_H","H",si)

    xx=extractHistoneCountsByTypeAndMutation(dd)
    tumorCounts=colSums(xx)
    xx=sweep(xx,2,tumorCounts,"/")
    colnames(xx)=paste0(colnames(xx)," (",tumorCounts,")")

    levels=13
    cls2=colorRampPalette(c("#EEEEEE","orange","red","black"))(levels)

    if(!interactive()) pdf(file=cc("heatMap_02_",coreHistone,".pdf"),height=16.5,width=12.75)

    {
        hh=plotHeatMap(xx)

        # Nsamps=6
        # topN=rownames(hh$carpet)[1:Nsamps]
        # rii=rowSums(xx[,topN])>0

        # hh=plotHeatMap(xx[rii,topN])
    }

    if(!interactive()) dev.off()

}

# na2zero<-function(x){ifelse(is.na(x),0,x)}
# xy=bind_rows(dx) %>%
#     unite(HistResidue,Histone,Mutated_Residue) %>%
#     spread(Cancer_Type_Display,Count) %>%
#     replace(is.na(.),0)

# xx=as.matrix(xy[,2:ncol(xy)])
# rownames(xx)=xy[[1]]

# levels=max(xx)+1
# cls2=colorRampPalette(c("#EEEEEE","orange","red","black"))(levels)

# hGroup=factor(sapply(strsplit(rownames(xx),"_"),"[",1))
# cls=c("red","green","blue","grey")

# if(!interactive()) pdf(file=cc("heatMap_02_","AllHistones",".pdf"),height=16.5,width=12.75)

# heatmap.2(xx,trace="none",col=cls2,dendrogram="none",margins =c(14,5),
#     main="All",breaks=seq(levels+1)-1.5,
#     keysize=0.75,lwid=c(1,3.5),lhei=c(1,8),density.info="none",labRow="",
#     Rowv=F,RowSideColors=cls[hGroup])


# y1=(rev(table(hGroup))/nrow(xx))/4
# text(.165,cumsum(rev(table(hGroup))/(1.1*nrow(xx)))-y1,rev(levels(hGroup)),xpd=T,cex=2,pos=2)


# if(!interactive()) dev.off()
