--------------------- ## POPHELPER ## ----------------------------
#pophelper to summarise faststructure data
# install dependencies
install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")
# install pophelper package from GitHub
remotes::install_github('royfrancis/pophelper')
# load library for use
library(pophelper)
library(viridis)

#read matrices
setwd("/path/to/sNMF_Q/matrices")

##Read Q files from sNMF analysis ----> clone corrected population

afiles <- c("pop.1.Q", "pop.2.Q",
            "pop.3.Q", "pop.4.Q",
            "pop.5.Q", "pop.6.Q",
            "pop.7.Q", "pop.8.Q",
            "pop.9.Q", "pop.10.Q")

afiles <- readQBasic(files = afiles)
afiles

# collate/tabulate a qlist
qlist_tab <- tabulateQ(afiles)

# summarise an output from tabulateQ()
qlist_sum <- summariseQ(qlist_tab)

## read metadata to annotate the plot
metadata <- read.csv("pop_metadata_str_clonecorrected.csv")
metadata$Year <- as.character(metadata$Year)

head(metadata)
twolabset <- metadata[,2:3]
threelabset <- metadata[,c(4,2,3)]
inds <- metadata[,5]
head(inds)
str(inds)


p1 <- plotQ(afiles[2],returnplot=T,exportplot=F,basesize=11)
print(p1$plot[[1]])

rownames(afiles[[1]]) <- inds
if(length(unique(sapply(afiles,nrow)))==1) afiles <- lapply(afiles,"rownames<-",inds)

# show row names of all runs and all samples
lapply(afiles, rownames)[1:2]

plotQ(afiles[1:10],imgoutput="join",
      clustercol = c('#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499', '#BBBBBB'),
      #clustercol = c('#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499', '#BBBBBB'),
      showindlab=T,useindlab=T,grplab=twolabset, grplabcol="black",
      subsetgrp=c("AR", "LA", "TX", "CU"),selgrp="Location",ordergrp=T,showlegend=T,legendrow = 1,
      legendkeysize = 5, legendtextsize = 5,
      height=1.9,indlabsize=NULL,indlabheight=0.1, grplabangle = 90, grplabheight = 3,  #indlabspacer=-1,
      indlabvjust = 0.8, indlabcol="black", splabcol="black", linecol = "black", pointcol = "black",
      barbordercolour="white",barbordersize=0,outputfilename="plotQ_cc_K1-K10_Plot2",imgtype="png",dpi = 600,width=22,
      exportpath=getwd())


plotQ(afiles[2:3],imgoutput="join",
      #clustercol = c("#000004FF", "#2D1160FF", "#721F81FF", "#B63679FF", "#F1605DFF", "#FEAF77FF", "#FCFDBFFF"),
      showindlab=T,useindlab=T,grplab=threelabset,
      subsetgrp=c("1993", "1999", "2000", "2001", "2002", "2004", "2005", "2008", "2009",
                  "2012", "2013", "2020", "2021", "2022"),selgrp="Year",ordergrp=T,showlegend=T,legendrow = 1,
      legendkeysize = 5, legendtextsize = 5, grplabangle = 90, grplabheight = 4,
      indlabsize=NULL,indlabheight=0.4,#indlabspacer=-1, #height=1.5,
      indlabvjust = 0.8,
      barbordercolour="white",barbordersize=0,outputfilename="plotQ_cc_K1-K3_year",imgtype="png",dpi = 500,width=15, height = 3,
      exportpath=getwd())

plotQ(afiles[2:3],imgoutput="join",
      clustercol = c('#332288', '#88CCEE', '#44AA99'),
      #clustercol = c("#000004FF", "#2D1160FF", "#721F81FF", "#B63679FF", "#F1605DFF", "#FEAF77FF", "#FCFDBFFF"),
      showindlab=T,useindlab=T,grplab=twolabset,
      subsetgrp=c("AR", "LA", "TX", "CU"),selgrp="Location",ordergrp=T,showlegend=T,legendrow = 1,
      legendkeysize = 5, legendtextsize = 5, grplabangle = 90, grplabheight = 4,
      indlabsize=NULL,indlabheight=0.4,#indlabspacer=-1, #height=1.5,
      indlabvjust = 0.8,
      barbordercolour="white",barbordersize=0,outputfilename="plotQ_cc_K1-K3_byLoc_v2",imgtype="pdf",dpi = 300,width=15, height = 3,
      exportpath=getwd())
