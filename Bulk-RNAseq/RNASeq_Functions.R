# RScript for Sourcing RNAseq Functions
# 7/20/2022
version <- 1.0

# Packages
library(ggplot2) # most figures
library(ggrepel) # gene labels on scatterplot
library(VennDiagram)
library(RColorBrewer) # color palettes
library(reshape2) # data manipulation
library(dplyr) # data manipulation
library(Hmisc) # capitalization/formatting
library(pheatmap) # heatmaps
library(matrixStats) # matrix input for heatmaps

# ----------
# FUNCTIONS
# ----------

# -------------------------------------
# Reading in Input Files from Pipeline
# -------------------------------------
read.inputs <- function(projectdir="IBD-IL11colon-RNAseq",basepath="/data/cbadger/RStudio/FRI/",inputdir="Inputs/", filenames=NULL, samplefile="",projectname=NULL) {
  # --------------------------
  # Formatting the Inputs Path
  # --------------------------
  if (substr(projectdir, nchar(projectdir), nchar(projectdir)) == "/") {
    project.dir <- projectdir
  } else {
    project.dir <- paste0(projectdir, "/")
  }
  if (substr(inputdir, nchar(inputdir), nchar(inputdir)) == "/") {
    input.dir <- inputdir
  } else {
    input.dir <- paste0(inputdir, "/")
  }
  if (substr(basepath, nchar(basepath), nchar(basepath)) == "/") {
    base.path <- basepath
  } else {
    base.path <- paste0(basepath, "/")
  }
  inputs.path <- paste0(base.path, project.dir, input.dir)
  print(paste0(inputs.path, " is the set inputs path."))
  
  # -----------------------------------------------
  # Reducing the project name for object assignment
  # -----------------------------------------------
  if (is.null(projectname)) {
    project <- projectdir
  } else {
    project <- projectname
  }
  project.red <- gsub("IBD", "", project, ignore.case=T)
  project.red <- gsub("RMMH", "", project.red, ignore.case=T)
  project.red <- gsub("[[:punct:]]", "", project.red)
  project.red <- gsub("RNAseq", "", project.red, ignore.case=T)
  print(paste0(project.red, " is the reduced project name for object assignment."))
  
  # ---------------------------------
  # Setting Filenames if not provided
  # ---------------------------------
  if (is.null(filenames)) {
    filenames=c(samplefile, "rnaseq-summary_GeneDESigs.txt", "rnaseq-summary_GeneDEMerged_anno.txt", "geneanno.txt", "gene.results.merged.tpm.filtered.auto.txt", "gene.results.merged.tpm.txt")
  } else {
    filenames=c(samplefile,filenames)
  }
  
  # ----------------------
  # Reading in Input Files
  # ----------------------
  
  # Sample Annotation File
  sample.index <- grep('samp', filenames, ignore.case=T)
  if (length(filenames[sample.index]) > 1) warning('More than one sample annotation file provided. Only reading in the first.')
  if (is.na(filenames[sample.index[1]])) {
    warning("No sample annotation file found. Skipping.")
  } else {
    sample.file <- filenames[sample.index[1]]
    sample.path <- paste0(inputs.path, sample.file)
    if (file.exists(sample.path)) {
      obj.name <- paste0(project.red,".sample.anno")
      sample.data <- read.delim(file=sample.path,quote="",header=T,sep="\t", fill=T)
      eval(call("<<-", as.name(obj.name), sample.data))
      print(paste0("Read in sample annotation file ",sample.file," as ",obj.name))
    } else {
      warning(paste0("Warning! File ",sample.file," was not found."))
    }
  }
  
  # Significance File: GeneDESigs.txt
  sigs.index <- grep('sigs', filenames, ignore.case=T)
  if (length(filenames[sigs.index]) > 1) warning('More than one sigs file provided. Only reading in the first.')
  if (is.na(filenames[sigs.index[1]])) {
    warning("No significance file found. Skipping.")
  } else {
    sigs.file <- filenames[sigs.index[1]]
    sigs.path <- paste0(inputs.path, sigs.file)
    if (file.exists(sigs.path)) {
      obj.name <- paste0(project.red,".sigs")
      sigs.data <- read.delim(file=sigs.path,fill=TRUE,header=TRUE)
      sigs.data[is.na(sigs.data)]<-0
      eval(call("<<-", as.name(obj.name), sigs.data))
      print(paste0("Read in significance file ",sigs.file," as ",obj.name))
    } else {
      warning(paste0("Warning! File ",sigs.file," was not found."))
    }
    
  }
  
  # Comprehensive File: GeneDEMerged_anno.txt
  merged.index <- grep('anno.*merg|merg.*anno', filenames, ignore.case=T)
  if (length(filenames[merged.index]) > 1) warning('More than one merged annotation file provided. Only reading in the first.')
  if (is.na(filenames[merged.index[1]])) {
    warning("No merged annotation file found. Skipping.")
  } else {
    merged.file <- filenames[merged.index[1]]
    merged.path <- paste0(inputs.path, merged.file)
    if (file.exists(merged.path)) {
      obj.name <- paste0(project.red,".merged.anno")
      merged.data <- read.delim(file=merged.path,quote="",header=T,sep="\t",row.names=1)
      eval(call("<<-", as.name(obj.name), merged.data))
      print(paste0("Read in merged annotation file ", merged.file, " as ", obj.name))
    } else {
      warning(paste0("Warning! File ",merged.file," was not found."))
    }
  }
  
  # Gene Name File: geneanno.txt
  anno.index <- grep('geneanno', filenames, ignore.case=T)
  if (length(filenames[anno.index]) > 1) warning('More than one gene annotation file provided. Only reading in the first.')
  if (is.na(filenames[anno.index[1]])) {
    warning("No gene annotation file found. Skipping.")
  } else {
    anno.file <- filenames[anno.index[1]]
    anno.path <- paste0(inputs.path, anno.file)
    if (file.exists(anno.path)) {
      obj.name <- paste0(project.red,".geneanno")
      anno.data <- read.delim(file=anno.path,quote="",header=T,sep="\t", fill=T,row.names=1)
      eval(call("<<-", as.name(obj.name), anno.data))
      print(paste0("Read in gene annotation file ", anno.file, " as ", obj.name))
    } else {
      warning(paste0("Warning! File ",anno.file," was not found."))
    }
  }
  
  # TPM File Filtered: Tpm.filtered.auto.txt
  tpm.filt.index <- grep('tpm.*filt|filt.*tpm', filenames, ignore.case=T)
  if (length(filenames[tpm.filt.index]) > 1) warning('More than one filtered TPM file provided. Only reading in the first.')
  if (is.na(filenames[tpm.filt.index[1]])) {
    warning("No filtered TPM file found. Skipping.")
  } else {
    tpm.filt.file <- filenames[tpm.filt.index[1]]
    tpm.filt.path <- paste0(inputs.path, tpm.filt.file)
    if (file.exists(tpm.filt.path)) {
      obj.name <- paste0(project.red,".tpm.filtered")
      tpm.filt.data <- read.delim(file=tpm.filt.path,quote="",header=T,sep="\t", row.names=1, fill=T)
      eval(call("<<-", as.name(obj.name), tpm.filt.data))
      print(paste0("Read in filtered TPM file ", tpm.filt.file, " as ", obj.name))
    } else {
      warning(paste0("Warning! File ",tpm.filt.file," was not found."))
    }
  }
  
  # TPM File Unfiltered: Tpm.txt
  tpm.index <- grep('tpm.txt', filenames, ignore.case=T)
  if (length(filenames[tpm.index]) > 1) warning('More than one unfiltered TPM file provided. Only reading in the first.')
  if (is.na(filenames[tpm.index[1]])) {
    warning("No unfiltered TPM file found. Skipping.")
  } else {
    tpm.file <- filenames[tpm.index[1]]
    tpm.path <- paste0(inputs.path, tpm.file)
    if (file.exists(tpm.path)) {
      obj.name <- paste0(project.red,".tpm.unfiltered")
      tpm.data <- read.delim(file=tpm.path,quote="",header=T,sep="\t", row.names=1, fill=T)
      eval(call("<<-", as.name(obj.name), tpm.data))
      print(paste0("Read in unfiltered TPM file ", tpm.file, " as ", obj.name))
    } else {
      warning(paste0("Warning! File ",tpm.file," was not found."))
    }
  }
}

# ----------------
# Venn Set Inputs
# ----------------

make.venn.sets <- function(sigs=NULL, allsets=F, upsets=F, downsets=F, dir="/data/cbadger/RStudio/FRI/IBD-IL11colon-RNAseq/Venn/", files=F) {
  # -------------------
  # Checking for Inputs
  # -------------------
  if (is.null(sigs)) {
    stop("No significance data provided.")
  }
  if (allsets==F && upsets==F && downsets==F) {
    stop("No sets selected. Try changing allsets, upsets, or downsets to T.")
  }
  # -------------------------
  # Checking Directory Format
  # -------------------------
  file.path <- ""
  if (substr(dir, nchar(dir), nchar(dir)) == "/") {
    file.path <- dir
  } else {
    file.path <- paste0(dir, "/")
  }
  # --------------------
  # Formatting Sigs File
  # --------------------
  project.sigs <- deparse(substitute(sigs)) # Sigs input object as string
  project.name <- strsplit(project.sigs,"[.]")[[1]][1] # Removing .sigs for project name
  
  # --------------
  # Creating Sets
  # --------------
  if (allsets==T) {
    print("Generating Total DEG Sets for all comparisons.")
    for (comparison in colnames(sigs[,-1])) {
      # Gene Set
      sigs.rows <- sigs[sigs[,comparison]!=0,]
      sigs.set <- paste(sigs.rows$Gene, sep="")
      # Set Name
      comp.temp <- gsub("[[:punct:]]", "", comparison)
      comp.name <- gsub("vs",".",comp.temp)
      set.name <- paste0(project.name,".",comp.name,".set")
      # Global Set
      print(paste0("Created set for ",comparison," as ",set.name,"."))
      eval(call("<<-", as.name(set.name), sigs.set))
      if (files==T) {
        # Generating a File
        file.name <- paste0(file.path,set.name,".txt")
        write.table(sigs.set,file=paste0(file.name), sep=",",quote=F, row.names=F, col.names=F)
        
      }
    }
  }
  if (upsets==T) {
    print("Generating Up DEG Sets for all comparisons.")
    for (comparison in colnames(sigs[,-1])) {
      # Gene Set
      sigs.rows.up = sigs[sigs[,comparison]==1,]
      sigs.up.set <- paste(sigs.rows.up$Gene, sep="")
      # Set Name
      comp.temp <- gsub("[[:punct:]]", "", comparison)
      comp.name <- gsub("vs",".",comp.temp)
      set.name <- paste0(project.name,".",comp.name,".up.set")
      # Global Set
      print(paste0("Created set for ",comparison," as ",set.name,"."))
      eval(call("<<-", as.name(set.name), sigs.up.set))
      if (files==T) {
        # Generating a File
        file.name <- paste0(file.path,set.name,".txt")
        write.table(sigs.up.set,file=paste0(file.name), sep=",",quote=F, row.names=F, col.names=F)
        
      }
    }
  }
  if (downsets==T) {
    print("Generating Down DEG Sets for all comparisons.")
    for (comparison in colnames(sigs[,-1])) {
      # Gene Set
      sigs.rows.down = sigs[sigs[,comparison]==-1,]
      sigs.down.set <- paste(sigs.rows.down$Gene, sep="")
      # Set Name
      comp.temp <- gsub("[[:punct:]]", "", comparison)
      comp.name <- gsub("vs",".",comp.temp)
      set.name <- paste0(project.name,".",comp.name,".down.set")
      # Global Set
      print(paste0("Created set for ",comparison," as ",set.name,"."))
      eval(call("<<-", as.name(set.name), sigs.down.set))
      if (files==T) {
        # Generating a File
        file.name <- paste0(file.path,set.name,".txt")
        write.table(sigs.down.set,file=file.name, sep=",",quote=F, row.names=F, col.names=F)
        
      }
    }
  }
}

# ------------------------
# For Generating Venn File
# ------------------------
venn.file <- function(set1, set2, name1, name2, col1="#cc5d43", col2="#8177cc", filepath="/data/cbadger/RStudio/test_venn.png") {
  if (length(set2) > length(set1)) {
    rotate = 180
  }
  else {
    rotate = 0
  }
  venn.diagram(
    x = list(set1, set2),
    category.names = c(name1 , name2),
    filename = filepath,
    scaled = FALSE,
    output=TRUE,
    # Output features
    imagetype="png",
    height = 480, 
    width = 480, 
    resolution = 300,
    
    # Circles
    lwd = 2,
    lty = 1,
    col=c(col1, col2),
    fill = c(alpha(col1,0.3), alpha(col2,0.3)),
    #fill = c(col1, col2),
    rotation.degree = rotate,
    
    # Numbers
    cex = 0.8,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.5,
    cat.pos = c(0,180),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
  )
}

# --------------------------
# Generating All 2-Way Venns
# --------------------------

# ----------------------
# For Drawing Venn in R
# ----------------------
venn.draw <- function(set1, set2, name1, name2, col1="#cc5d43", col2="#8177cc") {
  if (length(set2) > length(set1)) {
    rotate = 180
  }
  else {
    rotate = 0
  }
  grid.newpage()
  output <- venn.diagram(
    x = list(set1, set2),
    category.names = c(name1 , name2),
    filename = NULL,
    scaled = FALSE,
    output=TRUE,
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    
    # Circles
    lwd = 2,
    lty = 1,
    col=c(col1, col2),
    fill = c(alpha(col1,0.3), alpha(col2,0.3)),
    #fill = c(col1, col2),
    rotation.degree = rotate,
    
    # Numbers
    cex = 0.8,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.5,
    cat.pos = c(0,180),
    cat.dist = c(0.035, 0.035),
    cat.fontfamily = "sans"
  )
  grid.draw(output)
}

# ------------------------
# 3-way Venn Diagram File
# ------------------------

# Venn Diagram PNG File
venn3.file <- function(set1, set2, set3, name1="Set 1", name2="Set 2", name3="Set 3", cols=c("#B3E2CD", "#FDCDAC", "#CBD5E8"),filepath=NULL) {
  # File Name
  if (is.null(filepath)) {
    stop("No filepath provided.")
  } else if (substr(filepath, nchar(filepath), nchar(filepath)) == "/") {
    stop("Filepath must include the file name (.png)")
  } else if (dir.exists(filepath)) {
    stop("Filepath cannot be a directory. Please specify the full path and file name (.png)")
  } else {
    print(paste("Path and file name set as ", filepath))
  }
  # Circle Total Areas
  area.1 <- length(set1)
  area.2 <- length(set2)
  area.3 <- length(set3)
  # Intersections
  intersect.1.2 <- intersect(set1, set2)
  intersect.2.3 <- intersect(set2, set3)
  intersect.1.3 <- intersect(set1, set3)
  intersect.1.2.3 <- intersect(intersect.1.2, set3)
  # Intersections to Count
  n.1.2 <- length(intersect.1.2)
  n.2.3 <- length(intersect.2.3)
  n.1.3 <- length(intersect.1.3)
  n.1.2.3 <- length(intersect.1.2.3)
  # Venn Plot PNG
  png(filename=filepath)
  venn.plot <- draw.triple.venn(
    area1 = area.1,
    area2 = area.2,
    area3 = area.3,
    n12 = n.1.2,
    n23 = n.2.3,
    n13 = n.1.3,
    n123 = n.1.2.3,
    category = c(name1, name2, name3),
    cat.col = c("black","black","black"),
    
    # Circles
    euler.d = FALSE,
    scaled = FALSE,
    lwd = 2,
    lty = 'blank',
    fill = cols,
    
    # Numbers
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1.1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 20, 180),
    cat.dist = c(0.065, 0.065, 0.035),
    cat.fontfamily = "sans",
    rotation = 1
  )
  dev.off()
}

# --------------------------
# Drawing 3-way Venn Diagram
# --------------------------

venn3.draw <- function(set1, set2, set3, name1="Set 1", name2="Set 2", name3="Set 3", cols=c("#B3E2CD", "#FDCDAC", "#CBD5E8")) {
  # Circle Total Areas
  area.1 <- length(set1)
  area.2 <- length(set2)
  area.3 <- length(set3)
  # Intersections
  intersect.1.2 <- intersect(set1, set2)
  intersect.2.3 <- intersect(set2, set3)
  intersect.1.3 <- intersect(set1, set3)
  intersect.1.2.3 <- intersect(intersect.1.2, set3)
  # Intersections to Count
  n.1.2 <- length(intersect.1.2)
  n.2.3 <- length(intersect.2.3)
  n.1.3 <- length(intersect.1.3)
  n.1.2.3 <- length(intersect.1.2.3)
  # Venn Plot
  output <- draw.triple.venn(
    area1 = area.1,
    area2 = area.2,
    area3 = area.3,
    n12 = n.1.2,
    n23 = n.2.3,
    n13 = n.1.3,
    n123 = n.1.2.3,
    category = c(name1, name2, name3),
    cat.col = c("black","black","black"),
    euler.d=FALSE,
    scaled=FALSE,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = cols,
    
    # Numbers
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1.1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 20, 180),
    cat.dist = c(0.065, 0.065, 0.035),
    cat.fontfamily = "sans",
    rotation = 1
  )
  grid.draw(output)
  grid.newpage()
}

# -------------------
# Scatterplot Inputs
# -------------------
make.scatter.input <- function(data=NULL,project="",count=T,dir="",files=F) {
  # ---------------
  # Checking Inputs
  # ---------------
  if (is.null(data)) {
    stop("No data provided.")
  }
  # -------------------------
  # Checking Directory Format
  # -------------------------
  #file.path <- "" # Default path, for current directory
  #if (is.empty(dir)) {} # needs rapport package, find base-R method?
  #if (dir=="") {} # This should work
  #if (substr(dir, nchar(dir), nchar(dir)) == "/") {
  #file.path <- dir
  #} else {
  #file.path <- paste0(dir, "/")
  #}
  # ---------------
  # Picking Columns
  # ---------------
  data.logs.sigs <- data[,grepl("Log2|Signi",names(data), ignore.case=T)] # FC and Sigs
  data.logs <- data.logs.sigs[,!grepl("Signi",names(data.logs.sigs), ignore.case=T)] # FC
  data.sigs <- data.logs.sigs[,grepl("Signi",names(data.logs.sigs), ignore.case=T)] # Sigs
  data.gene.col <- grep("gene_name", names(data), value = TRUE)[1] # 1st geneID colname
  if (is.na(data.gene.col)) {
    stop("Input data needs to have a gene_name column.")
  }
  data.genes <- cbind(rownames(data), data[,data.gene.col]) # EnsID and GeneID
  # -----------------
  # Prepping for Loop
  # -----------------
  colcount <- length(colnames(data.sigs)) # Ex. 7
  colcount.adj <- length(colnames(data.sigs)) - 1 # Ex. 6
  # ------------------------
  # Generating Scatter Input
  # ------------------------
  for (i in 1:colcount.adj) { # Ex. 1 to 6
    ind <- i+1
    for (j in ind:colcount) { # Ex. 2 to 7
      # -------
      # Marks
      # -------
      data.sigset <- cbind(data.sigs[i],data.sigs[j])
      data.sigset[is.na(data.sigset)]<-0
      data.marks<-apply(data.sigset,1,function(x) {
        if(x[1]==0 & x[2]==0) {
          "NS"
        } else if(x[1]!=0 & x[2]==0) {
          "X Only"
        } else if(x[1]==0 & x[2]!=0) {
          "Y Only"
        } else if(x[1]!=0 & x[2]!=0 & x[1]==x[2]) {
          "X,Y Both"
        } else if (x[1]!=0 & x[2]!=0 & x[1]!=x[2]) {
          "X,Y Reverse"
        } else {
          "Others"
        }
      })
      # Dataset
      output.data <- cbind(data.logs[,i],data.logs[,j],data.marks) # Data and Mark
      rownames(output.data)<-rownames(data) # Ensembl IDs
      output.data <- output.data[output.data[,3]!="NS",] # Sig rows only
      rownames(output.data)<-make.names(data[rownames(output.data),data.gene.col],unique = T) # GeneID
      if (count) {
        # Counting Occurences of each Mark
        x.count <- length(which(output.data[,"data.marks"]=="X Only"))
        y.count <- length(which(output.data[,"data.marks"]=="Y Only"))
        xy.both.count <- length(which(output.data[,"data.marks"]=="X,Y Both"))
        xy.rev.count <- length(which(output.data[,"data.marks"]=="X,Y Reverse"))
        # Replacing Mark Values to include Count
        output.data[,"data.marks"] <- sub("X Only", paste0("X Only (",x.count,")"), output.data[,"data.marks"])
        output.data[,"data.marks"] <- sub("Y Only", paste0("Y Only (",y.count,")"), output.data[,"data.marks"])
        output.data[,"data.marks"] <- sub("X,Y Both", paste0("X,Y Both (",xy.both.count,")"), output.data[,"data.marks"])
        output.data[,"data.marks"] <- sub("X,Y Reverse", paste0("X,Y Reverse (",xy.rev.count,")"), output.data[,"data.marks"])
      }
      # Naming and Export
      comp1.temp <- gsub("[[:punct:]]", "", colnames(data.logs)[i])
      comp1.temp2 <- gsub("Log2|FC|fold|change|MLE|Group", "", comp1.temp, ignore.case=TRUE)
      comp1.name <- gsub("vs", ".", comp1.temp2, ignore.case=TRUE)
      comp2.temp <- gsub("[[:punct:]]", "", colnames(data.logs)[j])
      comp2.temp2 <- gsub("Log2|FC|fold|change|MLE|Group", "", comp2.temp, ignore.case=TRUE)
      comp2.name <- gsub("vs", ".", comp2.temp2, ignore.case=TRUE)
      output.name <- paste0(project,".",comp1.name,".vs.",comp2.name,".scatterset")
      # Global Set
      print(paste0("Created scatter set for ",comp1.name," and ",comp2.name," as ",output.name,"."))
      eval(call("<<-", as.name(output.name), output.data))
      # File Output
      #if (files==T) {
      # Generating a File
      #file.name <- paste0(file.path,set.name,".txt")
      #write.table(sigs.set,file=paste0(file.name), sep=",",quote=F, row.names=F, col.names=F)
    }
  }
}

# ----------------
# Scatterplot
# ----------------
scatter_plot<-function(xval,yval,colorby=NULL,shapeby=NULL,sizeby=NULL,labels=NULL,xlim=NULL,ylim=NULL,xlab="",ylab="",na.val=0,na.rm=T,main="",colors=NULL,regress=T,cor=T,show_labels=NULL,textrepel=T,repelcol=NULL,repelsize=4,maxoverlaps=500,zeroline=NULL,diagline=NULL,toplabel=NULL,repeldir="both",repelxlim=c(NA,NA),repelylim=c(NA,NA)) {
  
  #remove/replace NA values
  
  if(na.rm) {
    val.sel<-!is.na(xval) & !is.na(yval)
    xval<-xval[val.sel]
    yval<-yval[val.sel]
    
    if(!is.null(colorby)) {
      colorby<-colorby[val.sel]
    }
    
    if(!is.null(shapeby)) {
      shapeby<-shapeby[val.sel]
    }      
    
    if(!is.null(sizeby)) {
      sizeby<-as.numeric(sizeby[val.sel])
    }      
    
    if(!is.null(labels)) {
      labels<-labels[val.sel]
    }
    
  } else {
    xval[is.na(xval)]=na.val
    yval[is.na(xval)]=na.val
  }
  
  #convert x/y to numeric value
  xval<-as.numeric(xval)
  yval<-as.numeric(yval)
  
  #create df for plot
  #data<-data.frame(xval=as.numeric(xval),yval=as.numeric(yval))
  data<-data.frame(xval=as.numeric(xval),yval=as.numeric(yval),magnitude=abs(as.numeric(xval)*as.numeric(yval)),label=labels)
  
  #add color,shape, size inforamtion
  if(!is.null(colorby)) {
    data$colorby=colorby
  }
  if(!is.null(shapeby)) {
    data$shapeby=shapeby
  }
  if(!is.null(sizeby)) {
    data$sizeby=abs(as.numeric(sizeby))
  }
  if(!is.null(labels)) {
    #hide labels not in show_labels
    if(!is.null(show_labels)) {
      if(!is.null(toplabel)) {
        top_labels <- head(data[order(data$magnitude, decreasing = TRUE), ], toplabel)[,"label"]
        if(length(show_labels) > 1) { 
          show_labels <- show_labels[show_labels %in% top_labels]
        } else {
          show_labels <- top_labels
        }
      }
      labels[!(labels %in% show_labels)]=""
    }
    
    data$label=labels
  }   
  
  # Plot
  plot <- ggplot(data, aes(x = xval, y =yval,color=colorby,shape=shapeby,size=sizeby))
  
  # Zero Line
  if(!is.null(zeroline)){
    plot <- plot + geom_vline(xintercept=c(0), linetype = "dotted", color="red") + geom_hline(yintercept=c(0), linetype = "dotted", color="red")
  }
  # Diag Line
  if(!is.null(diagline)) {
    if(grepl("eg", diagline, ignore.case = TRUE)) {
      slope.sign <- -1
    } else if(grepl("os", diagline, ignore.case = TRUE)) {
      slope.sign <- 1
    } else {
      print("Diagonal Line not recognized. Use \"pos\" for slope 1 or \"neg\" for slope -1.")
    }
  }
  if(exists("slope.sign")) {
    plot <- plot + geom_abline(aes(slope=slope.sign,intercept=0,color="gray"), show.legend=FALSE, linetype="dashed", size=1)    
  }

  
  plot <- plot + geom_point()
  
  #choose color by names, need to be expanded by customized selection
  if(!is.null(colors)) {
    #colors now need to be named vector
    #names(colors)<-levels(as.factor(colorby))
    plot <- plot + scale_colour_manual(name ="colorby",values=colors)
  }
  
  #cusomized shape, need to be expanded by customized selection
  #to be implemented
  if(!is.null(shapeby)) {
    plot <- plot + scale_shape_discrete(name="shapeby")
  }
  
  #customized size, need to be expanded by customized selection
  if(!is.null(sizeby)) {
    plot <- plot + scale_size_continuous(name="sizeby")
  }
  
  #xlim & ylim    
  if(!is.null(xlim)) {
    plot <-plot + xlim(xlim)
  } else {
    xlim<-range(xval)
  }
  
  if(!is.null(ylim)) {
    plot <-plot + ylim(ylim)
  } else {
    ylim<-range(yval)
  }
  
  
  #title, background, xlab, ylab
  plot <- plot + ggtitle(label = main) +
    xlab(xlab) + # Change X-Axis label
    ylab(ylab) + # Change Y-Axis label
    theme_bw(base_size = 14) + # change overall theme
    #theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #plot.title = element_text(size = 10, face = "bold",hjust = 0.5),axis.title=element_text(size=8)) # change the legend
    theme(legend.position = "right",panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold",hjust = 0.5),axis.title=element_text(size=8)) + # change the legend
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
    theme(plot.title = element_text(size = rel(1.4)))
  
  #regression line
  if(regress) {
    #below line doesn't work well with multiple color/shape setting
    #plot <-plot +geom_smooth(data=data,method = "lm", se = F,aes(x = xval, y =yval),color="black")
    data.lm <- lm(yval ~ xval, data)
    plot<- plot + geom_abline(slope = coef(data.lm)[[2]], intercept = coef(data.lm)[[1]])    
  }
  
  #correlation
  if(cor) {
    pearsoncor<-round(cor(xval,yval,method="pearson"),3)
    spearmancor<-round(cor(xval,yval,method="spearman"),3)
    
    xcoord<-ycoord<-0
    
    if(pearsoncor>0) {
      xcoord<-(xlim[2]-xlim[1])*0.2+xlim[1]
      ycoord<--(ylim[2]-ylim[1])*0.1+ylim[2]
    } else {
      xcoord<--(xlim[2]-xlim[1])*0.2+xlim[2]
      ycoord<--(ylim[2]-ylim[1])*0.1+ylim[2]       
    }
    
    plot <-plot + annotate("text",x=xcoord,y=ycoord,label=paste("Pearson Cor:",pearsoncor,"\n","Spearman Cor:",spearmancor,sep=""),colour = "black", size = 5)
  }
  
  if(!is.null(show_labels)) {
    #data.show<-as.data.frame(cbind(xval,yval))
    #rownames(data.show)<-labels
    #data.show$label=labels
    #data.show<-data.show[show_labels,]
    #data.show<-data.show[!is.na(data.show[,1]),]
    
    if(textrepel) {
      if(!is.null(repelcol)) {
        plot<-plot+geom_text_repel(aes(label=label),
                                   size=repelsize,color=repelcol,
                                   hjust=0, vjust=0, 
                                   max.overlaps=maxoverlaps, direction = repeldir,
                                   min.segment.length = 0,
                                   xlim = repelxlim, ylim = repelylim)
      } else {
        plot<-plot+geom_text_repel(aes(label=label),
                                   size=repelsize, hjust=0, vjust=0, 
                                   max.overlaps=maxoverlaps, direction = repeldir,
                                   min.segment.length = 0,
                                   xlim = repelxlim, ylim = repelylim)
      }
    } else {
      #plot<-plot+geom_text(aes(label=label),color = "black",hjust=0, vjust=0)
      plot<-plot+geom_text(aes(label=label),color="black",hjust=0, vjust=0, nudge_x = 0.1, nudge_y = 0.1, size=4)
    }
    
  }
  
  plot
  
}

# ---------------
# Boxplot
# ---------------
box.draw <- function(geneanno=NULL, tpm=NULL, samples=NULL, gene="Il11", color=NULL, reorder=F, comp.order=NULL, title=gene, legend=T, ylims=NULL, transform=T, groups=NULL) {
  
  # Original default for color was "#8675ca"
  
  # Checking Inputs
  if (is.null(geneanno)) {
    stop("No Gene Annotation file provided. Use geneanno=filename")
  }
  if (is.null(tpm)) {
    stop("No TPM file provided. Use tpm=filename")
  }
  if (is.null(samples)) {
    stop("No Samples file provided. Use samples=filename")
  }
  
  # Converting Gene Name to Ensembl ID
  ensembl <- rownames(geneanno[geneanno$gene_name==gene,])
  if (length(ensembl) == 0) {
    stop("Gene name provided does not match any row in geneanno.")
  } else {
    print(paste0("Found Ensembl ID ",ensembl," for gene ", gene, "."))
  }
  
  # TPM Row for our Gene
  tpm.red <- tpm[rownames(tpm)==ensembl,]
  if (transform==T) {
    tpm.log <- log2(tpm.red+0.01)
  } else {
    tpm.log <- tpm.red
  }
  tpm.melt <- melt(tpm.log)
  # Assigning Groups from Sample Annotation
  tpm.melt$variable <- as.factor(samples$Group)
  colnames(tpm.melt) <- c("Group","TPM")
  
  # Group Selection, Group Order
  if (!is.null(groups)) {
    tpm.melt <- tpm.melt[tpm.melt$Group %in% groups,]
    tpm.melt$Group <- factor(tpm.melt$Group, levels=groups)
    tpm.melt <- tpm.melt[order(tpm.melt$Group),]
  }
  
  # Single Fill Color
  if (!is.null(color)) {
    plot <- ggplot(tpm.melt, aes(x=Group, y=TPM, fill=Group)) + geom_boxplot(width = 0.5, fill=color)
  } else {
    plot <- ggplot(tpm.melt, aes(x=Group, y=TPM, fill=Group)) + geom_boxplot(width = 0.5)
  }
  
  # Y Label
  if (transform==T) {
    plot <- plot + ylab("Log2(TPM+0.01)")
  } else {
    plot <- plot + ylab("TPM")
  }
  
  # Legend
  if (legend==F) {
    plot <- plot + theme(axis.text.x = element_text(size=12,angle=90), legend.position="none")
  } else {
    plot <- plot + theme(axis.text.x = element_text(size=12,angle=90))
  }
  
  # Y Limits
  if (!is.null(ylims)) {
    plot <- plot + coord_cartesian(ylim = ylims)
  }
  plot <- plot + ggtitle(toupper(title)) + geom_jitter(width=0.1, show.legend=FALSE)
  plot
}

# -----------
# Barplot
# -----------

# ----------------------------
# Barplots with Sortby = Comp1
# ----------------------------
barplot.sorted <- function(basepath="/data/cbadger/RStudio/FRI/RMMH-PGF2a-N582707-RNAseq/Gstools/Metabase/", in1="PGF2a_b1_48h_vs_Control_48h_gs-fisher.txt", in2="N582707_48h_vs_PGF2a_b2_48h_gs-fisher.txt", out1=NULL, out2=NULL, qvalue=0.05, top=10, db="") {
  
  library(Cairo)
  library(RColorBrewer)
  library(ggplot2)
  library(scales)
  
  # COMPARISON NAMES
  name1 <- gsub("_gs-fisher.txt", "", in1)
  name2 <- gsub("_gs-fisher.txt", "", in2)
  
  # AUTOMATIC OUTPUT NAMES
  if (is.null(out1)) {
    out1 <- paste0(name1,"_",db,"_barplot_bhp.png")
    out1 <- gsub("__", "_", out1)
  }
  if (is.null(out2)) {
    out2 <- paste0(name2,"_",db,"_barplot_bhp.png")
    out2 <- gsub("__", "_", out2)
  }
  
  # SIGNIFICANCE
  siglevel <- as.numeric(qvalue)
  
  # READING IN DATA
  data1<-read.table(paste0(basepath,in1),header=T,row.names=1,sep="\t",quote="",comment.char="")
  data2<-read.table(paste0(basepath,in2),header=T,row.names=1,sep="\t",quote="",comment.char="")
  # REMOVING FIRST ROW
  data1<-data1[-1,]
  data2<-data2[-1,]
  # ROW NAMES
  rownames(data1) <- gsub("REACTOME_", "", rownames(data1), ignore.case=TRUE) # Removing REACTOME Tag
  rownames(data2) <- gsub("REACTOME_", "", rownames(data2), ignore.case=TRUE)
  rownames(data1) <- gsub("GOBP_", "", rownames(data1), ignore.case=TRUE) # Removing GOBP Tag
  rownames(data2) <- gsub("GOBP_", "", rownames(data2), ignore.case=TRUE)
  if (tolower(db)=="reactome" | grepl("go",tolower(db), fixed=TRUE)) { # Rownames to Lowercase
    rownames(data1) <- capitalize(tolower(rownames(data1)))
    rownames(data2) <- capitalize(tolower(rownames(data2)))
  }
  rownames(data1) <- make.names(substr(rownames(data1),1,50),unique=T) # Reducing rowname length
  rownames(data2) <- make.names(substr(rownames(data2),1,50),unique=T)
  
  # VALUES 1 (in1)
  data1.z<-as.numeric(as.vector(unlist(data1[,4])))
  data1.p<-as.numeric(as.vector(unlist(data1[,5])))
  data1.p[is.na(data1.p)]<-1
  data1.p<- -log10(data1.p+1e-32) #conversion
  data1.q<-as.numeric(as.vector(unlist(data1[,7])))
  data1.q[is.na(data1.q)]<-1
  data1.q<- -log10(data1.q+1e-32) #conversion
  
  # VALUES 2 (in2)
  data2 <- data2[match(rownames(data1), rownames(data2)), ] # Using same comparisons order
  data2.z<-as.numeric(as.vector(unlist(data2[,4])))
  data2.p<-as.numeric(as.vector(unlist(data2[,5])))
  data2.p[is.na(data2.p)]<-1
  data2.p<- -log10(data2.p+1e-32) #conversion
  data2.q<-as.numeric(as.vector(unlist(data2[,7])))
  data2.q[is.na(data2.q)]<-1
  data2.q<- -log10(data2.q+1e-32) #conversion
  
  # HOW MANY PATHS SHOWN (min 5)
  topnum <- top
  topnum.q=max(min(topnum,sum(data1.q > -log10(siglevel))),5)
  # sum() = how many qvalues are above significance (Ex. 8)
  # min(topnum,sum) = lower of the two (Ex. min(10,8) = 8) incase there aren't 10 significant pathways
  # max(min,5) = showing at least 5 pathways, or more depending on significance/topnum
  
  # BARPLOT FUNCTION
  # Precaulculations
  p1=data1.q[1:topnum.q]
  p2=data2.q[1:topnum.q]
  pmax=max(p1,p2)
  # Function
  gs_bar<-function(z,p,gs,zname,pname,signame="p",name,siglevel=0.05) {
    
    # DATA FRAME
    gs<-make.names(gs,unique=T)
    
    data.df<-data.frame(factor(gs,levels=rev(gs)),		
                        as.numeric(unlist(z)),
                        as.numeric(unlist(p)))
    
    colnames(data.df)<-c("Gene Set",zname, pname)			  
    
    # FINALIZED INPUT NAMES
    zname.rev<-paste("`",zname,"`",sep="")
    pname.rev<-paste("`",pname,"`",sep="")
    
    plot1<-ggplot(data.df, aes_string(x="`Gene Set`", y=pname.rev, fill=zname.rev)) +geom_bar(stat='identity')+theme_classic() +scale_fill_gradient2(low = "blue",  mid="grey",high = "red", space = "Lab", limit = c(-4, 4),oob=squish) + coord_flip()+geom_hline(yintercept = -log10(siglevel), linetype="dashed",color="black") + annotate(geom="text", x=1, y=-log10(0.05), label=paste(signame,"<",siglevel,sep=""),color="black",size=4,hjust = 0) +theme(axis.text.x = element_text(angle = 0,size = 12 ),axis.text.y = element_text(angle = 0,size = 12 )) +ylab(paste0(pname,"\n",name)) +ggtitle(db) +ylim(0,pmax)
    
    print(plot1)  
  }
  # BARPLOT FILES
  # Comparison 1 Figure Names
  figure1=paste0(basepath,out1)
  # Comparison 1 Image (Auto Calculate Length)
  CairoPNG(file=figure1,res = 300,width = 3.5+floor(max(nchar(rownames(data1)[1:topnum.q])/10)),height = topnum.q*0.3,units = "in")
  # Comparison 2 Plotting
  gs_bar(z=data1.z[1:topnum.q],p=data1.q[1:topnum.q],gs=rownames(data1)[1:topnum.q],zname="ZScore",pname="-Log10BHP",siglevel=siglevel,signame="BHP",name=name1)
  dev.off()
  
  # Comparison 2 Figure Names
  figure2=paste0(basepath,out2)
  # Comparison 2 Image (Auto Calculate Length)
  CairoPNG(file=figure2,res = 300,width = 3.5+floor(max(nchar(rownames(data2)[1:topnum.q])/10)),height = topnum.q*0.3,units = "in")
  # Comparison 2 Plotting
  gs_bar(z=data2.z[1:topnum.q],p=data2.q[1:topnum.q],gs=rownames(data2)[1:topnum.q],zname="ZScore",pname="-Log10BHP",siglevel=siglevel,signame="BHP",name=name2)
  dev.off()
}

# -------------------
# TPM Gene Annotated
# -------------------
# Useful for Heatmap and Lineplot
tpm.anno <- function(tpm = NULL, geneanno = NULL) {
  geneanno.red <- geneanno[,c(1,2)] # Reduced gene annotation
  tpm.anno <- merge(geneanno.red, tpm, by = 0) # Annotating TPM
  tpm.anno <- tpm.anno[!duplicated(tpm.anno$gene_name), ] # Removing duplicate gene names
  rownames(tpm.anno) <- tpm.anno$gene_name # Genes as Rownames
  tpm.anno$Row.names <- NULL # Formatting
  tpm.anno$gene_name <- NULL # Formatting
  tpm.anno$description <- NULL # Formatting
  return(tpm.anno)
}

# -------------------
# Lineplot Function
# -------------------
lineplot <- function(tpm.data=NULL, sample.anno=NULL, gene=NULL, log2=T, main=NULL, ylims=NULL) {
  # NOTE: This function is written specifically for TPM input and Treatment and Time annotation
  
  # Checking Inputs
  if(is.null(tpm.data)) {
    stop("No TPM data provided (tpm.data).")
  }
  if(is.null(sample.anno)) {
    stop("No Sample Annotation provided (sample.anno)")
  }
  if(is.null(gene)) {
    stop("No Gene specified (gene)")
  }
  
  # Colors
  cols <- c("#b4943e","#727cce","#60a862","#c15ca5","#cb5a4c")
  
  # Gene TPM
  tpm <- tpm.data[gene,]
  # Log Transformation
  if(log2==T) {
    tpm <- log2(tpm+0.01)
    print("Log2 Transformed TPM Data.")
    is.na(tpm) <- sapply(tpm, is.infinite) # Only used if transformation is Log2(x+0)
    tpm <- tpm[,colSums(is.na(tpm))==0] # Only used if transformation is Log2(x+0)
  }
  # Formatting
  tpm.melt <- melt(tpm)
  # TPM for Selected Samples
  tpm.main <- merge(tpm.melt, sample.anno[,c("Sample", "Treatment", "Time")], by = 1)
  colnames(tpm.main)[1] <- "Sample"
  colnames(tpm.main)[2] <- "TPM"
  
  # Mean/Sd Calculation
    # Data is Gene TPM data frame
    # Varname is TPM
    # Groupnames is Treatment and Time
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  # Calculating Mean and Sds
  tpm.main.sd <- data_summary(tpm.main, varname="TPM", 
                              groupnames=c("Treatment", "Time"))
  # Ordering by Time
  tpm.final <- tpm.main.sd[order(tpm.main.sd$Time), ]
  # Convert Time to a factor variable
  tpm.final$Time <- as.factor(tpm.final$Time)
  
  # Line Plot
  p <- ggplot(tpm.final, aes(x=Time, y=TPM, group=Treatment, color=Treatment)) + 
    geom_line(size=1.5) +
    geom_errorbar(aes(ymin=TPM-sd, ymax=TPM+sd),
                  width=.1,position=position_dodge(0.05),
                  color="black") + 
    geom_point(size=3)
  
  # Axis Labels
  if(log2==T) {
    ylabel <- "Log2(TPM+0.01)"
  } else {
    ylabel <- "TPM"
  }
  if(!is.null(main)) {
    tlabel <- main
  } else {
    tlabel <- gene
  }
  p <- p + labs(title=tlabel, x="Time (hrs)", y = ylabel)
  
  # Axis Scale
  if(!is.null(ylims)) {
    if(length(ylims)!=2) {
      stop("Must specify upper and lower ylimits (Ex. c(-5,5)).")
    }
    p <- p + ylim(ylims)
  }
  
  # Theme, Formatting, Legend, Colors
  p <- p + theme_classic() + theme(axis.title = element_text(size = 26), plot.title = element_text(size = 22), axis.text = element_text(size = 18), legend.text = element_text(size = 18), legend.title = element_text(size = 18)) +
    scale_color_manual(values=cols[1:length(unique(tpm.final$Treatment))])
  return(p)
}

# -------------------
# Custom Scatterplot 
# -------------------
# Function for creating Jej 604 Indo Scatterplot
scatter_plot.2<-function(xval,yval,colorby=NULL,shapeby=NULL,sizeby=NULL,labels=NULL,bold=F,xlim=NULL,ylim=NULL,xlab="",ylab="",na.val=0,na.rm=T,main="",colors=NULL,shapes=NULL,regress=T,cor=T,show_labels1=NULL,show_labels2=NULL,textrepel=T,repelcol=NULL,repelsize=4,maxoverlaps=500,zeroline=NULL,diagline=NULL,toplabel=NULL,repeldir="both",repelxlim=c(NA,NA),repelylim=c(NA,NA)) {
  
  #remove/replace NA values
  
  if(na.rm) {
    val.sel<-!is.na(xval) & !is.na(yval)
    xval<-xval[val.sel]
    yval<-yval[val.sel]
    
    if(!is.null(colorby)) {
      colorby<-colorby[val.sel]
    }
    
    if(!is.null(shapeby)) {
      shapeby<-shapeby[val.sel]
    }      
    
    if(!is.null(sizeby)) {
      sizeby<-as.numeric(sizeby[val.sel])
    }      
    
    if(!is.null(labels)) {
      labels<-labels[val.sel]
    }
    
  } else {
    xval[is.na(xval)]=na.val
    yval[is.na(xval)]=na.val
  }
  
  #convert x/y to numeric value
  xval<-as.numeric(xval)
  yval<-as.numeric(yval)
  
  #create df for plot
  #data<-data.frame(xval=as.numeric(xval),yval=as.numeric(yval))
  data<-data.frame(xval=as.numeric(xval),yval=as.numeric(yval),magnitude=abs(as.numeric(xval)*as.numeric(yval)),label=labels)
  #data2<-data.frame(xval=data$xval[shapeby=="BothTissue"],yval=data$yval[shapeby=="BothTissue"],label=labels[shapeby=="BothTissue"])
  
  #add color,shape, size inforamtion
  if(!is.null(colorby)) {
    data$colorby=colorby
  }
  if(!is.null(shapeby)) {
    data$shapeby=shapeby
  }
  if(!is.null(sizeby)) {
    data$sizeby=abs(as.numeric(sizeby))
  }
  if(!is.null(labels)) {
    #hide labels not in show_labels
    if(!is.null(show_labels1)) {
      if(!is.null(show_labels2)) {
        labels1 <- labels
        labels2 <- labels
        labels[!(labels %in% c(show_labels1,show_labels2))]=""
        labels1[!(labels1 %in% show_labels1)]=""
        labels2[!(labels2 %in% show_labels2)]=""
        data$label1 = labels1
        data$label2 = labels2
      }
    }
    data$label = labels
  }   
  
  # Plot
  plot <- ggplot(data, aes(x = xval, y =yval,color=colorby,shape=shapeby,size=sizeby))
  
  # Zero Line
  if(!is.null(zeroline)){
    plot <- plot + geom_vline(xintercept=c(0), linetype = "dotted", color="red") + geom_hline(yintercept=c(0), linetype = "dotted", color="red")
  }
  # Diag Line
  if(!is.null(diagline)) {
    if(grepl("eg", diagline, ignore.case = TRUE)) {
      slope.sign <- -1
    } else if(grepl("os", diagline, ignore.case = TRUE)) {
      slope.sign <- 1
    } else {
      print("Diagonal Line not recognized. Use \"pos\" for slope 1 or \"neg\" for slope -1.")
    }
  }
  if(exists("slope.sign")) {
    plot <- plot + geom_abline(aes(slope=slope.sign,intercept=0,color="gray"), show.legend=FALSE, linetype="dashed", size=1)    
  }
  
  
  plot <- plot + geom_point()
  
  #choose color by names, need to be expanded by customized selection
  if(!is.null(colors)) {
    #colors now need to be named vector
    #names(colors)<-levels(as.factor(colorby))
    plot <- plot + scale_colour_manual(name ="colorby",values=colors)
  }
  
  #cusomized shape, need to be expanded by customized selection
  #to be implemented
  if(!is.null(shapeby)) {
    plot <- plot + scale_shape_manual(name="shapeby",values=shapes)
    #points(x=data2$xval,y=data2$yval, pch = "O")
    #points(x=data$xval[shapeby=="BothTissue"],y=data$yval[shapeby=="BothTissue"], pch = "O")
  }
  
  #customized size, need to be expanded by customized selection
  if(!is.null(sizeby)) {
    plot <- plot + scale_size_continuous(name="sizeby")
  }
  
  #xlim & ylim    
  if(!is.null(xlim)) {
    plot <-plot + xlim(xlim)
  } else {
    xlim<-range(xval)
  }
  
  if(!is.null(ylim)) {
    plot <-plot + ylim(ylim)
  } else {
    ylim<-range(yval)
  }
  
  
  #title, background, xlab, ylab
  plot <- plot + ggtitle(label = main) +
    xlab(xlab) + # Change X-Axis label
    ylab(ylab) + # Change Y-Axis label
    theme_bw(base_size = 14) + # change overall theme
    #theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #plot.title = element_text(size = 10, face = "bold",hjust = 0.5),axis.title=element_text(size=8)) # change the legend
    theme(legend.position = "right",panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold",hjust = 0.5),axis.title=element_text(size=8)) + # change the legend
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
    theme(plot.title = element_text(size = rel(1.4)))
  
  #regression line
  if(regress) {
    #below line doesn't work well with multiple color/shape setting
    #plot <-plot +geom_smooth(data=data,method = "lm", se = F,aes(x = xval, y =yval),color="black")
    data.lm <- lm(yval ~ xval, data)
    plot<- plot + geom_abline(slope = coef(data.lm)[[2]], intercept = coef(data.lm)[[1]])    
  }
  
  #correlation
  if(cor) {
    pearsoncor<-round(cor(xval,yval,method="pearson"),3)
    spearmancor<-round(cor(xval,yval,method="spearman"),3)
    
    xcoord<-ycoord<-0
    
    if(pearsoncor>0) {
      xcoord<-(xlim[2]-xlim[1])*0.2+xlim[1]
      ycoord<--(ylim[2]-ylim[1])*0.1+ylim[2]
    } else {
      xcoord<--(xlim[2]-xlim[1])*0.2+xlim[2]
      ycoord<--(ylim[2]-ylim[1])*0.1+ylim[2]       
    }
    
    plot <-plot + annotate("text",x=xcoord,y=ycoord,label=paste("Pearson Cor:",pearsoncor,"\n","Spearman Cor:",spearmancor,sep=""),colour = "black", size = 5)
  }
  
  if(!is.null(show_labels1) & !is.null(show_labels2)) {
    #data.show<-as.data.frame(cbind(xval,yval))
    #rownames(data.show)<-labels
    #data.show$label=labels
    #data.show<-data.show[show_labels,]
    #data.show<-data.show[!is.na(data.show[,1]),]
    if(textrepel) {
      if(!is.null(repelcol)) {
        # Not complete for label2
        plot<-plot+geom_text_repel(aes(label=label1),
                                   size=repelsize,color=repelcol,
                                   hjust=0, vjust=0, 
                                   max.overlaps=maxoverlaps, direction = repeldir,
                                   min.segment.length = 0,
                                   xlim = repelxlim, ylim = repelylim,
                                   fontface = "plain")
      } else {
        # Not complete for label2
        plot<-plot+geom_text_repel(aes(label=label1),
                                   size=repelsize, hjust=0, vjust=0, 
                                   max.overlaps=maxoverlaps, direction = repeldir,
                                   min.segment.length = 0,
                                   xlim = repelxlim, ylim = repelylim,
                                   fontface = "plain")
      }
    } else {
      plot<-plot+geom_text(aes(label=label1),color="black",hjust=0, vjust=0, 
                           nudge_x = 0.1, nudge_y = 0.1, size=4,
                           fontface = "bold")
      plot<-plot+geom_text(aes(label=label2),color="black",hjust=0, vjust=0,
                           nudge_x = 0.1, nudge_y = 0.1, size=4,
                           fontface = "plain")
    }
    
  }
  
  plot
  
}

# -----------------
# Volcano Function
# -----------------
volcano_helper <- function(data, anno, labels, fc, q, sig, xlab, ylab, fc_cutoff, q_cutoff, xlim, ylim) {
  
  # Default FC and Q (incase NA/blank input)
  if (is.na(fc_cutoff)) {
    fc_cutoff <- 2
  }
  if (is.na(q_cutoff)) {
    q_cutoff <- 0.05
  }
  fc_cutoff <- as.numeric(fc_cutoff)
  q_cutoff <- as.numeric(q_cutoff)
  
  # Removing NA Values
  data[fc][is.na(data[fc])] <- 0
  data[q][is.na(data[q])] <- 1
  
  # Recalculating Sigs for Cutoffs
  sig <- apply(data[,c(fc,q)],1,function(x) {
    if(x[2] < q_cutoff & x[1] >= fc_cutoff) {
      1
    } else if (x[2] < q_cutoff & x[1] <= -fc_cutoff) {
      -1
    } else {
      0
    }
  })
  
  # Gene Annotation
  data <- merge(data, anno, by = 0)
  data$gene_name <- make.unique(as.character(data$gene_name))
  gene <- data[,"gene_name"]
  
  # FC and Q Vectors
  fc <- data[[fc]]
  q <- data[[q]]
  
  return(enhanced_volcano_plot(gene = gene, labels = labels, fc = fc, q = q, sig = sig, xlab = xlab, ylab = ylab, main = "main",
                               fc_cutoff = fc_cutoff, q_cutoff = q_cutoff, xlim = xlim, ylim = ylim)
  )
}

# ----------------------
# Core Volcano Function
# ----------------------
enhanced_volcano_plot <- function(gene, fc, q, sig, labels = NULL,
                                  upcol = "red2", downcol = "blue2",
                                  xlab = "Log2FC", ylab = "-log10 P", main = "Volcano Plot",
                                  fc_cutoff = 0, q_cutoff = 0, xlim = "-10,10", ylim = "0,30") {
  # Gene Names, FC, and Q-values for DF
  fc <- as.numeric(unlist(fc))
  q <- as.numeric(unlist(q))
  
  # Remove NA values
  q[is.na(q)] <- 1
  sig[is.na(sig)] <- 0
  
  # GeneAnno Gene Names
  gene <- make.names(as.character(unlist(gene)), unique = T) # Unique = F for duplicates
  
  # Sig Counts
  up.number <- sum(sig == 1)
  down.number <- sum(sig == -1)
  total.number <- length(sig)
  
  # Setting xLim and yLim, Text -> Numeric
  if (length(xlim) == 1) {
    xlim <- as.numeric(unlist(strsplit(xlim,",")))
  }
  if (length(ylim) == 1) {
    ylim <- as.numeric(unlist(strsplit(ylim,",")))
  }
  
  # Default y-axis upper limit
  if (min(q) > 0) {
    max_pval <- -log10(min(q)) + 0.1
  } else {
    max_pval <- -log10(min(q[q > 0]) / 10) + 0.1 # changed here for pval=0
  }
  
  # Default yLim
  if (length(ylim) < 1) {
    ylim <- c(0, max_pval)
  }
  
  # Default x-axis upper/lower limit
  max_fc <- max(fc, na.rm = T)
  min_fc <- min(fc, na.rm = T)
  # Symmetric x-axis
  if (max_fc > abs(min_fc)) {
    min_fc <- -max_fc
  } else {
    max_fc <- -min_fc
  }
  
  # Default xLim
  if (length(xlim) < 1) { 
    xlim <- c(min_fc, max_fc)
  }
  
  # Filtered Data
  #newfc <- fc[abs(fc) <= max(xlim) & q >= 10^-max(ylim)]
  #newq <- q[abs(fc) <= max(xlim) & q >= 10^-max(ylim)]
  #newsig <- sig[abs(fc) <= max(xlim) & q >= 10^-max(ylim)]
  #newgene <- gene[abs(fc) <= max(xlim) & q >= 10^-max(ylim)]
  newfc <- fc[fc <= max(xlim) & fc >= min(xlim) & q >= 10^-max(ylim)]
  newq <- q[fc <= max(xlim) & fc >= min(xlim) & q >= 10^-max(ylim)]
  newsig <- sig[fc <= max(xlim) & fc >= min(xlim) & q >= 10^-max(ylim)]
  newgene <- gene[fc <= max(xlim) & fc >= min(xlim) & q >= 10^-max(ylim)]
  
  # DF for Volcano Plotting
  df <- data.frame(gene_name = as.character(newgene), lfc = newfc, q = newq, stringsAsFactors = F)
  rownames(df) <- newgene
  
  # Colors
  # Create a named vector of custom colors
  col_scheme <- c("Up" = upcol, "Down" = downcol, "N.S." = "grey")
  # Setting Colors
  cols <- rep(col_scheme["N.S."], length(newsig))
  cols[newsig == 1] <- col_scheme["Up"]
  cols[newsig == -1] <- col_scheme["Down"]
  names(cols)[cols == col_scheme["Up"]] <- "Up"
  names(cols)[cols == col_scheme["Down"]] <- "Down"
  names(cols)[cols == col_scheme["N.S."]] <- "N.S."
  
  # Labels: Default is Top 10 Up and down genes (based on FC)
  if (is.null(labels) == TRUE) {
    tmp.df <- data.frame(newgene, newfc, cols)
    rownames(tmp.df) <- rownames(df)
    tmp.df <- tmp.df[tmp.df$cols != col_scheme["N.S."], ]
    tmp.df <- tmp.df[order(tmp.df$newfc), ]
    labels <- c(rownames(tmp.df)[1:10], tail(rownames(tmp.df), n = 10))
  }
  # Render the Volcano Plot
  plt <- EnhancedVolcano(df,
                         x = "lfc", y = "q", lab = df$gene_name,
                         pCutoff = q_cutoff, FCcutoff = fc_cutoff,
                         gridlines.major = FALSE, gridlines.minor = FALSE,
                         drawConnectors = F, legendLabSize = 12,
                         cutoffLineCol = "red", colAlpha = 0.75,
                         cutoffLineType = "dashed", border = "full",
                         colCustom = cols, legendPosition = "right",
                         pointSize = 2, cutoffLineWidth = 0.4,
                         labFace = "plain", labSize = 4, subtitle = main,
                         # ylim = c(0, max_pval), xlim = c(min_fc, max_fc),
                         ylim = ylim, xlim = xlim,
                         axisLabSize = 12, captionLabSize = 12,
                         xlab = xlab, ylab = ylab, title = "",
                         caption = paste0("Total = ", total.number, " features"),
                         typeConnectors = "closed", legendIconSize = 2,
                         selectLab = labels, borderWidth = 1.5
  )
  # changed drawConnectors = T,labSize=3
  # Add numbers of up and down DE genes to the plot
  plt <- plt + geom_text(
    x = min(xlim), y = max(ylim), label = down.number,
    col = col_scheme["Down"], size = 5
  )
  plt <- plt + geom_text(
    x = max(xlim), y = max(ylim), label = up.number,
    col = col_scheme["Up"], size = 5
  )
  return(plt) # generate plot
}