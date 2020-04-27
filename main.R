addPathsToPedFile <- function(ped){

  colnames(ped) <- c("family_id","proband_id", "maternal_id", "paternal_id", "sex") 
  vcfs <- dir("./", ".final.vcf.gz$", recursive=T)
  bams <-  dir("./", "*bam$", recursive=T)
  vcfBamFiles <- do.call(rbind, lapply(1:nrow(ped), function(i){
    if (length(grep(paste0("/", ped$proband_id[i],"\\."),vcfs)) == 0){return (NULL)}
    if (length(grep(paste0("/", ped$proband_id[i],"\\."),bams)) == 0){return (NULL)}
    if (length(grep(paste0("/", ped$paternal_id[i],"\\."),vcfs)) == 0){return (NULL)}
    if (length(grep(paste0("/", ped$paternal_id[i],"\\."),bams)) == 0){return (NULL)}
    if (length(grep(paste0("/", ped$maternal_id[i],"\\."),vcfs)) == 0){return (NULL)}
    if (length(grep(paste0("/", ped$maternal_id[i],"\\."),bams)) == 0){return (NULL)}
    
    probandVCF = vcfs[grep(paste0("/", ped$proband_id[i],"\\."), vcfs)]
    x <- strsplit(probandVCF, "/")[[1]] 
    probandDir <- paste(x[1:(length(x)-1)],collapse="/")
    if (length(grep(probandDir, bams)) != 1){return (NULL)}
    
    paternalVCF = vcfs[grep(paste0("/", ped$paternal_id[i],"\\."), vcfs)]
    x <- strsplit(paternalVCF, "/")[[1]] 
    paternalDir <- paste(x[1:(length(x)-1)],collapse="/")
    if (length(grep(paternalDir, bams)) != 1){return (NULL)}
    
    maternalVCF = vcfs[grep(paste0("/", ped$maternal_id[i],"\\."), vcfs)]
    x <- strsplit(maternalVCF, "/")[[1]] 
    maternalDir <- paste(x[1:(length(x)-1)],collapse="/")
    if (length(grep(maternalDir, bams)) != 1){return (NULL)}
    
    data.frame(probandVCF = probandVCF, 
                probandBAM = bams[grep(probandDir, bams)],
                pateralVCF = paternalVCF,
                pateralBAM =bams[grep(paternalDir, bams)],
                maternalVCF =maternalVCF,
                maternalBAM = bams[grep(maternalDir, bams)])
  }))
  return( vcfBamFiles)
}





getPileup <- function(	data, sampleId, bamPath)
{
    library(parallel)
    pil <- mclapply(1:nrow(data), function(i){
      #  print(i)
        rr <- data[i,]
      
        ref <- rr$REF
        alt <- rr$ALT
        if (nchar(rr$REF)>1)
        {
            ref <- alt
            alt <- "dels"
        }else if (nchar(rr$ALT)>1)
        {
            alt <- "dups"
        }
        tmp <- tryCatch({
                                       
                                        pileup <- checkSide2(bamPath, rr)
                                        pRef <- pileup[ref]
                                        pAlt <- pileup[alt]
                                        (t(c(pRef,pAlt, pileup)))
                                        
                                    },
                                    error=function(e){
                                        print(e)
                                        (c(NA,NA))
                                    })
                                    
                        
                                
		colnames(tmp)<- paste0(sampleId, c("_REF","_ALT","_A","_C","_T","_G","_dels","_dups"))
                 return(tmp)
    }, mc.cores=10)
    return(do.call(rbind, pil))
 
}
#### check pileups


checkSide2 <- function(fl, rr)
{
    outFile <- paste( fl,".pileup.pos_",rr$POS, sep="")
    
    lines <- system (paste("samtools mpileup '",fl,"' -r ",paste0(rr$"#CHROM"), ":", rr$POS, "-",rr$POS , sep=""), intern=T)
    tt <- strsplit(lines[1], "\t")[[1]]
    pil <- tt[5]
    kk <- strsplit(pil,"")[[1]]
    aa <- length(which(tolower(kk)== "a"))
    cc <- length(which(tolower(kk)== "c"))
    tt <- length(which(tolower(kk)== "t"))
    gg <- length(which(tolower(kk)== "g"))
    dels <- length(which(tolower(kk)== "-"))
    dups <- length(which(tolower(kk)== "+"))
    res <- c(aa, cc,tt,gg, dels, dups)
    names(res) <- c("A", "C", "T", "G", "dels","dups")
    res
    
}

### retrieve the pileup information for all of the variants in the data object. This function utilizes the above checkSide2 function.

addPileups <- function(	data)
{
    library(parallel)
    pileups <- mclapply(1:nrow(data), function(i){
        print(i)
        rr <- data[i,]
      prBam <- rr$probandBAM
      p1Bam <- rr$pateralBAM
      p2Bam <- rr$maternalBAM
        
        if (rr$potDN == FALSE){
            tmp <- rep(NA,3*8)
            names(tmp)<- do.call(c, lapply(c("pr","p1","p2"), function(x){t(paste(x,c("_REF","_ALT","_A","_C","_T","_G","_dels","_dups"),sep=""))}))
            return(tmp);
        }
        ref <- rr$REF
        alt <- rr$ALT
        if (nchar(rr$REF)>1)
        {
            ref <- alt
            alt <- "dels"
        }else if (nchar(rr$ALT)>1)
        {
            alt <- "dups"
        }
        tmp <-do.call(cbind, lapply(1:3, function(j)
								{
                                    tryCatch({
                                        bam <- prBam
                                        if (j ==2){bam <- p1Bam}
                                        if (j == 3){bam<- p2Bam}
                                        pileup <- checkSide2(bam, rr)
                                        pRef <- pileup[ref]
                                        pAlt <- pileup[alt]
                                        return(t(c(pRef,pAlt, pileup)))
                                        
                                    },
                                    error=function(e){
                                        print(e)
                                        return(c(NA,NA))
                                    })
                                    
                                }))
                                
                                colnames(tmp)<- do.call(c, lapply(c("pr","p1","p2"), function(x){t(paste(x,c("_REF","_ALT","_A","_C","_T","_G","_dels","_dups"),sep=""))}))
                                return(tmp[])
    }, mc.cores=12)
	selectedPotMosaicWithPileups <- cbind(data,do.call(rbind, pileups))
	selectedPotMosaicWithPileups$pr_vRtR <- selectedPotMosaicWithPileups$pr_ALT /(selectedPotMosaicWithPileups$pr_REF + selectedPotMosaicWithPileups$pr_ALT)
	selectedPotMosaicWithPileups$p1_vRtR <- selectedPotMosaicWithPileups$p1_ALT /(selectedPotMosaicWithPileups$p1_REF + selectedPotMosaicWithPileups$p1_ALT)
	selectedPotMosaicWithPileups$p2_vRtR <- selectedPotMosaicWithPileups$p2_ALT /(selectedPotMosaicWithPileups$p2_REF + selectedPotMosaicWithPileups$p2_ALT)	
    return(selectedPotMosaicWithPileups)    
   
}


getPotMosaicFromVCFs <- function(pedWithPaths, nrofcores)
{	
	getAF <- function(vars){
		as.numeric(sapply(strsplit(vars[,10][[1]], ":"), function(x){x[3]}))
	}
		getDP <- function(vars){
		as.numeric(sapply(strsplit(vars[,10][[1]], ":"), function(x){x[4]}))
	}
	getKey <- function(vars){
		paste0(vars$"#CHROM", ":", vars$POS, "_", vars$REF, ">", vars$ALT)
	}

	getPotentialMosaics <- function(sel ,typeProband, typeParent1, typeParent2){
    		print("Reading VCF data - START")
		print("sel: -START")
		print(sel)
		print("sel: -END")
		
		print("colnames(sel): -START")
		print(colnames(sel))
		print("colnames(sel): -END")
		
		
		print("nrow(sel): -START")
		print(nrow(sel))
		print("nrow(sel): -END")
		
		print("typeProband: -START")
		print(typeProband)
		print("typeProband: -END")
		
		print("typeParent1: -START")
		print(typeParent1)
		print("typeParent1: -END")
		
		print("typeParent2: -START")
		print(typeParent2)
		print("typeParent2: -END")
		
		
		print(paste0("reading: ", sel[, typeProband]))
		probandVars  <- fread(paste0("zcat ",sel[, typeProband]), sep="\t",skip="#CHROM", header=T)    
    		print(paste0("reading: ", sel[, typeParent1]))
		parent1Vars <- fread(paste0("zcat ",sel[,typeParent1]), sep="\t",skip="#CHROM", header=T)
    		print(paste0("reading: ", sel[, typeParent2]))
		parent2Vars <-  fread(paste0("zcat ",sel[,typeParent2]), sep="\t",skip="#CHROM", header=T)		
		print("Reading VCF data - END")
		print(paste0("nrow(probandVars): ", nrow(probandVars )))
		print(paste0("ncol(probandVars): ", ncol(probandVars )))
		
		print("Filtering out non-PASS - START")
		#if(length(which(probandVars$FILTER == "PASS")) > 1){    probandVars <- probandVars[which(probandVars$FILTER == "PASS"),]}
		print("Filtering out non-PASS - END")
		print(paste0("nrow(probandVars): ", nrow(probandVars )))
		print(paste0("ncol(probandVars): ", ncol(probandVars )))
		      
		print("getAF - START")
		probandVars$VAF <- getAF (probandVars)
		parent1Vars$VAF <- getAF(parent1Vars)
		parent2Vars$VAF <- getAF(parent2Vars)
		print("getAF - END")		
		print(paste0("nrow(probandVars): ", nrow(probandVars )))
		print(paste0("ncol(probandVars): ", ncol(probandVars )))
		
		probandVars$key <- getKey(probandVars)
		parent1Vars$key <- getKey(parent1Vars)
		parent2Vars$key <- getKey(parent2Vars)
		
		print("getDP - START")
		probandVars$DP <- getDP (probandVars)
		parent1Vars$DP <- getDP(parent1Vars)
		parent2Vars$DP <- getDP(parent2Vars)
		print("getDP - END")
		print(paste0("nrow(probandVars): ", nrow(probandVars )))
		print(paste0("ncol(probandVars): ", ncol(probandVars )))
		
		print("DP and VAF filtering - START")		
		probandVars2 <- probandVars[which(probandVars$DP >= 20 & probandVars$VAF > 0.3 & probandVars$VAF < 0.7 ),]
		print("DP and VAF filtering - START")		
		print(paste0("nrow(probandVars): ", nrow(probandVars )))
		print(paste0("ncol(probandVars): ", ncol(probandVars )))
				      
    		parent1Vars$ToRemove <- parent1Vars$DP < 20 | parent1Vars$VAF > 0.1
		parent2Vars$ToRemove <- parent2Vars$DP < 20 | parent2Vars$VAF > 0.1
		print("Selecting variants absent in parents - START")		      
		selKeys <- setdiff(probandVars2$key, union(parent1Vars$key[which(parent1Vars$ToRemove)], parent2Vars$key[which(parent2Vars$ToRemove)]))
		print("Selecting variants absent in parents - END")		
		probandVars3 <- probandVars2[probandVars2$key %in% selKeys, ]
		print(paste0("nrow(probandVars): ", nrow(probandVars )))
		print(paste0("ncol(probandVars): ", ncol(probandVars )))
		      
		probandVars3
	}

	allPotMosaic <- mclapply((1:nrow(pedWithPaths)), function(i){
		print(paste0("getPotentialMosaic:", i))
		snps <- getPotentialMosaics(pedWithPaths[i,], "probandVCF", "pateralVCF", "maternalVCF")
		if (nrow(snps) ==0 ){return (NULL)}
		snps2 <- snps[which(nchar(snps$ALT) == 1 & nchar(snps$REF) == 1),]
		if (nrow(snps2) ==0 ){return (NULL)}
		df2 <-cbind(pedWithPaths[i,],snps2)
		df2$potDN <- TRUE
		df3 <- addPileups(df2)
		df3
	}, mc.cores=nrofcores)
	allPotMosaicDf <- rbindlist(allPotMosaic)
	allPotMosaicDf	
}

### MAIN ####
library(parallel)
library(data.table)
library(Rsamtools)
path_to_ped <- "" #set path to ped file
ped <- read.csv(path_to_ped, stringsAsFactors=F, header=F)

print("addPathsToPedFile() - START")
pedWithPaths <- addPathsToPedFile(ped)
print("addPathsToPedFile() - END")
print(paste0("pedWithPaths nrow: ", nrow(pedWithPaths)))
print("getPotentialMosaicFromVCF - START")

## second param is the number of cores
allPotMosaicDf <- getPotMosaicFromVCFs(pedWithPaths, 1)
print("getPotentialMosaicFromVCF - END")
print(paste0("nrow (allPotMosaicDf): ", nrow(allPotMosaicDf)))
print(paste0("ncol (allPotMosaicDf): ", ncol(allPotMosaicDf)))

# saving data:
dir.create("tomek")
write.table(allPotMosaicDf, paste0("tomek/", "allPotMosaicDf.tsv"), sep="\t", row.names=F)


