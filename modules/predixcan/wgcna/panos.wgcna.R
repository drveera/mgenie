args <- commandArgs(FALSE)

expr.file <- gsub("^-","",args[7])
output.name <- gsub("^-","",args[8])
output.name <- paste0(output.name,".",basename(expr.file))


library(coexpp)

datExpr <- readRDS(expr.file)

coexpression_network_pipeline <- function(
    NETWORK_EXPR_MAT = t(datExpr),
    EXPRESSION_DATA_NAME = output.name,

    RESULTS_DIR = ".",

    COEXPP_NUM_THREADS = 32, # Specify number of threads to run for network analysis
    SDOUT = 4, # Specify SD in absolute number for exclusion of subjects based on IAC

    CLUSTERING_CUT_METHOD = "tree", # choices are: "tree", "hybrid"

    NUM_HIGH_VAR_GENES_TO_RETAIN = NULL,

    REMOVE_SAMPLE_OUTLIERS = TRUE,
    BETA = NULL, # 6 # (hard threshold default for unsigned networks)
    R_SQUARED_CUT = 0.8, 

    TOM_PLOT.RAISE_POWER = 10

    ) {


   coexppSetThreads(COEXPP_NUM_THREADS)

   BASE_OUT_FILE_PREFIX = paste(RESULTS_DIR, "/", EXPRESSION_DATA_NAME, sep="")
   SAMPLE_CLUSTERING.FILE = paste(BASE_OUT_FILE_PREFIX, ".sample_clustering", ".pdf", sep="")
   SOFT_THRESHOLDING.FILE = paste(BASE_OUT_FILE_PREFIX, ".soft_thresholding.pdf", sep="")
   SCALE_FREE_PLOT.FILE = paste(BASE_OUT_FILE_PREFIX, ".scale_free_plot.pdf", sep="")
   GENE_CLUSTERING.FILE = paste(BASE_OUT_FILE_PREFIX, ".gene_clustering.pdf", sep="")
   TOM_CLUSTERING.FILE = paste(BASE_OUT_FILE_PREFIX, ".topological_overlap.png", sep="")
   MODULES.FILE = paste(BASE_OUT_FILE_PREFIX, ".clustering_modules.tsv", sep="")


    ## remove genes/samples woth missing data
    gene.set = goodSamplesGenes(NETWORK_EXPR_MAT, verbose=3)
    if (!gene.set$allOK) {
        if (sum(!gene.set$goodGenes) > 0)
            writeLines(paste("\nRemoving genes: ", paste(colnames(NETWORK_EXPR_MAT)[!gene.set$goodGenes], collapse = ", "), sep=""))
        if (sum(!gene.set$goodSamples) > 0)
            writeLines(paste("\nRemoving samples: ", paste(rownames(NETWORK_EXPR_MAT)[!gsg$goodSamples], collapse = ", "), sep=""))
        
        # Remove the "bad" genes and samples from the data:
        NETWORK_EXPR_MAT = NETWORK_EXPR_MAT[gene.set$goodSamples, gene.set$goodGenes]
    }

    # Check for sample outliers:
    IAC = cor(t(NETWORK_EXPR_MAT), use = "p")
    meanIAC = apply(IAC, 2, mean)
    sdCorr = sd(meanIAC)
    numbersd = (meanIAC - mean(meanIAC)) / sdCorr

    # Remove individuals that are outliers based on the SDOUT criterion and SET_LABELS selection
    if (REMOVE_SAMPLE_OUTLIERS) {
        removeSamples = rownames(NETWORK_EXPR_MAT)[abs(numbersd) > SDOUT]
        writeLines(paste("\nRemoving ", length(removeSamples), " outlier samples at IAC Absolute SD > ", sprintf("%.2f", SDOUT), " ", paste(removeSamples, collapse=", "), sep=""))
        NETWORK_EXPR_MAT = NETWORK_EXPR_MAT[setdiff(rownames(NETWORK_EXPR_MAT), removeSamples), ]
    }

   
    pdf(SAMPLE_CLUSTERING.FILE, width=10, height=5)
    par(mfrow = c(1,2))
    title = paste(EXPRESSION_DATA_NAME, ": ", "sample clustering to detect outliers", sep="")
    hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
    plot(numbersd)
    abline(h=SDOUT, col="red")
    abline(h=-SDOUT, col="red")
    dev.off()
    

 
    
    ################################################################################
    # Run WGCNA using coexpp:
    ################################################################################
    
    writeLines(paste("\nStarting ", EXPRESSION_DATA_NAME, " analysis on ", ncol(NETWORK_EXPR_MAT), " genes in ", nrow(NETWORK_EXPR_MAT), " samples.", sep=""))
    wgcnaCoexNet = coexpressionAnalysis(NETWORK_EXPR_MAT, beta=BETA, RsquaredCut=R_SQUARED_CUT, cut=CLUSTERING_CUT_METHOD)

    chosenBeta = wgcnaCoexNet$clusters@beta
    if (chosenBeta == 1) {
        stop("ERROR: BETA IS 1")
    }
    writeLines(paste("\nBETA IS: ", chosenBeta, sep=""))

    # Soft threshold plot:
    if (is.null(BETA)) {
        pdf(SOFT_THRESHOLDING.FILE)
        cex1 = 0.9
        cex2 = 1.2
        
        # Scale-free topology fit index as a function of the soft-thresholding power:
        softThreshStats = wgcnaCoexNet$clusters@sftStatistics
        softThresholds = softThreshStats[, "Power"]
        signedRsquaredScaleFreeTopo = - sign(softThreshStats[, "slope"]) * softThreshStats[, "SFT.R.sq"]
        title = paste(EXPRESSION_DATA_NAME, ": ", "Scale independence", sep="")
        plot(softThresholds, signedRsquaredScaleFreeTopo, xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=title)
        text(softThresholds, signedRsquaredScaleFreeTopo, labels=softThresholds, cex=cex1, col="red")        
        abline(h=R_SQUARED_CUT, col="red", lty="dotted")
        abline(v=chosenBeta, col="red", lty="dotted")

        meanConnectivity = softThreshStats[, "mean.k."]
        title = paste(EXPRESSION_DATA_NAME, ": ", "Mean connectivity", sep="")
        plot(softThresholds, meanConnectivity, xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main=title, log="y")
        text(softThresholds, meanConnectivity, labels=softThresholds, cex=cex1, col="red")
        meanConnAtBeta = softThreshStats[softThresholds == chosenBeta, "mean.k."]
        abline(h=meanConnAtBeta, col="red", lty="dotted")
        text(sort(softThresholds)[2], meanConnAtBeta, pos=3, labels=paste("|k|=", sprintf("%.1f", meanConnAtBeta), sep=""), cex=cex2, col="darkgray")
        abline(v=chosenBeta, col="red", lty="dotted")
        dev.off()
    }

    # Scale-free plot for selected beta:
    pdf(SCALE_FREE_PLOT.FILE)
    ScaleFreePlot = plotScaleFree(wgcnaCoexNet$clusters)
    print(ScaleFreePlot)
    dev.off()

    # Clustering plot for selected beta:
    pdf(GENE_CLUSTERING.FILE, width=25, height=5)
    GeneClusterPlot = plotClustering(wgcnaCoexNet$clusters, wgcnaCoexNet$geneModules, dendroLabels=FALSE, 
    main=paste(EXPRESSION_DATA_NAME, ": ", "Gene clustering", sep=""))
    print(GeneClusterPlot)
    dev.off()

    # Topological overlap matrix heatmap for selected beta:
    png(TOM_CLUSTERING.FILE, width=12, height=12, units="in", res=300)
    TOMHeatmapPlot = plotTOMHeatmap(wgcnaCoexNet$clusters, wgcnaCoexNet$geneModules, samplingThreshold=1000, plot.raise_power=TOM_PLOT.RAISE_POWER, 
    main=paste(EXPRESSION_DATA_NAME, ": ", "topological overlap matrix", sep=""))
    print(TOMHeatmapPlot)
    dev.off()

    # Output modules for the selected beta:
    GENES_TO_MODULES = data.frame(colnames(NETWORK_EXPR_MAT), wgcnaCoexNet$geneModules)
    colnames(GENES_TO_MODULES) = c("Gene", "Module")
    write.table(GENES_TO_MODULES, file=MODULES.FILE, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

coexpression_network_pipeline()
