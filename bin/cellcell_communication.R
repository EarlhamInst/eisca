#!/usr/bin/env Rscript

library(CellChat)
library(patchwork)
library(Seurat)
library(Matrix)
library(jsonlite)
options(stringsAsFactors = FALSE)



library(argparse)

parser <- ArgumentParser()
parser$add_argument("--outdir", help="Output directory")
parser$add_argument("--count", help="Input count matrix mtx file")
parser$add_argument("--metadata", help="Input metadata csv file")
parser$add_argument("--gids", help="Input gene ID list csv file")
parser$add_argument("--cids", help="Input cell ID list csv file")
parser$add_argument("--group", default="majority_voting", help="Specify the column for grouping the cells")
parser$add_argument("--n_actgrps", type="integer", default=12, help="Number of top active groups")
parser$add_argument("--n_pairLR", type="integer", default=100, help="Number of top LR pairs")
parser$add_argument("--normalize", action="store_true", help="Indicates whether to normalize the counts")
parser$add_argument("--db", default="human",  choices=c('human', 'mouse'), help="Specify the species of CellChatDB.")
parser$add_argument("--dbc", help="The categories of CellChatDB, e.g. Secreted Signaling")
# parser$add_argument("--dbv", help="The version of CellChatDB, e.g. v1")
parser$add_argument("--dbenps", action="store_true", help="Use all CellChatDB excepting 'Non-protein Signaling'")
parser$add_argument("--threads", type="integer", default=12, help="Number of threads for parallel runs")
parser$add_argument("--mean_method", default="triMean",  choices=c('triMean', 'truncatedMean'), 
        help="Specify the method for calculating the average gene expression per cell group")
parser$add_argument("--mincells", type="integer", default=10, help="The minimum number of cells required in each cell group")
parser$add_argument("--meta", default="auto", choices=c('auto', 'sample', 'group'), help="Specify a metadata column to define separate subsets of cells for analysis")
parser$add_argument("--scale", type="double", default=1, help="Adjust circular plot size by scaling plot's width and height")
parser$add_argument("--pdf", action="store_true", help="Whether to generate figure files in PDF format")
# parser$add_argument("--plotcoef", type="double", default=1.0, help="plot size equals the coef times default plot size")

args <- parser$parse_args()


dir.create(args$outdir, showWarnings = FALSE)

# create a cellchat object from counts and metadata
meta <- read.csv(args$metadata, header = TRUE, row.names = 1)
counts <- t(readMM(args$count))
genes <- readLines(args$gids)
cells <- readLines(args$cids)
rownames(counts) <- genes
colnames(counts) <- cells
if ("sample" %in% colnames(meta)) {
    colnames(meta)[colnames(meta) == "sample"] <- "samples"
    meta$samples <- as.factor(meta$samples)
}


# normalize the count data if input data is raw counts
if(args$normalize){
    library.size <- Matrix::colSums(counts)
    counts <- as(log1p(Matrix::t(Matrix::t(counts)/library.size) * 10000), "CsparseMatrix")
}

# set a column to define separate batches for analysis
batch <- "sample"
if (args$meta == "auto") {
  if ("group" %in% colnames(meta)) {
    batch <- "group"
  } else if ("plate" %in% colnames(meta)) {
    batch <- "plate"
  }
} else {
  batch <- args$meta
}     


# Set the ligand-receptor interaction database
if(args$db == "human"){
    CellChatDB <- CellChatDB.human
}else if(args$db == "mouse"){
    CellChatDB <- CellChatDB.mouse
}
CellChatDB.use <- CellChatDB
if(!is.null(args$dbc) && nchar(args$dbc) > 0){
    CellChatDB.use <- subsetDB(CellChatDB, search = args$dbc, key = "annotation")
}else if(args$dbenps){
    CellChatDB.use <- subsetDB(CellChatDB)
}


for(sid in unique(meta[[batch]])){
    sub_idx <- which(meta[[batch]] == sid)
    counts_s <- counts[, sub_idx]
    meta_s <- meta[sub_idx,]
    path_outdir_s <- paste0(args$outdir, '/', batch, '_', sid)
    dir.create(path_outdir_s, showWarnings = FALSE)

    cellchat <- createCellChat(object = counts_s, meta = meta_s, group.by = args$group)

    if(!is.null(CellChatDB.use)){cellchat@DB <- CellChatDB.use}

    # Preprocessing the expression data
    cellchat <- subsetData(cellchat)
    future::plan("multisession", workers = args$threads)
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    # Inference of cell-cell communication network
    cellchat <- computeCommunProb(cellchat, type = args$mean_method)
    cellchat <- filterCommunication(cellchat, min.cells = args$mincells)
    # Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)

    # Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)

    # Compute the network centrality scores
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    groupSize <- as.numeric(table(cellchat@idents))

    # Identify rows and columns with any non-zero interaction weights
    mat <- cellchat@net$count
    actgrp_idx <- which(rowSums(mat) > 0 & colSums(mat) > 0)
    mat_active <- mat[actgrp_idx, actgrp_idx]
    groupSize_active <- groupSize[actgrp_idx]
    grpnames_active <- rownames(mat_active)
    activity_score <- rowSums(mat_active) + colSums(mat_active)
    n_select <- min(args$n_actgrps, length(activity_score))
    top_idx <- order(activity_score, decreasing = TRUE)[1:n_select]
    mat_top <- mat_active[top_idx, top_idx]
    groupSize <- groupSize_active[top_idx]
    grpnames_top <- grpnames_active[top_idx]
    mat_weight_active <- cellchat@net$weight[actgrp_idx, actgrp_idx]
    mat_weight_top <- mat_weight_active[top_idx, top_idx]
    topgrp_idx <- actgrp_idx[top_idx]

    png(file=paste0(path_outdir_s, "/aggregated_network_all.png"), width=10*args$scale,height=10*args$scale, units="in", res=150)
    netVisual_circle(
        mat_top, 
        vertex.weight = groupSize, 
        weight.scale = T, 
        label.edge= F, 
        title.name = "Number of interactions",
        vertex.label.cex = 1.5
    )
    dev.off()
    png(file=paste0(path_outdir_s, "/aggregated_network_all_weights.png"), width=10*args$scale,height=10*args$scale, units="in", res=150)
    netVisual_circle(
        mat_weight_top, 
        vertex.weight = groupSize, 
        weight.scale = T, 
        label.edge= F, 
        title.name = "Interaction weights/strength",
        vertex.label.cex = 1.5
    )
    dev.off()

    # mat <- cellchat@net$weight
    mat <- mat_weight_top
    # grpnames <- rownames(mat)
    Ngrp <- nrow(mat)
    png(file=paste0(path_outdir_s, "/aggregated_network_groups.png"),width=9,height=3*ceiling(Ngrp/3), units = "in", res=100*ceiling(Ngrp/3))
    par(mfrow = c(ceiling(Ngrp/3),3), xpd=TRUE)
    for (i in 1:Ngrp) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(
            mat2, 
            vertex.weight = groupSize, 
            weight.scale = T, 
            edge.weight.max = max(mat), 
            title.name = rownames(mat)[i],
        )
    }
    dev.off()

    if(args$pdf){
        pdf(file=paste0(path_outdir_s, "/aggregated_network_all.pdf"), width=10*args$scale,height=10*args$scale)
        netVisual_circle(
            mat_top, 
            vertex.weight = groupSize, 
            weight.scale = T, 
            label.edge= F, 
            title.name = "Number of interactions",
            vertex.label.cex = 1.5     
        )
        dev.off()
        pdf(file=paste0(path_outdir_s, "/aggregated_network_all_weights.pdf"), width=10*args$scale,height=10*args$scale)
        netVisual_circle(
            mat_weight_top, 
            vertex.weight = groupSize, 
            weight.scale = T, 
            label.edge= F, 
            title.name = "Interaction weights/strength",
            vertex.label.cex = 1.5
        )
        dev.off()

        pdf(file=paste0(path_outdir_s, "/aggregated_network_groups.pdf"),width=9,height=3*ceiling(Ngrp/3))
        par(mfrow = c(ceiling(Ngrp/3),3), xpd=TRUE)
        for (i in 1:Ngrp) {
            mat2 <- matrix(0, nrow = Ngrp, ncol = ncol(mat), dimnames = dimnames(mat))
            mat2[i, ] <- mat[i, ]
            netVisual_circle(
                mat2, 
                vertex.weight = groupSize, 
                weight.scale = T, 
                edge.weight.max = max(mat), 
                title.name = rownames(mat)[i],
            )
        }
        dev.off()
    }

    # Visualization of cell-cell communication network for individual signalling pathways
    cellchat_subset <- cellchat
    cellchat_subset@idents <- cellchat@idents[topgrp_idx]
    cellchat_subset@net$weight <- cellchat@net$weight[topgrp_idx, topgrp_idx]
    cellchat_subset@net$count  <- cellchat@net$count[topgrp_idx, topgrp_idx]
    cellchat_subset <- computeCommunProbPathway(cellchat_subset)
    cellchat_subset <- netAnalysis_computeCentrality(cellchat_subset)
    for(pathway in cellchat_subset@netP$pathways){
        try({
        png(file=paste0(path_outdir_s, "/pathway_network_circle_", pathway, '.png'), width=10*args$scale,height=10*args$scale, units="in", res=150)
        par(mar=c(1, 3, 1, 3), cex=1.5)
        netVisual_aggregate(cellchat_subset, signaling = pathway, layout = "circle", signaling.name = pathway)
        dev.off()
        png(file=paste0(path_outdir_s, "/pathway_network_chord_", pathway, '.png'), width=10*args$scale,height=10*args$scale, units="in", res=150)
        netVisual_chord_cell(cellchat_subset, signaling = pathway, lab.cex=1)
        dev.off()
        png(file=paste0(path_outdir_s, "/pathway_network_heatmap_", pathway, '.png'), width=1.5*Ngrp,height=1*Ngrp, units="in", res=100)
        print(netVisual_heatmap(cellchat_subset, signaling = pathway, color.heatmap = "Reds", font.size = 12, font.size.title = 15))
        dev.off()
        }, silent = TRUE)
        if(args$pdf){
            try({
            pdf(file=paste0(path_outdir_s, "/pathway_network_circle_", pathway, '.pdf'), width=10*args$scale,height=10*args$scale)
            par(mar=c(1, 3, 1, 3), cex=1.5)
            netVisual_aggregate(cellchat_subset, signaling = pathway, layout = "circle", signaling.name = pathway)
            dev.off()
            pdf(file=paste0(path_outdir_s, "/pathway_network_chord_", pathway, '.pdf'), width=10*args$scale,height=10*args$scale)
            netVisual_chord_cell(cellchat_subset, signaling = pathway, lab.cex=1)
            dev.off()
            pdf(file=paste0(path_outdir_s, "/pathway_network_heatmap_", pathway, '.pdf'), width=1.5*Ngrp,height=1*Ngrp)
            print(netVisual_heatmap(cellchat_subset, signaling = pathway, color.heatmap = "Reds", font.size = 12, font.size.title = 15))
            dev.off()
            }, silent = TRUE)  
        }

        # Compute the contribution of each ligand-receptor pair to the overall signaling pathway
        pairLR <- extractEnrichedLR(cellchat_subset, signaling = pathway, geneLR.return = FALSE)
        png(file=paste0(path_outdir_s, "/pathway_network_contribution_", pathway, '.png'), width=6,height=1.5*nrow(pairLR), units="in", res=100)
        print(netAnalysis_contribution(cellchat_subset, signaling = pathway, font.size = 12, font.size.title = 15))
        dev.off()
        if(args$pdf){
            pdf(file=paste0(path_outdir_s, "/pathway_network_contribution_", pathway, '.pdf'), width=6,height=1.5*nrow(pairLR))
            print(netAnalysis_contribution(cellchat_subset, signaling = pathway, font.size = 12, font.size.title = 15))
            dev.off()
        }

        # visualize cell-cell communication mediated by a single ligand-receptor pair
        for(i in 1:nrow(pairLR)){
            LR <- pairLR[i,]
            try({
            png(file=paste0(path_outdir_s, "/pathway_network_LR_circle_", pathway, '.png'), width=10*args$scale,height=10*args$scale, units="in", res=200)
            netVisual_individual(cellchat_subset, signaling = pathway, pairLR.use = LR, layout = "circle")
            dev.off()
            png(file=paste0(path_outdir_s, "/pathway_network_LR_chord_", pathway, '.png'), width=10*args$scale,height=10*args$scale, units="in", res=200)
            netVisual_individual(cellchat_subset, signaling = pathway, pairLR.use = LR, layout = "chord")
            dev.off()
            }, silent = TRUE)
            if(args$pdf){
                try({
                pdf(file=paste0(path_outdir_s, "/pathway_network_LR_circle_", pathway, '.pdf'), width=10*args$scale,height=10*args$scale)
                netVisual_individual(cellchat_subset, signaling = pathway, pairLR.use = LR, layout = "circle")
                dev.off()
                pdf(file=paste0(path_outdir_s, "/pathway_network_LR_chord_", pathway, '.pdf'), width=10*args$scale,height=10*args$scale)
                netVisual_individual(cellchat_subset, signaling = pathway, pairLR.use = LR, layout = "chord")
                dev.off()
                }, silent = TRUE)
            }
        }

        # Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
        for(k in seq_along(topgrp_idx)){
            i = topgrp_idx[k]
            grpname = gsub("[/ ]", "_", grpnames_top[k])
            try({
                gg <- netVisual_bubble(cellchat_subset, sources.use = i, remove.isolate = FALSE)
                nx <- length(levels(cellchat_subset@idents))
                ny <- nrow(subsetCommunication(cellchat_subset, sources.use = i))
                ggsave(
                    filename=paste0(path_outdir_s, "/cellcell_LR_bubble_", grpname, '.png'), 
                    plot=gg, 
                    width=max(6, 0.3*nx),
                    height=max(6, 0.05*ny), 
                    units = 'in', 
                    dpi = 200, 
                    limitsize = FALSE)
                if(args$pdf){
                    ggsave(
                        filename=paste0(path_outdir_s, "/cellcell_LR_bubble_", grpname, '.pdf'), 
                        width=max(6, 0.3*nx),
                        height=max(6, 0.05*ny), 
                        limitsize = FALSE)
                }
            },silent = TRUE)
        }
        for(k in seq_along(topgrp_idx)){
            i = topgrp_idx[k]
            grpname = gsub("[/ ]", "_", grpnames_top[k])
            pairLR <- subsetCommunication(cellchat_subset, sources.use = i)
            pairLR <- pairLR[order(pairLR$prob, decreasing = TRUE), ]
            pairLR <- head(pairLR, args$n_pairLR)
            try({
            png(file=paste0(path_outdir_s, "/cellcell_LR_chord_", grpname, '.png'), width=12*args$scale,height=10*args$scale, units="in", res=150)
            netVisual_chord_gene(cellchat_subset, sources.use = i, pairLR.use = pairLR, lab.cex = 1.2, legend.pos.y = 30)
            dev.off()            
            if(args$pdf){
                pdf(file=paste0(path_outdir_s, "/cellcell_LR_chord_", grpname, '.pdf'), width=12*args$scale,height=10*args$scale)
                netVisual_chord_gene(cellchat_subset, sources.use = i, pairLR.use = pairLR, lab.cex = 1.2,legend.pos.y = 30)
                dev.off()
            }
            },silent = TRUE)
        }

        # Plot the signaling gene expression distribution using violin
        try({
        active_idents <- levels(cellchat_subset@idents)
        # png(file=paste0(path_outdir_s, "/pathway_genes_violin_", pathway, '.png'), width=1*Ngrp,height=10, units="in", res=150)
        # print(plotGeneExpression(cellchat, signaling = pathway, enriched.only = TRUE))
        png(file=paste0(path_outdir_s, "/pathway_genes_violin_", pathway, '.png'), width=max(8, 0.8*length(active_idents)), height=10, units="in", res=150)
        print(plotGeneExpression(cellchat, signaling = pathway, enriched.only = TRUE, group.by = args$group))
        dev.off()
        if(args$pdf){
            pdf(file=paste0(path_outdir_s, "/pathway_genes_violin_", pathway, '.pdf'), width=max(8, 0.8*length(active_idents)),height=10)
            print(plotGeneExpression(cellchat_subset, signaling = pathway, enriched.only = TRUE, group.by = args$group))
            dev.off()
        }
        }, silent = TRUE)

        # visualize the network centrality scores
        try({
        num_active <- length(levels(cellchat_subset@idents))
        # png(file=paste0(path_outdir_s, "/pathway_network_centrality_", pathway, '.png'), width=1*Ngrp,height=4, units="in", res=100)
        png(file=paste0(path_outdir_s, "/pathway_network_centrality_", pathway, '.png'), width=max(6, 0.8*num_active),height=4, units="in", res=100)
        netAnalysis_signalingRole_network(cellchat_subset, signaling = pathway, width=1.5*num_active,height=4, font.size = 10)
        dev.off()
        if(args$pdf){
            pdf(file=paste0(path_outdir_s, "/pathway_network_centrality_", pathway, '.pdf'), width=max(6, 0.8*num_active),height=4)
            netAnalysis_signalingRole_network(cellchat_subset, signaling = pathway, width=1.5*num_active,height=4, font.size = 10)
            dev.off()
        }
        }, silent = TRUE)
    }

    # Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
    heat_w <- unit(0.3*Ngrp, "in")
    heat_h <- unit(0.8*Ngrp, "in")
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = heat_w, height = heat_h)
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = heat_w, height = heat_h)
    png(file=paste0(path_outdir_s, "/heatmap_signaling_patterns.png"), units="in", width=0.5*Ngrp+4,height=1*Ngrp, res=150)
    print(ht1 + ht2)
    dev.off()
    if(args$pdf){
        pdf(file=paste0(path_outdir_s, "/heatmap_signaling_patterns.pdf"), width=0.5*Ngrp+4,height=1*Ngrp)
        print(ht1 + ht2)
        dev.off()
    }


    #  save dataframe consisting of all the inferred cell-cell communications at the level of ligands/receptors
    saveRDS(subsetCommunication(cellchat), file = paste0(path_outdir_s, "/inferred_cellcell_comm.rds"))
    # Save the CellChat object
    saveRDS(cellchat, file = paste0(path_outdir_s, "/cellchat.rds"))


}


# save analysis parameters into a json file
params <- list()
params[["--count"]] <- as.character(args$count)
params[["--metadata"]] <- as.character(args$metadata)
params[["--gids"]] <- as.character(args$gids)
params[["--cids"]] <- as.character(args$cids)
params[["--group"]] <- as.character(args$group)
params[["--n_actgrps"]] <- as.character(args$n_actgrps)
params[["--n_pairLR"]] <- as.character(args$n_pairLR)
if(args$normalize){params[["--normalize"]] <- as.character(args$normalize)}
params[["--db"]] <- as.character(args$db)
if (!is.null(args$dbc)) {params[["--dbc"]] <- as.character(args$dbc)}
if(args$dbenps){params[["--dbenps"]] <- as.character(args$dbenps)}
params[["--mean_method"]] <- as.character(args$mean_method)
params[["--mincells"]] <- as.character(args$mincells)
params[["--threads"]] <- as.character(args$threads)
params[["--scale"]] <- as.character(args$scale)
write_json(params, file.path(args$outdir, "parameters.json"), pretty = TRUE, auto_unbox = TRUE)



# yaml::write_yaml(
# list(
#     'MTX_TO_SEURAT'=list(
#         'Seurat' = paste(packageVersion('Seurat'), collapse='.')
#     )
# ),
# "versions.yml"
# )
