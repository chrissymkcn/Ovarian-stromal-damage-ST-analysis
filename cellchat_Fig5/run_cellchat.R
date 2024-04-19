## Load the required libraries
Sys.setenv(RETICULATE_PYTHON = "/home/chrissy1/.conda/envs/dcats/bin/python")
library(reticulate)
py_config()
options(future.globals.maxSize = 100000 * 1024^2)
workers = 30
memory.limit(size = 200000)  # Set to 2.5 GB (size in MB)

ptm = Sys.time()

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(anndata)
library(Matrix)
library(dplyr)
# options(future.globals.maxSize = +Inf)


# # Part I: Data input & processing and initialization of CellChat object
# CellChat requires four user inputs:

# * **Gene expression data of spots/cells**: genes should be in rows with rownames and cells in columns with colnames. Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1) is required as input for CellChat analysis. If user provides count data, we provide a `normalizeData` function to account for library size and then do log-transformed.

# * **User assigned cell labels**: a data frame (rows are cells with rownames) consisting of cell information, which will be used for defining cell groups.

# * **Spatial locations of spots/cells**: a data matrix in which each row gives the spatial locations/coordinates of each cell/spot. For 10X Visium, this information is in `tissue_positions.csv`.

# * **Scale factors and spot diameters of the full resolution images**: a list containing the scale factors and spot diameter for the full resolution images. scale.factors must contain an element named `spot.diameter`, which is the theoretical spot size (e.g., 10x Visium (spot.size = 65 microns)); and another element named `spot`, which is the number of pixels that span the diameter of a theoretical spot size in the original, full-resolution image.
# For 10X Visium, scale.factors are in the file `scalefactors_json.json`. `spot` is the `spot.size.fullres`.


#### For each sample, run the following code to generate the sliced cellchat object
#### special step to subset each adata into two parts
data_dir = '/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed'
cellchat_dir = paste0(data_dir, '/cellchat')
setwd(cellchat_dir)
adata_full = read_h5ad(paste0(data_dir, '/SCT_copy.h5ad'))
endo = read_h5ad('/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed/annot_celltype_ind_analysis/ENDO_SM(CLDN5)/SCT.h5ad')
# rownames(endo$obs) = endo$obs[['Unnamed: 0']]
# write_h5ad(endo, ('/home/chrissy1/spatial/stomics/ovary_froz/redo/seurat/bin50_processed/annot_celltype_ind_analysis/ENDO_SM(CLDN5)/SCT.h5ad'))
val = sapply(endo$obs$leiden, function(x) paste0('ENDO_SM_', x))
adata_full$obs$annot = adata_full$obs$annot %>% as.character()
adata_full$obs[endo$obs_names, 'annot'] = val
samples = unique(adata_full$obs$batch)

for (key in samples) {
# key = samples[[1]]
    adata = adata_full[adata_full$obs$batch == key, ]
    adata = adata[sample(rownames(adata), size = floor(nrow(adata))), ]
    data.input = adata$X %>% t(.) %>% as(.,'CsparseMatrix')
    colnames(data.input) = rownames(adata)
    rownames(data.input) = adata$var_names
    meta = adata$obs
    rownames(meta) = colnames(data.input)
    # load spatial imaging information
    # Spatial locations of spots from full (NOT high/low) resolution images are required
    spatial.locs = adata$obsm$spatial[, 1:2] %>% as.data.frame()
    rownames(spatial.locs) = colnames(data.input)   
    scale.factors = list(spot.diameter = 25, spot = 50, # these two information are required
                        ratio = 0.22, tol = 0.11,
                        fiducial = 100, hires = 0.7, lowres = 0.1 # these three information are not required
        )
    ## Create a CellChat object
    # **USERS can create a new CellChat object from a data matrix or Seurat.** If input is a Seurat object, the meta data in the object will be used by default and USER must provide `group.by` to define the cell groups. e.g, group.by = "ident" for the default cell identities in Seurat object.
    # **NB: If USERS load previously calculated CellChat object (version < 1.6.0), please update the object via `updateCellChat`**
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "annot", datatype = "spatial", coordinates = spatial.locs, spatial.factors = scale.factors)
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
    ## Set the ligand-receptor interaction database
    # Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse. CellChatDB v2 contains ~3,300 validated molecular interactions, including ~40% of secrete autocrine/paracrine signaling interactions, ~17% of extracellular matrix (ECM)-receptor interactions, ~13% of cell-cell contact interactions and ~30% non-protein signaling. Compared to CellChatDB v1, CellChatDB v2 adds more than 1000 protein and non-protein interactions such as metabolic and synaptic signaling. It should be noted that for molecules that are not directly related to genes measured in scRNA-seq, CellChat v2 estimates the expression of ligands and receptors using those molecules’ key mediators or enzymes for potential communication mediated by non-proteins. **Critically, synaptic signaling interactions can only be used when inferring neuron-neuron communication.**
    # Users can update CellChatDB by adding their own curated ligand-receptor pairs.Please check our tutorial on how to do it.
    CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
    # use a subset of CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
    # use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
    ## Preprocessing the expression data for cell-cell communication analysis
    # To infer the cell state-specific communications, we identify over-expressed ligands or receptors in one cell group and then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
    # We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. Specifically, a diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors. One might be concerned about the possible artifact introduced by this diffusion process, however, it will only introduce very weak communications. USERS can also skip this step and set `raw.use = TRUE` in the function `computeCommunProb()`.
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    # future::plan("multisession", workers = workers)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
    # cellchat <- projectData(cellchat, PPI.mouse)
    # # Part II: Inference of cell-cell communication network
    remove(adata)
    remove(data.input)
    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05,
                            distance.use = TRUE, interaction.range = 250, 
                            scale.distance = 1.2, contact.dependent = FALSE, contact.range = NULL)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    # ## Infer the cell-cell communication at a signaling pathway level
    # CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.
    # NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.
    cellchat <- computeCommunProbPathway(cellchat)
    ## Calculate the aggregated cell-cell communication network
    # We can calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. USER can also calculate the aggregated network among a subset of cell groups by setting `sources.use` and `targets.use`.
    cellchat <- aggregateNet(cellchat)
    # cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
    saveRDS(cellchat, file = paste0(key, "_endo_subbed.rds"))
}

