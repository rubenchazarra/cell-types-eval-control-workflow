#Generate k-fold cell indexes of input SCE object. 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
  make_option(
    c("-i", "--input-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input SCE object in rds format'
  ),
  make_option(
    c("-k", "--k-folds-number"),
    action = "store",
    default = 5,
    type = 'integer',
    help = 'Number of groups to split the input sce object'
  ), 
  make_option(
    c("-o", "--output-cell-labels"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the output cell labels vectors in rds format.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "output_cell_labels"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input file does not exist.")
sce <- readRDS(opt$input_sce_object)
#if object has no colnames, append them
if(is.null(colnames(sce))) colnames(sce) <- sce$Barcode
#create Folds with cell indexes
#suppressPackageStartupMessages(require(caret))
suppressPackageStartupMessages(require(caret))
#list of folds with cell indexes
print(sce)
cell_index_list <- createFolds(y = colnames(sce), k = opt$k_folds_number)
#list of folds with cell names
cell_labels_list <- lapply(cell_index_list, function(x) colnames(sce)[x])

#save SCE object with size Factors
print(cell_labels_list)
saveRDS(cell_labels_list$Fold1, file=opt$output_cell_labels)
