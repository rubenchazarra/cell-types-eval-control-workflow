#!/usr/bin/env Rscript 

#This script is to perform subsetting of dataset based on the cel labels generated in the k-fold cv previous process. 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
 make_option(
    c("-i", "--input-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to the SCE object to perform the test/train split."
  ), 
  make_option(
    c("-c", "--input-cell-labels"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input cell labels to perform test/train data subseting'
  ),
 make_option(
    c("-q", "--test-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to the test sce split where cell identity will be inferred based on the model generated over train sce split."
  ),
 make_option(
    c("-t", "--train-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to the train sce split where the the model for inferring cell identity will be generated."
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "input_cell_labels", "test_sce_object", "train_sce_object"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input SCE file does not exist.")
sce <- readRDS(opt$input_sce_object)
#if object has no colnames, append them
if(is.null(colnames(sce))) colnames(sce) <- sce$Barcode

#read input cell labels 
if(!file.exists(opt$input_cell_labels)) stop("Input cell labels file does not exist.")
test_cell_labels <- readRDS(opt$input_cell_labels)

#split data into test and train
test_sce <- sce[, test_cell_labels]
train_sce <- sce[, !(colnames(sce) %in% test_cell_labels)]

print(test_sce)
print(train_sce)
