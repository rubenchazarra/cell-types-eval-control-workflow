#!/usr/bin/env nextflow 

// Meta-workflow that controls execution of steps in cell type 
// prediction tools evaluation framework

// TODO: parse input data in a standardised format 
// NB: need to specify input data in workflow-wide config and provide them as 
// argument to sub-workflows 

// TODO: should we add cross-validation indices generation??
// TODO: add control statements on which tools need to be run ?


// after input data are obtained, run individual workflows 


// run garnett 
//process run_garnett_workflow {} 


// run scmap-cell 
process run_scmap_cell_workflow {
    publishDir "${params.tool_outputs_dir}", mode: 'copy'
    conda 'envs/nextflow.yaml'

    output: 
        file("scmap-cell_output.txt") into SCMAP_CELL_OUTPUT

    """
    RESULTS_DIR=\$PWD
    SUBDIR="scmap_cell"
    mkdir -p $WORK_DIR/\$SUBDIR     

    nextflow run $SCMAP_BASE_DIR/main.nf\
                        -work-dir $WORK_DIR/\$SUBDIR\
                        --results_dir \$RESULTS_DIR\
                        --projection_method ${params.scmap_cell.projection_method}
    """
}


// run scmap-cluster 
process run_scmap_cluster_workflow {
    publishDir "${params.tool_outputs_dir}", mode: 'copy'
    conda 'envs/nextflow.yaml'

    output:
        file("scmap-cluster_output.txt") into SCMAP_CLUST_OUTPUT


    """
    RESULTS_DIR=\$PWD
    SUBDIR="scmap_clust"
    mkdir -p $WORK_DIR/\$SUBDIR 

    nextflow run $SCMAP_BASE_DIR/main.nf\
                        -work-dir $WORK_DIR/\$SUBDIR\
                        --results_dir \$RESULTS_DIR\
                        --projection_method ${params.scmap_cluster.projection_method}

    """
}

// run scpred 
process run_scpred_workflow {
    publishDir "${params.tool_outputs_dir}", mode: 'copy'
    conda 'envs/nextflow.yaml'

    output:
        file("scpred_output.txt") into SCPRED_OUTPUT


    """
    RESULTS_DIR=\$PWD
    SUBDIR="scpred"
    mkdir -p $WORK_DIR/\$SUBDIR 

    nextflow run $SCPRED_BASE_DIR/main.nf\
                        -work-dir $WORK_DIR/\$SUBDIR\
                        --results_dir \$RESULTS_DIR\
                        --method ${params.scpred.method}\
                        --training_10x_dir ${params.scpred.training_10x_dir}\
                        --prediction_10x_dir ${params.scpred.training_10x_dir}\
                        --metadata_file ${params.scpred.metadata_file}

    """

}

// run predicted labels analysis
TOOL_OUTPUTS_DIR = Channel.fromPath(params.tool_outputs_dir)
REF_LAB_FILE = Channel.fromPath(params.label_analysis.ref_labels_file)
process run_label_analysis {
    conda 'envs/nextflow.yaml' 
    publishDir "${params.label_analysis_output}", mode: 'copy'

    input:
        file(tool_outputs_dir) from TOOL_OUTPUTS_DIR
        file(ref_lab_file) from REF_LAB_FILE

    output:
        file("${params.cell_anno_table}") into CELL_ANNO_TABLE
         file("${params.label_analysis.tool_table_pvals}") into TOOL_PERF_PVALS


    """
    RESULTS_DIR=\$PWD
    SUBDIR="label_analysis"
    mkdir -p $WORK_DIR/\$SUBDIR 

    nextflow run $LABEL_ANALYSIS_BASE_DIR/main.nf\
                        -work-dir $WORK_DIR/\$SUBDIR\
                         --results_dir \$RESULTS_DIR\
                         --input_dir ${tool_outputs_dir}\
                         --ref_labels_file ${ref_lab_file}\
    """

}
