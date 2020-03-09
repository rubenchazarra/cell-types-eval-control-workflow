#!/usr/bin/env nextflow 

// Workflow to perform k-fold cross validation over the imput dataset

//cross validation
if(params.run == "True"){
	process generate_folds{
		
		publishDir "${params.output_dir}", mode: 'copy'
        	conda "${baseDir}/envs/r-caret.yaml"

		errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
        	maxRetries 5
        	memory { 16.GB * task.attempt }

        	output:
        	    file("output_cell_indexes.rds") into K_FOLD_CELL_LABELS

        	"""
        	caret-create-folds.R\
        	        --input-barcodes-tsv ${params.barcodes}\
        	        --k-folds-number ${params.generate_folds.k_folds_num}\
        	        --output-cell-indexes ${params.generate_folds.output_cell_indexes}
        	"""
	}
      
	process split_train_test{
		
		publishDir "${params.output_dir}", mode: 'copy'
		conda "${baseDir}/envs/r-caret.yaml"
		
		errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
                maxRetries 5
                memory { 16.GB * task.attempt }

                input:
		    file(output_cell_indexes) from K_FOLD_CELL_LABELS
		
		output:
                    //file("*.rds") into TEST_TRAIN_OBJECTS 
                    file("test.*") into TEST_OBJECTS 
                    file("train.*") into TRAIN_OBJECTS 

                """
                split-train-test-data.R\
                        --input-matrix ${params.matrix}\
                        --input-cell-indexes ${output_cell_indexes}\
			--input-barcodes-tsv ${params.barcodes}\
                        --test-name ${params.split_train_test.test_name}\
                        --train-name ${params.split_train_test.train_name}
                """
        }
	

	//Channel for features.tsv
	GENES_CH = Channel.value(params.genes)
	
	//Generate tuples  as [Fold_i, [test_file, train_file]] to later group tuples by fold number to send to the methods
	TEST_OBJECTS_TUPLE = TEST_OBJECTS.flatten().map{f -> tuple("${f.baseName}".split("\\.")[1], f) }.groupTuple()
	TRAIN_OBJECTS_TUPLE = TRAIN_OBJECTS.flatten().map{f -> tuple("${f.baseName}".split("\\.")[1], f) }.groupTuple()

//Get class of object	
	//TRAIN_OBJECTS_TUPLE.subscribe{ println it.getClass().getSimpleName() }

	
	
	//Add gene file to the cross validations [NOTE: At the cost of losing the list name t -> t[1]. I don't know how to do it otherwise
	TEST_DATA = TEST_OBJECTS_TUPLE.map{t -> t[1]}.combine(GENES_CH)
	TRAIN_DATA = TRAIN_OBJECTS_TUPLE.map{t -> t[1]}.combine(GENES_CH)  
	
	//CURRENT SITUATION!// SO FAR WE HAVE A CHANNEL OF LISTS, EACH LIST CONTAINING: barcodes_i.tsv, matrix_i.tsv, features.tsv. 
	// WE'D LIKE TO SAVE EACH OF THIS LISTS TO A DIR, RENAMING THEIR ELEMENTS.
	
	TEST_DATA.subscribe{ println it } 
	TRAIN_DATA.subscribe{ println it }

	TEST_DATA.map{
		query_dir="$baseDir/mamaso/query_10x_dir/"
		file(it[0]).copyTo(query_dir + 'barcodes.tsv')
		file(it[1]).copyTo(query_dir + 'matrix.mtx')
		file(it[2]).copyTo(query_dir + 'features.tsv')
	}	
	
	TRAIN_DATA.map{
		ref_dir="$baseDir/mamaso/reference_10x_dir/"
		file(it[0]).copyTo(ref_dir + 'barcodes.tsv')
		file(it[1]).copyTo(ref_dir + 'matrix.mtx')
		file(it[2]).copyTo(ref_dir + 'features.tsv')
	}	



	OK! AFTER TRYING LITTERALLY EVERYTHING. THE BEST OPTION I CAN THINK OF IS TO CREATE A MULTIINPUT CHANNEL. INPUT_1: TUPLE_CH, INPUT_2: GENE_CH	
	//process distribute_test_data{
	//	
	//	input:
	//	val x from TEST_DATA
	//	
	//	script:
	//	"""
	//	echo $x
	//	"""
	//}
}

