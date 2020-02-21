#!/usr/bin/env bash

cd $EVAL_WORKFLOW_ROOT

# Clone or update the cell-types-eval-workflows repo containing submodules for individual pipelines
if [ ! -d 'cell-types-eval-workflows' ]; then
    git clone --recursive https://github.com/ebi-gene-expression-group/cell-types-eval-workflows
fi

pushd cell-types-eval-workflows > /dev/null
git checkout $CONTROL_WORKFLOW_BRANCH > /dev/null
git pull > /dev/null
git submodule update > /dev/null
popd > /dev/null


# run nextflow command 
nextflow run $CONTROL_WORKFLOW_ROOT/main.nf -profile cluster 