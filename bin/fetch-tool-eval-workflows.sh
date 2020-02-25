#!/usr/bin/env bash
EVAL_WORKFLOW_ROOT="$PWD"
EVAL_WORKFLOWS="$PWD/cell-types-eval-workflows"


# Clone or update the cell-types-eval-workflows repo containing submodules for individual pipelines
if [ ! -d 'cell-types-eval-workflows' ]; then
    git clone --recursive https://github.com/ebi-gene-expression-group/cell-types-eval-workflows $EVAL_WORKFLOWS
fi

pushd $EVAL_WORKFLOWS > /dev/null
git checkout develop > /dev/null
git pull origin develop > /dev/null
git submodule update --recursive --remote > /dev/null
popd > /dev/null
