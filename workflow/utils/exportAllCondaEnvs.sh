#!/bin/bash


source ~/.bashrc
# modified from https://github.com/conda/conda/issues/5165

EXPORT_DIR="$HOME/condaEnvExport/"
mkdir -p $EXPORT_DIR
ENVS=$(conda env list | grep '^\w' | cut -d' ' -f1)
for env in $ENVS; do
	echo "Exporting $env..."
    conda activate $env
#    conda env export |grep -vP "^prefix: " > $EXPORT_DIR/$env.with-build.yml
    conda env export --no-builds |grep -vP "^prefix: " > $EXPORT_DIR/$env.yml
    echo "Done."
done
echo "Cleaning up conda install..."
conda clean --all --yes --quiet
echo "Done"
