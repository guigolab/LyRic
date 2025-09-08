#!/bin/bash -i

# remove cluster logs from previous runs so they don't accumulate:
\mv LyRic/qsub_logs LyRic/qsub_logs_bkp; mkdir LyRic/qsub_logs; rm -rf LyRic/qsub_logs_bkp/ &

# activate snakemake conda environment: 
conda activate snakemake; 

# launch LyRic workflow in cluster/DRMAA mode:
snakemake -p --reason --latency-wait 100 --use-conda --configfile LyRic/config_EXAMPLE.json -s LyRic/master.smk -j 500 --jobname {rulename}.{jobid}._{config["PROJECT_NAME"]}_.sh --cluster-config LyRic/cluster_config.json --max-jobs-per-second 5 --drmaa " -V -q {cluster.queue} -l disk={cluster.disk} -l virtual_free={cluster.virtual_free} -l h_rt={cluster.h_rt}  -o {cluster.out} -e {cluster.err} {cluster.threads} -P {cluster.project}" --rerun-incomplete --keep-going  --show-failed-logs

#generate DAG to visualize workflow as SVG image:
snakemake -p --configfile LyRic/config_EXAMPLE.json -s LyRic/master.smk --forceall --rulegraph > EXAMPLE_dag.dot
cat EXAMPLE_dag.dot| dot -Tsvg > EXAMPLE_dag.svg
