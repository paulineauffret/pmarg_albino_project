#!/usr/bin/env bash
. ~/snakemake/5.4.0/env.sh

snakemake $@ --verbose --cluster-config cluster.yml --configfile cluster.yml -p --jobs 10 --cluster "qsub -S '/bin/bash' -N {rule} -m n -q {cluster.queue} -l ncpus={cluster.ncpus} -l walltime={cluster.walltime} -l mem={cluster.memory}" --latency-wait 420

. ~/snakemake/5.4.0/delenv.sh
