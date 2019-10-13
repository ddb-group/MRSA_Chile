Some commands to run the pipelines

Dry run:

```
snakemake -n -r -p -j 4 --local-cores 4 -w 90 --max-jobs-per-second 8 --cluster-config cluster.json --cluster "qsub -N {cluster.name} -v PATH -S /bin/sh -m beas -M juan@ugalde.bio -q externos.q -cwd -o {cluster.output} -e {cluster.error}"
```
Official run:
```
snakemake -j 3 --local-cores 4 -w 90 --max-jobs-per-second 8 --cluster-config cluster.json --cluster "qsub -N {cluster.name} -v PATH -S /bin/sh -m beas -M juan@ugalde.bio -q externos.q -cwd -o {cluster.output} -e {cluster.error}"

```

```
snakemake -j 4 --local-cores 4 -w 90 --max-jobs-per-second 8 --cluster-config cluster.json --cluster "qsub -N {cluster.name} -v PATH -S /bin/sh -m beas -M juan@ugalde.bio -q icim.q -cwd -o {cluster.output} -e {cluster.error}" --use-conda

```