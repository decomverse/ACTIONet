Benchmark datasets can be obtained from: s3://insitro-tardis/benchmark_datasets/

* CellSIUS dataset
```
aws s3 cp --recursive s3://insitro-tardis/benchmark_datasets/CellSIUS_Wegmann19 .
```
** for evaluating clustering approaches (reduction --> network construction -->  visualization --> clustering )

* fetal_liver_100k+
```
aws s3 cp --recursive s3://insitro-tardis/benchmark_datasets/fetal_liver_100k+ .
```
** For testing scalability

* PBMC_Granja19
```
aws s3 cp --recursive s3://insitro-tardis/benchmark_datasets/PBMC_Granja19 .
```
** Pre-annotated dataset for computing marker genes

* Monaco
```
aws s3 cp --recursive s3://insitro-tardis/benchmark_datasets/Monaco .
```
** Bulk RNA-seq datasets
