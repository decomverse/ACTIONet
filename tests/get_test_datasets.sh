#!/bin/bash
#copies benchmark datasets from S3 (insitro root)

aws s3 cp --recursive s3://insitro-tardis/benchmark_datasets .
