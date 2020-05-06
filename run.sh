#!/bin/bash
if [[ "${COMP_S3_PROBLEM_PATH}" == *".xz" ]];
then
    aws s3 cp s3://${S3_BKT}/${COMP_S3_PROBLEM_PATH} supervised-scripts/task.cnf.xz
    unxz supervised-scripts/task.cnf.xz
else
    aws s3 cp s3://${S3_BKT}/${COMP_S3_PROBLEM_PATH} supervised-scripts/task.cnf
fi

painless/painless -c=63 -solver=maple -wkr-strat=4 -split-heur=2 -copy-mode=2 -shr-strat=1 -lbd-limit=4  -init-act=1 -init-pol=1 -gaussian supervised-scripts/task.cnf
