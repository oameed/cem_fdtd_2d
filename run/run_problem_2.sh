#! /bin/bash

source activate lnxpy

cd     run

python fdtd_2d_ps.py -p 'problem_2' -v '02' -t 800
