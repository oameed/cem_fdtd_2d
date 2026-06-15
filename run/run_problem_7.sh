#! /bin/bash

source activate lnxpy

cd     run

python fdtd_2d_tfsf.py -p 'problem_3' -v '07'
