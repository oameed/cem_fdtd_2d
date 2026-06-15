#! /bin/bash

source activate lnxpy

cd     run

python fdtd_2d_tfsf.py -p 'problem_4' -v '08'
