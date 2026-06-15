#! /bin/bash

source activate lnxpy

cd     run

python fdtd_1d.py -p 'problem_13' -v '13' -t 1500
