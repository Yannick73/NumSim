#!/bin/bash

for i in 24 20 16 12 8 6 4 2 1
do
        echo
        echo "SIMULATE WITH $i RANK and async"
        echo
        mpirun -n $i -oversubscribe /home/yannick/Git/numeric-simulation/Ex2/build/numsim_parallel input/scenario3_6.txt
        wait
done
