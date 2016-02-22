#!/bin/bash

declare -a arr=(65536)
make randmst
for i in "${arr[@]}"
do
    ./randmst 0 $i 5 0 >> results.txt
    ./randmst 0 $i 5 2 >> results.txt
    ./randmst 0 $i 5 3 >> results.txt
    ./randmst 0 $i 5 4 >> results.txt
done
