#!/bin/bash

declare -a arr=(16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536)
rm results.txt
make randmst
for i in "${arr[@]}"
do
    ./randmst 0 $i 5 0 >> results.txt
    ./randmst 0 $i 5 2 >> results.txt
    ./randmst 0 $i 5 3 >> results.txt
    ./randmst 0 $i 5 4 >> results.txt
done