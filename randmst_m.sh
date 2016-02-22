#!/bin/bash

echo "finding max edge weight in mst"
make randmst
./randmst 0 1024 250 0 >> max.txt
./randmst 0 1024 250 2 >> max.txt
./randmst 0 1024 250 3 >> max.txt
./randmst 0 1024 250 4 >> max.txt


