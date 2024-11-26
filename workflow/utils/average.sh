#!/bin/bash
#prints average of a series of numbers given as input (one column, one number per line)
awk '{f += $1}END{ printf("%.9f\n",f/NR)}' |tail -n1
