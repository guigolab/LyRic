#!/bin/bash

awk -v no=$2 '{if(no=="") print $1; else print $no}' $1 > /tmp/coucoufromsarah.txt
echo 'l<-scan("/tmp/coucoufromsarah.txt");summary(l)' | R --vanilla --slave
rm /tmp/coucoufromsarah.txt
