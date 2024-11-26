#!/bin/bash
#takes BED as input
perl -lane 'if(($F[3]=~/Donor/ && $F[5] eq "+") || ($F[3]=~/Acceptor/ && $F[5] eq "-")){$F[1]=$F[1]+1; $F[2]=$F[2]+1} elsif(($F[3]=~/Donor/ && $F[5] eq "-") || ($F[3]=~/Acceptor/ && $F[5] eq "+")){$F[1]=$F[1]-1; $F[2]=$F[2]-1} print join("\t",@F)'
