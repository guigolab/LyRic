#!/bin/bash
dir=$1
dir=$(echo $dir | sed 's/\///g')
echo $dir
tar -cf $dir.tgz $dir;
echo XXdoneXX
