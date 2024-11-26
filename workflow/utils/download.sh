#!/bin/bash

set -e

user=$1
password=$2
host=$3
suffixUrl=$4

echoerr "### Downloading..."
echoerr "###   Source: $host/$suffixUrl"
echoerr "###   Destination: $host/$suffixUrl" 
#install -D /dev/null $host/$suffixUrl; 
#lftp -u $user,$password -e "pget -c -n5 $host/$suffixUrl -o .$suffixUrl; quit "; 
#remove file if exists. this is safer.
rm -rf $host/$suffixUrl
rm -rf $host/$suffixUrl.md5
wget  --no-verbose --continue --mirror --header="X-Auth-Challenge: true"  --tries=10 --user=$user --password=$password $host/$suffixUrl;
#wget --auth-no-challenge --user=$user --password=$password --tries=50 $host/$suffixUrl 
cd `dirname $host/$suffixUrl`;
md5sum  `basename $host/$suffixUrl` > `basename $host/$suffixUrl`.md5
echoerr \"$host/$suffixUrl :: XXdoneXX\"; 
