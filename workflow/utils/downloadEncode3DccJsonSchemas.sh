#!/bin/bash

#set -e
source $HOME/.bashrc
cd $ENCODE3_JSONSCHEMAS_DIR
echo "###############################"
printf "now: "; now 
for baseUrl in `cat urls.list | skipcomments`; do
echo $baseUrl
for json_object in `cat jsonSchemas.list`; do
url="https://${baseUrl}/profiles/${json_object}"
\mv $baseUrl/${json_object} $baseUrl/${json_object}.bkp
processEncode3DccJsonObject.pl $url GET $baseUrl/${json_object}
done
done
printf "now: "; now
echo "All done!"
exit 0;
