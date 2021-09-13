#!/bin/bash
source ~jlagarde/.bashrc
cd $ENCODE3_COLLECTIONS_DIR
set -e;
#get credentials
encodeUser=`awk '{print $1}' $ENCODE3_CREDENTIALS_FILE`
encodePass=`awk '{print $2}' $ENCODE3_CREDENTIALS_FILE`


echo "###############################" >&2
printf "now on $HOSTNAME: " >&2
now >&2
for baseUrl in `cat urls.list |skipcomments`; do
echo $baseUrl >&2
cd $ENCODE3_COLLECTIONS_DIR
echo "GETting collections from $baseUrl" >&2
rm -f download.sh
# for collection in `cat collections.list|skipcomments |cut -f1`; do
cat collections.list|skipcomments | while read collection term; do
  #url="https://${encodeUser}:${encodePass}@${baseUrl}/$collection"
#  url="https://${baseUrl}/$collection"
url="https://${baseUrl}/search/?type=$term"
touch $baseUrl/${collection}.json
\cp -p $baseUrl/${collection}.json $baseUrl/${collection}.json.bkp
echo -e "\nURL: ${url}&frame=object&limit=all&format=json" >&2
#processEncode3DccJsonObject.pl $url GET $baseUrl/${collection}.json
echo "wget --no-verbose --auth-no-challenge --user=$encodeUser --password=$encodePass --tries=50 \"${url}&frame=object&limit=all&format=json\" -O $baseUrl/${collection}.json ; cat $baseUrl/${collection}.json |json_xs > tmp; \mv tmp $baseUrl/${collection}.json" >> download.sh
#echo "Flattening json collections" >&2
#encode3jsonToIndexFile.pl $baseUrl/${collection}.json uuid > $baseUrl/`basename $baseUrl/${collection}.json .json`.deep.json.flat.txt;
 done

bash download.sh
echo "done downloading" >&2
echo "Merging all collections into one (indexed by uuid and accession)">&2
 makeFullEncode3JsonCollection.pl $baseUrl | json_xs> $baseUrl/all_objects_indexed_by_uuid_accession_alias.json;
echo "Done" >&2


done
printf "now: " >&2
now >&2
echo "All done!">&2
exit 0;
