#!/bin/bash

set -e

# this function must be called like: 
# source $function
# if [ "$(check_exit)" != "0" ]; then echo "ERROR" >&2; exit 1; fi


function check_exit {
## check the pipe status to ensure that piped command
## succeeded
local BAK=(${PIPESTATUS[@]})
for s in "${BAK[@]}"; do
if [ "$s" != "0" ]; then
echo "bash pipe failed: ${BAK[@]}" 1>&2;
echo $s;
break;
fi;
done
echo 0
}
