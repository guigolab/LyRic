#!/bin/bash

yearmonth=$(date +%Y%m)
jobID=$1

qacct -j $jobID
if [ $? != 0 ]; then
echo "Trying /var/user_accounting/uge_accounting_$yearmonth"
qacct -j $jobID -f /var/user_accounting/uge_accounting_$yearmonth
fi
