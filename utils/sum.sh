awk 'BEGIN { SUM=0 } { SUM += ($1);} END { print SUM }'
