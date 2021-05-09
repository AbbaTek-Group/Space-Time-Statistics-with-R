 wc -l OBS_NODE_mod.OUT | awk '{print $1}' > tmpfile
read val < tmpfile
val2=6012
if [ "$val" -lt "$val2" ]; then
  cp OBS_NODE_blank.OUT OBS_NODE_mod.OUT 
fi

rm tmpfile
