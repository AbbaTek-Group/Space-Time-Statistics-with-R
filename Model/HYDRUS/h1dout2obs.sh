awk 'NR<12' OBS_NODE.OUT | awk '{print $0}' > OBS_NODE_mod.OUT
awk 'NR>12' OBS_NODE.OUT | awk '{if ($1 % 1.0 == 0.) print $0}' >> OBS_NODE_mod.OUT
