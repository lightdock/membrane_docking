#!/bin/bash

################
# You may change these variables according to your needs
SWARMS=400
GLOWWORMS=200
STEPS=100
CORES=24
################

# Setup
if test -f "restraints.list"; then
    lightdock3_setup.py receptor.pdb ligand.pdb ${SWARMS} ${GLOWWORMS} --noxt --noh -rst restraints.list
else
    lightdock3_setup.py receptor.pdb ligand.pdb ${SWARMS} ${GLOWWORMS} --noxt --noh
fi

# Simulation
lightdock3.py setup.json ${STEPS} -s fastdfire -c ${CORES}

# Generate predictions in PDB format
s=`ls -d swarm_* | wc -l`
swarms=$((s-1))
for i in $(seq 0 $swarms)
  do
    cd swarm_${i}; lgd_generate_conformations.py ../receptor.pdb ../ligand.pdb  gso_${STEPS}.out ${GLOWWORMS}; cd ..;
  done

# Cluster per swarm
for i in $(seq 0 $swarms)
  do
    cd swarm_${i}; lgd_cluster_bsas.py gso_${STEPS}.out; cd ..;
  done

# Generate ranking
lgd_rank.py ${swarms} ${STEPS}

echo "Done."
