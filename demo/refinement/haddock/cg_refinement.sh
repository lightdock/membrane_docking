#!/bin/bash

rm -rf run1/haddock.out

sed -i "s/structures_0=0/structures_0=100/g" run1/run.cns
sed -i "s/structures_1=0/structures_1=100/g" run1/run.cns
sed -i "s/anastruc_1=0/anastruc_1=100/g" run1/run.cns
sed -i "s/waterrefine=0/waterrefine=100/g" run1/run.cns
sed -i "s/rotate180_it0=true/rotate180_it0=false/g" run1/run.cns
sed -i "s/crossdock=true/crossdock=false/g" run1/run.cns
sed -i "s/randorien=true/randorien=false/g" run1/run.cns
sed -i "s/rigidmini=true/rigidmini=false/g" run1/run.cns
sed -i "s/rigidtrans=true/rigidtrans=false/g" run1/run.cns
sed -i "s/ntrials=5/ntrials=0/g" run1/run.cns
sed -i "s/initiosteps=500/initiosteps=0/g" run1/run.cns
sed -i "s/cool1_steps=500/cool1_steps=0/g" run1/run.cns
sed -i "s/cool2_steps=1000/cool2_steps=0/g" run1/run.cns
sed -i "s/cool3_steps=1000/cool3_steps=0/g" run1/run.cns
sed -i "s/dielec_0=rdie/dielec_0=cdie/g" run1/run.cns
sed -i "s/dielec_1=rdie/dielec_1=cdie/g" run1/run.cns
