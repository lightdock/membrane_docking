#!/bin/bash

sed -i "s/structures_0=1000/structures_0=0/g" run1/run.cns
sed -i "s/structures_1=200/structures_1=0/g" run1/run.cns
sed -i "s/anastruc_1=200/anastruc_1=0/g" run1/run.cns
sed -i "s/waterrefine=200/waterrefine=0/g" run1/run.cns
