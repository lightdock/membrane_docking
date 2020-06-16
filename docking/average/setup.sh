#!/bin/bash

for i in `cat ../../pdbs.list`;do mkdir $i; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; scp -r brianj@alcazar:/home/jroel/Light-memb/benchmark/lightdock/membrane-AVE/${i}/setup.json .; cd ..; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; scp -r brianj@alcazar:/home/jroel/Light-memb/benchmark/lightdock/membrane-AVE/${i}/lightdock_*.pdb* .; cd ..; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; scp -r brianj@alcazar:/home/jroel/Light-memb/benchmark/lightdock/membrane-AVE/${i}/ligand.pdb .; cd ..; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; scp -r brianj@alcazar:/home/jroel/Light-memb/benchmark/lightdock/membrane-AVE/${i}/receptor_membrane.pdb .; cd ..; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; scp -r brianj@alcazar:/home/jroel/Light-memb/benchmark/lightdock/membrane-AVE/${i}/restraints.list .; cd ..; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; scp -r brianj@alcazar:/home/jroel/Light-memb/benchmark/lightdock/membrane-AVE/${i}/clustered/lgd_*.list .; cd ..; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; head -100 lgd_clustered_rank.list | awk '{print $1}' > top100.list; cd ..; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; mkdir top100; cd top100;for j in `cat ../top100.list`;do echo $j; scp -r brianj@alcazar:/home/jroel/Light-memb/benchmark/lightdock/membrane-AVE/${i}/clustered/${j} .; done; cd ..; cd ..; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; mv top100.list top100; tar zcf top100.tgz top100; cd ..; done
for i in `cat ../../pdbs.list`;do echo $i; cd $i; rm -rf top100; chmod 0664 *; cd ..; done
