# Integrative Modeling of Membrane-associated Protein Assemblies: Demo

## 1. Introduction

This is an exhaustive guide for installing and running the software necessary for reproducing the results of the manuscript:

**Integrative Modeling of Membrane-associated Protein Assemblies**<br>
Jorge Roel-Touris, Brian Jiménez-García and Alexandre M.J.J. Bonvin<br>
*bioRxiv* 2020.07.20.211987; doi: [https://doi.org/10.1101/2020.07.20.211987](https://doi.org/10.1101/2020.07.20.211987)

Reproducing all the results of the different simulations would require of HPC capabilities in order to accomplish them in a moderate short time. To overcome this problem, we have set a full demo using the smallest complex (*3x29*) in our [dataset and the default scenario](../docking/lightdock/membrane/3x29) which will demonstrate the **reproducibility** and **repeatability** of our protocol.


## 2. Installation

### 2.1. LightDock

Lightdock software is compatible and it has been tested with the followings Operating Systems:

- **macOS**: El Capitan, Sierra, High Sierra, Mojave, Catalina.
- **GNU/Linux**: Ubuntu 16+, Debian Stretch+, Scientific Linux 6+, CentOS 6+.

**Microsoft Windows is not supported**, despite many parts of the protocol might be able to run. Please use it at your own risk.

#### 2.1.1. Python

Lightdock software is written in the version 3 of the Python programming language. Python3 supported versions by LightDock are 3.6, 3.7, 3.8 and 3.9.

Current versions of macOS and GNU/Linux ship by default Python3 interpreter. To check your current version, please type in the command line terminal:

```bash
python -V
```

or 

```bash
python3 -V
```

If Python3 is not installed in your OS, there are different options to install it:

- You may try download and installing it from the [official source](https://www.python.org/downloads/).
- In macOS, there are alternatives using [Homebrew](https://brew.sh/) or [MacPorts](https://www.macports.org/).
- In GNU/LINUX, you may use the system package manager to install it. For example, in Ubuntu with `sudo apt-get update && apt-get install python3`.

#### 2.1.2. *pip* and virtual environments

***pip*** is the reference Python package manager. It's used to install and update Python packages. *pip* is usually installed by default if the Python interpreter is available. We will need it to install LightDock and its dependencies. Please follow the [official guide](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) if *pip* is not installed in your OS.

#### 2.1.3. Installing LightDock via *pip*

We will use *pip* and *venv* to create a "virtual" isolated Python installation and install packages into that virtual installation. Then, we will install LightDock in that virtual environment.

* First, we create and activate the virtual environment:

```bash
cd ~
python3 -m venv lightdock-env
source lightdock-env/bin/activate
```
* Finally, we will install `numpy` Python package and then the `lightdock` Python package:

```bash
cd lightdock-env
pip3 install numpy
pip3 install lightdock
```

If installation was successful, you should see a similar output to this if you try to execute `lightdock_3_setup.py` command:

```bash
(lightdock-env) user@machine ~/lightdock-env # lightdock3_setup.py 
usage: lightdock_setup [-h] [--seed_points STARTING_POINTS_SEED]
                       [-ft ftdock_file] [--noxt] [--noh] [--verbose_parser]
                       [-anm] [--seed_anm ANM_SEED] [-anm_rec ANM_REC]
                       [-anm_lig ANM_LIG] [-rst restraints] [-membrane]
                       receptor_pdb_file ligand_pdb_file swarms glowworms
lightdock_setup: error: the following arguments are required: receptor_pdb_file, ligand_pdb_file, swarms, glowworms
```

### 2.2. HADDOCK

Please visit the [tutorial describing the use of a local version of HADDOCK2.4](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-local-tutorial/) for an exhaustive guide on how to install and how to use the local version of the HADDOCK2.4 software.


## 3. Docking

For the *3x29* demo, we will download the input structures and a bash script to easily run the simulation and analysis pipeline:

```bash
mkdir demo
cd demo
curl -O https://raw.githubusercontent.com/lightdock/membrane_docking/master/docking/lightdock/membrane/3x29/run.sh
curl -O https://raw.githubusercontent.com/lightdock/membrane_docking/master/demo/docking/lightdock/receptor_membrane.pdb
curl -O https://raw.githubusercontent.com/lightdock/membrane_docking/master/demo/docking/lightdock/ligand.pdb
chmod u+x run.sh
```

Please note that the structures are simulation-ready and have been prepared according to the [**Membrane-associated protein docking**](https://lightdock.org/tutorials/membrane) tutorial.

The `run.sh` script might be tunned in the main variables section. The values shown here are the ones used for all the results in the manuscript:

```bash
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
    lightdock3_setup.py receptor_membrane.pdb ligand.pdb ${SWARMS} ${GLOWWORMS} --noxt --noh -membrane -rst restraints.list
else
    lightdock3_setup.py receptor_membrane.pdb ligand.pdb ${SWARMS} ${GLOWWORMS} --noxt --noh -membrane
fi

# Simulation
lightdock3.py setup.json ${STEPS} -s fastdfire -c ${CORES}

# Generate predictions in PDB format
s=`ls -d swarm_* | wc -l`
swarms=$((s-1))
for i in $(seq 0 $swarms)
  do
    cd swarm_${i}; lgd_generate_conformations.py ../receptor_membrane.pdb ../ligand.pdb  gso_${STEPS}.out ${GLOWWORMS}; cd ..;
  done

# Cluster per swarm
for i in $(seq 0 $swarms)
  do
    cd swarm_${i}; lgd_cluster_bsas.py gso_${STEPS}.out; cd ..;
  done

# Generate ranking
lgd_rank.py ${swarms} ${STEPS}

echo "Done."
```

For running a simulation, it is as easy as running this script:

```bash
time ./run.sh
```

**Before running this demo simulation, please see next section "3.1. Timing" to get an expected completion time**.

### 3.1. Timing

Here two different simulations are presented for the same demo target *3x29*. Please note that only one CPU core is used for clustering for the sake of simplicity of the `run.sh` script, despite this step could be parallelized too using the script included in LightDock called `ant_thony.py` (see its usage [here](https://lightdock.org/tutorials/membrane)).

#### 3.1.1. Reduced simulation of the *3x29* complex

This is a reduced simulation with very few number of glowworms and steps. This should not be considered as a biologically relevant simulation, but as a simple example to get the feeling of a complete and fast simulation.

**Variables:**

- CPU: 8 cores Intel(R) i7 CPU X980 @ 3.33GHz
- Swarms: 400, after setup number of swarms to simulate is 72.
- Glowworms: 20
- Steps: 10

|   Step     |      Time     | Cores  |
|------------|:-------------:|:------:|
| setup      |  15.598s      | 1      |
| simulation |  11.799s      | 8      |
| clustering |  3m10.348s    | 1      |


#### 3.1.2. Manuscript simulation of the *3x29* complex

The following configuration is exactly the same as the results reported in the manuscript:

**Variables:**

- CPU: 8 cores Intel(R) i7 CPU X980 @ 3.33GHz
- Swarms: 400, after setup number of swarms to simulate is 72.
- Glowworms: 200
- Steps: 100

|   Step     |      Time     | Cores  |
|------------|:-------------:|:------:|
| setup      |  16.685s      | 1      |
| simulation |  74m53.112s   | 8      |
| clustering |  28m27.955s   | 1      |


### 3.2. Results

Precalculated results for your convenience are available [here](docking/lightdock/3x29/).

The file [rank\_by\_scoring.list](docking/lightdock/3x29/rank_by_scoring.list) is a list of the top clustered simulation predicted models ranked by scoring function (the higher the scoring term the better).

Inside the [clustered](docking/lightdock/3x29/clustered) folder there is a file containing the models included in `rank_by_scoring.list`, but analyzed compared to the reference crystal structure: [lgd\_clustered\_rank.list](docking/lightdock/3x29/clustered/lgd_clustered_rank.list). Here you can see the top 5 models according to LightDock from `lgd_clustered_rank.list` file:

```
#model             Fnatt      i-RMSD      l-RMSD     score
swarm_60_174.pdb  0.186667    4.465       12.592     32.364
swarm_64_171.pdb  0.64        1.643        5.028     30.785
swarm_42_162.pdb  0.133333    4.317       13.809     27.664
swarm_58_157.pdb  0.666667    1.293        3.071     26.081
swarm_33_164.pdb  0.426667    1.774        5.108     26.013
```


## 4. Refinement
