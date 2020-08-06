# Integrative Modeling of Membrane-associated Protein Assemblies: Demo

## 1. Introduction

This is an exhaustive guide for installing and running the software necessary for reproducing the results of the manuscript:

**Integrative Modeling of Membrane-associated Protein Assemblies**<br>
Jorge Roel-Touris, Brian Jiménez-García and Alexandre M.J.J. Bonvin<br>
*bioRxiv* 2020.07.20.211987; doi: [https://doi.org/10.1101/2020.07.20.211987](https://doi.org/10.1101/2020.07.20.211987)

Reproducing all the results of the different simulations would require of HPC capabilities in order to accomplish them in a moderate amount of time. Here, we have set a full demo using *3x29* as a test case in order to demonstrate the **reproducibility** and **repeatability** of our protocol.


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

***pip*** is the reference Python package manager and it is used to install and update Python packages. *pip* is usually installed by default if the Python interpreter is available. In this case, LightDock and its dependencies will be installed via *pip*, so please follow the [official guide](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) if *pip* is not installed in your OS.

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

Here two different simulations are presented for the same demo target *3x29*. Please note that only one CPU core is used for clustering, despite this step could be parallelized too using the script included in LightDock called `ant_thony.py` (see its usage [here](https://lightdock.org/tutorials/membrane)).

#### 3.1.1. Reduced simulation of the *3x29* complex

This is a reduced simulation with very few number of glowworms and steps. **It is intended for testing the workflow and settings, but not for a real production run.**

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

The following configuration is exactly the same as the results reported in the [manuscript](https://www.biorxiv.org/content/10.1101/2020.07.20.211987v1):

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

Inside the [clustered](docking/lightdock/3x29/clustered) folder there is a file containing the models included in `rank_by_scoring.list`, but analyzed compared to the reference crystal structure: [lgd\_clustered\_rank.list](docking/lightdock/3x29/clustered/lgd_clustered_rank.list). Here you can see the top 5 models of the LightDock docking simulation: `lgd_clustered_rank.list` file:

```
#model                Fnat      i-RMSD      l-RMSD     score
swarm_60_174.pdb    0.186667    4.465       12.592     32.364
swarm_64_171.pdb    0.64        1.643        5.028     30.785
swarm_42_162.pdb    0.133333    4.317       13.809     27.664
swarm_58_157.pdb    0.666667    1.293        3.071     26.081
swarm_33_164.pdb    0.426667    1.774        5.108     26.013
```

`Fnat` stands for *fraction of native contacts*, `i-RMSD` for interface-RMSD, `l-RMSD` for ligand-RMSD and `score` is the LightDock scoring (using DFIRE scoring function).

Please see "*Metrics for the evaluation of model quality and success rate*" in the "*Materials and Methods*" section of the manuscript for a description of the different metrics.


## 4. Refinement

The top 100 docked models are refined using HADDOCK to remove potential clashes at the interface. This is a *two-step* procedure: 

- (1) The generation of the all-atom and coarse-grained topologies.
- (2) The coarse-grained refinement of the generated models.

In [here](refinement/haddock/) you can find all the data needed to perform the refinement of the top100 docked models for the *3x29* example. In `run.param` file, you have to define a number of variables as:

- `HADDOCK_DIR`: The directory of your HADDOCK local instalation.
- `N_COMP`: The number of components.
- `PDB_FILE1`: The path to your first *all-atom* receptor PDB file.
- `PDB_LIST1`: A single column file listing all of your *all-atom* receptor PDB files.
- `CGPDB_FILE1`: The path to your first *coarse-grained* receptor PDB file.
- `CGPDB_LIST1`: A single column file listing all of your *coarse-grained* receptor PDB files.
- `PROT_SEGID_1`; The segid record to be used during the simulation of your first (receptor) component.
- `PDB_FILE2`: The path to your first *all-atom* ligand PDB file.
- `PDB_LIST2`: A single column file listing all of your *all-atom* ligand PDB files.
- `CGPDB_FILE2`: The path to your first *coarse-grained* ligand PDB file.
- `CGPDB_LIST2`: A single column file listing all of your *coarse-grained* ligand PDB files.
- `PROT_SEGID_2`: The segid record to be used during the simulation of your second (ligand) component.
- `CGTOAA_TBL`: A HADDOCK-like restraints file of the CG to AA mapping.
- `PROJECT_DIR`: The directory of your run.
- `RUN_NUMBER`: A run number.

`receptor` and `ligand`folders contain the top 100 receptor and ligand PDB files respectively without the bead bilayer. The bead bilayer is not required for the coarse-grained refinement since the proteins have been already docked.

### 4.1. Generation of topologies

For this task, please download all the content of the [refinement](refinement/haddock/) folder as:

```bash
cd ~
mkdir refinement
cd refinement
curl -O https://raw.githubusercontent.com/lightdock/membrane_docking/master/demo/refinement/haddock/*
```

First, you need to execute HADDOCK once as:

```bash
haddock2.4
```

Then, in order to edit the `run.cns` to generate the topologies, execute the `generate_toppar.sh` as:

```bash
./generate_toppar.sh
```

And execute HADDOCK again as:

```bash
cd run1
haddock2.4 &> haddock.out &
```

Check that the **all-atom** and **coarse-grained** topologies have been generated without issues in `begin-aa` and `begin` folders respectively.

### 4.2. Coarse-grained refinement

For this, we need to modify a handful of parameters within the HADDOCK parameter file `run.cns` including:

- `rotate180_it0=false` (to skip sampling 180° complementary interfaces)
- `crossdock=false` (to refine receptor – ligand from the structures provided)
- `rigidmini=false` (to skip it0 stage)
- `randorien=false` (to skip it0 stage)
- `rigidtrans=false` (to skip it0 stage)
- `ntrials=1` (to skip it0 stage)
- `structures_0=100` (for it0 stage)
- `structures_1=100` (for it1 stage; must always be ≤ than structures_0)
- `anastruc_1=100` (for analysis purposes at it1 stage)
- `waterrefine=100` (for itw stage; this is the number of final output models)
- `initiosteps=0` (to skip it1 stage)
- `cool1_steps=0` (to skip it1 stage)
- `cool2_steps=0` (to skip it1 stage)
- `cool3_steps=0` (to skip it1 stage)
- `dielec_0=cdie` (to switch a constant dieletric constant when CG is used)
- `dielec_1=cdie` (to switch a constant dieletric constant when CG is used)

For your convenience, you can just execute the `cg_refinement.sh` script as:

```bash
./cg_refinement.sh
```

And execute HADDOCK again:

```bash
cd run1
haddock2.4 &> haddock.out &
```

Once the simulation is done, you will find the refined models under `run1/structures/it1/water` with `run1/structures/it1/water/file.list` containing the ranking according to HADDOCK score.
