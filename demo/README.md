# Integrative Modeling of Membrane-associated Protein Assemblies: Demo

## 1. Introduction

This is an exhaustive guide for installing and running the software necessary for reproducing the results of the manuscript:

**Integrative Modeling of Membrane-associated Protein Assemblies**<br>
Jorge Roel-Touris, Brian Jiménez-García and Alexandre M.J.J. Bonvin<br>
*bioRxiv* 2020.07.20.211987; doi: [https://doi.org/10.1101/2020.07.20.211987](https://doi.org/10.1101/2020.07.20.211987)

Reproducing all the results of the different simulations would require of HPC capabilities in order to accomplish them in a moderate short time. To overcome this problem, we have set a full demo using the smallest complex (*3x29*) in our [dataset and the default scenario](../docking/lightdock/membrane/3x29) which will demonstrate the **reproducibility **and **repeatability** of our protocol.


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
pip install numpy
pip install lightdock
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
