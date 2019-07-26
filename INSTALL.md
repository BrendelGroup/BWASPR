# BWASPR Installation and Setup

## Installation as a singularity container

Assuming _git_ and  _singularity_ are installed on your system, you can get the
BWASPR code from GitHub and the container from the
[Singularity Hub](https://www.singularity-hub.org/collections/3301) as follows:

```bash
git clone https://github.com/BrendelGroup/BWASPR.git
cd BWASPR
singularity pull --name bwaspr.simg shub://BrendelGroup/BWASPR
```

For a gentle introduction to singularity, see our group
[handbook article](https://github.com/BrendelGroup/bghandbook/blob/master/doc/06.2-Howto-Singularity-run.md).


## Optional: System-wide Installation

BWASPR use via the singularity container is highly recommended, with no known
drawbacks.
However, if desired, you can of course install all the required R packages
individually on your computer system.
The singularity [recipe file](./Singularity) in this repository should serve as
a guide to perform such an installation.
The `bwaspr.simg` container was built on the 
[current long-term supported Ubuntu 18.04 distribution](https://www.ubuntu.com/download/desktop)
and thus the instructions apply to that particular Linux version.
For different Linux distributions, you will have to install the equivalent
packages using your distribution's package manager.


## Finally

proceed to the [demo/README](./demo/README) document, and we'll tell you how to
execute sample workflows (or, equally easy, your very own data analyses).
