# NorfolkIslandHybridization
 Scripting for Norfolk Island hybridization project


## Setting up conda environment for the first time
```
mamba create -n NorfolkIslandHybridisation
mamba activate NorfolkIslandHybridisation
```

## Install required programmes
```
mamba install -c bioconda fastp
mamba install -c bioconda bwa-mem2
```

## Directly compiled programmes
```
#Make 'bin' for compiled programmes
mkdir bin
cd bin

#Download and compile DSuite
git clone https://github.com/millanek/Dsuite.git
cd Dsuite
make
cd ../

```
