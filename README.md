# HfStudiesALICE3
Repository for HF performance studies with ACTS for ALICE3

# install ACTS with aliBuild
```bash
mkdir -p alice
cd alice
aliBuild init O2@dev
```
Currently a small modification is needed in `Detectors/Upgrades/ALICE3/TRK/reconstruction/CMakeLists.txt`, i.e. lines `21` and `50` should be commented out.
```
git clone https://gitlab.cern.ch/hepmc/HepMC3.git
aliBuild build O2 --defaults o2-acts
```

# enter the environment
```bash
alienv enter O2/latest ACTS/latest
```

# run the simulation
```bash
./run_sim.sh config.yml
```
examples of config files can be found in the `config` directory

