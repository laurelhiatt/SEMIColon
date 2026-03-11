# SEMIColon
Scripts and data for SEMIColon: Somatic Exploration of Mosaicism in Colon

# Contents

We have a scripts folder and a data folder.

Within the scripts folder, there are (planned) subfolders specific to different processes

# Cleaning Branch To-Do's
- [ ] Containerize software
  - [ ] Create conda environments in container
  - [ ] Create environments in container for respective envmodules
- [ ] Hardcode paths
  - [ ] Remove hardcoded paths
  - [ ] Replace hardcoded paths with relative paths, input parameters, etc.
- [ ] Replace redundant copy+paste rules with reused rules using `use rule ... as ... with:`
- [ ] Determine need of symlinks for inSTRbility, MuSiCal.
- [ ] Determine how `figures`, `filtering`, `mutations`, and `pks` relate to the pipeline.
  - [ ] Include `filtering` and `pks` Snakefiles as modules in main Snakefile?
- [ ] 