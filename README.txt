My PC environment
- Windows 10 Home
- WSL2 (windows subsystems on linux)
- ubuntu LTS20.04
- R 4.0.5
- rstudio server 1.3.1093
- cmdstanr 0.3.0.9000

(2) download the Rotterdam Eye Study data from http://www.rodrep.com/

(3) create a project at Rstudio (e.g. 'BiacedDie') and put files as following.

BiacedDie/
├── 01_enumerate_response_pattern.R
├── 02_estimates_Henson_FOS.R
├── 03_estimate_RODREP_FOS.R
├── BiacedDie.Rproj
├── biaced_die_model.stan
├── explore_OPI_package_(optional).R
├── full_threshold_2.R
├── original_data
│   ├── Bryan2013.csv
│   ├── CHANGELOG.txt
│   ├── DESCR.txt
│   ├── Erler2014.csv
│   ├── FORMAT.txt
│   ├── LICENSE.txt
│   ├── Patients.csv
│   ├── README.txt
│   ├── VFPoints.csv
│   └── VisualFields.csv
└── results (empty folder)


(5) execute .R files

01_enumerate_response_pattern.R
02_estimates_Henson_FOS.R
03_estimate_RODREP_FOS.R

All codes are under MIT licence except for full_threshold_2.R (GPL licence)
