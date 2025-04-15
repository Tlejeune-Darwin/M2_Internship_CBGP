# M2_Internship_CBGP
 This project aims to link the two disciplines of population dynamics and population genetics within the same model. This link will ultimately make it possible to study biological processes that can bias estimates between the effective size and the actual size of the populations studied.

### LINUX USE ###

1. Download the "Linux" branch 
2. Set the directory :
    Documents/GitHub/M2_Internship_CBGP/Project_Model
3. Get the executables working with : 
    3.1. chmod +x slim
    3.2. chmod +x Ne2x
    3.3. chmod +x *.sh
4. Launch the main script : python3 Scripts/Python/main_linux.py

### PYTHON CONFIGURATION ###

Install dependencies via pip:
    pip install tskit pyslim msprime numpy matplotlib pandas
or use the requirements.txt in the "Project_Model" directory : 
    pip install -r requirements.txt

### SLIM & NEESTIMATOR INSTALLATION ###

There is no real installation for this part. Just make sure that the executables are in the right directory (Bin) and are working (see point 3. of the LINUX USE part)

### R USE ### 

Depending on the Rstudio version (the version is recommanded when you launch Rstudio even if you don't have the right version)
1. Download and install R : 
    sudo apt update && sudo apt upgrade -y
    sudo apt install -y software-properties-common dirmngr gnupg apt-transport-https ca-certificates
    sudo gpg --keyserver keyserver.ubuntu.com --recv-key 'E298A3A825C0D65DFD57CBB651716619E084DAB9'
    sudo gpg --export 'E298A3A825C0D65DFD57CBB651716619E084DAB9' | sudo tee /etc/apt/trusted.gpg.d/cran.gpg > /dev/null
    sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
    sudo apt update
    sudo apt install -y r-base
2. Download and install Rstudio : 
    wget https://download1.rstudio.org/electron/jammy/amd64/rstudio-2024.04.1-748-amd64.deb
    sudo apt install ./rstudio-2024.04.1-748-amd64.deb
3. Install main compilation libraries : 
    sudo apt install -y \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libgit2-dev \
        libglpk-dev \
        libgmp-dev
4. Install tidyverse package directly from Rstudio

### SIMULATIONS FOLDER ###

All simulations are saved to `~/Bureau/simulations/` or `~/Desktop/simulations/` depending on your system language.

Each simulation creates:
- A folder `sim_YYYYMMDD_HHMMSS_local/`
- A config file used by SLiM (`slim_config.txt`)
- A NeEstimator info file (`info`, `option`)
- Outputs: `.gen`, `.trees`, `summary.txt`, and optional logs
- All summaries are aggregated into: `summary_table.csv` (It's the part exploited by the R script)
- Other files are created during the simulation but are deleted at the end to gain space (see "run_simulation_linux.py" (section 14.1))

### TROUBLESHOOTING ###

Every errors I encountered would be present here (if I remember to add them)
- Use "python" instead of "python3" for Windows
- Watch out for OneDrive directories
- 