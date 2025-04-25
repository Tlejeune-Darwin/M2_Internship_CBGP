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

### PYTHON INSTALLATION/CONFIGURATION ###

1. In the command console, install python : 
    sudo apt update
    sudo apt install python3 python3-pip python3-venv
2. Install dependencies via pip:
    pip install tskit pyslim msprime numpy matplotlib pandas
3. It then can be used by typing "python3" before the script you want to run
4. Python needs some packages that need to be installed beforehand to run a simulation without error : 
    numpy
    pandas
    tskit
    pyslim
    You can use the "requirements.txt" provided on the GitHub to install those packages quickly : 
        pip install requirements.txt

### SLIM INSTALLATION ###

1. Download / Installation : 
    Go to the website "Messerlab" or the GitHub of the same name for downloading the corresponding files.
    Please follow the installation section in the manual made by Benjamin C. Haller & Philipp W. Messer for SLiM.
    The manual is a really great source of informations for every script you would like to make and can explain how to SLiM for Windows/Mac/Linux.
2. Make sure your executable is activated by moving to your executable directory and typing : 
    chmod +x "Executable_Name" (Usually "slim")

### NEESTIMATOR ###

1. Download : 
    Please go to the GitHub associated (Named "NeEstimator2.X") to download the .zip file.
    From there, depending of your OS, you will have to search a way to install the executable.
2. For Linux/Ubuntu : 
    2.1. first install base tools for C++ projects :
        sudo apt update
        sudo apt install build-essential cmake g++
    2.2. Then download the from the GitHub : 
        git clone https://github.com/bunop/NeEstimator2.X
        cd NeEstimator
    2.3. Eventually, you will have to compile with cmake : 
        mkdir build
        cd build
        cmake ..
        make
3. Make sure your executable is activated by moving to your executable directory and typing : 
    chmod +x "Executable_Name" (Usually "Ne2-L" or "Ne2x")

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
        libgmp-dev \
4. Install tidyverse package directly from Rstudio.
5. The abcrf package needs more dependencies so please use :
    install.packages("abcrf", dependencies = T)

### SIMULATIONS FOLDER ###

All simulations are saved to `~/Bureau/simulations/` or `~/Desktop/simulations/` depending on your system language.

Each simulation creates:
- A folder `sim_000000*/`
- A config file used by SLiM (`slim_config.txt`)
- A NeEstimator info file (`info`, `option`)
- Outputs: `.gen`, `.trees`, `summary.txt`, and optional logs
- All summaries are aggregated into: `summary_table.csv` (It's the part exploited by the R script)
- Other files are created during the simulation but are deleted at the end to gain space (see "run_simulation_linux.py" (section 14.1))

### TROUBLESHOOTING ###

Every errors I encountered would be present here (if I remember to add them)

### CONTENT ###

The main GitHub branch contains multiple folders (One for Windows and the second for Linux). Both contains mainly the same files but there are still a few changes :

1. Config folder
    Contains two files ("info" & "option"). Essential to the functionning of NeEstimator.
    Those files are used by NeEstimator to assess the different parameters of the output files created at the end of the resuming stat calculations.
    The "info" file helps the executable to run the calculation according to the methods you want to see the results. It also provides essential informations for some calculations (like the generation of the first and last sample, important for the temporal methods)
    The "option" file is used by NeEstimator to create the output files. Most of those files are created but deleted afterwards to reduce space in the folder. But they can be reactivated by checking the end of "run_simulation_linux.py"
    Please read the comments in the "run_simulation_linux.py" script at the beginning (Section #3) for more information of the content of each file.
2. Script folder
    - labels
    The less important script there. It is only used to give "pretty" names for the output files at the end.
    Instead of reading outputs like "num_loci" you will see "Num of loci". Only a Quality of Life Update.
    - main_linux
    This part is really important. It is the one you run from the console (with "python3 main_linux.py"). Once you press enter, you will have two requests appear consecutively. It demands what the name of the folder you want the simulations you are about to create will be. If you put a folder name already used, it simply will add the new simulations in the current folder without deleting the old version. The second request demands the number of simulations you want to run. Again, if there are already simulations in the folder, the process will just continue from the latest name of simulation (e.g. if the folder contains simulations up to "sim_000010", the next simulation it will create will be named "sim_000011").
    - run_simulation_linux
    This script is the core of the entire simulation workflow. It handles the generation of input files,
    the execution of SLiM, tree sequence processing, mutation simulation, formatting for NeEstimator,
    and summarizing results. It is structured into clearly defined sections:

    Section 1: Imports
        Loads all Python libraries and utility functions (e.g., labels.py) necessary for the simulation and analysis.

    Section 2: Initialization and Paths
        Sets up folders, timestamps, global config loading, and determines the path to save simulations.

    Section 3: Config File Generation
        Builds input files required for SLiM (slim_config.txt) and NeEstimator (info, option).

    Section 4: Run SLiM Simulation
        Executes SLiM using the .slim script and logs the output.

    Section 5: Tree Sequence Processing
        Loads .trees file from SLiM and filters/simplifies it to retain relevant individuals.

    Section 6: Recapitation
        Uses msprime to recapitate the tree sequence based on the specified effective size.

    Section 7: Simulate Mutations
        Simulates mutations using a stepwise mutation model (SMM) to produce allelic variability.

    Section 8: Format Data for NeEstimator
        Formats genotype data in the GENEPOP format (.gen) for use with NeEstimator.

    Section 9: Run NeEstimator
        Launches the NeEstimator binary with prepared input files and collects the output.

    Section 10: Extract NeEstimator Results
        Parses the simulation_dataNe.txt file to extract LD, HE, Coancestry, and Temporal estimates.

    Section 11: Genetic Diversity Summaries
        Parses the simulation_dataLoc.txt file to compute observed/expected heterozygosity, 
        allele counts, and allelic variances.

    Section 12: Write Summary Report
        Compiles all outputs into a human-readable summary (summary.txt) and formats per-locus tables.

    Section 13: Append to Global Summary CSV
        Adds the simulation results to the master file (summary_table.csv) used for ABC analysis.

    Section 14: Cleanup Temporary Files
        Deletes intermediate simulation files to save disk space and avoid clutter.

    Each section is marked with a header like:
        # ---___---___--- 6. Recapitation ---___---___---___--- #
    to help with navigation and readability in the script.

### LICENSES ###

This project relies on several open-source tools and libraries. Below is a description of each of them, including authorship and licensing details, to ensure transparency, traceability, and ethical compliance.

1. Python
   - Python is an interpreted, high-level programming language created by the Python Software Foundation.
   - It is widely used for scientific computing, data processing, and scripting.
   - License: Python Software Foundation License (PSF), compatible with GPL.
   - Website: https://www.python.org

2. Python Packages

   - NumPy: Developed by Travis Oliphant and the NumPy community. It provides support for numerical arrays and mathematical operations.
     License: BSD 3-Clause.

   - pandas: Created by Wes McKinney, pandas offers high-performance data structures and analysis tools.
     License: BSD 3-Clause.

   - matplotlib: A plotting library originally developed by John D. Hunter and maintained by a large community.
     License: PSF/BSD-compatible.

   - msprime: Created by Jerome Kelleher and collaborators, msprime is used for efficient simulation of genealogies and genetic data.
     License: MIT.

   - pyslim: Developed by Ben Haller and Jerome Kelleher to link SLiM simulations with tree sequence processing.
     License: MIT.

   - tskit: Developed by the tskit-dev team, it underlies both msprime and pyslim for tree sequence handling.
     License: MIT.

3. SLiM (Selection on Linked Mutations)

   - SLiM is a flexible and powerful framework for forward-time population genetic simulations.
   - It was developed by Benjamin C. Haller (Institute for Systems Biology) and Philipp W. Messer (Cornell University).
   - License: GNU General Public License v3 (GPL-3.0).
   - Website: https://messerlab.org/slim/
   - Source: https://github.com/MesserLab/SLiM

4. NeEstimator v2.1

   - NeEstimator is a software tool for estimating effective population size from genetic data.
   - Original authors include Paul A. H. Waples and Tim Do.
   - The original version is distributed under a restrictive academic license; some forks are under GPL or BSD—check the specific repository used.
   - Note: This project uses a fork compiled with cmake/make for Linux.

5. R

   - R is a language and environment for statistical computing and graphics, maintained by the R Core Team.
   - License: GNU General Public License v2 or later.
   - Website: https://www.r-project.org/

6. R Packages

   - abcrf: Developed by François Lemaire, Jean-Michel Marin, Arnaud Estoup, and collaborators, this package performs Random Forest-based Approximate Bayesian Computation.
     License: GPL-2.

   - abc: Created by François Csilléry, Olivier François, and Michael G.B. Blum to implement Approximate Bayesian Computation methods.
     License: GPL-2.

   - randomForest: By Andy Liaw and Matthew Wiener, implements the Breiman and Cutler's Random Forest algorithm.
     License: GPL-2.

   - Rcpp: Developed by Dirk Eddelbuettel and Romain François, it facilitates the integration of C++ code in R.
     License: GPL-2.

   - tidyverse: A meta-package maintained by Hadley Wickham and the RStudio team, composed of several data manipulation packages.
     License: MIT/GPL (depending on component package).

This file is intended to support the reproducibility and compliance of scientific software environments used during the M2 internship project at CBGP.

