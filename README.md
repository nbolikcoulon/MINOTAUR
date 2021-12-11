# MINOTAUR

1. Requirements
MINOTAUR has to be run using Python 3. The following libraries need to be installed: corner, emcee, imageio, numpy, matplotlib, multiprocessing, PIL, scipy, tkinter. MINOTAUR is called using the following command in a terminal: python run.py. Alternatively, if the C-script are already compiled, you can run directly MINOTAUR by  writing  in a terminal: python MINOTAUR.py Path_to_folder, with Path_to_folder the path to the folder containing the C-scripts.

2. Preparing for a run
The Demo folder gives an example of C-files containing expressions of relaxation rates and a relaxation matrix for a 13C1H2H2 methyl group with a vicinal deuterium, as well input files. The spectral density function is written using the extended model-free formalism and the carbon-13 CSA is among the fitted parameters. We will refer to the given files in the following.

2.1. C-scripts
In this section, we will detail how to implement the expressions for the relaxation rates and relaxation matrix and make it readable by MINOTAUR. When running the run.py script, the user will be prompted to select for a folder containing the C-scripts. This folder should contain the following files:
- setup.py
- Parameters.py
- Rates.h
- Rates.c
- _Rates.c
- itemize

setup.py: This script will be used to compile the C-functions. The part related to the relaxation matrix (names _RelaxMat) should not be changed, but all the others are related to the relaxation rates that can be added to the intensity decay and can (have to) be changed  only for the name of the function which must follow the synthax: _RateNamecalculation, where RateName is a user-defined name for the considered relaxation rate. 

Parameters.py: This script is largely user-based and care has to be taken when editing it. The table RelaxationRates contains the list of RateName. Make sure to use the same as in the setup.py script, and throughout all the subsequent scripts. The variable name PositionAuto is the position of the operator of interest in the basis where the relaxation matrix is written. The dictionary Names contains the variable names that will be determined during the analysis. Make sure to keep the same keys (OrderParam, CorrTimes, others). If one key has no associated variable, simply leave an empty table: []. The variables names do not have to match the ones used in the subsequent scripts, but will be used for labelling purposes. The ImportFunc function contains the list of the C-functions that are defined (see below). The same nomenclature as in the setup.py script has to be used. The Cons function sets constrains on the variables for the MCMC script. The constrains can be edited as wished by the user, but we recommend using the same bounds for the parameter lnf (see the emcee documentations for further informations).

Rates.h: This file introduces the functions that will be used. It starts with the relaxation matrix which, in C language is defined as a struct with dimensions which has to be set accordingly (in our example, the relaxation matrix is a 3x3 matrix). Then, the functions defined subsequently have to be written following the same synthax as used in the given example.

Rates.c: This script is user-based only. It does not necessarily contains definitions of the spectral density functions, as these will not be used in MINOTAUR. It must however contain the definition of all the relaxation rates defined in the setup.py, parameters and Rates.h scripts.

_Rates.c: This scripts is meant to transform the python inputs into C-readable inputs for the functions, and transform the C-outputs into python readable arrays. For each of the functions defined in the Rates.h script, the following parts must be present in this script:
- a docstring: an explanation of what the associated functions does
- the line: static PyObject *RateNamecalculation_RateNamecalculation(PyObject *self, PyObject *args);
- an additional element in module_methods written as: {"RateNamecalculation", RateNamecalculation _RateNamecalculation, METH_VARARGS, RateNamecalculation_docstring}
- a PyMODINIT\_FUNC module. Make sure to change its name as in the struct part.
- a static PyObject with names matching the ones from the previous scripts. A simple copy-paste from the already given code works (except for the relaxation matrix), but changing the name of the function being called is still necessary.

In all the above items, the ones associated to the relaxation matrix do not have to be changed. 

After choosing the folder containing the C-scripts, the user can start MINOTAUR. This will automatically compile the scripts. Then a second GUI opens to allow the user to load the data sets.


2.2. Inputs for MINOTAUR
The GUI is devided into different parts:
- the fitting parameters part: the user has to set the global tumbling correlation time, the type of shuttling, the number of steps in the MCMC and the number of chains (see emcee  documentations for more information). The default value for the number of chains is twice the number of free parameters, but this can be increased. In addition, the first 500 steps of the MCMC will be discarded when reporting the median and standard deviations or the first 20% if less than 500 steps are calculated. The user can always calculate the distribution afterwards as the whole trajectories are saved.
- the limits for the initial parameters guess. This part will only be used to generate the initial random parameters. Make sure the given ranges also comply with the Cons function from the Parameters.py script.
- PDB: a 4-letter PBD ID can be provided. It will generate scripts to color the structure in PyMOL according to parameters value if more than one residue is analyzed.
- files: this is the part where all the file are loaded. We review below how to format them.

A number of files need to be loaded, and follow the logic explained here. Each columns in each file must be separated by a tab:

- field calibration: two-column file, the first one with the height from the center of the magnet, the second the associated magnetic field. The height can be given in any units as long as it corresponds to the same unit as provided below for the experimental setup.
- experimental setup: contains one line for each relaxometry experiment. All delays are given in milliseconds. Column 1 is the experiment number. Column 2 is the shuttling height. Column 3 is the delay after shuttling to low field. Column 4 is the delay after shuttling to high field. Column 5 is the waiting time before shuttling to low field. Column 6 is the waiting time after shuttling to low field. Column 7 is the delay before shuttling back to high field. Column 8 is the delay after shuttling back to high field. All the subsequent columns are the relaxation delay times.
- the relaxometry intensity folder: it contains one file for each experiments mentioned in the experimental setup file. Each of these files contains one line for each considered amino acid. The lines contain first the amino acid number, and then the intensity associated to the first relaxation delay, the corresponding experimental error, then the  intensity for the second relaxation delay and the error, etc... If some data are missing, replace the values by NA.
- other input file: the C-function are writen so that they can accept as input the magnetic field, the global tumbling correlation time, the different parameters of the model, and other amino acid-dependent variables contained in this file. It contains N+1 columns for N other parameters, the first column being the amino acid number. In  the  example provided here, the other inputs are the position of the vicinal deuterium, which in effect are set to infinity.
- high field rate files. There should be one file for each type of rate (which has to be specified) and each magnetic field (which also has to be specified, in Tesla). The file contains 3 columns: the amino acid number, the relaxation rate and its experimental error. Upon  loading a file, a 'Add high field rates' button will appear, allowing the user to add all its high-fields datasets. If some data are missing (example of the R2 for residue 44), replace the values by NA.

Loading all the files can be time-consuming. A 'save current parameters' allows the user to generate a file that is readbale in the GUI by clicking the 'load previous parameters' button and which will automatically fill the GUI as previously saved. Such file is also generated by MINOTAUR after starting a run. This file can also be generated manually.
