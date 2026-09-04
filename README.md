# Project Title

pNMS: Python package for virtual experiments on the neuromuscular system at multiple scales

## Description

The pNMS was developed to enable the automatic construction and parallel simulation of biophysical models for heterogeneous neuromuscular cells at the levels of the motoneuron, muscle-tendon fiber, and motor unit under various physiological conditions. The software package provides a flexible, efficient application programming interface (API) with functions that can be used in the terminal or integrated into application scripts.
```
Project/
├────── pNMS/
│       ├────── __init__.py
│       ├────── application_programming.py
│       ├────── simulation_engine.py
│       ├────── ppserver.sh
│       └────── license.txt
├────── parameters/
│       ├────── Iaxon/     
│       ├────── Isoma/
│       ├────── Isyn/
│       ├────── MF_Parameters/
│       ├────── MN_Parameters/
│       ├────── RN Dpath/
│       └────── Xm/
├────── examples/
│       ├────── article.pdf     
│       └────── scripts.zip
├────── pNMS_API.pdf
├────── README.md
└────── LICENSE.md      
```
## Getting Started

### Dependencies

Python (3.12.7), Pandas (2.2.2), Scipy (1.13.1), Matplotlib (3.9.2), PyQt5 (12.13.0), pp (1.6.4.4)

### Installing

The pNMS can be installed by cloning it from GitHub and moving it to the folder site-packages in the Python library folder under the folder where the Python interpreter is located. The folder parameters, including user-defined parameter files, should be in the same folder as the application program.

### Using API functions

[1] How to use the pNMS API 
#### import API library from the pNMS folder
```
import pNMS.application_programming as cf
```
#### use API functions 
```
cf.function_name(…)
```

[2] pNMS API functions
#### Creation of a homogeneous population model 
```
createPool
```
#### Transition to a heterogeneous population model 
```
setParameters 
plot_params
```
#### Setup of simulation conditions 
```
setInitialValues 
setSimulTimes
```
#### Setup of input conditions 
```
genNeuronInputSignals 
genSynConSignals 
genSpikeSignals 
genMuscleLengthSignals 
importNeuronInputSignals 
importSynConSignals 
importSpikeSignals 
importMuscleLengthSignals 
setNeuronInputSignals 
setSynConSignals 
setSpikeSignals 
setMuscleLengthSignals 
plotNeuronInputSignal 
plotSynConSignal 
plotSpikeSignal 
plotMuscleLengthSignal
```
#### Setup of parallel computing environment 
```
setComputeNode
```
#### Execution of parallel simulation 
```
runSimulation
```
#### Display and saving of simulation results 
```
plotSimulResult 
saveSimulationResults 
plotImportData
```
Please see the user manual for the pNMS API (pNMS_API.pdf) for details.

### Executing program

#### Multicore computer environment 
 
The number of cores to use is specified by an integer or ‘autodetect’ to use all available cores on the local computer. An example code is as follows:
```
node_list = []
cf.setComputeNode(node_list) 
cf.runSimulation('autodetect', None)
```
The application code can be run in any Python development environment or command window.

#### High-performance computing environment 
 
The venv module is first installed in the shared folder on the management node. Within the venv folder, the relevant Python version and libraries are installed to run simulations using the pNMS software. The pNMS folder is then moved to the Python library folder (e.g., /venv/lib/python2.6/site-packages) within the venv folder. The parallel environment is specified by the number of cores on the management node and the names of the participating computational nodes (e.g., mupool-c01 and mupool-c02). An example code is as follows: 
```
node_list = ['mupool-c01', 'mupool-c02'] 
cf.setComputeNode(node_list)
cf.runSimulation(0, ‘ALL’) # with no core in the management node and all computational nodes
```
To run the application code on the management node, the venv module needs to be activated first, as shown in the following example.
```
[root@mupool-m 2024]# source venv/bin/activate
```

### Usage examples

Please refer to the folder example for the previously published application scripts.

## Help

Please get in touch with the corresponding author (Hojeong Kim) for any advice or help with problems or issues on the pNMS.

## Authors

Contributors: Jungwhan Gwak, Yeongjae Kim, Minjung Kim, and Hojeong Kim
Corresponding author: Hojeong Kim
Contact: hojeong.kim03@gmail.com

## Version History

Version 3.0
  - pNMS for Python 2.7
  
Version 4.0
  - pNMS for Python 3.12

## License

This project is licensed under the GNU version 3 License - see the LICENSE.md file for details

## Citation

When you use the pNMS in your research, please cite the following article: 
Kim H (2026). Automating population construction and parallel simulation of biophysical models for neuromuscular cells: An inverse approach. PLoS Comput Biol 22(4): e1014184. https://doi.org/10.1371/journal.pcbi.1014184

## Acknowledgments

The authors thank Hyejin Choi and Sohyun Yang for their assistance in the early versions of pNMS. This project was funded by the Ministry of Science and ICT (DGIST R&D Program) and the National Research Foundation of Korea (2021R1F1A1062265).
