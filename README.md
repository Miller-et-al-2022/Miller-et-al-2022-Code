# Miller-et-al-2022-Code

The code contained here accompanies Miller et al. "In vitro experiments and kinetic models of pollen hydration show that MSL8 is not a simple tension-gated osmoregulator" Current Biology, 2022.

The code for fitting the parameters includes files for fitting epsmin, kinact, and non-linear elastic deformation parameters (Pc and kp). Each set has two parts: the main script (e.g. NLE_fit.m) and an associated function (e.g. NLE_fit_funct.m). All parameter fitting codes use iterative searching to find the parameter value(s) that minimize the error when compared to the experimental data. 

The main code for running simulation is "BasicModel.m". This code defines the parameters and assumptions of the basic model for running both the WT and msl8 simulation. 

Note that the hydration experimental data is also included for comparison to the output (WT_Data.csv; msl8_Data.csv). These datasets are called for comparison when parameter fitting and the graphing at the end of the simulation code. 

Please refer to the manuscript for required changes to simulate the model variations, including parameter values and assumptions.

Contact:
Elizabeth Haswell, Washington University in St. Louis
ehaswell@wustl.edu 
