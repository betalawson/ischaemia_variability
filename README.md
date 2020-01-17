## Variability for Re-entry Induction in Ischaemic and Heterogeneous Tissue

This repository stores the code used for Figure generation, as well as the tools used for  
determining steady state conditions in 1D Fibre, processing of the simulation data,  
emulation, main effect calculation, and so on.

All data is included, except for the raw simulation data, due to filesize reasons.  
Full data for the visualised simulations ( #597, #2613 ) are included  
Data is stored in the following MATLAB data files:  

**data_without_sims.mat:**               Metric data in 1D and 2D, as well as classifications of simulations.  
                                     1D data is repeated because multiple 2D sims were run for each set of parameters  
**data1Dmat.mat:**                       Only the 1D data (without repeats)  
**main_effect_data.mat:**                Calculated main effects. To understand format, check FIGURE_MainEffectPlots  
**main_effect_data_shortrange.mat:**     Calculated main effects over a shorter range.  

#####Additional raw data is included in the .txt files:

**param_SS_1000configs.txt:**         Steady state data obtained for the 1000 parameter configurations  
**bifurcation1D-table.txt:**          Information from bifurcating fibre simulations

#####Metrics stored are:

**1D:** [ APD, CV, WL, V_rest, V_peak, APA, block_level, delay1, delay2, delay3, delay4, block_susceptibility ]  
(note that delays are post subtraction of delay1)

**2D:** [ wave_disorder, blockRatio, blockCount, APDmean, APDstd, WLmean, WLstd, DF, F_std, multiActRatio, multiActCount ]

#####Classifications ("flags") are:

**-1:**       Non-propagating  
** 0:**       Propagates, no type of re-entry  
** 1:**       At least 0.1% of sites were activated more than once  
** 2:**       At least 2.5% of sites were activated more than once, but re-entry doesn't persist  
** 3:**       As per 2, but re-entry persisted for the full t = 1000ms

----------------------------------------------------------------

#####Key functions are:

**createFibreProblem:**      Creates the 1D Fibre problem in a format used by the simulator  
**findSteadyState:**         Finds steady state configuration (and checks if excitation propagates at all) for a given problem (including parameters)  
**createRunInfo:**           Selects parameters from the specified parameter space, and finds the associated steady state values (used as initial condition)

**calculateMetrics:**        Calculates the metrics for 2D simulation data

**generate1DEmulatorData:**  Takes the 1D data and builds partitioned emulators for the metrics of interest

**FIGURE_[fig_desc]:**       Files that read in the saved data and create the figures