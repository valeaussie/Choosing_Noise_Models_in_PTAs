# Choosing_Noise_Models_in_PTAs

To run the simulations use the python file "noise_injections.py" in the folder "sims_software".
You can run the simulations using Libstempo (https://github.com/vallis/libstempo/blob/master/demo/libstempo-demo.ipynb) only, or first create a dataset with PTASimulate (https://bitbucket.org/psrsoft/ptasimulate/src/master/) then inject the different noise sources with Libstempo.
Use the file "WN_PTASim.imp" (modified as needed) as the imput file for PTASimulate.

The analysis is done with Enterprise https://github.com/nanograv/enterprise. 
Different models are found in the folder "posteriors".
For example the file WN_TN_CRN.py uses a model containing white noise (WN), timing noise (TN) and a common red noise (CRN).
You can create your own file with the combinations of noise in the model that you need for your analysis.
It is important to use the same ephemeris for simulations and parameter estimations. PTASimulate uses DE421.

Once you have the outputs from the Enterprise analysis, you will need to run "makeCommonCorner.py" to burn and analyse the chains. To make the corner plots use "corner_chainplotter.py". This will make plots like in figure 1, 2, 3 and 4.

To make the mahalanobis distances first run "mahalanobis.py" then "plot_distribution.py" for the plot in the appendix, or "plot_mahalanobis.py" for the plot in Figure 5 of the paper.

DATA

The folder "sims_sanity_checks" contains the simulations for the Figure in the appendix as well as the chains.
The misspecified model compares chains 3 and 4 and the correctly specified model compares chains 4 and 5. They are in the folder "sims_sanity_check/sim_0/chains/commonNoise_pl_nocorr_freegam_DE421" for each simulation (in this case sim_0)

The folder "" contains the simulations for Figure 5 of the paper

The folder "" contains the simulations for Figure 1, 2, 3 and 4 of the papers

