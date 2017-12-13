# BLCMforMI
R code for the simulation studies carried out in the paper (Vidotto, Vermunt, Van Deun) - "Bayesian Latent Class Models for The Multiple Imputation of Categorical Data"



This repository includes files, datasets, results, and R packages used for the paper 'Bayesian Latent Class 

models for the Multiple Imputation of categorical data' (Vidotto, Vermunt, Van Deun). 

In order to run the code for the frequentist Latent Class imputation models, LatentGOLD (version 5.1) needs 

to be installed in the machine.

In the simulation studies, packages such as "doParallel" and "foreach" were used in order to run the 

simulations in parallel.

Before performing the Experiments, please install the R packages DPMMimpute and MMLCimpute (files *.tar.gz 

in the main folder). Steps for installation are:


-Make sure the Rgui file is in the Environment Variables of the system
-In Windows, open Command Prompt and go to the folder where the .gz files are located
-type Rcmd build --binary NAMEPACKAGE
-type Rcmd INSTALL NAMEPACKAGE_1.0.tar.gz


where NAMEPACKAGE must be replaced by DPMMimpute or by MMLCimpute.


IMPORTANT: the package MMLCimpute is an unofficial (and incomplete) version of an R package that will be 

released soon by the authors. Thus, it should only be used with the simulation studies of this repository. 

Furthermore, the package DPMMimpute was created for the only purpose of evaluating the results of the 

simulation study, and should not be evaluated in terms of efficiency. 


In the folder, you will find both the generated data (.RData) used for the simulations and the R code (.R) 

used to generate it. If you want to simulate data from scratch, we suggest to save the .RData file at the 

end of each step, in such a way to be able to perform the simulations in the next steps with the code 

generated in previous steps. 
