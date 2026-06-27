# Code simulating ALMA and ALMA+AtLAST observations, run a quality assessment on the images and compare them

# Folder structure and dependencies
Install Analysis Utils
Download DataComb and se t the path

# Folder structure

Folder trees:

(created by the user)

- /workingdirectory (of your choice)
  - /models : folder with ALL the models
  - /simulations : folder with all the simulated visibilities and SD images are stored 
  - /DataComb : folder with all the DataComb scripts

 
# How to run the ALPACA pipeline
  - 1.- Move to your working folder
  
  - 2.- Outside CASA, configure the DataComb file in your working directory 
	./path2DC/configure
	
  - 3.-  Adapts the script if needed:
  	set_parameters.py: if you need ot modify the names and files used in teh analysis
  	simulator.py: if you need to modify the integration time ratio using different configurations
  	sim_DCpars.py: if you need to modify the clean parameters and the data combination method used
  	sim_IQA_scripts.py: if you need to modify the QA parameters
  
  - 4.-  Start CASA
  -  Run the script
```
      CASA> exec(open("run_simulator_and_assessment.py").read(),globals())
```

