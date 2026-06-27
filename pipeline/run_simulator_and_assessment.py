#Before opening casa configure datacomb with the command ./DataComb/configure in the working directory

#====================================================================================
# run_simulator_and_assessment: ALMA + AtLAST simulator script
#====================================================================================

execfile("sim_scripts.py")

execfile("set_parameters.py", globals())          # user modifiable parameters (templates exist)

execfile("simulator.py",  globals())          	  # call the datacomb routines

print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
print("Starting the imaging procedure and the data combination")
print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

execfile("DC_locals.py",globals())

for i in range (np.size(sim_folder)):
	d=i
	s=i
	if d==2:
		s=1
		
	print(d,s)
	execfile('sim_DCpars.py', globals())          # user modifiable parameters (templates exist)
	execfile("DataComb/DC_run.py",  globals())         


execfile('IQA_script.py')

execfile('sim_IQA_script.py')
