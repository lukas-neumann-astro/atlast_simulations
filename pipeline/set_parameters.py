#====================================================================================
# set_parameters.py: Template of the parameters file
#====================================================================================
# All user-defined parameters are included (and can be modified) in this file

path2model='./models/'				#path to the model folder
model_files=['skymodel-b']			#the model to simulate

path2simulations='./simulations'			#path to the simulations folder
path2DCfolder=path2simulations+'/'+model_files[0]+'/'
sim_folder=['ALMA','C43-1+C43-4+AtLAST','C43-4+AtLAST']

alma_cycle='7'					#ALMA cycle to consider		
		
alma12m_array_configurations=['4','1']		#ALMA 12m array configuration to considered, ordered from the extended to the compact one the possibilities are 10, 9+6, 8+5, 7+4, 6+3+7m, 5+2+7m,
						# 4+1+7m, 3+7m, 2+7m, 1+7m
aca7m= True					#ACA-7m array to be considered or not

mapsize_12m='3arcmin'				#Area to observe with the 12m array: default ('') is the model image size
mapsize_7m='3.5arcmin'				#Area to observe with the 7m array: default ('') is the model image size. It should be larger than the area observed with the 12m array

SD_name=['TP','AtLAST']				#SD name
SD_diameter=['12','50']				#SD diameter in meters used to calculate the SD resolution

# the extended 12m configuration integration time per pointing is 30 sec observed twice in 2 consecutive days starting at different HA. 


#Modify the parameters of datacomb script for the data imaging and combination (see the parameters from the DC_pars.py)

thesteps=[1,2,3,5]#,2,3,5]  			#steps of the combination to execute

t_imsize      = [800,800] 
t_cell        = '0.21arcsec'   
Niters=100000					#Number of cleaning iterations

mydeconvolver='HB'

