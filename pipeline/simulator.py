# This script generates the simulated visibilities to use

#Load scripts

from astropy import constants as const 
import astropy.units as u 
import numpy as np 
     
myskymodel_fits = [path2model+s for s in model_files]

os.system('rm -r '+path2simulations)
os.system('mkdir '+path2simulations)

for i in range (0,np.size(myskymodel_fits)):
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print("Model selected: " + str(model_files[i]))
	print("Running now the simulator")
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

	importfits(fitsimage=str(myskymodel_fits[i])+'.fits', imagename=str(model_files[i]), overwrite=True)
	myskymodel=str(model_files[i])
    
	#path2DCfolder=path2simulations+'/'+myskymodel+'/'
	
	global path2DCfolder
	os.system('mkdir '+path2DCfolder)

    	#Extract the parameters from the model header
    
	model_params=imhead(imagename=myskymodel)  
	ra=round(model_params['refval'][0]*180/3.14159265359)
	dec=round(model_params['refval'][1]*180/3.14159265359)
    
	freq=round(model_params['refval'][2]/1e9)
    
	myphasecenter='J2000 '+str(ra)+'deg '+str(dec)+'deg'
    	#print(myphasecenter)
	model_size=model_params['shape']
	incr_size=model_params['incr'] 
	pix_size=incr_size[1]*(3600*180)/3.14159265359	#pixel size from rad to arcsec
	model_size_arcmin=pix_size*model_size[0]/60		#total size of the model in arcmin
	
	tmpdir=['ALMA','C43-'+alma12m_array_configurations[1]+'+C43-'+alma12m_array_configurations[0]+'+AtLAST','C43-'+alma12m_array_configurations[0]+'+AtLAST']
	path_2_tmp_sim=[path2DCfolder+t for t in tmpdir]
	for t in range (np.size(path_2_tmp_sim)):
		os.system('mkdir '+path_2_tmp_sim[t])	
	
	if alma12m_array_configurations==['']:
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
		print('ALMA 12m configurations not provided. Skipping ALMA 12 simulations.')
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    
    		
	else:
		alma12m_configurations=['alma.cycle'+alma_cycle+'.'+c+'.cfg' for c in alma12m_array_configurations]		#configuration of the 12m array
    	
		if alma12m_array_configurations[0]=='10':
    			print('---No combination aviable with C43-10-----')
		if alma12m_array_configurations[0]=='9':
    			print('The selcted extended configuration of the ALMA 12m array is C43-'+str(alma12m_array_configurations[0])+' in Cycle '+str(alma_cycle))
    			print('Integration time per pointing ratios are C43-'+str(alma12m_array_configurations[0])+': C43-'+str(alma12m_array_configurations[1])+' = 1 : 0.21')
		if alma12m_array_configurations[0]=='8':
    			print('The selcted extended configuration of the ALMA 12m array is C43-'+str(alma12m_array_configurations[0])+' in Cycle '+str(alma_cycle))
    			print('Integration time per pointing ratios are C43-'+str(alma12m_array_configurations[0])+': C43-'+str(alma12m_array_configurations[1])+' = 1 : 0.22')
		if alma12m_array_configurations[0]=='7':
    			print('The selcted extended configuration of the ALMA 12m array is C43-'+str(alma12m_array_configurations[0])+' in Cycle '+str(alma_cycle))
    			print('Integration time per pointing ratios are C43-'+str(alma12m_array_configurations[0])+': C43-'+str(alma12m_array_configurations[1])+' = 1 : 0.23')
		if alma12m_array_configurations[0]=='6':
    			print('The selcted extended configuration of the ALMA 12m array is C43-'+str(alma12m_array_configurations[0])+' in Cycle '+str(alma_cycle))
    			print('Integration time per pointing ratios are C43-'+str(alma12m_array_configurations[0])+': C43-'+str(alma12m_array_configurations[1])+' : ACA-7m = 1 : 0.25 : 0.6')
		if alma12m_array_configurations[0]=='5':
    			print('The selcted extended configuration of the ALMA 12m array is C43-'+str(alma12m_array_configurations[0])+' in Cycle '+str(alma_cycle))
    			print('Integration time per pointing ratios are C43-'+str(alma12m_array_configurations[0])+': C43-'+str(alma12m_array_configurations[1])+' : ACA-7m = 1 : 0.26 : 1.21')
		if alma12m_array_configurations[0]=='4':
    			print('The selcted extended configuration of the ALMA 12m array is C43-'+str(alma12m_array_configurations[0])+' in Cycle '+str(alma_cycle))
    			print('Integration time per pointing ratios are C43-'+str(alma12m_array_configurations[0])+': C43-'+str(alma12m_array_configurations[1])+' : ACA-7m = 1 : 0.34 : 2.4')
		if alma12m_array_configurations[0]=='3':
    			print('The selcted extended configuration of the ALMA 12m array is C43-'+str(alma12m_array_configurations[0])+' in Cycle '+str(alma_cycle))
    			print('Integration time per pointing ratios are C43-'+str(alma12m_array_configurations[0])+' : ACA-7m = 1 : 2.4')
		if alma12m_array_configurations[0]=='2':
    			print('The selcted extended configuration of the ALMA 12m array is C43-'+str(alma12m_array_configurations[0])+' in Cycle '+str(alma_cycle))
    			print('Integration time per pointing ratios are C43-'+str(alma12m_array_configurations[0])+' : ACA-7m = 1 : 4.7')
		if alma12m_array_configurations[0]=='1':
    			print('The selcted extended configuration of the ALMA 12m array is C43-'+str(alma12m_array_configurations[0])+' in Cycle '+str(alma_cycle))
    			print('Integration time per pointing ratios are C43-'+str(alma12m_array_configurations[0])+' : ACA-7m = 1 : 7')
    
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
		print('These simulator targets ALMA 12m C43-4 + C-43-1 + ACA-7m.')
		print("The visibilities are simulated following the time ratio reported in the ALMA THB,\nC43-4 : C-43-1 : ACA-7m = 1 : 0.34 : 2.4.")
		print('The C43-4 12m configuration integration time per pointing is 30 sec observed twice in 2 consecutive days,\nstarting at different HA. Total time per pointing: 2min.')
		print('The C43-1 12m configuration integration time per pointing is 20 sec observed twice in 1 day.\nTotal time per pointing: 40sec.')
		print('The ACA-7m integration time per pointing is 30 sec observed twice in 5 consecutive days starting at different HA.\nTotal time per pointing: 5min.')
		print('If you consider differnt configurations, follow the time ratio you selected printed above and modified the file\n simulator.py accordingly.')
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    
    # the extended and compact configurations of the 12m and the 7m used in the simulations: C43-4 as the extended configuration and the C43-1 as the compact one
    	# the time ratio calculated according to the ALMA THB is C43-4 : C43-1 : 7m : TP = 1 : 0.34 : 2.4 : 4.0
    
    #Simulating the extended 12m array visibilities (NOISE-FREE)
    
		for c in range (np.size(alma12m_array_configurations)):
			if c==0:
			#print(c)
				for d in range (0,2):
				#print(d)
    	#We calculated the observing time for the extended configuration of the 12m array (C43-4) as 2min per pointing: 30 sec of integration time per pointing observed twice per 2 days for 
    	#a total of 2min per pointing. The number of pointing is 67 (mapsize=3arcmin). Total time for the C43-4 conf. is 134 min = 2.23h 
    
					myproject_12m_ext=str(myskymodel)+'_12m_C43-'+alma12m_array_configurations[c]+'_day'+str(d)
    
					simobserve(project=myproject_12m_ext, skymodel=myskymodel, direction = myphasecenter, obsmode = "int", antennalist=alma12m_configurations[c],
				 mapsize=mapsize_12m, thermalnoise = '',setpointings=True, maptype='hex',refdate='2018/10/0'+str(d+1),hourangle=str(d)+'h',pointingspacing='nyquist',
				 integration='30s', totaltime='2')
				 
					os.system('cp -r '+myproject_12m_ext+'/*.ms '+path_2_tmp_sim[0])	
					os.system('cp -r '+myproject_12m_ext+'/*.ms '+path_2_tmp_sim[1])
					os.system('cp -r '+myproject_12m_ext+'/*.ms '+path_2_tmp_sim[2])	
				 
    
 	#We calculated the observing time for the compact configuration of the 12m array (C43-1) as 0.34x extended time=40s: 20 sec of integration time per pointing observed twice per 1 days 
 	#for a total of 0.68min per pointing. The number of pointing is 67 (mapsize=3arcmin). Total time for the C43-1 conf. is 46 min
		
			if c==1:
				#print(c)
				myproject_12m_comp=str(myskymodel)+'_12m_C43-'+alma12m_array_configurations[c]
     
				simobserve(project=myproject_12m_comp, skymodel=myskymodel, direction = myphasecenter, obsmode = "int", antennalist = alma12m_configurations[c], mapsize=mapsize_12m,
			 thermalnoise = '',setpointings=True, maptype='hex',refdate='2018/10/02',hourangle='0h',pointingspacing='nyquist', integration='20s', totaltime='2')
			 
				os.system('cp -r '+myproject_12m_comp+'/*.ms '+path_2_tmp_sim[0])	
				os.system('cp -r '+myproject_12m_comp+'/*.ms '+path_2_tmp_sim[1])


	if aca7m==True:
		aca7m_configurations='aca.cycle'+alma_cycle+'.cfg'		#configuration of the 7m array
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
		print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
		print('ACA-7m simulations included.')
 
		#Simulating the compact 7m array visibilities (NOISE-FREE)
		for d in range (0,5):
			myproject_7m=str(myskymodel)+'_7m_day'+str(d)

    			#We calculated the observing time for the 7m array as 2.4x extended time=4.8min: 30 sec of integration time per pointing observed twice per 5 days for a total of 
    			#5min per pointing. The number of pointing is 27. Total time for the compact conf. is 3.8h

			simobserve(project=myproject_7m, skymodel=myskymodel, direction = myphasecenter, obsmode = "int", antennalist =aca7m_configurations, mapsize=mapsize_7m, 
    			thermalnoise = '',setpointings=True, maptype='hex',refdate='2018/10/0'+str(d+1),hourangle=str(d-2)+'h',pointingspacing='nyquist', integration='30s',
    			totaltime='2',graphics='both')
    			
			os.system('cp -r '+myskymodel+'_7m_day'+str(d)+'/*.ms '+path_2_tmp_sim[0])	
    
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	for s in range (0,np.size(SD_name)):
    	
		res_radians=1.22*((const.c/(freq*u.GHz)).to(u.m))/(int(SD_diameter[s])*u.m)
		res=res_radians.value*(3600*180)/3.14159265359
		print('The '+SD_name[s]+' has a diameter of '+str(SD_diameter[s])+'m. Its final resolution at the model frequency '+str(freq)+'GHz is '+str(res)+' arcsec.')
		print('SImulating the '+SD_name[s]+'image by convolving the model to the SD resolution calculated above.')
		print('--------------------------------------------------------------------------------------------------------------------------------------')
    
		myproject_sd=str(myskymodel)+'_'+SD_name[s]
		os.system('rm -r '+ myproject_sd)
		os.system('mkdir '+ myproject_sd)

		imsmooth(imagename= myskymodel,outfile= myproject_sd,kernel='gauss',major=str(res)+'arcsec',minor=str(res)+'arcsec', pa='0deg',targetres=True,overwrite=True)
		hdval = 'Jy/beam'
		imhead(imagename=myproject_sd,mode='put',hdkey='BUNIT',hdvalue=hdval)
		
		if s==0:
			os.system('cp -r '+myproject_sd+' '+path_2_tmp_sim[0])	
		else:
			os.system('cp -r '+myproject_sd+' '+path_2_tmp_sim[1])	
			os.system('cp -r '+myproject_sd+' '+path_2_tmp_sim[2])		

	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('All the simulations produced and located in the simulations directory.')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++') 
	
	os.system('rm -r '+myskymodel+'*')	
		
