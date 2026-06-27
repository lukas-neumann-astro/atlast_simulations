#Running a QA where the reference image is the model convolved to the final resolution o the simulated images

execfile('IQA_script.py')

import os

#cleaning_threshold = 0.163 
beam_final = "1.5"  # in arcsec

datacombs_labels = ["ALMA12m-C43-4+C43-1+7m_TCLEAN","ALMA12m-C43-4+C43-1+7m+TP_FEATHER","ALMA12m-C43-4+C43-1+7m+TP_MACF","ALMA12m-C43-4+C43-1_TCLEAN","ALMA12m-C43-4+C43-1+AtLAST_FEATHER","ALMA12m-C43-4+C43-1+AtLAST_MACF","ALMA12m-C43-4_TCLEAN","ALMA12m-C43-4+AtLAST_FEATHER","ALMA12m-C43-4+AtLAST_MACF"]

clean_labels = ["ALMA12m-C43-4+C43-1+7m_TCLEAN","ALMA12m-C43-4+C43-1_TCLEAN","ALMA12m-C43-4_TCLEAN"]

feather_labels = ["ALMA12m-C43-4+C43-1+7m+TP_FEATHER","ALMA12m-C43-4+C43-1+AtLAST_FEATHER","ALMA12m-C43-4+AtLAST_FEATHER"]

macf_labels = ["ALMA12m-C43-4+C43-1+7m+TP_MACF","ALMA12m-C43-4+C43-1+AtLAST_MACF","ALMA12m-C43-4+AtLAST_MACF"]

images=['tclean','feather_f1.0','hybrid_f']
image_names=[path2DCfolder+ d+'/datacomb/INTimage.mfs_nt1_'+mydeconvolver+'_SD-INT-AM_nIA_n'+str(Niters)+'.' for d in sim_folder]
#AtLAST_ext_comp='../AtLAST+ext+comp/AtLAST+ext+comp.mfs_nt1_HB_AM_nIA_n100000.'
#AtLAST_ext='../AtLAST+ext/AtLAST+ext.mfs_nt1_HB_AM_nIA_n100000.'
target_image=[]

macf_image=[]
feather_image=[]
clean_image=[]

for i in range (0,np.size(images)):
	clean_image_pb=image_names[i]+images[0]+'.image.pbcor'
	feather_image_pb=image_names[i]+images[1]+'.image.pbcor'
	macf_image_pb=image_names[i]+images[2]+'.image.pbcor'
	macf_image.append(macf_image_pb)
	feather_image.append(feather_image_pb)
	clean_image.append(clean_image_pb)
	for j in range (0,np.size(images)):
		pb_image=image_names[i]+images[j]+'.image.pbcor'
		target_image.append(pb_image)

datacombs_image = [s for s in target_image]

print(datacombs_image)
print('------')
print(clean_image)
print('------')
print(feather_image)
print('------')
print(macf_image)
print('------')

# Plus one PB image (e.g. TCLEAN)
pb_image = image_names[0]+"tclean.pb"

path2model='./models/'				#path to the model folder
model_files=['skymodel-b']
model=model_files[0]
skymodel=path2model+model+'.fits'

########################################################################################
## Comparison with Skymodel
if (beam_final != "-1"):
	# Import skymodel into CASA + convolve it into a circular beam
	importfits(fitsimage=skymodel, imagename=skymodel+".image", overwrite=True)
	os.system("rm -rf *_conv"+str(beam_final))
	os.system("rm -rf *_conv"+str(beam_final)+".fits")
	get_convo(convo_file=skymodel+".image",beam_final=str(beam_final))
	#imregrid(imagename=skymodel+".image_conv"+str(beam_final),template= target_image[0],output=skymodel+".image_conv"+str(beam_final)+'_imregrid',overwrite=True)
	os.system('rm -r '+ skymodel+".image_conv"+str(beam_final)+'_imtrans')
	imtrans(imagename=skymodel+".image_conv"+str(beam_final),outfile=skymodel+".image_conv"+str(beam_final)+'_imtrans',order='0132')
	skymodel_image = skymodel+".image_conv"+str(beam_final)+'_imtrans'
	exportfits(imagename=skymodel_image,fitsimage=skymodel_image+".fits",overwrite=True)

	##
	# IQAs
	get_IQA(ref_image=skymodel_image,target_image=datacombs_image, masking_RMS=0.35, target_beam_index=0)

	#Get plots 
	for i in datacombs_image:
		show_Apar_map(ref_image=skymodel_image,#+"_masked",
			target_image=i,
			channel=0,
			save=True,
			titlename=i,
			plotname=i+"_Apar")
		#show_Fidelity_map(ref_image=skymodel_image,#+"_masked",
			#target_image=i,
			#channel=0,
			#save=True,
			#plotname=i+"_Fidelity")

	for i in datacombs_image:
		Compare_Apar_signal(ref_image=skymodel_image,#+"_masked",
			target_image=[i],
			save=True,
			plotname=i+"_Apar_vs_signal")

	Compare_Apar(ref_image=skymodel_image,#+"_masked",
			target_image=datacombs_image,
			save=True,
			plotname="ALL_Comparison_Apar",
			labelname=datacombs_labels)
			
	Compare_Apar(ref_image=skymodel_image,#+"_masked",
			target_image=clean_image,
			save=True,
			plotname="Clean_Comparison_Apar",
			labelname=clean_labels)
	Compare_Apar(ref_image=skymodel_image,#+"_masked",
			target_image=feather_image,
			save=True,
			plotname="Feather_Comparison_Apar",
			labelname=feather_labels)
	Compare_Apar(ref_image=skymodel_image,#+"_masked",
			target_image=macf_image,
			save=True,
			plotname="MACF_Comparison_Apar",
			labelname=macf_labels)

	#Compare_Fidelity(ref_image=skymodel_image,#+"_masked",
			#target_image=datacombs_image,
			#save=True,
			#plotname="ALL_Comparison_Fidelity",
			#labelname=datacombs_labels)


	genmultisps(fitsimages=[s + "_convo2ref.fits" for s in datacombs_image], 
			save=True, 
			plotname='ALL_PS', 
			labelname=datacombs_labels,
		        titlename='')
		        
	exportfits(imagename=i+'_convo2ref_Apar',fitsimage=i+'_convo2ref_Apar'+".fits",overwrite=True)	        
	
	get_IQA_nomask(ref_image=skymodel_image,target_image=i, masking_RMS=0, target_beam_index=0)

