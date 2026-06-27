##-------------------------------------------------------
"""
EMERGE PIPELINE TOOLS
ALMA data reduction pipeline for the EMERGE project

More info: https://emerge.univie.ac.at/
==============================
Authors: A.Hacar (Univ. Vienna), D. Petry (ESO), F.Bonanomi (Univ. Vienna), A.Socci (Univ. Vienna)

"""

# Comments:
# - In general:
#	datadir : folder where data (.MS or FITS) are stored
#	workdir : folder where the products (after applying a function) will be stored/copied


###############################################
## Libraries
###############################################

# Dependencies: analysisUtils

# General
import os
import sys
import numpy as np
from importlib import reload 
import re

from casashell import *
from casatasks import *
from casatools import image as iatools
import casatools as ctools
ia = iatools()

# CASA tasks
import casatasks as ct

# Astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
import astropy.units as u
import astropy.constants as c

# DC scripts
import IQA_script as iqa



## To-do
# Improve GOUS selection

# Check consistencia before concat. Transform if necessary

# Implement imaging of self-cal.TBD

# Add option to create mosaics in GALACTIC
# imregrid(imagename="input.image", output="output.image", template="B1950")



#######################################################
# Import EMERGE line catalogs
#######################################################

from emerge_catalogs import *


#######################################################
# Section 0: Environment + general tools
#######################################################

def gen_foldertree(scriptdir,DCdir,workdir="./"):
	'''

	gen_foldertree (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Generates a standard EMERGE folder tree

	Arguments
	----------
	  workdir : str 
		(Default = ./)
		Working directory
	  path2scripts : str
		Absolute path to 

	Outputs
	----------
	  Generate a series of subfolders: ["scripts_tmp","vis_raw", "vis_tmp", "FITS_raw", "FITS_tmp", "imaging_tmp", "reports", "products"]

	'''

	print("==================================================================")
	print(" Starting gen_foldertree() ... ")
	# Starting
	print(" Geneating EMERGE folder tree...")

	# list of folders
	folder_tree = ["scripts_tmp","vis_raw", "vis_tmp", "FITS_raw", "FITS_tmp", "imaging_tmp", "reports", "products"]

	for myfolders in folder_tree:
		# Does output directory exist?
		if not os.path.exists(str(workdir) + str(myfolders)):
			os.makedirs(str(workdir) + str(myfolders))
			print(str(str(workdir) + str(myfolders)) + " directory created successfully")
			# Copy scripts if path exist
			if (myfolders == "scripts_tmp"):
				os.system("cp -r " + str(scriptdir) + "emerge_DCpars_template.py ./scripts_tmp/.")
				os.system("cp -r " + str(scriptdir) + "emerge_run_assess_ms_public_template.py ./scripts_tmp/.")
				os.system("cp -r " + str(DCdir) + " ./scripts_tmp/DataComb")
		else:
    			print(str(str(workdir) + str(myfolders)) + " directory already exists. Continuing.")

	# Done
	print(" gen_foldertree() ... Done")
	print("==================================================================")


def rm_foldertree(workdir="./",rmproducts=False):
	'''

	rm_foldertree (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Removed standard EMERGE folder tree

	Arguments
	----------
	  workdir : str 
		(Default = ./)
		Working directory
	  rmproducts : boolean
		(Default = False)
		Also remove /products & /reports? 

	Outputs
	----------
	  Remove ALL EMERGE subfolders: ["scripts_tmp","vis_raw", "vis_tmp", "FITS_raw", "FITS_tmp", "imaging_tmp", "reports", "products"]

	'''

	print("==================================================================")
	print(" Starting rm_foldertree() ... ")
	# Starting
	print(" Removing EMERGE folder tree...")

	# list of folders
	folder_tree = ["scripts_tmp","vis_raw", "vis_tmp", "FITS_raw", "FITS_tmp", "imaging_tmp"]#, "reports", "products"]

	for myfolders in folder_tree:
		# Does output directory exist?
		if os.path.exists(str(workdir) + str(myfolders)):
			print(" Removing " + str(workdir) + str(myfolders))
			os.system("rm -rf " + str(workdir) + str(myfolders))
	if (rmproducts):
		print(" Also removing products!")
		os.system("rm -rf " + str(workdir) + "reports")
		os.system("rm -rf " + str(workdir) + "products")

	# Done
	print(" rm_foldertree() ... Done")
	print("==================================================================")


def listofMSfiles(datadir="./",listexpr=".cal",ends=False):
	'''

	listofMSfiles (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Generates a list of .MS files such as path2MS/listexpr

	Arguments
	----------
	  path2MS : str 
		(Default = ./)
		Path to .MS file to be examined
	  listexpr: str
		(Default = *.cal)
		Condition
	  ends : boolean
		(Default = False)
		If endswith is required (True), otherwise taken at any position
	Outputs
	----------
	  Returns a list with all folder names within path2MS 

	'''
	# list of ALL contents in pathtoMS
	entries = os.listdir(datadir)

	# Select only those with pathtoMS/listexpr
	if (ends):
		myfolders = [entry for entry in entries if os.path.isdir(os.path.join(datadir, entry)) and entry.endswith(listexpr)]
	else:
		myfolders = [entry for entry in entries if os.path.isdir(os.path.join(datadir, entry)) and listexpr in entry]
	
	# Return a list of folders
	return myfolders


def listofFITS(datadir="./",listexpr=".fits"):
	'''

	listofFITS (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Generates a list of .MS files such as path2MS/listexpr

	Arguments
	----------
	  path2MS : str 
		(Default = ./)
		Path to .MS file to be examined
	  listexpr: str
		(Default = *.fits)
		Ending condition

	Outputs
	----------
	  Returns a list with all folder names within path2FITS files 

	'''
	# list of ALL contents in datadir
	entries = os.listdir(datadir)

	# Select only those with pathtodata/listexpr
	myfiles = [entry for entry in entries if os.path.isfile(os.path.join(datadir, entry)) and listexpr in entry]
	# Add folder
	myfiles = [str(datadir) + filename for filename in myfiles]

	# Return a list of folders
	return myfiles



def do_copyMSfromArchive(datadir="./archive/",listexpr="*",workdir="./vis_raw/."):
	'''

	do_copyMSfromArchive (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Generates a list of .MS files such as path2MS/listexpr

	Arguments
	----------
	  datadir : str 
		(Default = ./archive)
		Path to data folder. Usually path to EMERGE .MS local archive
	  listexpr: str
		(Default = *)
		Ending condition
	  workdir : str 
		(Default = ./vis_raw)
		Path to local folder with _only_ relevant .MS files

	Outputs
	----------
	  Copy data

	'''
	# Copy all files (usually .MS files/folders) from path2Archive/listexpr* to path2rawvis
	os.system("cp -r " + str(datadir) + str(listexpr) + " " + str(workdir)) 


##-----------------------------------------------------

def getlocalDB(targetcoord,targetfreq_Hz,mindV,maxdV,R0):
	'''

	do_getlocalDB (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function to create a local database with project information

	Arguments
	----------
	  targetcoord: list
		Target coordinates in degrees
	  targetfreq_Hz: float OR str
		If float [in Hz] then uses it as skyfreq
		if str then uses name in emerge_lincat_full to get skyfreq
	  mindV,maxdV: float
		min/max channelwidth in ms
	  R0: float
		Searching radius in deg

	Outputs
	----------
	  dataframe with all projects satisfying the above selection criteria

	'''
	import alminer

	# Run an ALminer query with the general searching parameters defined in emerge_setip.py
	myquery = alminer.conesearch(ra=targetcoord[0],dec=targetcoord[1],search_radius=R0*60.,public=True)

	skyfreq_GHz=targetfreq_Hz/1E9 # Line frequency in GHz

	print(" Selecitng data with Sky freq: " + str(skyfreq_GHz) + " GHz")

	localDB = myquery[(myquery["min_freq_GHz"] < skyfreq_GHz) & 
                       (myquery["max_freq_GHz"] > skyfreq_GHz) &
                       (myquery['vel_res_kms'] > (mindV/1E3)) &
                       (myquery['vel_res_kms'] < (maxdV/1E3))
                      ]

	# Return a dataframe with all projects in ASA satisfying the selection criteria
	return localDB


	# Create local database localDB
	#restfreq = emerge.emerge_linecata_full[targetline]*1E9	# Line frequency in GHz
	#locaDB = emerge.createDB(alldata=myquery,freq_GHz=restfreq,
	#	min_freqres_kms=mindV/1E3,max_freqres_kms=maxdV/1E3,FoV=0.,radius=R0,mytracer=targetline) 



def find_GOUS(targetcoord,targetline,targetfreq_Hz,mindV,maxdV,R0,get_mosaics=True):
	"""
	find_GOUS

	Indentifies which 12m datasets fulfill a series of selection criteria based on freq, resolution, etc... 

	Arguments
	----------
	  targetcoord: list
		target coordinates in degrees (e.g. [234.4,-5.67])
	  targetline: str
		Line name in emerge catalog
	  targetfreq_Hz: float
		Frequency of interest in Hz (default = 93.173E9 - N2H+ (1-0))
	  min/maxDV: float
		minumum/maximum frequency resolution required in m/s
	  get_mosaics: logical
		(default = True)
		Get only mosaics?
	  R0: float
		(default = 0.1 deg)
		Search radius in deg for other complementary datasets 

        
	"""
	import pandas as pd

	# Definitions
	ALMA_available = False
	ACA_available = False
	TP_available = False

	# Info
	print("===============================================")
	print(" Serching criteria:")
	print(" Freqeuncy = " + str(targetfreq_Hz/1E9) + " GHz")
	print(" Resolution = " + str(mindV/1E3) + " - " + str(maxdV/1E3) + " km/s")
	print("-----------------------------------------------")
	print(" Results:")

	# Run localDB
	selected = getlocalDB(targetcoord,targetfreq_Hz,mindV,maxdV,R0)	
	
	# Calculate baselines
	selected['shortbaseline_m']=get_shortbaseline(freq_GHz=selected['frequency'],MRS_arcsec=selected['spatial_scale_max']) 
	selected['longbaseline_m']=get_longbaseline(freq_GHz=selected['frequency'],res_arcsec=selected['s_resolution'])

	# Create a GOUS dictionary
	GOUSdict = {}
	
	# Selection per GOUS
	for i in np.unique(selected['group_ous_uid']):
		# Define 
		ALMA12mdata = pd.DataFrame([]); ALMA12mdata["member_ous_uid"] = []
		ACAdata = pd.DataFrame([]); ACAdata["member_ous_uid"] = []
		TPdata=pd.DataFrame([]); TPdata["member_ous_uid"] = []

		# Select all entries from the same GOUS
		ALMAGOUS = selected[selected['group_ous_uid'] == i]
		for j in np.unique(ALMAGOUS['target_name']):
			#proj_id = ALMA12mdata[ALMA12mdata['target_name'] == i[0]]['project_code'].values[0]
			print("========================================================================================")
			target = j
			proj_id = np.unique(ALMAGOUS['project_code'])
			print(" Target = " + str(target))
			print(" Proj. id = " + str(proj_id))
			print(" PI = " + str(ALMAGOUS['pi_name'].values[0]))
			print(" Title = " + str(ALMAGOUS['obs_title'].values[0]))
			print(" GOUS = " + str(np.unique(ALMAGOUS["group_ous_uid"])))
			print(" MOUS = " + str(np.unique(ALMAGOUS["member_ous_uid"])))
			#print(" ASDM = " + str(np.unique(ALMAGOUS["asdm_uid"])))

			# Separate ALMA (12m) data from Zero-spacing data (ACA+TP)
			ALMA12mdata = ALMAGOUS[(ALMAGOUS['shortbaseline_m'] >12.)]
			if (np.shape(ALMA12mdata)[0] != 0):
				print("--------------------------------------")
				print("  ALMA 12m data available")
				print(" phasecenter " + str(ALMA12mdata["RAJ2000"].values[0]) + " " + str(ALMA12mdata["DEJ2000"].values[0]))

				# info about 12m data
				print("--------------------------------------")
				print(" 12m - longest baseline = " + str(np.round(np.min(ALMA12mdata['longbaseline_m']),1)) + " m")
				print(" 12m - spatial resolution = " + str(np.round(np.min(ALMA12mdata['ang_res_arcsec']),1)) + " arcsec")
				print(" 12m - shortest baseline = " + str(np.round(np.min(ALMA12mdata['shortbaseline_m']),1)) + " m")
				print(" 12m - MRS = " + str(np.round(np.min(ALMA12mdata['LAS_arcsec']),1)) + " arcsec")
				print(" 12m - max vel. resolution = " + str(np.round(np.min(ALMA12mdata['vel_res_kms']),3)) + " km/s")
				print(" GOUS = " + str(np.unique(ALMA12mdata["group_ous_uid"])))
				print(" MOUS = " + str(np.unique(ALMA12mdata["member_ous_uid"])))
				#print(" ASDM = " + str(np.unique(ALMA12mdata["asdm_uid"])))
			#else:
			#	ALMA12mdata["member_ous_uid"] = ["None"]

			# Zero-Spacing available?
			zs_data = ALMAGOUS[ALMAGOUS['shortbaseline_m'] <= 12. ]
			if (np.shape(zs_data)[0] != 0):
				ACAdata = zs_data[zs_data['longbaseline_m'] >= 7.]
				TPdata = zs_data[zs_data['longbaseline_m'] < 7.]
				if (np.shape(ACAdata)[0] != 0):
					ACA_available = True
					print("--------------------------------------")
					print(" ACA data available")
					print(" phasecenter " + str(ACAdata["RAJ2000"].values[0]) + " " + str(ACAdata["DEJ2000"].values[0]))
					print("--------------------------------------")
					# info about 7m data
					print(" 7m - longest baseline = " + str(np.round(np.min(ACAdata['longbaseline_m']),1)) + " m")
					print(" 7m - spatial resolution = " + str(np.round(np.min(ACAdata['ang_res_arcsec']),1)) + " arcsec")
					print(" 7m - shortest baseline = " + str(np.round(np.min(ACAdata['shortbaseline_m']),1)) + " m")
					print(" 7m - MRS = " + str(np.round(np.min(ACAdata['LAS_arcsec']),1)) + " arcsec")
					print(" 7m - max vel. resolution = " + str(np.round(np.min(ACAdata['vel_res_kms']),3)) + " km/s")
					print(" GOUS = " + str(np.unique(ACAdata["group_ous_uid"])))
					print(" MOUS = " + str(np.unique(ACAdata["member_ous_uid"])))
					#print(" ASDM = " + str(np.unique(ACAdata["asdm_uid"])))
					#print(" 7m - MOUS = " + str(np.unique(ACAdata["member_ous_uid"])))
				else:
					print(" ACA data NOT available")
					#ACAdata["member_ous_uid"] = ["None"]

				if (np.shape(TPdata)[0] != 0):
					TP_available = True
					print("--------------------------------------")
					print(" TP data available")
					print(" phasecenter " + str(TPdata["RAJ2000"].values[0]) + " " + str(TPdata["DEJ2000"].values[0]))
					print("--------------------------------------")
					print(" TP - max vel. resolution = " + str(np.round(np.min(TPdata['vel_res_kms']),3)) + " km/s")
					print(" GOUS = " + str(np.unique(TPdata["group_ous_uid"])))
					print(" MOUS = " + str(np.unique(TPdata["member_ous_uid"])))
					#print(" ASDM = " + str(np.unique(TPdata["asdm_uid"])))
					#print(" TP - MOUS = " + str(np.unique(TPdata["member_ous_uid"])))
				else:
					print(" TP data NOT available in project.")
					# Trying if there are TP in any other
					TPdata = selected[selected['longbaseline_m'] < 7.]
					if (np.shape(TPdata)[0] != 0):
						print(" Taking TP from GOUS " +str(np.unique(TPdata["group_ous_uid"])))
					#else:
						#TPdata["member_ous_uid"] = ["None"]

			else:
				print(" WARNING: NO short-spacing data available!")
				#ACAdata["member_ous_uid"] = ["None"]
				#TPdata["member_ous_uid"] = ["None"]

                            
		# Return GOUS dictionary
		# If ACA and TP data available add them to the GOUS dictionary
		#if (ACA_available & TP_available):
		GOUSdict[str(i)] = {
			"ProjID" : str(proj_id),
			"fieldname" : str(j),
			"GOUS" : str(i),
			"MOUS" :{
				"MOUS_12m" : np.unique(ALMA12mdata["member_ous_uid"]),
				"MOUS_7m" : np.unique(ACAdata["member_ous_uid"]),
				"MOUS_TP" : np.unique(TPdata["member_ous_uid"])
				}
			}
	# Return a GOUS dictionary
	return GOUSdict



def get_shortbaseline(freq_GHz,MRS_arcsec):
	"""
	  get_shortbaseline

	  - Estimates shortest baseline based on the available metadata
	  - L_min[m] = 0.93*wavelength[m] / MRS[rad]
          	(see formulas in Handbook: https://almascience.eso.org/documents-and-tools/cycle8/alma-technical-handbook)
	  - Note that this is a bit reverse engeneering although this is the only way to do it from the information 
		provided by the archive

	  freq_GHz: float
		Central frequency in GHz
	  MRS_arcsec: float
		Maximum recoverable scale in arcsec
	"""
	wavelength_m=3.E8/(freq_GHz*1E9)
	MRS_rad = MRS_arcsec/206265.
	shortbaseline=0.93*wavelength_m/MRS_rad
	return shortbaseline

def get_longbaseline(freq_GHz,res_arcsec):
	"""
	get_longbaseline
	    
	  - Estimates longest baseline based on the available metadata
	  - L_max[m] = 0.574*wavelength[m] / resoution[rad]
		(see formulas in Handbook: https://almascience.eso.org/documents-and-tools/cycle8/alma-technical-handbook)
	- Note that this is a bit reverse engeneering although this is the only way to do it from the information
		provided by the archive
          
	freq_GHz: float
		Central frequency in GHz
	res_arcsec: float
		resolution in arcsec
	"""
	wavelength_m=3.E8/(freq_GHz*1E9)
	res_rad = res_arcsec/206265.
	longbaseline=0.574*wavelength_m/res_rad
	return longbaseline

def getTPfiles(allGOUS,skyfreq_Hz,mindV,maxdV):
	"""
	getTPfiles
	    
	Function to download TP FITS files from JVO archive
          
	allGOUS: GOUS dictionary (see aslo emerge.find_GOUS)

	"""
	# Add pandas module
	import pandas as pd

	# Download JVO FITS archive metadata, if needed
	if (os.path.exists("./alma-meta.psv.gz")):
	#os.system("wget -O alma-meta.psv.gz https://www.dropbox.com/scl/fi/mcpgjpem3jp4vw1efr9al/alma-meta.psv.gz?rlkey=vx6na5dzponyorz6cwpv4c30d&st=loul6nxx&dl=0")
		os.system("gunzip alma-meta.psv.gz")
        # read JVO metadata file
	print(" Reading data from JVO metadata file")
	JVOfiles = pd.read_csv("alma-meta.psv", sep = "|", skiprows=83)
	JVOfiles = JVOfiles[JVOfiles.iloc[:,12] != "GHz"]
	JVOfiles.iloc[:,7] = JVOfiles.iloc[:,7].astype(str)
	JVOfiles.iloc[:,12] = JVOfiles.iloc[:,12].astype(float)
	JVOfiles.iloc[:,13] = JVOfiles.iloc[:,13].astype(float)
	JVOfiles.iloc[:,42] = JVOfiles.iloc[:,42].astype(float)
	JVOfiles.iloc[:,54] = JVOfiles.iloc[:,54].astype(str)

	# Loop over all MOUSid within allGOUS
	FITScounter = 0
	minfV_Hz = float(mindV/3E8*skyfreq_Hz)
	maxfV_Hz = float(maxdV/3E8*skyfreq_Hz)
	print(" Selection: " + str(skyfreq_Hz) + " - " + str(minfV_Hz) + "-" + str(maxfV_Hz))

	for mm in allGOUS.keys():
		allMOUS = allGOUS[str(mm)]["MOUS"]["MOUS_TP"]; print(allMOUS)
		for myMOUS in allMOUS:
			print("---" + myMOUS + " ---" + allMOUS)
			if (len(myMOUS) != 0):
				# Adapt MOUS id to JVO format
				myMOUS2 = myMOUS.split("//")[1].replace("//","___").replace("/","_")
				# Identify linsk and names of the FITS files
				fits2download = JVOfiles[(JVOfiles.iloc[:,54]==str(myMOUS2)) & (JVOfiles.iloc[:,12] >= skyfreq_Hz/1E9) & (JVOfiles.iloc[:,13] <= skyfreq_Hz/1E9) & (JVOfiles.iloc[:,7] == "CUBE")  & (JVOfiles.iloc[:,42] <= maxfV_Hz) & (JVOfiles.iloc[:,42] >= minfV_Hz)]
				#fits2download = JVOfiles[(JVOfiles.iloc[:,54]==str(myMOUS2))  & (JVOfiles.iloc[:,7] == "CUBE") ]
				#fits2download = JVOfiles[(JVOfiles.iloc[:,12] >= skyfreq_Hz/1E9) & (JVOfiles.iloc[:,13] <= skyfreq_Hz/1E9)]
				# Download all FITS files
				for FITSfiles,FITSnames in zip(fits2download.iloc[:,59],fits2download.iloc[:,55]):
					os.system("wget -O ./FITS_raw/" + str(FITSnames) +  " " + str(FITSfiles))
					FITScounter += 1
	# Remiving weigth files 
	os.system("rm -rf ./FITS_raw/*weight*")

	# Return number of FITS files
	return FITScounter



#######################################################
# Section #1: Pipeline tools
#######################################################

def getMSinfo(myMSfile):
	'''

	getMSinfo (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function to extract the relevant information for an .MS file

	Arguments
	----------
	  myMSfile: string
		Path to .MS file to be examined

	Outputs
	----------
	  Dictionary including all relevant info for eahc of the science targets
	  inside the .MS file
	  Keywords:
	  [source_name]
		[phasecenter] with RAdeg,DEdeg
		[mosaic] with RAmindeg, RAmaxdeg, DEmindeg, DEmaxdeg
		[Source_Vel_ms-1] in ms-1
		[SPWinfo_LSRK_Hz] with one entry per SPS including [mean,min,max,chanwidth] in Hz plus [Nchans]
		[Antenna] information with [Diameter_m]
		[Baseline_info] baseline information including [antenna_diameter_m, L80_beamsize_arcsec, L05_MRS_arcsec, min_baeline_m, max_baeline_m]
		[Project_info] with info such as [ProjID, MOUS]

	'''

	print("=====================================================")
	print(" Starting getMSinfo...")
	print(" MS file:" + str(myMSfile))

	# Use analysisUtils to derive all propreties from MS file

	# Creat a (empty) dictionary to store the results
	MS_info = {}
	
	# Create a list of Science targets within the MSfile
	MS_Sciencetargets = np.unique(au.getScienceTargets(myMSfile))
	
	# Loop over all targets in the MS file
	for mytarget in MS_Sciencetargets:
		# Phase center
		# au.getPhaseCenterForField outputs = [Nfield, xmax, xmin, ymax, ymin] in arcsec
		print(" Science Target = " + str(mytarget))
		MS_mosaic = au.plotmosaic(myMSfile,doplot=False,sourceid=str(mytarget))
		# au.getRADec = [RA,DEC] in rad
		MS_phasecenter = au.getRADecForField(myMSfile,field=mytarget,forcePositiveRA=False)
		if (MS_phasecenter[0][0] < 0):
			# For some reason, some MS files have RA as negative value -> turn it positive
			MS_phasecenter[0][0] = 2.*np.pi+MS_phasecenter[0][0]

		# Trick to get the phase center right
		phasecenter_fieldID = MS_mosaic[0]	# get phasecenter field id
		mytb = ctools.table()
		mytb.open(myMSfile+'/FIELD')
		myphasedir = mytb.getcol('PHASE_DIR')
		MS_phasecenter[0][0] = myphasedir[0][0][phasecenter_fieldID]
		MS_phasecenter[1][0] = myphasedir[1][0][phasecenter_fieldID]
		mytb.close()

		if (MS_phasecenter[0][0] < 0):
			# For some reason, some MS files have RA as negative value -> turn it positive
			MS_phasecenter[0][0] = 2.*np.pi+MS_phasecenter[0][0]

		# SPWs information (for Science only fields)
		# au.spwToLSRK = [SPW, minfreq, maxfreq, chanwidth] one entry per SPW (in LSRK frame)
		MS_SPWsinfo_LSKR = au.spwToLSRK(myMSfile,field=mytarget,units='Hz')
		# Number of channels per SPW
		MS_SPWSinfo_Nchans =  au.getNChans(myMSfile,spw=MS_SPWsinfo_LSKR.keys())
		# Add Nchans per SPW to MS_SPWsinfo_LSKR info
		for i,j in zip(MS_SPWsinfo_LSKR.keys(),MS_SPWSinfo_Nchans):
			MS_SPWsinfo_LSKR[i]['Nchans'] = j	

		# source Velocity
		MS_SourveVel = au.radialVelocity(myMSfile,source=str(mytarget))

		# Proj. information
		MS_proj = au.projectCodeFromDataset(myMSfile)
		
		# Get Baseline information
		MS_baselines = au.getBaselineLengths(myMSfile)#,field=str(mytarget))
		MS_baselines = [sublist[1] for sublist in MS_baselines]

		# Get Antenna info
		MS_AntDiam = au.pickDishDiameter(myMSfile)

		# Get GOUS information
		MS_GOUS = extractgous(myMSfile)
		if (type(MS_GOUS) == bool):
			MS_GOUS = "XXXX"
		else:
			MS_GOUS = MS_GOUS[6:]	# skip "uid___"

		#
		# Create a dictionary with all relevant values
		MS_info[mytarget] = {
			# Phase center section
			'phasecenter':{'RAdeg': MS_phasecenter[0][0]*180./np.pi,
				'DEdeg': MS_phasecenter[1][0]*180./np.pi},
			# Mosaic parameters section
			'mosaic':{'RAmindeg': MS_phasecenter[0][0]*180./np.pi+MS_mosaic[2]/3600.,
				'RAmaxdeg': MS_phasecenter[0][0]*180./np.pi+MS_mosaic[1]/3600.,
				'DEmindeg': MS_phasecenter[1][0]*180./np.pi+MS_mosaic[4]/3600.,
				'DEmaxdeg': MS_phasecenter[1][0]*180./np.pi+MS_mosaic[3]/3600.},
			# Velocity
			'Source_Vel_ms-1': np.unique(MS_SourveVel),
			# SPW info section (in Hz) + Nchans
			'SPWinfo_LSRK_Hz': MS_SPWsinfo_LSKR,
			# Baseline information/estimates:
			'Baseline_info:':{
				# Get dish diameter
				'antenna_diameter_m': au.pickDishDiameter(myMSfile),
				# Estimate beamsize, MRS
				'L80_beamsize_arcsec': au.estimateSynthesizedBeam(myMSfile),	# L80%
				'L05_MRS_arcsec': au.estimateMRS(myMSfile),			# L05%
				'min_baeline_m': min(MS_baselines),
				'max_baeline_m': max(MS_baselines)},
			# Antenna
			'Antenna':{
			'Diameter_m': MS_AntDiam},
			# Project
			'Project_info':{
				'ProjID': MS_proj[0],
				'GOUS': MS_GOUS,
				'MOUS': MS_proj[2]
				}
		}

	#
	print(" ... Done")
	print("=====================================================")
	# Returns dictionary
	return MS_info


def extractgous(myms):
    """

	extract GOUS (EMERGE)
	
	Author: D. Petry + EMERGE team

    	Return GOUS UID for given CalMS-processed MS
    	as a string.
    	Returns False in case of problems.
    """
    mytb = ctools.table()
    try:
        mytb.open(myms+'/HISTORY')
        mymessage = mytb.getcell('MESSAGE',mytb.nrows()-1).split(' ')
        mytb.close()
    except:
        print('Error accessing MS ', myms)
        mytb.close()
        return False
    
    if mymessage[0] != 'GOUS_UID':
        print('Could not find GOUS UID.')
        return False

    if len(mymessage)!=3:
        print('Unexpected format: ', mymessage)
        return False

    return mymessage[2]


def extractgous2(myms):
    """
    DP (ESO)
    Return GOUS UID for given CalMS-processed MS
    as a string.
    Returns False in case of problems.
    """
    mytb = tbtool()
    rval = True
    try:
        mytb.open(myms+'/HISTORY')
        mymessage = mytb.getcell('MESSAGE',mytb.nrows()-1).split(' ')
        mytb.close()
    except:
        print('Error accessing MS ', myms)
        mytb.close()
        rval = False

    if rval:
        if mymessage[0] != 'GOUS_UID':
            print('Could not find GOUS UID.')
            rval = False

    if rval:
        if len(mymessage)!=3:
            print('Unexpected format: ', mymessage)
            rval = False

    if rval: # extraction from the HISTORY table has worked
        rval = mymessage[2]
    else: # try via Archive web service
        if type(myms) != str:
            return False

        import os
        
        if myms[-1]=='/':
            myms = myms[:-1]
        
        myasdm = os.path.basename(myms).split('.')[0].replace('___','://').replace('_','/')
        print('Trying webservice on ASDM UID ', myasdm)

        try:
            os.system('curl -L -d "LANG=ADQL&RESPONSEFORMAT=csv&QUERY=select top 1 group_ous_uid from ivoa.obscore WHERE asdm_uid = \''+myasdm+'\' " https://almascience.eso.org/tap/sync > asdmInfo.out 2>/dev/null')

            with open('asdmInfo.out') as f:
                mykey = f.readline()[:-1]
                if mykey == 'group_ous_uid':
                    rval = f.readline()[:-1]
                    if rval=='':
                        rval = False
                else:
                    print('Unexpected return value ', mykey)
                    rval = False

        except:
            print('Error while querying Archive ...')
            rval = False

    return rval




def getMSinfo_short(myMSfile):
	'''

	getMSinfo_short (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function to extract the relevant information for an .MS file

	Arguments
	----------
	  myMSfile: string
		Path to .MS file to be examined

	Outputs
	----------
	  Dictionary including all relevant info for eahc of the science targets
	  inside the .MS file
	  Keywords:
	  [source_name]
		[phasecenter] with RAdeg,DEdeg
		[mosaic] with RAmindeg, RAmaxdeg, DEmindeg, DEmaxdeg
		[Source_Vel_ms-1] in ms-1
		[SPWinfo_LSRK_Hz] with one entry per SPS including [mean,min,max,chanwidth] in Hz plus [Nchans]
		[Antenna] information with [Diameter_m]
		[Baseline_info] baseline information including [antenna_diameter_m, L80_beamsize_arcsec, L05_MRS_arcsec, min_baeline_m, max_baeline_m]
		[Project_info] with info such as [ProjID, MOUS]

	'''

	print("=====================================================")
	print(" Starting getMSinfo_short...")
	print(" MS file:" + str(myMSfile))

	# Use analysisUtils to derive all propreties from MS file

	# Creat a (empty) dictionary to store the results
	MS_info = {}
	
	# Create a list of Science targets within the MSfile
	MS_Sciencetargets = np.unique(au.getScienceTargets(myMSfile))
	
	# Loop over all targets in the MS file
	for mytarget in MS_Sciencetargets:
		# Phase center
		# au.getPhaseCenterForField outputs = [Nfield, xmax, xmin, ymax, ymin] in arcsec
		print(" Science Target = " + str(mytarget))
		MS_mosaic = au.plotmosaic(myMSfile,doplot=False,sourceid=str(mytarget))
		# au.getRADec = [RA,DEC] in rad
		MS_phasecenter = au.getRADecForField(myMSfile,field=mytarget,forcePositiveRA=False)
		if (MS_phasecenter[0][0] < 0):
			# For some reason, some MS files have RA as negative value -> turn it positive
			MS_phasecenter[0][0] = 2.*np.pi+MS_phasecenter[0][0]

		# Trick to get the phase center right
		phasecenter_fieldID = MS_mosaic[0]	# get phasecenter field id
		mytb = ctools.table()
		mytb.open(myMSfile+'/FIELD')
		myphasedir = mytb.getcol('PHASE_DIR')
		MS_phasecenter[0][0] = myphasedir[0][0][phasecenter_fieldID]
		MS_phasecenter[1][0] = myphasedir[1][0][phasecenter_fieldID]
		mytb.close()

		if (MS_phasecenter[0][0] < 0):
			# For some reason, some MS files have RA as negative value -> turn it positive
			MS_phasecenter[0][0] = 2.*np.pi+MS_phasecenter[0][0]
		
		# SPWs information (for Science only fields)
		# au.spwToLSRK = [SPW, minfreq, maxfreq, chanwidth] one entry per SPW (in LSRK frame)
		MS_SPWsinfo_LSKR = au.spwToLSRK(myMSfile,field=mytarget,units='Hz')
		# Number of channels per SPW
		MS_SPWSinfo_Nchans =  au.getNChans(myMSfile,spw=MS_SPWsinfo_LSKR.keys())
		# Add Nchans per SPW to MS_SPWsinfo_LSKR info
		for i,j in zip(MS_SPWsinfo_LSKR.keys(),MS_SPWSinfo_Nchans):
			MS_SPWsinfo_LSKR[i]['Nchans'] = j	

		# source Velocity
		MS_SourveVel = au.radialVelocity(myMSfile,source=str(mytarget))

		#
		# Create a dictionary with all relevant values
		MS_info[mytarget] = {
			# Phase center section
			'phasecenter':{'RAdeg': MS_phasecenter[0][0]*180./np.pi,
				'DEdeg': MS_phasecenter[1][0]*180./np.pi},
			# Mosaic parameters section
			'mosaic':{'RAmindeg': MS_phasecenter[0][0]*180./np.pi+MS_mosaic[2]/3600.,
				'RAmaxdeg': MS_phasecenter[0][0]*180./np.pi+MS_mosaic[1]/3600.,
				'DEmindeg': MS_phasecenter[1][0]*180./np.pi+MS_mosaic[4]/3600.,
				'DEmaxdeg': MS_phasecenter[1][0]*180./np.pi+MS_mosaic[3]/3600.},
			# Velocity
			'Source_Vel_ms-1': np.unique(MS_SourveVel),
			# SPW info section (in Hz) + Nchans
			'SPWinfo_LSRK_Hz': MS_SPWsinfo_LSKR
		}

	#
	print(" ... Done")
	print("=====================================================")
	# Returns dictionary
	return MS_info



def islineinMS(myMSfile,myrestfreq_Hz):
	'''

	islineinMS (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function determine wither a line (freq) is covered in the .MS file

	Arguments
	----------
	  myMSfile : string 
		Path to .MS file to be examined
	  myrestfreq_Hz: float
		Line rest frequency [Hz]

	Outputs:
	  Returns a True/False boolean

	'''

	# Defaults
	lineinMS = False
	targetFieldname = 'None'
	lineSPW = []
	targetFieldname = []
	chanresSPW = 1E9

	# Get info from MS
	MS_info = getMSinfo(myMSfile)

	# loop over all sources
	for sou in MS_info.keys():
		myfreq = getskyfreq(myrestfreq_Hz,MS_info[sou]['Source_Vel_ms-1'][0])
		# loop over all SPWs
		for spws in MS_info[sou]['SPWinfo_LSRK_Hz'].keys():
			if ((myfreq > MS_info[sou]['SPWinfo_LSRK_Hz'][spws]['min']) & (myfreq < MS_info[sou]['SPWinfo_LSRK_Hz'][spws]['max'])):
				# If inside the SPW then sets boolean to True
				lineinMS = True
				# 
				targetFieldname.append(sou)
				if (MS_info[sou]['SPWinfo_LSRK_Hz'][spws]['chanwidth'] < chanresSPW):
					chanresSPW = MS_info[sou]['SPWinfo_LSRK_Hz'][spws]['chanwidth']
					lineSPW.append(spws)
					print(" Line is included in SPW = " +str(spws))

	# Returns True/False boolean	
	lineSPW = np.unique(lineSPW)
	targetFieldname = np.unique(targetFieldname)
	return lineinMS, targetFieldname, lineSPW , chanresSPW



def istargetinMS(myMSfile,mycoords,R0,units='deg'):
	'''

	istargetinMS (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function determine wither a line (freq) is covered in the .MS file

	Arguments
	----------
	  myMSfile : str
		Path to .MS file to be examined
	  mycoords : list or str 
		Target coordinates [RA,DE] if in degrees
		e.g., [REdeg,DEdeg] = [4.232545,-32.9452343]
		     target coordinates ['RADEC string'] if in absolute coordinates
		      [RAabs,DEabs]] = ['12:56:11.16658 -05.47.21.5246']
	  units : string
		(default = 'deg')
		options: 'deg' or 'abs'

	Outputs
	----------
	  Returns a True/False boolean + the corresponding field name

	'''

	# Boolean set to False
	targetinMS = False
	# Field name associated to target
	targetFieldname = 'None'

	# Convert from Absolute to degrees if necessary
	if (units == 'abs'):
		mycoords = au.radec2deg(mycoords)

	# Get info from MS
	MS_info = getMSinfo(myMSfile)
	
	# loop over all sources
	targetFieldname = []
	for sou in MS_info.keys():
		# Calculate distance to mosaic
		xoff = mycoords[0]-MS_info[sou]['phasecenter']['RAdeg']
		yoff = mycoords[1]-MS_info[sou]['phasecenter']['DEdeg']
		Dmosaic = np.sqrt(xoff**2.+yoff**2.)
		if (Dmosaic <= R0):
			# If inside the mosaic bounduaries then sets boolean to True
			targetinMS = True
			print(" Target is included in MS file")
			targetFieldname.append(sou)
#		# Is wihtin the mosaic?
#		if ((mycoords[0] > MS_info[sou]['mosaic']['RAmindeg']) &
#			(mycoords[0] < MS_info[sou]['mosaic']['RAmaxdeg']) &
#			(mycoords[1] > MS_info[sou]['mosaic']['DEmindeg']) &
#			(mycoords[1] < MS_info[sou]['mosaic']['DEmaxdeg']) ):
			# If inside the mosaic bounduaries then sets boolean to True
#3			targetinMS = True
#			print(" Target is included in MS file")
#			# 
#			targetFieldname = sou

	targetFieldname = np.unique(targetFieldname)
	# Returns True/False boolean	
	return targetinMS, targetFieldname




def istarANDlineinMS(myMSfile,myrestfreq,mycoords,R0,units='deg'):
	'''

	istarANDlineinMS (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function determine wither a line (freq) AND the position are covered in the .MS file

	Arguments
	----------
	  myMSfile : string 
		Path to .MS file to be examined
	  myrestfreq: float
		Line rest frequency [Hz]
	  mycoords : list or str 
		Target coordinates [RA,DE] if in degrees
		e.g., [REdeg,DEdeg] = [4.232545,-32.9452343]
		     target coordinates ['RADEC string'] if in absolute coordinates
		      [RAabs,DEabs]] = ['12:56:11.16658 -05.47.21.5246']
	  R0 : float
		Search radius in [deg]

	Outputs:
	  Returns a True/False boolean

	'''

	# Defaults
	lineinMS = False
	targetFieldname = 'None'
	lineSPW = []
	targetFieldname = []

	# Boolean set to False
	targetinMS = False

	# Convert from Absolute to degrees if necessary
	if (units == 'abs'):
		mycoords = au.radec2deg(mycoords)

	# Get info from MS
	MS_info = getMSinfo(myMSfile)
	
	# loop over all sources
	for sou in MS_info.keys():
		chanresSPW = 1E9
		# For whatever reason sometimes this keyword is a float, sometimes a list...
		if (type(MS_info[sou]['Source_Vel_ms-1']) == float):
			# If float
			myfreq = getskyfreq(myrestfreq,MS_info[sou]['Source_Vel_ms-1'])
		else:
			# If list, take the first element
			myfreq = getskyfreq(myrestfreq,MS_info[sou]['Source_Vel_ms-1'][0])
		
		# loop over all SPWs
		#print(" (0/3) Checking MS...")
		for spws in MS_info[sou]['SPWinfo_LSRK_Hz'].keys():
			if ((myfreq > MS_info[sou]['SPWinfo_LSRK_Hz'][spws]['min']) & (myfreq < MS_info[sou]['SPWinfo_LSRK_Hz'][spws]['max'])):
				# If inside the SPW then sets boolean to True
				lineinMS = True
				#print(" (1/3) Line frequency covered by SPW" + str(spws) + " in the MS file")
				# 
				if (MS_info[sou]['SPWinfo_LSRK_Hz'][spws]['chanwidth'] < chanresSPW):
					#print(" (2/3) MS file contains right channel width")
					chanresSPW = MS_info[sou]['SPWinfo_LSRK_Hz'][spws]['chanwidth']
					

					# Calculate distance to mosaic
					xoff = mycoords[0]-MS_info[sou]['phasecenter']['RAdeg']
					yoff = mycoords[1]-MS_info[sou]['phasecenter']['DEdeg']
					Dmosaic = np.sqrt(xoff**2.+yoff**2.)
					if (Dmosaic <= R0):
						# If inside the mosaic bounduaries then sets boolean to True
						targetinMS = True
						lineSPW.append(spws)
						targetFieldname.append(sou)
						print(" Selected source = " +str(sou))
						print(" Selected SPW = " +str(spws))

					# Is wihtin the mosaic?
#					if ((mycoords[0] > MS_info[sou]['mosaic']['RAmindeg']) &
#						(mycoords[0] < MS_info[sou]['mosaic']['RAmaxdeg']) &
#						(mycoords[1] > MS_info[sou]['mosaic']['DEmindeg']) &
#						(mycoords[1] < MS_info[sou]['mosaic']['DEmaxdeg']) ):
#						# If inside the mosaic bounduaries then sets boolean to True
#						targetinMS = True
#						#print(" (3/3) Target position is included in MS file")
#						# 
#						targetFieldname = sou

	# Returns True/False boolean	
	lineSPW = np.unique(lineSPW)
	targetFieldname = np.unique(targetFieldname)
	return targetinMS, targetFieldname, lineSPW , chanresSPW



def getskyfreq(myrestfreq,myVLSR):
	'''

	getskyfreq (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function determine skyfreq

	Arguments:
	  myrestfreq : float
		Molecule rest frequency [Hz]
	  myVLSR : float
		Source VLSR [m/s]

	Outputs:
	  Returns skyfreq in Hz

	'''
	# Calculate skyfreq based on rest freq and VLSR
	skyfreq = myrestfreq*(1.-myVLSR/c.c.value)
	#
	return skyfreq




def genContSubMask(myMSfile,mysource,linecat=emerge_linecat,linetol_LSRK_ms=15.E3,targetvel_LSRK_ms=-1E9):
	'''

	genContSubMask (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function determine skyfreq

	Arguments
	----------
	  myMSfile : string 
		Path to .MS file to be examined
	  mysource : string
		Source name
	  linecat: dictionary (optional)
		(Default = emerge_linecat)
		Line catalog used for line identification
	  linetol_LSRK_ms: float (optional)
		(Default = 15E3 m/s)
		Line tolereance in LSRK [m/s] to be used as line window
	  targetvel_LSRK_ms: float (optional)
		(Default = -1E9 -> will use .MS source velocity)
		Target LRSK velocity [m/s]

	Outputs
	----------
	  Returns a string with the expected SPW mask to be used in ContSub
	  E.g., '25:109.747054~109.770492GHz;109.789230~109.794692GHz'

	'''

	print("=====================================================")
	print(" Starting genContSubMask...")
	print("-----------------------------------------------------")
	print(" Visbility file: " +str(myMSfile))

	# Defaults
	myContSubmask = []	

	# Get info from myMSfile
	MS_info = getMSinfo(myMSfile)
	
	# select target of interest
	mytarget = MS_info[mysource]

	# Create a dictionary with all SPWs
	SPW_lines = {key:[] for key in MS_info[mysource]['SPWinfo_LSRK_Hz'].keys()}

	# Add values of the band edges (in Hz)
	for myspws in MS_info[mysource]['SPWinfo_LSRK_Hz'].keys():
		SPW_lines[myspws].append(MS_info[mysource]['SPWinfo_LSRK_Hz'][myspws]['min'])
		SPW_lines[myspws].append(MS_info[mysource]['SPWinfo_LSRK_Hz'][myspws]['max'])
	
	# Get min/max frequency covered by observations; it will speed up line identification later on
	SPW_lines_minfreq = min(min(SPW_lines.values()))
	SPW_lines_maxfreq = max(max(SPW_lines.values()))  

	# Add extra key for lines inside bands
	SPW_lines["lines"]=[]

	# loop over all lines in dictionary
	for mylin in linecat.keys():

		# Get line frequency from catalog
		linefreq_res_Hz = linecat[mylin]*1E9
		#print(" Line rest freq. [Hz] = " + str(linefreq_res_Hz))

		# If target velocity not specified, then use MS velocity
		if (targetvel_LSRK_ms == -1E9):
			if (type(mytarget['Source_Vel_ms-1']) == float):
				# If float
				target_velo = mytarget['Source_Vel_ms-1']
			else:
				# If list, take the first element
				target_velo = mytarget['Source_Vel_ms-1'][0]

		else:
			target_velo = targetvel_LSRK_ms
		
		linefreq_sky_Hz = getskyfreq(linefreq_res_Hz,target_velo)
		#print(" Line Rest freq. [MHz] = " + str(linefreq_res_Hz/1E6))
		#print(" Target velocity [km/s] = " + str(target_velo/1E3))
		#print(" Line sky freq. [MHz] = " + str(linefreq_sky_Hz/1E6))

		# Is line within range?
		if ((linefreq_sky_Hz > SPW_lines_minfreq) & (linefreq_sky_Hz < SPW_lines_maxfreq)):

			# Is this line/freq. in the MS file?
			isinMS, souinMS, lineSPW, linechanres = islineinMS(myMSfile,linefreq_res_Hz)

			# If line is observed in mysource then identify the SPWs
			if (mysource in souinMS):
				for myspws in lineSPW:
					# Identify the line and frequencies
					print(" Line = " + str(mylin) + " with freq.= " + str(linefreq_sky_Hz/1E9) + "GHz included in the SPW" + str(lineSPW) )
					# Store lines in variable (in Hz)
					SPW_lines[myspws].append(linefreq_sky_Hz)

					# Add line name identified inside
					SPW_lines["lines"].append(mylin)

	# generate a Contsub mask for the lines in the SPWs
	mykeys = list(SPW_lines.keys())
	# loop over all keys but last one ('lines' = string)
	for myspws in mykeys[:-1]:

		# If any line was added to the BW edges then generate the windows
		if (len(SPW_lines[myspws]) > 2):

			# Order freqs.
			myfreqs = np.sort(SPW_lines[myspws])

			# calculate tolerance based on input value and current freq.
			linetol_Hz = linetol_LSRK_ms/c.c.value*myfreqs[0]
			print("line tolerance = " +str(linetol_LSRK_ms) + " m/s OR "+ str(linetol_Hz) + "Hz")

			# Add windows without overlapping frequencies
			# This loop will go through all frequencies and add a small window around each
			# of them with Delta(freq)=+/-linetol_Hz. If overlapping with the previous/next
			# it wont be added
			for k in np.arange(0,len(myfreqs)):
				if (k == 0):
					myfreqs_final = [myfreqs[0]+linetol_Hz]
				if ((k > 0) & (k < (len(myfreqs)-1))):
					lowfreq =  myfreqs[k]-linetol_Hz
					highfreq =  myfreqs[k]+linetol_Hz
					if ((lowfreq > myfreqs_final[k-1]) & (lowfreq > (myfreqs[k-1]+linetol_Hz))):
						myfreqs_final.append(lowfreq)
					if (highfreq < (myfreqs[k+1]-linetol_Hz)):
						myfreqs_final.append(highfreq)
				if (k == (len(myfreqs)-1)):
					myfreqs_final.append(myfreqs[k]-linetol_Hz)

			# Convert them into GHz and make them a list of strings	
			myfreqs_final = [str(num/1E9) for num in myfreqs_final]
			myfreqs_final = ['~'.join(pair) for pair in zip(myfreqs_final[::2], myfreqs_final[1::2])]
			myfreqs_final = [num+"GHz" for num in myfreqs_final]
			print(" SPW: " + str(myspws) + ", (LSRK)= " + str(myfreqs_final))
			
			# Generate mask
			# If first time, create a new one
			mymask = au.LSRKRangesToTopo(myMSfile, myfreqs_final, spw=myspws, field=mysource)
			mymask = mymask[next(iter(mymask))]
			print(" SPW: " + str(myspws) + ", (TOPO)= "+str(mymask))
			if (myContSubmask == []):	
				myContSubmask = mymask
			# Otherwise append new SPWS window mask
			else:
				myContSubmask = myContSubmask + "," + mymask

	# Final mask
	print("-----------------------------------------------------")
	print(" Final mask (TOPO) = "+str(myContSubmask))
	print(" Masked lines : " + str(SPW_lines['lines']))
	print(" genContSubMask() ... Done")
	print("=====================================================")
	# Return mask
	return myContSubmask


def do_contsub(myMSfile,mysource,workdir="./vis_tmp",mycat="",linetol_LSRK_ms=15.E3,targetvel_LSRK_ms=-1E9,myfitorder=1):
	'''

	do_contsub (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function to carry out uvcontsub in a give .MS file

	Arguments
	----------
	  myMSfile : string 
		Path to .MS file to be examined
	  mysource : string
		Source name
	  workdir : string (optional)
		(Default = ./)
		Path to output/working folder
	  mycat: dictionary (optional)
		(Default = emerge_linecat)
		Line catalog used for line identification
	  linetol: float (optional)
		(Default = 15E3 m/s)
		Line tolereance in LSRK [m/s] to be used as line window
	  targetvel_LSRK_ms: float (optional)
		(Default = -1E9 -> will use .MS source velocity)
		Target LRSK velocity [m/s]
	  myfitorder: integer (optional)
		Fit order for baseline

	Outputs
	----------
	  Generates a continuum subtracted visibility file (.contsub): path2output+MOUS+".contsub"

	'''

	print("=====================================================")
	print(" Starting do_contsub() ...")

	# Does output directory exist?
	print(" Output directory = " + str(workdir))
	if not os.path.exists(workdir):
    		os.makedirs(workdir)
    		print(str(workdir) + " directory created successfully")
	else:
    		print(str(workdir) + " directory already exists. Continuing.")


	# Get info from myMSfile
	MS_info = getMSinfo(myMSfile)

	# line catalog
	if (mycat == ""):
		linecat = emerge_linecat
	else:
		linecat = mycat

	# Identify relevant SPWS
	myspws = [sublist for sublist in MS_info[mysource]['SPWinfo_LSRK_Hz']]
	myspws = ",".join(map(str,myspws))

	# Create uvcontmask
	uvcontmask = genContSubMask(myMSfile,mysource,linecat,linetol_LSRK_ms)#,targetvel_LSRK_ms)
	print(" Final mask (TOPO) = "+str(uvcontmask))
	print("-----------------------------------------------------")

	# execute uvcont substraction
	print(" Executing uvcontsub()...")
	print(" Target MS = " + str(myMSfile))
	print(" SPWs = " + str(myspws))

	# Define output path + name
	if not workdir.endswith("/"):
		workdir += "/"

	# Create output path
	split_myMSfile = myMSfile.split("/")
	##myoutputvis = workdir+str(split_myMSfile[-1])+".contsub"
	# NOTE: somehow the uvcont later on misses the Project information, we then store it here
	##---- myoutputvis = workdir+str(split_myMSfile[-1])+ "_Proj" + str(MS_info[mysource]['Project_info']['ProjID'])+ "_GOUS" + str(MS_info[mysource]['Project_info']['GOUS']) + "_Sou" + str(mysource) + "_Dim" + str(int(MS_info[mysource]['Antenna']['Diameter_m'])) + "m" + ".contsub"
	myoutputvis = workdir+str(split_myMSfile[-1])+ "_Proj" + str(MS_info[mysource]['Project_info']['ProjID']) + "_Sou" + str(mysource) + "_Dim" + str(int(MS_info[mysource]['Antenna']['Diameter_m'])) + "m" + ".contsub"

	# Execute uvcontsub depending on the CASA version
	casaversion = casalog.version()[13:18]
	if ((casaversion >= "6.5.3")):
		ct.uvcontsub(vis=myMSfile, outputvis=myoutputvis, spw=myspws, fitspec=uvcontmask, fitorder=myfitorder)
	else:
		ct.uvcontsub(vis=myMSfile, spw=myspws, fitspw=uvcontmask, fitorder=myfitorder,combine='spw')
		# Rename outputfile to make it easier to track
		os.system("mv " + str(myMSfile) + ".contsub " + str(myoutputvis))

	# Return file name with output visibilities
	print(" New uvcontsub visbilities stored in : " + myoutputvis)
	return myoutputvis
	print("-----------------------------------------------------")

	print(" do_consub() ... Done")
	print("=====================================================")


def do_split_singelspw(myMSfile,mysource,myrestfreq_Hz,mindV,maxdV,workdir="./vis_tmp"):
	'''

	do_split_singelspw (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function to split a single spw from a give .MS file (usually after uvcontsub)

	Arguments
	----------
	  myMSfile : string 
		Path to .MS file to be examined
	  mysource : string
		Source name
	  myrestfreq_Hz: float
		Line rest frequency [Hz]
	  workdir : string (optional)
		(Default = ./)
		Path to output/working folder

	Outputs
	----------
	  Generates a split file (.split): path2output+MOUS+".split"

	'''

	print("=====================================================")
	print(" Starting do_split() ...")

	# Get info from myMSfile
	MS_info = getMSinfo(myMSfile)
	
	# Check if line is in .MS file
	lineinMS, targetFieldname, lineSPW , chanresSPW = islineinMS(myMSfile=myMSfile,myrestfreq_Hz=myrestfreq_Hz)

	# If so, then split this SPW
	print("-----------------------------------------------------")
	print(" Executing split()...")

	# Define output path + name
	if not workdir.endswith("/"):
		workdir += "/"

	split_myMSfile = myMSfile.split("/")

	# Carry out split
	myoutputvis_list = []
	for myspws in lineSPW:
		# Get MS info again
		MSinfo = getMSinfo(myMSfile)
		chanwidth_info = MSinfo[mysource]["SPWinfo_LSRK_Hz"][myspws]['chanwidth']
		meanfreq_info = MSinfo[mysource]["SPWinfo_LSRK_Hz"][myspws]['mean']
		velres_info = c.c.value*chanwidth_info/meanfreq_info

		# split only those windows that satisfy the resolution criteria
		if ((velres_info >= mindV) & (velres_info <= maxdV)):
			#myoutputvis = workdir+split_myMSfile[-1]+".spw"+str(myspws)
			myoutputvis = workdir+split_myMSfile[-1]+"_spw"+str(myspws)+"_f0"+str(np.round(myrestfreq_Hz/1E9,5))+"GHz"
			ct.split(vis=myMSfile,outputvis=myoutputvis,field=mysource,spw=str(myspws),datacolumn='all')
			myoutputvis_list.append(myoutputvis)

	# Remove previous MSfile to get only the split one
	os.system("rm -rf " + myMSfile)

	# Return file name with output visibilities
	print(" New split visbilities stored in : " + str(myoutputvis_list))
	return myoutputvis_list
	print("-----------------------------------------------------")

	print(" do_split() ... Done")
	print("=====================================================")



def do_regridspw(myMSfile,myrestfreq_Hz,mysourceVel_ms=0.0,mychanwidth_ms=0.1,myNchan=300):
	'''

	do_regridspw (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Regrid the corresponding .MS file (after uvcontsub & split)

	Arguments
	----------
	  myMSfile : str 
		Path to .MS file to be examined
	  myrestfreq_Hz: float
		Line rest frequency [Hz]
	  mysourceVel_ms: float
		(Default = 0.0 [m/s])
		Central velocitity of the spectrum [m/s], usually the source velocity
	  mychanwidth: float
		Channel width in [m/s]
	  myNchan: int
		Number of channels

	Outputs
	----------
	  Regrid target .MS file

	'''

	# SPW bandwidth [Hz] (symmetric wrt center)
	myBW_ms = myNchan*mychanwidth_ms 

	# Copy .MS file into a new one
	os.system("cp -r " + str(myMSfile) + " " + str(myMSfile) + "_reframed")

	# Regrid SPW
	ms.open(str(myMSfile) + "_reframed",nomodify=False)
	ms.regridspw(outframe='LSRK', mode='vrad', restfreq= myrestfreq_Hz, interpolation='LINEAR', center = mysourceVel_ms, bandwidth= myBW_ms, chanwidth=- mychanwidth_ms, hanning=True)
	
	ms.close()





#######################################################
# Section #2: Mosaic tools
#######################################################

##-------------------------------------------------------
# General Functions

# --- Coordinate transformation
def coord_hours2deg(mycoords="10:00:00 +30:00:00"):
	c = SkyCoord(mycoords, frame=ICRS, unit=(u.hourangle, u.deg))
	return c.ra.deg,c.dec.deg

def coord_deg2hours(mycoords="150.0 30.0"):
	c = SkyCoord(mycoords, frame=ICRS, unit=(u.deg, u.deg))
	return c.ra.to_string(unit=u.hour,sep=":",pad=True), c.dec.to_string(unit=u.degree,sep=":",pad=True)



# --- find all files within a folder given an expression and return the list
def read_allfiles(listexp="*",areimages=False):
	"""
	read_allfiles (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Generate a list of files following listexp criteria

	Arguments
	----------
	  listexp : string 
		(Default = * == all files)
		Criteria for ls command 
	  areimages : boolean
		False = FITS, True = CASA images (default False)
	
	Output
	----------
	  Returns a list of files

	"""
	# Generate a list of files
	if (areimages):
		os.system("ls -d " + str(listexp) + " > tmp.txt")
	else:
		os.system("ls " + str(listexp) + " > tmp.txt")
	with open('tmp.txt', 'r') as file:
    		file_content = file.readlines()

	listfiles = [line.rstrip('\n') for line in file_content]
	os.system("rm -rf tmp*")
	#
	return listfiles


# --- select relevant FITS
def select_FITS(mycoords=[0.0,0.0],units="deg",R0=0.5,freq0=100.E9,dVmin=0.0E3,dVmax=1.0E3,datadir=".",workdir=".",do_copy=True):
	"""

	select_FITS (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Identify FITS files that satisfy position + freq. selection criteria

	Arguments
	----------
	  mycoords : list of floats
		RA,DEC coordinates [deg]
	  units : str
		(Defaults = deg)
		RA,DEC units	
	  R0 : float
		(Default = 0.5 deg)
		search radius [deg] 
	  freq0 : float 
		Frequency [Hz]
	  dVmin : float (optional) 
		(Default = 0.0E3 m/s)
		minimum velocity resolution [km/s] 
	  dVmax : float (optional) 
		(Default = 1.0E3 m/s)
		maximum velocity resolution [km/s] 
	  datadir : string (optional) 
		(Default = ".")
		path to FITS data 
	  workdir : string (optional) 
		(Default = ".")
		working and writing directory  
	  do_copy : boolean
		(Default = True)
		Copy data from datadir to workdir

	Outputs
	----------
	  Copy those FITS files satisfying the previous selection criteria
	  found in the datadir folder into the workdir folder

	"""

	print("==============================================================")
	print(" select_FITS(): select FITS files overlaping with criteria")
	print("==============================================================")

	# Convert from Absolute to degrees if necessary
	if (units == 'abs'):
		mycoords = au.radec2deg(mycoords)

	# Check paths to datadir and workdir
	print(" FITS directory = " + str(datadir))
	if not os.path.exists(datadir):
		print(str(datadir) + " does not exist. Leaving script.")
		sys.exit()

	print(" Working directory = " + str(workdir))
	if not os.path.exists(workdir):
    		os.makedirs(workdir)
    		print(str(workdir) + " directory created successfully")
	else:
    		print(str(workdir) + " directory already exists. Continuing.")

	
	# Generate a list of FITS files from datadir
	os.system("rm -rf tmp.txt")
	os.system("ls " + str(datadir) + "/*.fits > tmp.txt")
	with open('tmp.txt', 'r') as file:
    		file_content = file.readlines()

	allFITS = [line.rstrip('\n') for line in file_content]

	# Loop over all files
	nfiles = 0
	for myFITS in allFITS: 
		
		#print(" Checking " + str(myFITS))

		# 0.- Get FITS information
		hdul = fits.open(myFITS)
		header = hdul[0].header
		# Identify axis
		for i in (np.arange(0,header['NAXIS'])+1):
			if (header["CTYPE"+str(i)] != ''):
				if (header["CTYPE"+str(i)][0] == 'R'):
					axis_RA = i
					RA_shape = header["NAXIS"+str(i)]
					RA_refpixel = header["CRPIX"+str(i)]
					RA_value = header["CRVAL"+str(i)]
					RA_incre = header["CDELT"+str(i)]
					#print("RA = " + str(RA_refpixel) + ", " + str(RA_value) + ", " + str(RA_incre))
				if (header["CTYPE"+str(i)][0] == 'D'):
					axis_DE = i
					DE_shape = header["NAXIS"+str(i)]
					DE_refpixel = header["CRPIX"+str(i)]
					DE_value = header["CRVAL"+str(i)]
					DE_incre = header["CDELT"+str(i)]
					#print("DE = " + str(DE_refpixel) + ", " + str(DE_value) + ", " + str(DE_incre))
				if (header["CTYPE"+str(i)][0] == 'F'):
					axis_freq = i
					freq_shape = header["NAXIS"+str(i)]
					freq_refpixel = header["CRPIX"+str(i)]
					freq_value = header["CRVAL"+str(i)]
					freq_incre = header["CDELT"+str(i)]
					freq_unit = header["CUNIT"+str(i)]
					#print("FREQ = " + str(freq_refpixel) + ", " + str(freq_value) + ", " + str(freq_incre))
					if (freq_unit[0:2] != "Hz" ):
						print(" WARNING: freq. units are in " + str(freq_unit))
					# is spectral axis in freq?
					spectralaxis_infreq = True
					vel_res_ms = abs(freq_incre/freq_value*3E8)	# in m/s
					

				if (header["CTYPE"+str(i)][0] == 'V'):
					axis_vel = i
					vel_shape = header["NAXIS"+str(i)]
					vel_refpixel = header["CRPIX"+str(i)]
					vel_value = header["CRVAL"+str(i)]
					vel_incre = header["CDELT"+str(i)]
					#vel_unit = header["CUNIT"+str(i)]
					restfreq = header["RESTFREQ"]
					#print("VEL = " + str(vel_refpixel) + ", " + str(vel_value) + ", " + str(vel_incre))
					# is spectral axis in freq?
					spectralaxis_infreq = False
					vel_res_ms = abs(vel_incre)	# in m/s

		

		# 1.- Selection criteria 
		# a) selection by velocity resolution
		if ((vel_res_ms >= dVmin) and (vel_res_ms <= dVmax)):
			# b) frequency inside BW
			# b.1. if working in freq. units
			if (spectralaxis_infreq):
				if (freq_incre > 0):
					freqmin = freq_value+(0.0-freq_refpixel)*freq_incre		
					freqmax = freq_value+(freq_shape-freq_refpixel)*freq_incre
				else:
					freqmax = freq_value+(0.0-freq_refpixel)*freq_incre		
					freqmin = freq_value+(freq_shape-freq_refpixel)*freq_incre
			# b.2. if working in vel. units
			else:
				if (vel_incre > 0):
					velmin = vel_value+(0.0-vel_refpixel)*vel_incre		
					velmax = vel_value+(vel_shape-vel_refpixel)*vel_incre
				else:
					velmax = vel_value+(0.0-vel_refpixel)*vel_incre		
					velmin = vel_value+(vel_shape-vel_refpixel)*vel_incre
				# Get values in freq. units
				freqmin = getskyfreq(restfreq,velmax)
				freqmax = getskyfreq(restfreq,velmin)

			#print(" freq = " + str(freqmin/1E9) + " - " + str(freqmax/1E9))

			if ((freq0 < freqmax) and (freq0 > freqmin)):
				# c) selection by position (phase center)
				RA_center = RA_value + (RA_shape/2.-RA_refpixel)*RA_incre
				DE_center = DE_value + (DE_shape/2.-DE_refpixel)*DE_incre
				Dist = np.sqrt((RA_center-mycoords[0])**2.+(DE_center-mycoords[1])**2.)
				#print(" R = " + str(RA_center) + " ," + str(DE_center)+ " ," + str(Dist))
				if (Dist <= R0):
					# then write into workdir
					nfiles = nfiles+1
					print(" --------------------------------------------")
					print(" Short spacing FITS file found: " + str(myFITS))
					print(" (RA,DEC) = " + str(RA_center) + " ," + str(DE_center))
					print(" (freq_min, freq_max) = (" + str(freqmin/1E9) + " ," + str(freqmax/1E9) + ") GHz")
					print(" DV = " + str(vel_res_ms/1E3) + " km/s")
					if (do_copy):
						print(" Copying file into " + str(workdir) + " directory")
						os.system("cp -r " + str(myFITS) + " " + str(workdir) + "/.")

	# Clean up
	os.system("rm -rf tmp*")

	if (nfiles == 0):
		print(" --------------------------------------------")
		print(" No files found that overlaps with the current criteria")
		return 0

	if (nfiles != 0):
		print(" --------------------------------------------")
		if (do_copy):
			print(" Files copied = " + str(nfiles))
			filelist = read_allfiles(listexp=str(workdir) + "/*.fits",areimages=False)
			return len(filelist)
		else:
			print(" Files found = " + str(nfiles))
			return nfiles

	print(" ... DONE")
	print("==========================================================")


# --- Extract lines from cubes
def extract_alllines(myimage,output_name="",linecat=emerge_linecat_full,VLSR_ms=0.0E3,Nchans=100, doFITS=False):
	"""
	
	extract_alllines (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function to extract ALL individual lines from spectral cubes in FITS format

	Arguments
	----------
	  myimage : string
		Image to be explored
	  output_name : string
		Name of the output file (if "" then = myimage)
	  linecat : string (optional) 
		(default = emerge_linecat)
		Line catalog 
	  VLSR_ms : float (optional) 
		(Default = 0.0E3)
		Source velocity in m/s 
	  Nchans : integer (optional) 
		(Default = 100)
		Number of channels to be extracted 
	  doFITS : boolean
		(Default = False)
		Create also FITS files?

	Outputs
	----------
	  If a line in linecat is included in the FITS file then
	  it creates all possible sub-cubes of Nchans centred at the line freq. (corrected by VLSR_ms)
	  labeled as "outputfile_cNchans_dVReskms-1"

	"""

	print("==============================================================")
	print(" extract_alllines(): extract a list of lines from cubes")
	print("==============================================================")

	# Info
	print("==========================================================")
	print(" Extracting line cubes from image: " + str(myimage) )
	# get axes information
	ia.open(myimage)
	axes = ia.coordsys().names()
	ia.done()
	# get frequencies
	freqaxis = axes.index("Frequency")
	refpixel = imhead(myimage,mode='get',hdkey='crpix'+str(freqaxis+1))
	refvalue = imhead(myimage,mode='get',hdkey='crval'+str(freqaxis+1))['value']
	refincre = imhead(myimage,mode='get',hdkey='cdelt'+str(freqaxis+1))['value']
	refshape = imhead(myimage,mode='get',hdkey='shape')
	# max-min freqs in the image
	if (refincre > 0):
		freqmin = refvalue+(0.0-refpixel)*refincre		
		freqmax = refvalue+(refshape[freqaxis]-refpixel)*refincre
	else:
		freqmax = refvalue+(0.0-refpixel)*refincre		
		freqmin = refvalue+(refshape[freqaxis]-refpixel)*refincre

	print(" Freq. range = [" + str(np.round(freqmin/1E9,6)) + " , " + str(np.round(freqmax/1E9,6)) + "] GHz")

	# Loop over linecat
	for mylin in np.arange(0,len(linecat)):
		linname = list(linecat.items())[mylin][0]
		linrestfreq = list(linecat.items())[mylin][1]*1E9
		skyfreq = getskyfreq(linrestfreq,VLSR_ms)
		
		# print(" Line: " + str(linname) + " - " + str(linfreq))
		# if freq. within band then create a subimage (freqs in GHz)
		if ((skyfreq > freqmin) and (skyfreq < freqmax)):
			print(" Creating line map: " + str(linname) + " - " + str(linrestfreq))
			print(" Sky Freq. = " + str(skyfreq/1E9) + " GHz")
			os.system("rm -rf " + str(output_name) + "_" + str(linname))
			# Get channels
			chanmin = np.round(skyfreq-Nchans/2.*abs(refincre),3)
			chanmax = np.round(skyfreq+Nchans/2.*abs(refincre),3)
			chansexp = "range=[" + str(float(chanmin/1E9)) + "GHz,"+ str(float(chanmax/1E9)) +"GHz]"
			print(" Extracting channels = " + chansexp)
			# Naming
			if (output_name == ""):
				output_file = myimage + "_" + str(linname) + "_c" + str(Nchans) + "_dV" + str(float(np.round(abs(refincre/(skyfreq)*3E5),3))) + "kms-1"
			else:
				output_file = output_name + "_" + str(linname) + "_c" + str(Nchans) + "_dV" + str(float(np.round(abs(refincre/(skyfreq)*3E5),3))) + "kms-1"
			
			# 
			imsubimage(imagename=myimage,outfile=output_file,chans=(chansexp),overwrite=True)
			# Set new Reference Freq.
			imhead(output_file,mode="put", hdkey="restfreq", hdvalue=str(linrestfreq)+"Hz")
			# Info
			print(" Filename: " + str(output_file))

			# FITS files
			if (doFITS):
				exportfits(imagename=str(output_file),fitsimage=str(output_file)+".fits")
	
	print(" ... DONE")
	print("==========================================================")



# --- Extract lines from cubes
def extract_targetline(myimage,myline,output_name="",linecat=emerge_linecat,VLSR_ms=0.0E3,Nchans=100, doFITS=False):
	"""
	
	extract_targetline (EMERGE)
	(similar to extract_lines() but just for one)
	
	Author: A.Hacar + EMERGE team
	
	Function to extract an individual line from spectral cubes in FITS format

	Arguments
	----------
	  myimage : str
		Image to be explored
	  myline : str
		Line to be extracted (must be part of linecat)
	  output_name : str (optional)
		Name of the output file (if "" then = myimage)
	  linecat : str (optional) 
		(default = emerge_linecat)
		Line catalog 
	  VLSR_ms : float (optional) 
		(Default = 0.0E3)
		Source velocity in m/s 
	  Nchans : integer (optional) 
		(Default = 100)
		Number of channels to be extracted 

	Outputs
	----------
	  If a line in linecat is included in the FITS file then
	  it creates all possible sub-cubes of Nchans centred at the line freq. (corrected by VLSR_ms)
	  labeled as "outputfile_cNchans_dVReskms-1"

	"""

	print("==============================================================")
	print(" extract_targetline(): extract one single target line from a cube")
	print("==============================================================")

	# Info
	print("==========================================================")
	print(" Extracting line cubes from image: " + str(myimage) )
	# get axes information
	ia.open(myimage)
	axes = ia.coordsys().names()
	ia.done()
	# get frequencies
	freqaxis = axes.index("Frequency")
	refpixel = imhead(myimage,mode='get',hdkey='crpix'+str(freqaxis+1))
	refvalue = imhead(myimage,mode='get',hdkey='crval'+str(freqaxis+1))['value']
	refincre = imhead(myimage,mode='get',hdkey='cdelt'+str(freqaxis+1))['value']
	refshape = imhead(myimage,mode='get',hdkey='shape')
	# max-min freqs in the image
	if (refincre > 0):
		freqmin = refvalue+(0.0-refpixel)*refincre		
		freqmax = refvalue+(refshape[freqaxis]-refpixel)*refincre
	else:
		freqmax = refvalue+(0.0-refpixel)*refincre		
		freqmin = refvalue+(refshape[freqaxis]-refpixel)*refincre

	print(" Freq. range = [" + str(np.round(freqmin/1E9,6)) + " , " + str(np.round(freqmax/1E9,6)) + "] GHz")

	# Loop over linecat
	linname = myline
	linfreq = linecat[myline]*1E9
	skyfreq = getskyfreq(linfreq,VLSR_ms)
	# print(" Line: " + str(linname) + " - " + str(linfreq))
	# if freq. within band then create a subimage (freqs in GHz)
	if ((skyfreq > freqmin) and (skyfreq < freqmax)):
		print(" Creating line map: " + str(linname) + " - " + str(linfreq))
		os.system("rm -rf " + str(output_name) + "_" + str(linname))
		# Get channels
		chanmin = np.round(skyfreq-Nchans/2.*abs(refincre),3)
		chanmax = np.round(skyfreq+Nchans/2.*abs(refincre),3)
		chansexp = "range=[" + str(float(chanmin)) + "Hz,"+ str(float(chanmax)) +"Hz]"
		print("Linefreq: " + str(linfreq/1E9) + " GHz; VLSR: " + str(VLSR_ms/1E3) + " km/s; Skyfreq: " + str(skyfreq/1E9) + " GHz")
		print("Channel selection: " + chansexp)
		# Naming
		if (output_name == ""):
			output_file = myimage
		else:
			output_file = output_name 
		output_file = output_file + "_" + str(linname) + "_Nchan" + str(Nchans) + "_dV" + str(float(np.round(abs(refincre/(skyfreq)*3E5),3)))  + "kms-1_f0" + str(np.round(linfreq/1E9,5))+"GHz"
		# 
		imsubimage(imagename=myimage,outfile=output_file,chans=(chansexp),overwrite=True)
		# Set new Reference Freq.
		imhead(output_file,mode="put", hdkey="restfreq", hdvalue=str(linfreq)+"Hz")
		# Info
		print(" Filename: " + str(output_file))
	
	print(" ... extract_targetline() DONE")
	print("==========================================================")




def extract_targetline_vres(myimage,myline,output_name="",linecat=emerge_linecat,VLSR_ms=0.0E3,target_dv=0.0,Nchans=100, doFITS=False):
	"""
	
	extract_targetline (EMERGE)
	(similar to extract_lines() but just for one)
	
	Author: A.Hacar + EMERGE team
	
	Function to extract an individual line from spectral cubes in FITS format

	Arguments
	----------
	  myimage : str
		Image to be explored
	  myline : str
		Line to be extracted (must be part of linecat)
	  output_name : str (optional)
		Name of the output file (if "" then = myimage)
	  linecat : str (optional) 
		(default = emerge_linecat)
		Line catalog 
	  VLSR_ms : float (optional) 
		(Default = 0.0E3)
		Source velocity in m/s 
	  target_dv : float (optional)
		(Default = 0.0E3)
		Target spectral resolution in m/s
		If target_dv > image freq.resolution then myimage will be smoothed out in frequency
	  Nchans : integer (optional) 
		(Default = 100)
		Number of channels to be extracted 

	Outputs
	----------
	  If a line in linecat is included in the FITS file then
	  it creates all possible sub-cubes of Nchans centred at the line freq. (corrected by VLSR_ms)
	  labeled as "outputfile_cNchans_dVReskms-1"

	"""

	print("==============================================================")
	print(" extract_targetline_vres(): extract one single target line from a cube")
	print("==============================================================")

	# Info
	print("==========================================================")
	print(" Extracting line cubes from image: " + str(myimage) )
	# get axes information
	ia.open(myimage)
	axes = ia.coordsys().names()
	ia.done()
	# get frequencies
	freqaxis = axes.index("Frequency")
	refpixel = imhead(myimage,mode='get',hdkey='crpix'+str(freqaxis+1))
	refvalue = imhead(myimage,mode='get',hdkey='crval'+str(freqaxis+1))['value']
	refincre = imhead(myimage,mode='get',hdkey='cdelt'+str(freqaxis+1))['value']
	refshape = imhead(myimage,mode='get',hdkey='shape')
	# max-min freqs in the image
	if (refincre > 0):
		freqmin = refvalue+(0.0-refpixel)*refincre		
		freqmax = refvalue+(refshape[freqaxis]-refpixel)*refincre
	else:
		freqmax = refvalue+(0.0-refpixel)*refincre		
		freqmin = refvalue+(refshape[freqaxis]-refpixel)*refincre

	print(" Freq. range = [" + str(np.round(freqmin/1E9,6)) + " , " + str(np.round(freqmax/1E9,6)) + "] GHz")

	#

	# Loop over linecat
	linname = myline
	linfreq = linecat[myline]*1E9
	skyfreq = getskyfreq(linfreq,VLSR_ms)

	# Compare target_dv with actual freq. resolution
	# If current dv smaller than target_dv, then smooth spectral axis to target_dv
	dvincre_ms = abs(refincre)/skyfreq*c.c.value	# in m/s; Velocity resolution of the data
	target_nu_Hz = target_dv*c.c.value/skyfreq		# in Hz; target spectral resolution

	if (dvincre < target_dv):
		# Create a tmp file
		os.system("cp -r " + str(myimage) +  " tmp_newchan.image")
		# Modify its spectral resolution
		imhead("tmp_newchan.image",hdkey="cdelt" + str(myidx+1),mode='put',hdvalue=str(target_nu_Hz))
		# Regrid original image with tmp file (with new dv) as template
		imregrid(imagename = myimage, template="tmp_newchan.image", output = "tmp_newchan2.image")
		# Replace original image
		os.system("rm -rf " + str(myimage))
		os.system("rm -rf tmp_newchan.image")
		os.system("mv tmp_newchan2.image " + str(myimage))
		# And read again refincre
		refincre = imhead(myimage,mode='get',hdkey='cdelt'+str(freqaxis+1))['value']

	# print(" Line: " + str(linname) + " - " + str(linfreq))
	# if freq. within band then create a subimage (freqs in GHz)
	if ((skyfreq > freqmin/1E9) and (skyfreq < freqmax/1E9)):
		print(" Creating line map: " + str(linname) + " - " + str(linfreq))
		os.system("rm -rf " + str(output_name) + "_" + str(linname))
		# Get channels
		chanmin = np.round(skyfreq-Nchans/2.*abs(refincre/1E9),6)
		chanmax = np.round(skyfreq+Nchans/2.*abs(refincre/1E9),6)
		chansexp = "range=[" + str(chanmin) + "GHz,"+ str(chanmax) +"GHz]"
		# Naming
		if (output_name == ""):
			output_name = myimage + "_" + str(linname) + "_Nchan" + str(Nchans) + "_dV" + str(np.round(abs(refincre/(skyfreq*1E9)*3E5),3)) + "kms-1_f0" + str(np.round(linfreq,5))+"GHz"
		else:
			output_file = output_name + "_" + str(linname) + "_Nchan" + str(Nchans) + "_dV" + str(np.round(abs(refincre/(skyfreq*1E9)*3E5),3))  + "kms-1_f0" + str(np.round(linfreq,5))+"GHz"
		# 
		imsubimage(imagename=myimage,outfile=output_file,chans=(chansexp),overwrite=True)
		# Set new Reference Freq.
		imhead(output_file,mode="put", hdkey="restfreq", hdvalue=str(linfreq*1E9)+"Hz")
		# Info
		print(" Filename: " + str(output_file))
	
	print(" ... DONE")
	print("==========================================================")

	
##-------------------------------------------------------
# --- Specific Mosaicing functions

def remove_masked_fits(myFITS_file):
	hdul = fits.open(myFITS_file, mode='update')
	# Access the masked data
	data = hdul[0].data
	mymask = np.isnan(data) 
	# Set the masked data to zero
	data[mymask] = 0.0
	# Save the modified FITS file
	hdul.flush()

def make_mask0_zeros(myimage):
	# myimage = image
	ia.open(myimage)
	ia.replacemaskedpixels(pixels='0.',mask=myimage+':mask0')
	ia.close()
	makemask(mode="delete",inpimage=myimage,inpmask=myimage+':mask0')

def remove_mask0(myimage):
	# myimage = image
	makemask(mode="delete",inpimage=myimage,inpmask=myimage+':mask0')


def do_mosaic(fits_files, fits_files_ref = 0, output_mosaic="mymosaic", workdir="./FITS_tmp/", do_baseorder_before = -1, do_baseorder_after = -1, do_edgechan = 30, do_sigweig = True, do_linweig = False, debugging_files = False):
	"""
	do_mosaic (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Create a TP mosaic from a list of FITS files (data products from the ASA)
	
	Arguments
	----------
	  fits_files : string
		List of FITS files to be combined into a single mosaic
	  fits_files_ref : integer (optional)
		(default = 0)
		File ID from fits_files used as template
	  output_mosaic : string
		(Default = "mymosaic")
		output mosaic name
	  workdir : str
		(Default = "./FITS_tmp/"
		working dir for creating the mosaic
	  do_baseorder_before : integer (optional)
		(default = -1 - no subtraction)
		Baseline order to be subtracted after combination 
	  do_baseorder_after : integer (optional)
		(default = -1 - no subtraction) 
		Baseline order to be subtracted after combination
	  do_edgechan : integer (optional)
		(default = 30; if <=0 then no removal)
		Remove edge channels 
	  do_sigweig : boolean
		Apply sigma weights (default = True) (PREFERRED)
		I(x,y) = sum(I_i(x,y) * w_i(x,y)) / sum(w_i(x,y)) where w_i(x,y) = 1/RMS  where w_i = 1/RMS
	  do_linweig : boolean
		Apply linear weights (default = False)
		I(x,y) = sum(I_i(x,y))/N(x,y) where N(x,y) = sum(i) = number of measurements
	  debugging_files : boolean (optional)
		(default = True)
		Delete intermediate products after combination for debugging?

	Outputs
	----------
	  Combine mosaic: output_mosaic + keywords depending on the options
	
	Keywords
	----------
	  sigmaweighting - I = sum(I_i * w_i) / sum(w_i) where w_i = 1/RMS
	  linearweighting - I = sum(I_i)/N where N = sum(i) = number of measurements
	  smo - 1/3 beam smoothing to homogeneize map
	  bl - baseline subtraction
	  clear - removing edge channels

	"""
	print("==============================================================")
	print(" do_mosaic(): Create mosaic")
	print("==============================================================")


	# Output mosaic name
	os.system('rm -rf ' + str(workdir) + output_mosaic + '*')

	## --------------------------------------
	print("-----------------------------------------------------")
	print(" File information")

	list_image_edges = np.zeros((len(fits_files),4))
	list_phasecenter = np.zeros((len(fits_files),2))
	list_beam = np.zeros((len(fits_files),3))
	list_freq0 = np.zeros((len(fits_files),1))
	list_freqdelta = np.zeros((len(fits_files),1))
	list_Nchan = np.zeros((len(fits_files),1))
	consistency = True

	# Retrieving information
	for num in np.arange(0,len(fits_files)):
		hdul = fits.open(fits_files[num])
		header = hdul[0].header
		# Beams
		list_beam[num,0] = header['BMAJ']*3600
		list_beam[num,1] = header['BMIN']*3600
		list_beam[num,2] = header['BPA']
		# Phase Center
		list_phasecenter[num,0] = header['CRVAL1']
		list_phasecenter[num,1] = header['CRVAL2']
		# Map edges
		list_image_edges[num,0] = header['CRVAL1']+(header['NAXIS1']-header['CRPIX1'])*abs(header['CDELT1']) # highest RA
		list_image_edges[num,1] = header['CRVAL1']+(-header['CRPIX1'])*abs(header['CDELT1'])		# lowest RA
		list_image_edges[num,2] = header['CRVAL2']+(header['NAXIS2']-header['CRPIX2'])*abs(header['CDELT2'])	# highest DE
		list_image_edges[num,3] = header['CRVAL2']+(-header['CRPIX2'])*abs(header['CDELT1'])		# lowest DE
		#
		if (num==0):
			list_axis = [[header['CTYPE1'],header['CTYPE2'],header['CTYPE3'],header['CTYPE4']]]
		else:
			list_axis.append([header['CTYPE1'],header['CTYPE2'],header['CTYPE3'],header['CTYPE4']])
		list_freq0[num,0] = header['RESTFRQ']/1E9
		list_freqdelta[num,0] = header['CDELT3']/1E9
		list_Nchan[num,0] = header['NAXIS3']


	# Comaparison with reference
	for num in np.arange(0,len(fits_files)):
		#
		print("------------------------------------")
		if (num == fits_files_ref):
			print(" ---> Reference image")
		print(" Tile ID: "+str(num))
		print(" FITS: " + fits_files[num])
		print(" axis = " + str(list_axis[num]))
		print(" freq0 = " + str(np.round(list_freq0[num,0],6)) + " GHz")
		print(" freq_delta = " + str(np.round(list_freqdelta[num,0]*1E6,3)) + " kHz")
		if ((list_freq0[num,0]-list_freq0[fits_files_ref,0])/list_freqdelta[fits_files_ref,0] > 1.0):
			print("Tile ID " + str(num) + " ---> WARNING: Different central frequency (>1 channel) than referecence")
			#consistency = False
		if (list_freqdelta[num,0] < list_freqdelta[fits_files_ref,0]):
			print("Tile ID " + str(num) + " ---> WARNING: Spectral resolution smaller than referecence")
		print(" Nchan = " + str(list_Nchan[num,0]) + " channels")
		if (list_Nchan[num,0] != list_Nchan[fits_files_ref,0]):
			print("Tile ID " + str(num) + " ---> WARNING: different number of channels than referecence")
			#consistency = False
		print(" Beam = (" + str(np.round(list_beam[num,0],2)) + "," + str(np.round(list_beam[num,1],2)) + ") arcsec")
		print(" Phase center = (" + str(list_phasecenter[num,0]) + "," + str(list_phasecenter[num,1]) + ")")

	## --------------------------------------
	# Check consistency
	print("------------------------------------")
	print(" Checking mosaic consistency")
	if (consistency == True):
		print(" Files are consistent. Continuing.")
	else:
		print(" Consistency issues: please check above")
		sys.exit()


	## --------------------------------------
	# Mosaic definitions
	print("------------------------------------")
	print(" Defining Final Mosaic parameters")
	# Defining final cube dimentions and beam
	beam_final = np.max(list_beam[:,0:1])
	beam_diff = (np.max(list_beam[:,0:1])-np.min(list_beam[:,0:1]))/beam_final
	print(" Final beam = (" + str(np.round(beam_final,2)) + ") arcsec")

	if (beam_diff/beam_final > 0.1):
		print(" ---> WARNING: large beam differences between (>10%) images. Check image list")

	# Definining final mosaic parameters
	max_RA = np.max(list_image_edges[:,0])
	min_RA = np.min(list_image_edges[:,1])
	max_DE = np.max(list_image_edges[:,2])
	min_DE = np.min(list_image_edges[:,3])
	new_phasecenter_RA = (max_RA+min_RA)/2.
	new_phasecenter_DE = (max_DE+min_DE)/2.
	print(" Mosaic Phase Center = (" + str(new_phasecenter_RA) + "," + str(new_phasecenter_DE) + ")")

	## Frequency
	freq0_REF = list_freq0[fits_files_ref,0]
	freqdelta_REF = list_freqdelta[fits_files_ref,0]

	## --------------------------------------
	## --- Reference image
	print("-----------------------------------------------------")
	print("--- Creating Mosaic Template")
	ref_fits = fits_files[fits_files_ref]

	# Import FITS
	os.system('rm -rf' + str(workdir) + ' ref_image*')
	os.system('cp ' + str(ref_fits) + ' ' + str(workdir) + 'ref_fits_tmp')
	
	# Modify header to make reference image larger
	hdul = fits.open(str(workdir) + 'ref_fits_tmp')
	data = hdul[0].data*0.
	header = hdul[0].header
	orig_naxis1 = header['NAXIS1']
	orig_naxis2 = header['NAXIS2']
	orig_pixelsize1 = header['CDELT1']
	orig_pixelsize2 = header['CDELT2']
	

	# New diamentions
	oversampling = 1.	# Oversampling factor in the template image; helps to get better transitions
	new_pixelsize1 = orig_pixelsize1/oversampling
	new_pixelsize2 = orig_pixelsize2/oversampling
	new_naxis1 = int(1.1*(max_RA-min_RA)/abs(new_pixelsize1))	# 110% of the max-min size in RA
	new_naxis2 = int(1.1*(max_DE-min_DE)/abs(new_pixelsize2))	# 110% of the max-min size in DEC
	resized_data = np.zeros((header['NAXIS4'], header['NAXIS3'], new_naxis2, new_naxis1 ))
	header['CRVAL1'] = new_phasecenter_RA
	header['CRVAL2'] = new_phasecenter_DE
	header['NAXIS1'] = new_naxis1
	header['NAXIS2'] = new_naxis2
	header['CRPIX1'] = new_naxis1/2.
	header['CRPIX2'] = new_naxis2/2.
	header['CDELT1'] = new_pixelsize1
	header['CDELT2'] = new_pixelsize2
	header['BMAJ'] = beam_final/3600.
	header['BMIN'] = beam_final/3600.
	header['BPA'] = 0.0

	new_hdu = fits.PrimaryHDU(resized_data, header)
	os.system('rm -rf ' + str(workdir) + 'ref_fits')
	new_hdu.writeto(str(workdir) + 'ref_fits')

	# Creatre cubes of zeros as reference
	importfits(fitsimage= str(workdir) + 'ref_fits', imagename=str(workdir) + 'ref_image')
	#imhead('ref_image',mode="put", hdkey="restfreq", hdvalue=str(freq0_REF*1E9)+"Hz")

	## --------------------------------------
	## --- Regrid tiles and calculate weights
	print("-----------------------------------------------------")
	print(" Regridding individual tiles into template")
	os.system('rm -rf ' + str(workdir) + 'tile*')
	os.system('rm -rf ' + str(workdir) + '*_regrid*')
	os.system('rm -rf ' + str(workdir) + 'tmp*')

	for num in np.arange(0,len(fits_files)):
		print("-------------------------")
		print("--- Tile"+str(num))
		if (num == fits_files_ref):
			print(" ---> Reference field")
		os.system('cp ' + fits_files[num] + ' ' + str(workdir) + 'fits_tmp')
		os.system('rm -rf ' + str(workdir) + 'tmp*')
		importfits(fitsimage= str(workdir) + 'fits_tmp', imagename= str(workdir) +'tmp_image',overwrite=True)

		# Generate imputs for immath (combination)
		if (num==0):
			list_images = [str(workdir) +'tile'+ str(num) +'_regrid']
			list_images_weights = [str(workdir) +'tile'+ str(num) +'_regrid_ones']
			list_images_math = '(IM'+str(num)
			list_images_sigmas = [str(workdir) +'tile'+ str(num) +'_regrid_sigmas']
		else:
			list_images.append(str(workdir) +'tile'+ str(num) +'_regrid')	# list of files for immath
			list_images_weights.append(str(workdir) +'tile'+ str(num) +'_regrid_ones')	# list of files for immath
			list_images_math = list_images_math+"+IM"+str(num)	# expression for immath
			list_images_sigmas.append(str(workdir) +'tile'+ str(num) +'_regrid_sigmas')	# list of files for immathh
		
		# smooth 
		imres = imhead(str(workdir) +'tmp_image',hdkey='CDELT3',mode='get')["value"]
		smofactor = int(np.round(freqdelta_REF*1E9/imres))
		if (smofactor > 1):	
			print("--- Image freq.res. different than REF: applying smoothing factor = " + str(smofactor))
			os.system("mv " + str(workdir) + "tmp_image " + str(workdir) + "tmp_image1")
			specsmooth(imagename=str(workdir) + "tmp_image1", outfile= str(workdir) + "tmp_image", width=smofactor, axis=2, overwrite=True)
			os.system("rm -rf " + str(workdir) + "tmp_image1")	

		# recenter frequency if necessary
		imfreq = imhead(str(workdir) + 'tmp_image',hdkey='restfreq',mode='get')["value"]
		print("---> Image Freq0 = " + str(imfreq/1E9) + " GHz")
		if (imfreq != (freq0_REF*1E9)):
			print("--- Image freq. different than REF: Recentering")
			# Set new Reference Freq.
			os.system("mv " + str(workdir) + "tmp_image " + str(workdir) + "tmp_image1")
			imhead(str(workdir) + 'tmp_image1',mode="put", hdkey="restfreq", hdvalue=str(freq0_REF*1E9)+"Hz")
			imfreq = imhead(str(workdir) + 'tmp_image1',hdkey='restfreq',mode='get')["value"]
			os.system("cp -r " + str(workdir) + "tmp_image1 " + str(workdir) + "tmp_image")
			print("---> New Freq0 = " + str(imfreq/1E9) + " GHz")
			#os.system("cp -r tmp_image1 kk.image")

		#os.system("cp -r tmp_image kk" + str(num) +".image")		

		# Convolution into final beam (if needed)
		print("--- Smoothing image tile"+str(num))
		imbeam = imhead(str(workdir) +'tmp_image',hdkey='beammajor',mode='get')
		if (abs(1.-imbeam['value']/beam_final) > 0.05):
			print("--- Convolving image into final resolution/beam")
			imsmooth(imagename= str(workdir) +'tmp_image', outfile= str(workdir) +'tmp_image2',beam= {"major": str(beam_final) + "arcsec", "minor": str(beam_final) + "arcsec", "pa": "0deg"},targetres=True)
		else:
			print("--- Beam difference < 5%. No convolution needed")
			os.system("mv " + str(workdir) + "tmp_image " + str(workdir) + "tmp_image2")


		# Regrid into the ref_image frame
		print("--- Reprojecting image tile"+str(num))
		imregrid(imagename= str(workdir) +'tmp_image2', template= str(workdir) +'ref_image', output= str(workdir) +'tile'+ str(num) +'_regrid',interpolation='nearest')
		os.system('rm -rf ' + str(workdir) + 'tmp_image*')  # Clean up temporary files

		# Baseline subtraction (usually = 0 to remove platforming)
		if (do_baseorder_before != -1):
			print(" --- Baseline substraction")
			os.system('mv ' + str(workdir) +'tile'+ str(num) +'_regrid ' + str(workdir) + 'tmp_image3')
			imshape = imhead(str(workdir) +'tmp_image3',mode='get',hdkey='shape')
			Nchan = imshape[2]
			Nchan_13 = int(imshape[2]/3.)
			Nchan_23 = int(imshape[2]*2/3.)
			imcontsub(imagename= str(workdir) +'tmp_image3', linefile= str(workdir) +'tile'+ str(num) +'_regrid', chans='20~'+str(Nchan_13)+', '+str(Nchan_23)+'~'+str(Nchan-20), fitorder=do_baseorder_before, stokes="I") 
			os.system('rm -rf ' + str(workdir) + 'tmp_image*')

		# Get an image of 1's where there is emission
		immath(imagename= str(workdir) + 'tile'+ str(num) +'_regrid', expr='IM0/IM0', outfile= str(workdir) + 'tile'+ str(num) +'_regrid_ones')

		# Make all masked voxels = 0 and/or remove mask
		print("--- Remove masked data or convert them into zeros")
		# For whatever reason mask only works in the same folder, not with paths
		# Workaround: move the file to ., then mask, then move back
		os.system("mv " + str(workdir) +'tile'+ str(num) +'_regrid_ones .')
		make_mask0_zeros('tile'+ str(num) +'_regrid_ones')
		os.system('mv tile'+ str(num) +'_regrid_ones ' + str(workdir) +'tile'+ str(num) +'_regrid_ones')
		os.system("mv " + str(workdir) +'tile'+ str(num) +'_regrid .')
		remove_mask0('tile'+ str(num) +'_regrid')
		os.system('mv tile'+ str(num) +'_regrid ' + str(workdir) +'tile'+ str(num) +'_regrid')

		# Weigths
		print("--- Calculating weights for tile"+str(num))
		# Calculate RMS
		imshape = imhead(str(workdir) +'tile'+ str(num) +'_regrid',mode='get',hdkey='shape')
		Nchan = imshape[2]
		Nchan_13 = int(imshape[2]/3.)
		Nchan_23 = int(imshape[2]*2/3.)
		tile_sigma = imstat(str(workdir) +'tile'+ str(num) +'_regrid',chans='10~'+str(Nchan_13)+', '+str(Nchan_23)+'~'+str(Nchan-10))["rms"][0]
		imunits = imhead(str(workdir) +'tile'+ str(num) +'_regrid',mode='get',hdkey='bunit')
		print("---> tile RMS = " + str(tile_sigma) + " " + str(imunits))
		# Create RMs and RMS+Signal images
		os.system('rm -rf ' + str(workdir) + 'tile' + str(num) +'*sigma*')
		immath(imagename=[str(workdir) + 'tile'+ str(num) +'_regrid'], expr='(IM0)*'+str(tile_sigma), outfile= str(workdir) +'tile'+ str(num) +'_regrid_sigmaXsignal')
		immath(imagename=[str(workdir) +'tile'+ str(num) +'_regrid_ones'], expr='(IM0)*'+str(tile_sigma), outfile= str(workdir) +'tile'+ str(num) +'_regrid_sigmas')


		
		# Make all masked voxels = 0 in all images
		#print("--- Remove masked data or convert them into zeros")
		#remove_mask0('tile'+ str(num) +'_regrid')
		#remove_mask0('tile'+ str(num) +'_regrid_sigmaXsignal')
		#make_mask0_zeros('tile'+ str(num) +'_regrid_ones')
		#immath(imagename=['tile'+ str(num) +'_regrid_ones'], expr='(IM0)*'+str(tile_sigma), outfile='tile'+ str(num) +'_regrid_sigmas')
		#make_mask0_zeros('tile'+ str(num) +'_regrid_sigmas')




	# closing expression for immath
	list_images_math = list_images_math+")"	


	## --------------------------------------
	## --- Combine tiles into large Mosaic
	addfilename = ""
	os.system('rm -rf ' + str(workdir) + output_mosaic + '_linearweighting'+str(addfilename)+'_*')
	os.system('rm -rf ' + str(workdir) + output_mosaic + '_sigmaweighting'+str(addfilename)+'_*')

	print(" -------------------------------------------------------- ")
	print(" Combining individual tiles into a large mosaic")

	# a) Weights given RMS: value = sum(I_i*sigma_i)/sum(sigma_i) = DEFAULT
	if (do_sigweig):
		print("--- Creating mosaic using Weighted Combination = " + str(workdir) + output_mosaic + '_sigmaweighting')
		os.system('rm -rf ' + str(workdir) + output_mosaic + '_sigmatot')
		immath(imagename=list_images_sigmas, expr=list_images_math, outfile= str(workdir) + output_mosaic + '_sigmatot')

		# Combined using weighted mean (i.e. I_comb = sum(I_i*w_i)/sum(w_i))
		os.system('rm -rf ' + str(workdir) + output_mosaic + '_sigmaweighting')
		newlist = [s + "_sigmaXsignal" for s in list_images]
		newlist.append( str(workdir) + output_mosaic + '_sigmatot')
		math_expr = list_images_math+"/IM"+ str(num+1)
		immath(imagename=newlist, expr=math_expr, outfile=str(workdir) + output_mosaic + '_sigmaweighting')

		# Remove Nan's and Infs
		os.system("mv " + str(workdir) +output_mosaic + '_sigmaweighting ./tmp_sigmaweighting')
		ia.open('tmp_sigmaweighting')
		ia.calcmask('tmp_sigmaweighting < 10000.',name="masknan")
		#ia.replacemaskedpixels('0.0',mask=output_mosaic + '_sigmaweighting:mask12')
		ia.close()
		os.system("mv tmp_sigmaweighting " + str(workdir) +output_mosaic + '_sigmaweighting')

	# b) Linear weights
	if (do_linweig):
		print("--- Creating mosaic using Linear Combination = " + str(workdir) + output_mosaic + '_linearweighting')
		os.system('rm -rf ' + str(workdir) + output_mosaic + '_weights')
		immath(imagename=list_images_weights, expr=list_images_math, outfile= str(workdir) + output_mosaic + '_weights')

		# Combine using normal/linear mean (i.e. I_comb = sum(I_i))
		os.system('rm -rf ' + str(workdir) + output_mosaic + '_linearweighting')
		newlist = [s for s in list_images]
		newlist.append(str(workdir) + output_mosaic + '_weights')
		math_expr = list_images_math+"/IM"+ str(num+1)
		immath(imagename=newlist, expr=math_expr, outfile= str(workdir) + output_mosaic + '_linearweighting')
		
		# Remove Nan's and Infs
		os.system("mv " + str(workdir) +output_mosaic + '_linearweighting ./tmp_linearweighting ')
		ia.open('tmp_linearweighting')
		ia.calcmask('tmp_linearweighting < 10000.',name="masknan")
		#ia.replacemaskedpixels('0.0',mask=output_mosaic + '_linearweighting:mask12')
		ia.close()
		os.system("mv tmp_linearweighting " + str(workdir) +output_mosaic + '_linearweighting')



	# --- Smooth out final image to improve overlap
	print(" -------------------------------------------------------- ")
	print(" Small convolution by 1/3 beam")
	if (do_linweig):
		imsmooth(imagename= str(workdir) +output_mosaic + '_linearweighting'+str(addfilename),
			outfile= str(workdir) +output_mosaic + '_linearweighting'+str(addfilename)+'_smo',
			beam= {"major": str(beam_final/3.) + "arcsec", "minor": str(beam_final/3.) + "arcsec",
			"pa": "0deg"},targetres=False)
	if (do_sigweig):
		imsmooth(imagename= str(workdir) + output_mosaic + '_sigmaweighting'+str(addfilename),
		outfile= str(workdir) +output_mosaic + '_sigmaweighting'+str(addfilename)+'_smo',
		beam= {"major": str(beam_final/3.) + "arcsec", "minor": str(beam_final/3.) + "arcsec", 
		"pa": "0deg"},targetres=False)
	#add name
	addfilename = addfilename+'_smo'


	print(" -------------------------------------------------------- ")
	print(" Baseline substraction")
	if (do_baseorder_after != -1):
		print("--- Applying baseline subtraction of order = " + str(do_baseorder))
		# get number of channels
		imshape = imhead(str(workdir) + output_mosaic + '_linearweighting'+str(addfilename),mode='get',hdkey='shape')
		Nchan = imshape[2]
		Nchan_13 = int(imshape[2]/3.)
		Nchan_23 = int(imshape[2]*2/3.)
		if (do_linweig):
			imcontsub(imagename= str(workdir) + output_mosaic + '_linearweighting'+str(addfilename), 
				linefile= str(workdir) +output_mosaic + '_linearweighting'+str(addfilename)+'_b'+str(do_baseorder), 
				chans='50~'+str(Nchan_13)+', '+str(Nchan_23)+'~'+str(Nchan-50), fitorder=do_baseorder_after, stokes="I") 
		if (do_sigweig):
			imcontsub(imagename= str(workdir) + output_mosaic + '_sigmaweighting'+str(addfilename), 
				linefile=str(workdir) + output_mosaic + '_sigmaweighting'+str(addfilename)+'_b'+str(do_baseorder),
				chans='50~'+str(Nchan_13)+', '+str(Nchan_23)+'~'+str(Nchan-50), fitorder=do_baseorder, stokes="I") 
		# Add name
		addfilename = addfilename+"_b"+str(do_baseorder_after)
	else:
		print("--- No baseline subtraction applied")


	print(" -------------------------------------------------------- ")
	print(" Rebinning into Nyquist sampling")
	# Rebinning into Nyquist
	oversampling = beam_final/(abs(new_pixelsize1)*3600.)	# Calculate sampling (=beam/pixel) ratio
	rebin_factor = int(np.floor(oversampling/2.0))	# At least a factor of 4 pixels per beam
	if (rebin_factor >= 2):
		print("--- Original sampling >> Nyquist")
		print("--- Applying rebin factor = " + str(rebin_factor))
		if (do_linweig):
			imrebin(imagename= str(workdir) +output_mosaic + '_linearweighting'+str(addfilename),factor=[rebin_factor,rebin_factor,1,1],outfile= str(workdir) +output_mosaic + '_linearweighting'+str(addfilename)+'_rebin')
		if (do_sigweig):
			imrebin(imagename= str(workdir) +output_mosaic + '_sigmaweighting'+str(addfilename),factor=[rebin_factor,rebin_factor,1,1],outfile= str(workdir) + output_mosaic + '_sigmaweighting'+str(addfilename)+'_rebin')
		addfilename = addfilename+"_rebin"
	else:
		print("--- Original sampling ~ Nyquist")
		print("--- No rebining applied")


	print(" -------------------------------------------------------- ")
	print(" Removing channels at the end of the BW")
	if (do_edgechan > 0):
		if (do_linweig):
			imsubimage(imagename= str(workdir) + output_mosaic + '_linearweighting'+str(addfilename),outfile=str(workdir) + output_mosaic + '_linearweighting'+str(addfilename)+"_clean",chans=str(do_edgechan)+'~'+str(Nchan-do_edgechan))
		if (do_sigweig):
			imsubimage(imagename= str(workdir) + output_mosaic + '_sigmaweighting'+str(addfilename),outfile= str(workdir) + output_mosaic + '_sigmaweighting'+str(addfilename)+"_clean",chans=str(do_edgechan)+'~'+str(Nchan-do_edgechan))
		addfilename = addfilename+"_clean"


	# --- Cleaning up
	print(" -------------------------------------------------------- ")
	print(" Cleaning files")
	if(debugging_files):
		print(" Debugging mode ON = all intermediate products are available")
	else:
		print(" Cleaning up tmp & intermediate files/images")
		# Remove individual tiles
		os.system("rm -rf " + str(workdir) + "tmp*")
		os.system("rm -rf " + str(workdir) + "tile*")
		os.system("rm -rf " + str(workdir) + "ref_image*")
		# Remove intermediate products
		if (do_linweig):
			os.system("mv " + str(workdir) + output_mosaic + "_linearweighting"+str(addfilename) + " " + str(workdir) + "tmp_final1.image")
		if (do_sigweig):
			os.system("mv " + str(workdir) + output_mosaic + "_sigmaweighting"+str(addfilename) + " " + str(workdir) + "tmp_final2.image")
		os.system("rm -rf " + str(workdir) + output_mosaic + '_*')
		os.system("rm -rf " + str(workdir) + output_mosaic + '_*')
		if (do_linweig):
			os.system("mv " + str(workdir) + "tmp_final1.image " + str(workdir) + output_mosaic + "_linearweighting"+str(addfilename))
		if (do_sigweig):
			os.system("mv " + str(workdir) + "tmp_final2.image " + str(workdir) + output_mosaic + "_sigmaweighting"+str(addfilename))
		

	# Print a message indicating the completion of the task
	print(" -------------------------------------------------------- ")
	print(" FITS files combined successfully ")
	if (do_linweig):
		outfilename = str(workdir) + output_mosaic + '_linearweighting'+addfilename
		print(" Linear Mosaic = {}".format(outfilename))
		imbeam = imhead(outfilename,hdkey='beammajor',mode='get')
		print(" Final mosaic resolution = " + str(imbeam) + "arcsec")
		
	if (do_sigweig):
		outfilename = str(workdir) + output_mosaic + '_sigmaweighting'+addfilename
		print(" Weighted Mosaic = {}".format(outfilename))
		imbeam = imhead(outfilename,hdkey='beammajor',mode='get')
		
	print(" Final mosaic resolution = " + str(imbeam) + "arcsec")
	return outfilename

	print(" -------------------------------------------------------- ")
	print(" ... DONE")
	print("==========================================================")


def do_fullmosaic(myworkdir="./FITS_raw/",mosaic_name="FULLmosaic",coordsys="ABS"):
	"""
	full_mosaic (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Wrapper of do_mosaic
	
	Arguments
	"""
	# 
	os.system("rm -rf " +str(mosaic_name))
	# Collect all FITS files
	FITSfiles = listofFITS(myworkdir,listexpr=".fits")
	print("here" + str(FITSfiles))

	if (((coordsys == "ABS") | (coordsys == "GAL")) & (len(FITSfiles)>0)):
		# If conversion to GALACTIC coordintes
		if (coordsys=="GAL"):
			for ima in FITSfiles:
				importfits(fitsimage= str(ima),imagename=str(ima)+".image",overwrite=True)
				imregrid(imagename=str(ima)+".image",output=str(ima)+"_GAL.image",template='GALACTIC')
				exportfits(imagename=str(ima)+"_GAL.image",fitsimage=str(ima)+"_GAL")
				os.system("rm -rf "+ str(myworkdir) +str(ima)+"_GAL.image")
			# update list of files
			FITSfiles = listofFITS(myworkdir,listexpr="_GAL")
			print("here2" + str(FITSfiles))

		# Run mosaic
		Finalimage = do_mosaic(FITSfiles, output_mosaic = mosaic_name,workdir=myworkdir)
		print("here" + str(Finalimage))

		# cleaning up
		os.system("rm -rf *tmp ref* *last *_GAL *fits.image")

	#return
	return Finalimage


#######################################################
# Section 3: Pipeline
#######################################################

##----------------------------------------------------
## Step 0: Fetch data


def do_fetchdata(mycoords,R0,myline,mindV=0.0E3,maxdV=1.0E3,VLSR_ms=-1E6,datadir="./archive/",workdir="./",do_copy=False,do_IntOnly=False,MSext=".cal"):
	'''

	do_fetchdata (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Prepares data (visibilities + FITS files) for pipeline work
	The function will look for a line and will try to in

	Arguments
	----------
	  mycoords : list or str 
		Target coordinates [RA,DE] if in degrees
		e.g., [REdeg,DEdeg] = [4.232545,-32.9452343]
		     target coordinates ['RADEC string'] if in absolute coordinates
		      [RAabs,DEabs]] = ['12:56:11.16658 -05.47.21.5246']
	  R0 : float
		(Default = 0.5 deg)
		Seach radius in deg
	  myline : str or float
		If "str" the line ID will be identified in emerge.linecat (e.g. "N2H+10").
		If else "float" the line will be taken as frequency in Hz
		Note that do_fetchdata will try to infer the VLSR from MS files and apply it
		in order to calculte the skyfreq
	  dVmin : float (optional) 
		(Default = 0.0E3 m/s)
		minimum velocity resolution [km/s] 
	  dVmax : float (optional) 
		(Default = 1.0E3 m/s)
		maximum velocity resolution [km/s] 
	  VLSR_ms : float (optional)
		(Default = -1E6 m/s)
		If VLSR_ms Default, then it will look for source VLSR inside MS files
	  datadir : str
		(Default = ".")
		path to FITS data 
	  workdir : str
		(Default = ".")
		working and writing directory tree. See also gen_foldertree().
	  do_copy : boolean
		(Default = False)
		Copy data from datadir to workdir
	  do_IntOnly : boolean
		(Default = False)
		If False it will only consider datasets that contain short-spacing information.
		If True it will continue even if there is no short-spacing information

	Outputs
	----------
	  Copy data into corresponding folders (see EMERGE folder structure)

	'''
	print("==================================================================")
	print(" Starting do_fetchdata()...")

	#---------------------------------
	# Environment

	# Generate environment if needed
	#gen_foldertree(workdir=workdir)

	#---------------------------------
	# Frequency
	# If a string try to identify it in emerge.linecat
	if (type(myline) == str):
		myrestfreq_Hz = emerge_linecat[myline]*1E9		# in Hz
	# If numeric, it will attempt the given frequency
	if (type(myline) == float):
		myrestfreq_Hz = myline

	print(" Line frequency = " + str(myrestfreq_Hz/1E9) + " GHz")

	#---------------------------------
	# Short spacing data
	SDdata = 0
	TPdata = 0
	#ShortSpaData = SDdata+TPdata

	# Attempt to get velocity information from visibility files
	# Not the most elegant solution but we need some VLSR information for the FITS files selection
	allMS = listofMSfiles(datadir=datadir+"visibilities/")#,listexpr=MSext)
	VLSRfromMS = False

	# If VLRS_ms != 0.0 --> source VLSR provided by user 
	if (VLSR_ms != -1E6):
		print(" Using user VLSR = " + str(VLSR_ms) + " m/s to calculate skyfreq")
		myfreq = getskyfreq(myrestfreq_Hz,myVLSR=VLSR_ms)

	# If VLSR_ms == Default, then try to get it from MS information
	if (VLSR_ms == -1E6):
		# Explore MS files to look for a source velocity (only once)
		for myMS in allMS:
			# Print file name
			print(" Visibility file: " + str(datadir) + "visibilities/" +str(myMS))
			# Check if the MS file contains the right coordinates and frequency
			if (not(VLSRfromMS)):
				isinMS, mytargetID = istargetinMS(datadir + "visibilities/" + str(myMS),mycoords=mycoords,R0=R0,units='deg')
				if (isinMS):
					VLSRfromMS = True
					tmp = getMSinfo(datadir + "visibilities/" + str(myMS))
					print(" Extracting VLSR information from MS file : " + datadir + "visibilities/" + str(myMS))
					# If multiple sources satify criteria, take the first one
					if (isinstance(mytargetID,np.ndarray)):
						mytargetID = mytargetID[0]
					# For whatever reason sometimes this keyword is a float, sometimes a list...
					if (type(tmp[mytargetID]['Source_Vel_ms-1']) == float):
						# If float
						myVLSRfromMS = tmp[mytargetID]['Source_Vel_ms-1']
					else:
						# If list, take the first element
						myVLSRfromMS = tmp[mytargetID]['Source_Vel_ms-1'][0]
					print(" Using VLSR = " + str(myVLSRfromMS) + " m/s to calculate skyfreq")
					myfreq = getskyfreq(myrestfreq_Hz,myVLSRfromMS)
		if (not(VLSRfromMS)):
			print(" No VLSR found in .MS files, using VLSR = 0.0 m/s to calculate skyfreq instead")
			print(" WARNING: Note that in this case it is unlikely that we will find any .MS file that satisfies the searching criteria")
			myfreq = getskyfreq(myrestfreq_Hz,myVLSR=0.0)


	# Getting short-spacing data	
	# (1) Check if SD data exist, assumed to be from a big SD telescope (preferred)
	if (os.path.isdir(str(datadir)+"shortspacing_SD")):
		SDdata = select_FITS(mycoords,R0=R0,freq0=myfreq,dVmin=mindV,dVmax=maxdV,datadir=str(datadir)+"shortspacing_SD",workdir=workdir+"/FITS_raw",do_copy=do_copy)

	# (2) If no SD data, then try local ALMA-TP data
	if ((SDdata == 0) & (os.path.isdir(str(datadir)+"shortspacing_TP"))):
		# Check TP data exist
		TPdata = select_FITS(mycoords,R0=R0,freq0=myfreq,dVmin=mindV,dVmax=maxdV,datadir=str(datadir)+"shortspacing_TP",workdir=workdir+"/FITS_raw",do_copy=do_copy)

	# (3) If no local data, then Try to download data from JVO database
	if ((SDdata == 0) & (TPdata == 0)):
		try:
			import alminer
			allGOUS = find_GOUS(targetcoord=mycoords,targetline=myline,targetfreq_Hz=myfreq,maxdV=maxdV,mindV=mindV,R0=R0)
			TPdata = getTPfiles(allGOUS=allGOUS,skyfreq_Hz=myfreq,maxdV=maxdV,mindV=mindV); print(" Total of FITS files found =  " + str(TPdata))
		except:
			print("No Archival data found")
	# Are there short spacing data?
	if ((SDdata == 0) & (TPdata == 0) & (not(do_IntOnly))):
		# If no short-spacing information available, then exit
		print(" Neither SD nor ALMA-TP data available. No short-spacing information available. Exit.")
		print(" Note: If wanted to work with only Interferometric data, set do_IntOnly=True")
		sys.exit()
	if ((SDdata == 0) & (TPdata == 0) & (do_IntOnly)):
		print(" No short-spacing information available. Continuing with only interferometric data.")
	
	# Final number of short-spacing datasets
	ShortSpaData = SDdata+TPdata

	#---------------------------------
	# Interferometric data
	print(" Looking for MS files in archive...")
	Intdata = 0

	# Find all .cal MS-files in emerge_archive
	allMS = listofMSfiles(datadir=datadir+"visibilities/",listexpr=".ms")
	contsubms = listofMSfiles(datadir=datadir+"visibilities/",listexpr="_target.ms")
	selfcalms = listofMSfiles(datadir=datadir+"visibilities/",listexpr="_targets.ms")


	# Copy files from Archive to local vis_raw folder
	#
	# (README from calMS)
	# The most usual MS names:
	# *.ms.split.cal: the complete calibrated data for the EB in the DATA column.
	# *_target.ms: the result of the continuum subtraction stage of the imaging pipeline (Cycle 8)
	#	      for the science target(s) only. The CORRECTED_DATA column contains the
	#	      continuum-subtracted visibilities ready for line imaging, while the DATA column
	#	      contains the visibilities with continuum as in *.ms.split.cal.
	# *_targets.ms: like *.ms.split.cal but for the science target(s) only (Cycle >8).
	#	       From Cycle 10 onwards, these MSs can contain a CORRECTED_DATA column with the
	#	       self-calibrated visibilities if self-calibration was successful for the given target.
	# *_targets_line.ms: the result of the continuum subtraction stage of the imaging pipeline (Cycle >8)
	#	            for the science target(s) only: The DATA column contains the
	#	            continuum-subtracted visibilities ready for line imaging.
	#	            From Cycle 10 onwards, these MSs can contain a CORRECTED_DATA column with the
	#	            self-calibrated visibilities if self-calibration was successful for the given target.
	#
	for myMS in allMS:
		# Check if the MS file contains the right coordinates and frequency
		isinMS, mytargetID, mylineSPW , mychanresSPW = istarANDlineinMS(datadir + "visibilities/" + str(myMS),myrestfreq_Hz,mycoords,R0=R0,units='deg')

		#if so, then copy the MS file into ./vis_raw folder
		if (isinMS):
			print("=================================================")
			print(" Copying " + str(myMS) + " into "+ str(workdir) + "vis_raw folder")
			if (do_copy & (myMS.endswith("_targets_line.ms") == False)):
				if (myMS.endswith(".ms") | myMS.endswith(".ms.split.cal")) and not (myMS.endswith("_target.ms") | myMS.endswith("_targets.ms")):	# aka contsub only
					do_copyMSfromArchive(datadir=datadir+"visibilities/",
					listexpr=myMS,workdir=workdir+"vis_raw/")
					Intdata += 1
				if (myMS.endswith("_target.ms") | myMS.endswith("_targets.ms")):
					newMS = myMS[:myMS.index("_target")]
					split(vis= datadir + "visibilities/" + str(myMS), 
					outputvis = "./vis_raw/" + str(newMS) + ".ms", datacolumn="data")
					Intdata += 1
				# We ignore *_targets_line.ms files
			print(" ... DONE")
			print("=================================================")

	#---------------------------------
	# Return results: number of datasets for SD, TP, and Int data
	return ShortSpaData, Intdata

	print(" do_fetchdata() ... DONE")
	print("==================================================================")

##----------------------------------------------------
## Step 1: Generate .MS file

def do_extractvis(mycoords,R0,myline,mindV,maxdV,datadir="./vis_raw/",workdir="./",linetol_ms=15.E3,debugging=False,MSext=".ms",doassessMS=True,path2assessMS="./"):
	'''

	do_extractvis (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Extracts relevant visibilities given position, freq., and spectral criteria

	Arguments
	----------
	  mycoords : list or str 
		Target coordinates [RA,DE] if in degrees
		e.g., [REdeg,DEdeg] = [4.232545,-32.9452343]
		     target coordinates ['RADEC string'] if in absolute coordinates
		      [RAabs,DEabs]] = ['12:56:11.16658 -05.47.21.5246']
	  R0 : float
		(Default = 0.5 deg)
		Seach radius in deg
	  myline : str or float
		If "str" the line ID will be identified in emerge.linecat (e.g. "N2H+10").
		If else "float" the line will be taken as frequency in Hz
	  dVmin : float (optional) 
		(Default = 0.0E3 m/s)
		minimum velocity resolution [km/s] 
	  dVmax : float (optional) 
		(Default = 1.0E3 m/s)
		maximum velocity resolution [km/s] 
	  datadir : str
		(Default = ".")
		path to FITS data 
	  workdir : str
		(Default = ".")
		working and writing directory tree. See also gen_foldertree().
	  do_copy : boolean
		(Default = False)
		Copy data from datadir to workdir
	  do_IntOnly : boolean
		(Default = False)
		If False it will only consider datasets that contain short-spacing information.
		If True it will continue even if there is no short-spacing information
	  linetol : float
		(Default = 15E3 m/s) 
		Window tolerance around the lines to be masked for contsub.
	  debugging : boolean
		(Default = False)
		If False, then delete all intermediate files (mostly visibilities)
		If True, keep them for code debugging
	  MSext : str
		(Default = ".ms")
		type of file considered including "MSext" in their extension (at any positon)
	  doassessMS : boolean
		(Default = True)
		Run assess_ms scripts onto individual GOUS
	Outputs
	----------
	  Copy data into corresponding folders (see EMERGE folder structure)

	  Returns name of the visibility .ms file containing all visibilities

	'''
	print("==================================================================")
	print(" Starting do_extractvis()...")
	#---------------------------------
	# Environment

	# Generate environment if needed
	#gen_foldertree(workdir=workdir)

	#---------------------------------
	# Frequency
	# If a string try to identify it in emerge.linecat
	if (type(myline) == str):
		myrestfreq_Hz = emerge_linecat[myline]*1E9		# in Hz
	# If numeric, it will attempt the given frequency
	if (type(myline) == float):
		myrestfreq_Hz = myline

	print(" Line frequency = " + str(myrestfreq_Hz/1E9) + " GHz")

	#---------------------------------
	newvis = []

	# Find all .cal MS-files in emerge_archive
	allMS = listofMSfiles(datadir=datadir,listexpr=MSext)

	# Identify MS files
	for myMS in allMS:
		# Check if the MS file contains the right coordinates and frequency
		isinMS, mytargetID, mylineSPW , mychanresSPW = istarANDlineinMS(datadir + str(myMS),myrestfreq_Hz,mycoords,R0=R0,units='deg')

		# Calculate velocity resolution
		velres = mychanresSPW/myrestfreq_Hz*c.c.value

		# if satisfy criteria, then carry out contsub and split
		if (isinMS & (velres >= mindV) & (velres <= maxdV)):
			print(" Processing  ./vis_tmp/" + str(myMS))
			for mytargetID_list in mytargetID:
				newMScontsub = do_contsub(datadir + str(myMS),mysource=mytargetID_list,workdir="./vis_tmp/",mycat=emerge_linecat_full,linetol_LSRK_ms=linetol_ms)
				newMScontsubspw = do_split_singelspw(newMScontsub,mysource=mytargetID_list,mindV=mindV,maxdV=maxdV,myrestfreq_Hz=myrestfreq_Hz)
				if (type(newMScontsubspw)==str):
					newvis.append(newMScontsubspw)
				if (type(newMScontsubspw)==list):
					newvis = newvis+newMScontsubspw
			# FAST??
			#newMSspw = do_split_singelspw(datadir + str(myMS),mysource=mytargetID,myrestfreq_Hz=myrestfreq_Hz)
			#newMScontsub = do_contsub(newMSspw,mysource=mytargetID,workdir="./vis_tmp/",mycat=emerge_linecat,linetol_LSRK_ms=linetol_ms)
			#newvis.append(newMScontsub)

	#---------------------------------
	# Identify relevants GOUS
	newGOUS = []

	# Within a GOUS there is a unique combination of projID + sourceID
	pattern = re.compile(r'_Proj(.*?)_Dim')
	GOUS = [pattern.search(s).group(1) for s in newvis if pattern.search(s)]
	

	for myGOUS in np.unique(GOUS):

		# Indentify the number of files within the same GOUS
		myMS2concat = listofMSfiles(datadir="./vis_tmp/",listexpr=myGOUS)
		myMS2concat = ["./vis_tmp/" + item for item in myMS2concat]

		# Find frequency information
		pattern = re.compile(r'_f0(.*?)GHz')
		myf0 = [pattern.search(s).group(1) for s in myMS2concat if pattern.search(s)]

		# Find arrays
		whicharray = au.mixedDishDiameters
		pattern = re.compile(r'Dim(.*?).contsub')
		myarrays = [pattern.search(s).group(1) for s in myMS2concat if pattern.search(s)]
		myarrays = np.unique(myarrays)
		if (len(myarrays)==1):
			myarrays = myarrays[0]
		else:
			myarrays = "+".join(myarrays)


		# Generate concatenated visibilities
		myfinalMS = "GOUS" + str(myGOUS) + "_ALMA" + str(myarrays) + "_" + str(myf0[0]) + "GHz_allvis.ms"
		ct.concat(vis=myMS2concat, concatvis="./vis_tmp/" + myfinalMS)
		newGOUS.append("./vis_tmp/" + myfinalMS)

		# Regrid
		#myGOUSinfo = getMSinfo_short("./vis_tmp/" + myfinalMS)
		#myGOUS_sou = myGOUSinfo.keys()
		#myGOUS_chanwidth_Hz = np.max([myGOUSinfo[myGOUS_sou]['SPWinfo_LSRK_Hz'][spws]["chanwidth"] for spws in myGOUSinfo[myGOUS_sou]['SPWinfo_LSRK_Hz']])
		#myGOUS_chanwidth_ms = myGOUS_chanwidth_Hz/myrestfreq_Hz*c.c.value

		# -----------------------------------
		# Get ASSESSMS per GOUS
		if (doassessMS):
			if (os.path.isdir("./"+str(myGOUS)+"_tmp")): 
				# Remove folder if exist
				os.system("rm -rf "+str(myGOUS)+"_tmp")
			# Create folder for symlinks
			os.system("mkdir "+str(myGOUS)+"_tmp")
			mypath = os.getcwd()+"/"
			myMS_symlin = 0
			for myMS in myMS2concat:
				myMS_symlin += 1 
				os.system("ln -s " + str(mypath) + str(myMS) + " ./" + str(myGOUS)+"_tmp/" + str(myGOUS) + "_MS" + str(myMS_symlin))

			# Write variables into DCpar template
			# Open the input file, read its contents, and replace the keyword
			with open("./scripts_tmp/emerge_run_assess_ms_public_template.py", "r") as file:
			    file_contents = file.read()  # Read the entire file content
			
			# Replacements
			replacements = {
				"path2assessMS": str(path2assessMS),
				"TMPFOLDER": str(myGOUS)+"_tmp"
			}
			# Perform all replacements
			for keyword, replacement in replacements.items():
			    file_contents = file_contents.replace(keyword, replacement)
			

			# Write the updated contents to the output file
			with open("./scripts_tmp/emerge_run_assess_ms_public.py", "w") as file:
			    file.write(file_contents)

			# Execute assess_ms script
			exec(open("./scripts_tmp/emerge_run_assess_ms_public.py").read())

			# Cleaning up
			os.system("rm -rf " + str(myGOUS)+"_tmp")

	# If no debugging mode, then clean up
	if (debugging == False):
		os.system("rm -rf ./vis_raw")
		os.system("rm -rf ./vis_tmp/*split*")

	#---------------------------------
	# Return list of calibrated MS files 
	return newGOUS

	print(" do_extractvis() ... DONE")
	print("==================================================================")


##----------------------------------------------------
## Step 2: Prepare SD data

def do_shortspacing(myline,Nchans,VLSR0=0.0,BW=0.0):
	'''

	do_shortspacing (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Generate SD data (if available)

	Arguments
	----------
	  myline : str or float
		If "str" the line ID will be identified in emerge.linecat (e.g. "N2H+10").
		If else "float" the line will be taken as frequency in Hz
	  Nchans : integer
		Number of channels
	  VLSR0 : float (optional) 
		(Default = 0.0E3 m/s)
		Source velocity [km/s] 
	  BW	: float (optional)
		(Default = 0.0 Hz)
		Total bandwidth [in Hz]; if BW != 0.0 then Nchans adapted to match BW. Otherwise Nchans is taken as input
	Outputs
	----------
	  Creates the final SD image for DC with Nchan channels.
	  If more than one FITS file, then it combines it into a SD mosaic

	  This file is stored as ./FITS_tmp/shortspacingXXXX where XXXX are the spectral configuration of the SD

	'''
	print("==================================================================")
	print(" Starting do_shortspacing()...")

	# Find how many SD FITS data are available
	SPdata = read_allfiles(listexp="./FITS_raw/*",areimages=False)

	# Final short-spacing map
	SDmap = "SDmosaic"

	# If no SP data, then exit module
	if (len(SPdata) == 0):
		return None

	# If only one FITS (either SD or TP) then just use it directly
	if (len(SPdata) == 1):
		SDmap = "./FITS_tmp/" + str(SDmap)
		importfits(fitsimage=str(SPdata[0]),imagename=SDmap)

	# If more than one FITS (TP) data, then do a mosaic
	if (len(SPdata) > 1):
		SDmap = do_mosaic(fits_files=SPdata, fits_files_ref = 0, output_mosaic=SDmap, do_baseorder_before = -1, do_baseorder_after = -1, do_sigweig = True, do_linweig = False, debugging_files = False)

	os.system(" rm -fr ./FITS_tmp/*tmp* ./FITS_tmp/*fits")

	# ----------------------------------------
	# get axes information
	print(" File name =======" + str(SDmap))
	ia.open(SDmap)
	axes = ia.coordsys().names()
	ia.done()

	# get frequency axis
	freqaxis = axes.index("Frequency")
	refincre = imhead(str(SDmap),mode='get',hdkey='cdelt'+str(freqaxis+1))['value']
	# Recalculate Nchans
	BWini = np.abs(refincre)*Nchans
	if (BWini <= BW):
		Nchans = int(BW/abs(refincre))

	# ----------------------------------------
	# Extract target line

	extract_targetline(str(SDmap),myline=myline,output_name="./FITS_tmp/shortspacing",linecat=emerge_linecat,VLSR_ms=VLSR0,Nchans=Nchans,doFITS=True)
#	# Clean-up
#	print("-------- " + str(SDmap))
#	os.system("rm -rf " + str(SDmap))


	print(" do_shortspacing() ... DONE")
	print("==================================================================")


def do_compare_chanwidth(intvis,shortspacing,myline):
	'''

	do_compare_chanwidth() (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Compares channel widths of the INT vs SD data

	Arguments
	----------
	  intvis : str 
		Path to interferometric visibilities.
	  shortspacing : str 
		Path to SD image
	  myline : str or float
		If "str" the line ID will be identified in emerge.linecat (e.g. "N2H+10").
		If else "float" the line will be taken as frequency in Hz 

	Outputs
	----------
	  Compares the channelwidths of the INT vs SD data

	  Returns dictionary with channel values (max,INT,SD) [0,1,2] and whether SD needs to be resampled False/True [3]

	'''
	print("==================================================================")
	print(" Starting do_compare_chanwidth()...")

	#---------------------------------
	# Frequency
	# If a string try to identify it in emerge.linecat
	if (type(myline) == str):
		myrestfreq_Hz = emerge_linecat[myline]*1E9		# in Hz
	# If numeric, it will attempt the given frequency
	if (type(myline) == float):
		myrestfreq_Hz = myline

	print(" Line frequency = " + str(myrestfreq_Hz/1E9) + " GHz")

	#---------------------------------
	# Interferometric visibilities
	myMSinfo = getMSinfo_short(intvis)
	myMSinfo_sou = list(myMSinfo.keys())[0]
	myMSinfo_chanwidth_Hz = np.max([myMSinfo[myMSinfo_sou]['SPWinfo_LSRK_Hz'][spws]["chanwidth"] for spws in myMSinfo[myMSinfo_sou]['SPWinfo_LSRK_Hz']])
	myMSinfo_chanwidth_ms = myMSinfo_chanwidth_Hz/myrestfreq_Hz*c.c.value

	print(" Interometrid data: chanwidth = " + str(myMSinfo_chanwidth_ms) + " m/s")

	#---------------------------------
	# Short spacing information
	shortspacing_chanwidth_Hz = imhead(shortspacing,mode='get',hdkey='cdelt3')['value']
	shortspacing_chanwidth_ms = shortspacing_chanwidth_Hz/myrestfreq_Hz*c.c.value

	print(" Short-spacing data: chanwidth = " + str(shortspacing_chanwidth_ms) + " m/s")

	# Does SD needs resampling?
	resamp = False
	if (np.abs(np.round(shortspacing_chanwidth_ms)) < np.abs(np.round(myMSinfo_chanwidth_ms))):
		resamp = True

	#---------------------------------
	chan_values = {'chanmax_ms' : np.max([np.abs(myMSinfo_chanwidth_ms),np.abs(shortspacing_chanwidth_ms)]),
			'vis_chan_ms' : myMSinfo_chanwidth_ms,
			'short_chan_ms' : shortspacing_chanwidth_ms,
			'SD_needs_resamp' : resamp}

	# Return a dictionary of channel values values
	return chan_values

	print(" do_compare_chanwidth() ... DONE")
	print("==================================================================")



def do_fetchSDdata(mycoords,R0,myline,VLSR_ms=0.0E3,mindV=0.0E3,maxdV=1.0E3,datadir="./archive/",workdir="./",do_copy=False):
	'''

	do_fetchSDdata (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Prepares data (FITS files) for pipeline work
	The function will look for a line and will try to in

	Arguments
	----------
	  mycoords : list or str 
		Target coordinates [RA,DE] if in degrees
		e.g., [REdeg,DEdeg] = [4.232545,-32.9452343]
		     target coordinates ['RADEC string'] if in absolute coordinates
		      [RAabs,DEabs]] = ['12:56:11.16658 -05.47.21.5246']
	  R0 : float
		(Default = 0.5 deg)
		Seach radius in deg
	  myline : str or float
		If "str" the line ID will be identified in emerge.linecat (e.g. "N2H+10").
		If else "float" the line will be taken as frequency in Hz
	  VLSR : float  (optional) 
		(Default = 0.0E3 m/s)
		source velocity in [m/s]; If known use it, otherwise set to 0.0m/s
	  dVmin : float (optional) 
		(Default = 0.0E3 m/s)
		minimum velocity resolution [m/s] 
	  dVmax : float (optional) 
		(Default = 1.0E3 m/s)
		maximum velocity resolution [m/s] 
	  datadir : str
		(Default = ".")
		path to FITS data 
	  workdir : str
		(Default = ".")
		working and writing directory tree. See also gen_foldertree().
	  do_copy : boolean
		(Default = False)
		Copy data from datadir to workdir

	Outputs
	----------
	  Copy data into corresponding folders (see EMERGE folder structure)

	'''
	print("==================================================================")
	print(" Starting do_fetchSDdata()...")

	#---------------------------------
	# Environment

	# Generate environment if needed
	#gen_foldertree(workdir=workdir)

	#---------------------------------
	# Frequency
	# If a string try to identify it in emerge.linecat
	if (type(myline) == str):
		myrestfreq_Hz = emerge_linecat_full[myline]*1E9		# in Hz
	# If numeric, it will attempt the given frequency
	if (type(myline) == float):
		myrestfreq_Hz = myline

	print(" Line frequency = " + str(myrestfreq_Hz/1E9) + " GHz")

	#---------------------------------
	# Short spacing data
	SDdata = 0
	TPdata = 0
	#ShortSpaData = SDdata+TPdata

	# Get skyfrequency after VLSR correction
	myfreq = getskyfreq(myrestfreq_Hz,myVLSR=VLSR_ms)
	print(" Source velocity = " + str(VLSR_ms/1E3) + " km/s")
	print(" Sky frequency = " + str(myfreq/1E9) + " GHz")	
	
	# Getting short-spacing data	
	# (1) Check if SD data exist, assumed to be from a big SD telescope (preferred)
	if (os.path.isdir(str(datadir)+"shortspacing_SD")):
		SDdata = select_FITS(mycoords,R0=R0,freq0=myfreq,dVmin=mindV,dVmax=maxdV,datadir=str(datadir)+"shortspacing_SD",workdir=workdir+"/FITS_raw",do_copy=do_copy)

	# (2) If no SD data, then try local ALMA-TP data
	if ((SDdata == 0) & (os.path.isdir(str(datadir)+"shortspacing_TP"))):
		# Check TP data exist
		TPdata = select_FITS(mycoords,R0=R0,freq0=myfreq,dVmin=mindV,dVmax=maxdV,datadir=str(datadir)+"shortspacing_TP",workdir=workdir+"/FITS_raw",do_copy=do_copy)

	# (3) If no local data, then Try to download data from JVO database
	if ((SDdata == 0) & (TPdata == 0)):
		try:
			import alminer
			allGOUS = find_GOUS(targetcoord=mycoords,targetline=myline,targetfreq_Hz=myfreq,maxdV=maxdV,mindV=mindV,R0=R0)
			TPdata = getTPfiles(allGOUS=allGOUS,skyfreq_Hz=myfreq,maxdV=maxdV,mindV=mindV); print(" Total of FITS files found =  " + str(TPdata))
		except:
			print("No Archival data found")

	# Are there short spacing data?
	if ((SDdata == 0) & (TPdata == 0) ):
		# If no short-spacing information available, then exit
		print(" Neither SD nor ALMA-TP data available. No short-spacing information available. Exit.")
		sys.exit()
	
	# Final number of short-spacing datasets
	ShortSpaData = SDdata+TPdata

	#---------------------------------
	# Return results: number of datasets for SD, TP, and Int data
	return ShortSpaData

	print(" do_fetchSDdata() ... DONE")
	print("==================================================================")





def do_prepareSD(myvisfile,myline,Nchans,VLSR_ms=-1E6):
	'''

	do_prepareSD (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Compares channel widths of the INT vs SD data

	Arguments
	----------
	  myvisfile : str 
		Path to target visibilities.
	  myline : str or float
		If "str" the line ID will be identified in emerge.linecat (e.g. "N2H+10").
		If else "float" the line will be taken as frequency in Hz 
	  Nchans : int
		Number of channels
	  VLSR_ms : float
		(Dafult = -1E6 m/s)
		Source velocity. Leave it as Default if unknown or dubious

	Outputs
	----------
	  Prepares the SD for DC. After producing the SD mosaic (using do_shortspacing), it compares the channelwidths
	  of the INT vs SD data (see also do_compare_chanwidth) and resamples it if the resolution of the INT data is worse,

	  Returns the path to the SD data ready to be DC with the INT dataset.

	'''
	print("==================================================================")
	print(" Starting do_prepareSD()...")

	# Compare with visibilities
	allvisfile_info = getMSinfo_short(myvisfile[0])
	mysource = list(allvisfile_info.keys())[0]

	# Evaluate Bandwidth
	chanwidth = allvisfile_info[mysource]['SPWinfo_LSRK_Hz'][0]['chanwidth'] # Hz; channel width
	BW_target = Nchans*abs(chanwidth)

	# Generate SD mosaic
	if (VLSR_ms == -1E6):
		do_shortspacing(myline,Nchans,VLSR0=allvisfile_info[mysource]['Source_Vel_ms-1'],BW=BW_target)
	else:
		do_shortspacing(myline,Nchans,VLSR0=VLSR_ms,BW=BW_target)

	# SD image (variable) for DC script
	#sdimage_input = "./FITS_tmp/" + read_allfiles("./FITS_tmp/")[0]
	sdimage_input = listofMSfiles("./FITS_tmp/",listexpr="shortspacing")
	sdimage_input = "./FITS_tmp/" + str(sdimage_input[0])

	# Compare channels in INT vs SD data
	chans_ms = do_compare_chanwidth(myvisfile[0],sdimage_input,myline)

	# Rebin SD data?
	if (chans_ms['SD_needs_resamp']):

		# Identify 'Frequency' axis
		im_info = imhead(sdimage_input,mode="summary")
		myidx = list(im_info['axisnames']).index('Frequency')

		# Get channel info and update tmp_newchan file
		os.system("cp -r " + str(sdimage_input) +  " ./FITS_tmp/tmp_newchan.image")
		im_info = imhead("./FITS_tmp/tmp_newchan.image",mode="summary")
		myidx = list(im_info['axisnames']).index('Frequency')
		myrestfreq_Hz = imhead("./FITS_tmp/tmp_newchan.image",hdkey="restfreq",mode='get')['value']
		if (chans_ms["short_chan_ms"] < 0.0):		
			chan_sign = -1. 
		else:
			chan_sign = 1.
		chans_Hz = chan_sign*chans_ms["chanmax_ms"] * myrestfreq_Hz / c.c.value
		imhead("./FITS_tmp/tmp_newchan.image",hdkey="cdelt" + str(myidx+1),mode='put',hdvalue=str(chans_Hz))

		# Resample into new chanwidth
		pattern = re.compile(r'(.*?)dV')
		myf1 = pattern.search(sdimage_input)[0]	# 1st part of the filename	
		pattern = re.compile(r'km(.*)')
		myf2 = pattern.search(sdimage_input)[0]	# 2nd part of the filename
		newchanwidth = chans_ms['chanmax_ms']
		newfilename = myf1 + str(np.round(abs(chans_ms['chanmax_ms']/1E3),3)) + myf2
		imregrid(imagename = sdimage_input, template="./FITS_tmp/tmp_newchan.image", output = newfilename, overwrite = True)

		# Cleaning up
		os.system("rm -rf ./FITS_tmp/tmp*")
		
		# Final steps
		sdimage_input = newfilename

	# Return
	#os.system("cp -r " + str(sdimage_input) + " ./products/.")

	return sdimage_input

	print(" do_prepareSD() ... DONE")
	print("==================================================================")


##----------------------------------------------------
## Step 4: Imaging

	
# Imaging: Using a wrapper of DC scripts (by Plunkett et al 2023)
def do_lineimaging(myINTvis,mySD,Nchans,myline,Niters,VLSR_ms=-1E6,DCsteps=[1,2,3,5,8],DCdeconv="HB"):
	'''

	do_lineimaging (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Wrapper of the DC scripts (see Plunkett+2023) for line imaging and data combination

	Arguments
	----------
	  myINTvis: str 
		Path to target visibilities.
	  mySD: str
		Path to SD data
	  Nchans : int
		Number of channels used for deconvolution
		Note that the channel width will be automatically selected and will maximize
		the spectral resolution of the data within the selection parameters chosen before
	  myline : str or float
		Reference frequency for deconvolution. The sky frequency is selected after correction by
		the source velocity according to the MS file.
		If "str" the line ID will be identified in emerge.linecat (e.g. "N2H+10").
		If else "float" the line will be taken as frequency in Hz 
	  Niters : int
		Number of iterations for decolvolution
	  VLSR_ms : float (optional)
		(Default = -1E6 m/s)
		If Default (recommended), it will look for source VLSR in MS files
	  DCsteps : list
		Steps for DCscripts (see definitions in emerge_DCpars_template.py)
	  DCdeconv : str
		(Default =  "HB")
		Deconvolution method: HB = hogbom, MS = multiscale
		Note: we notice MS produces instabilities in several automatic reductions
		MS might be better for extended structures but requires fine tunning not suitable for 
		our automatic pipeline.

	How it works
	----------
		Our deconvolution process takes advantage of the DCscripts. 
		Most parameters are then taken as standard from these scripts:
			- weighting, masking, ...
		Specific parameters are automatically calculated by the script:
			- skyfreq, phasecenter, resolution, pixel size...
		Only a reduced number of parameters is then needed as imput (see above)

		Steps:
		(1) Calculate specific parameters 
		(2) Create a DCpar file using a predefined template
		(3) Run DC scripts
		(4) Copy relevant results into /product folder

	Outputs
	----------
	  

	'''
	print("==================================================================")
	print(" Starting do_lineimaging()...")

	## --------------------------------------
	## Step 1: Calcualte specific parameters
	##
	## Sampling parameters (hard coded)
	target_oversamp = 4		# Oversampling factor wrt beamsize_final; Note that the final image will be smoothed and Nyquist rebinned
					# (Default = 4.)
	target_smoothfactor = 1.1	# Smoothing factor wrt beamsize in .MS file; 
					# (Default = 1.1; i.e. 20% larger than expected beamsize)

	# Retrieve information from .MS file
	myINTvis_info = getMSinfo(str(myINTvis))
	mysource = list(myINTvis_info.keys())[0]

	# Estimate expected beamsize, chanwidth...
	beamsize_expected = au.estimateSynthesizedBeam(str(myINTvis))
	beamsize_final = np.round(target_smoothfactor*beamsize_expected,decimals=1)

	# Which arrays are included? (to be added to final image name)
	Intarrays = myINTvis.split("/")[-1]
	Intarrays = re.search("ALMA(.*?)m_",Intarrays)[0][:-1] + "+TP"

	# Spectral paramters
	target_dv_ms = do_compare_chanwidth(str(myINTvis),mySD,myline)	# spectral resolution, taken as the worst between SD and INT
	target_dv = target_dv_ms['chanmax_ms']/1E3
	if (VLSR_ms == -1E6):
		target_v0 = myINTvis_info[mysource]['Source_Vel_ms-1']/1E3	# km/s; reference velocity, taken from .MS file
	else:
		target_v0 = VLSR_ms

	target_vmin = float(np.round(target_v0-Nchans/2.*target_dv_ms['chanmax_ms']/1E3,1))	# km/s; starting velocity for tclean. Assumed as V0-dv*Nchan/2

	# Mask for DC to estimate noise minmaskchan:maxmaskchan
	minmaskchan = 5
	maxmaskchan = int(Nchans/5.)		

	# Mosaic parameters
	target_phasecenter = SkyCoord(ra=myINTvis_info[mysource]['phasecenter']["RAdeg"]*u.degree, dec=myINTvis_info[mysource]['phasecenter']["DEdeg"]*u.degree)
	target_pixsize = beamsize_final/target_oversamp			# arcsec; pixel size
	target_mosaic = myINTvis_info[mysource]['mosaic']		#
	dRA = abs(target_mosaic['RAmaxdeg']-target_mosaic['RAmindeg'])	# deg; Mosaic dimention in RA
	dDE = abs(target_mosaic['DEmaxdeg']-target_mosaic['DEmindeg'])	# deg; Mosaic dimention in RA
	target_RAnpix = int(1.05*dRA*3600./target_pixsize)			# pixels; map size in DE for tclean wrt pixsize
	target_DEnpix = int(1.05*dDE*3600./target_pixsize)			# pixels; map size in DE for tclean wrt pixsize
	
	# Check if they are odd numbers. Odd numbers must be avoided to prevent issues in FFTs.
	if target_RAnpix % 2 == 0:
		pass
	else:
		target_RAnpix = target_RAnpix+1
	if target_DEnpix % 2 == 0:
		pass
	else:
		target_DEnpix = target_DEnpix+1

	# Calculate effective SD dish size for SDint
	# SDint looks for a SD dishsize such as halfpb=Quantity(1.22*C::c/freq/dishDiam_p, "rad")
	SDinfo = imhead(mySD,mode="list")
	SD_beam_arcsec = SDinfo["beammajor"]["value"]  # in arcsec
	SD_restfreq_Hz = SDinfo["restfreq"]
	effdishdim = float(1.22*c.c.value/SD_restfreq_Hz/SD_beam_arcsec*3600.*180./np.pi)

	## --------------------------------------
	## Step 2: Create DCpars file
	##
	# Write variables into DCpar template
	# Open the input file, read its contents, and replace the keyword
	with open("./scripts_tmp/emerge_DCpars_template.py", "r") as file:
	    file_contents = file.read()  # Read the entire file content
	
	# Replacements
	replacements = {
		"mySD_tmp": str(mySD),
		"myINTvis_tmp": str(myINTvis),
		"mysource_tmp": str(mysource) + "_" + str(myline) + "_" + str(Intarrays),
		"mytarget_phasecenter": str(target_phasecenter.to_string("hmsdms")),
		"mytarget_RAnpix": str(target_RAnpix),
		"mytarget_DEnpix": str(target_DEnpix),
		"mytarget_pixsize": str(target_pixsize),
		"mytarget_vmin": str(target_vmin),
		"mytarget_dv": str(target_dv),
		"myNchans": str(Nchans),
		"myminmaskchan": str(minmaskchan),
		"mymaxmaskchan": str(maxmaskchan),
		"mytarget_freq": str(emerge_linecat[myline]),
		"myNiters": str(Niters),
		"myDCsteps": str(DCsteps),
		"myDCdeconv": str(DCdeconv),
		"mySDintdishdim" : str(effdishdim)
	}

	# Perform all replacements
	for keyword, replacement in replacements.items():
	    file_contents = file_contents.replace(keyword, replacement)
	

	# Write the updated contents to the output file
	with open("./scripts_tmp/emerge_DCpars.py", "w") as file:
	    file.write(file_contents)

	## --------------------------------------
	## Step 3: Execute DC_run
	##
	# DC scripts wrapper
	##exec(open("../DataComb/templates/DC_pars_Template.py").read(),globals())
	exec(open("./scripts_tmp/emerge_DCpars.py").read(),globals())
	exec(open("./scripts_tmp/DataComb/DC_run.py").read(),globals())


	print(" ... do_lineimaging() DONE")
	print("==================================================================")
	



##----------------------------------------------------
## Step 5: Copying data into ./products and cleaning up

# Smooth and regrid
def do_smooreb(myimage,mytargetbeam,target_oversamp=8):
	'''

	do_smooreb (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function to smooth image into a target resolution and resample it 

	Arguments
	----------
	  myimage : str 
		Image to be convolved

	  mytargetbeam : float 
		Target beamsize in arcsec

	  target_oversamp : int 
		(Default = 8)
		Pixel obsersampling of the beam

	Outputs
	----------
	  Creates a Nyquist sample image with a final resolution of mytargetbeam

	  Returns the name of the final image and the corresponding pbimage (in case is needed for further pbcut)

	'''
	print("==================================================================")
	print(" Starting do_smooreb()...")
	print(" Target image = " +str(myimage) )

	# Nyquist ratio
	Nyquist = 2					# i.e. every 2 beams
	factor_sampling = int(target_oversamp/2)	# resampling factor for imrebin (see below)

	# Smooth image to target beamsize
	imsmooth(imagename=myimage,outfile=myimage+".smoo"+str(mytargetbeam)+'arcsec',
		beam = {'major': str(mytargetbeam)+'arcsec', 'minor': str(mytargetbeam)+'arcsec', 'pa': '0deg'},targetres=True)
	# Rebin to Nyquist
	imrebin(imagename=myimage+".smoo"+str(mytargetbeam)+'arcsec',outfile=myimage+".smoo"+str(mytargetbeam)+'arcsec_Nyquist',factor=[factor_sampling,factor_sampling,1,1])

	# Relevant file names: final convolved image + pbimage
	image_pbcor = myimage						# pbimage
	image_final = myimage+".smoo"+str(mytargetbeam)+'arcsec_Nyquist'	# final image

	# Return relevant images filenames
	return image_final, image_pbcor

	print(" Smoothed image = " + str(image_final))
	print(" ... do_smooreb() DONE")
	print("==================================================================")



# Copy data/images to /product folder
def do_getproducts(doconvo=True,doQAs=True,pbcut=0.5,getCASAima=False,debugging=False):
	'''

	do_copyproducts (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function to get all FITS products into ./products folder

	Arguments
	----------
	  doconvo : boolean
		Do product convolution into circular beam?
	  doQAs : boolean
		Run IQAs from DCsripts?
	  pbcut : float [0,1]
		(Default = 0.5)
		Primary beam cut to be applied to all images
	  getCASAima : boolean
		(default = False)
		Get also CASA images in addition to FITS? 
	  debugging : boolean
		Copy list of images to /product folder

	Outputs
	----------
	  Copies images to /product folder

	'''
	print("==================================================================")
	print(" Starting do_getproducts()...")

	# Copy SD regridded image
	SDlist =  listofMSfiles("./imaging_tmp/",listexpr="SD_ro-rg_INTpar.image",ends=True)
	if (getCASAima):
		os.system("cp -r ./imaging_tmp/" + str(SDlist[0]) + " ./products/.")
	exportfits(imagename= "./imaging_tmp/" + str(SDlist[0]),fitsimage="./products/"+str(SDlist[0]) + ".fits")

	# Find the list of .pbcor files created after DC
	os.system("rm -rf ./imaging_tmp/*INTpar_template*") 			# Remove dirty images
	os.system("rm -rf ./imaging_tmp/*hybrid_f.*")
	DClist = listofMSfiles("./imaging_tmp/",listexpr=".image.pbcor",ends=True)
	PBima = listofMSfiles("./imaging_tmp/",listexpr=".pb",ends=True)  # PB[0] image, same for all

	for DCimage in DClist:
		# Get PBcut image (*.image.pbcorXX)
		impbcor(imagename= "./imaging_tmp/" + str(DCimage)[:-12] + ".image",pbimage= "./imaging_tmp/" + str(DCimage[:-12]) + ".pb", outfile="./imaging_tmp/" + str(DCimage)[:-12] + ".pbcor" + str(pbcut),cutoff=pbcut)
		# Copy CASA image (*.image.pbcorX)
		if (getCASAima):
			os.system("cp -r ./imaging_tmp/" + str(DCimage)[:-12] + ".pbcor" + str(pbcut) + " ./products/.")
		exportfits(imagename= "./imaging_tmp/" + str(DCimage)[:-12] + ".pbcor" + str(pbcut),fitsimage="./products/"+str(DCimage[:-12])+ ".pbcor" + str(pbcut) + ".fits")
		
		# Get moment 0 of them
		immoments("./imaging_tmp/" + str(DCimage)[:-12] + ".pbcor" + str(pbcut),outfile="./imaging_tmp/" + str(DCimage)[:-12] + ".pbcor" + str(pbcut) + ".mom0",moments=0)
		exportfits(imagename= "./imaging_tmp/" + str(DCimage)[:-12] + ".pbcor" + str(pbcut) + ".mom0",fitsimage="./products/" + str(DCimage)[:-12] + ".pbcor" + str(pbcut) + ".mom0" + ".fits")

		# Copy also .pb image
		if (getCASAima):
			os.system("cp -r ./imaging_tmp/" + str(DCimage[:-12]) + ".pb  ./products/.")
		exportfits(imagename="./imaging_tmp/" + str(DCimage[:-12])+".pb",fitsimage="./products/" + str(DCimage[:-12])+".pb.fits")

#		if (getCASAima):
#			os.system("cp -r ./imaging_tmp/" + str(PBima[0]) + ".pb  ./products/.")
#		exportfits(imagename=PBima[0],fitsimage="./products/" + str(DCimage[:-12])+".pb.fits")

		# Copy PSF (if exists; not always, e.g. hybryd_f1.0 , same as hybrid_f)
		if(os.path.exists("./imaging_tmp/" + str(DCimage[:-12]) + ".psf")):
			if (getCASAima):
				os.system("cp -r ./imaging_tmp/" + str(DCimage[:-12]) + ".psf  ./products/.")
			exportfits(imagename="./imaging_tmp/" + str(DCimage[:-12])+".psf",fitsimage="./products/" + str(DCimage[:-12])+".psf.fits")

		
	# Run data convolution?
	if (doconvo):
		# Convolution into a circular beam with beamsize = major beam
		myimage_smoo_final_all = []	# list with all file names after smoothing
		for DCimage in DClist:
			# full name
			myimage = "./imaging_tmp/" + DCimage[:-12] + ".pbcor" + str(pbcut)
			# Find major beam
			myimage_head = imhead(str(myimage))
			mytargetbeam = myimage_head['restoringbeam']['major']
			# Smooth to circular beam and resample image
			DCimage_smoo_final, DCimage_final = do_smooreb(myimage,np.round(1.1*mytargetbeam['value'],2),target_oversamp=8)
			# Copy smoothed image to ./products
			if (getCASAima):
				os.system("cp -r " + str(DCimage_smoo_final) + " ./products/.")
			DCimage_smoo_final_name = DCimage_smoo_final.split("/")[-1]
			exportfits(imagename=str(DCimage_smoo_final),fitsimage="./products/" + str(DCimage_smoo_final_name)+".fits")

			# Get moment 0 of them
			immoments(str(DCimage_smoo_final),outfile="./imaging_tmp/" + str(DCimage_smoo_final_name)+".mom0",moments=0)
			if (getCASAima):
				os.system("cp -r " + str(DCimage_smoo_final_name)+".mom0" + " " + "./products/" + str(DCimage_smoo_final_name)+".mom0")
			exportfits(imagename="./imaging_tmp/" + str(DCimage_smoo_final_name)+".mom0",fitsimage="./products/" + str(DCimage_smoo_final_name)+".mom0.fits")
			#imview(raster="./products/" + str(DCimage_smoo_final_name)+".mom0",out="./products/" + str(DCimage_smoo_final_name)+".mom0.png")

			# Append file names
			#myimage_smoo_final_all.append(myimage_smoo_final)

	if (doQAs):
		# Get file lists
		# SD image (reference)
		SDimage = listofMSfiles("./imaging_tmp/",listexpr="SD_ro-rg_INTpar.image",ends=True)
		SDimage = "./imaging_tmp/" + str(SDimage[0])

		if (doconvo):
			# Smoothed products (IQA targets)
			DCimage_final_all = listofMSfiles("./imaging_tmp/",listexpr="Nyquist",ends=True)
			DCimage_final_all = ["./imaging_tmp/" + s for s in DCimage_final_all]
		else:
			# Normal products (IQA targets)
			DCimage_final_all = listofMSfiles("./imaging_tmp/",listexpr="image.pbcor"+str(pbcut),ends=True)
			DCimage_final_all = ["./imaging_tmp/" + s for s in DCimage_final_all]

		# Get PB images
		DCimage_final_allPBs = listofMSfiles("./imaging_tmp/",listexpr="pb",ends=True)
		DCimage_final_allPBs = ["./imaging_tmp/" + s for s in DCimage_final_allPBs]

		# Get labels for plots
		#pattern1 = re.compile('products/(.*?).cube')
		#DCnames = [pattern1.search(s).group(1) for s in DCimage_final_all if pattern1.search(s)]
		#pattern2 = re.compile('_nIA_(.*?).image.pbcor')
		#DCtypes = [pattern2.search(s).group(1) for s in DCimage_final_all if pattern2.search(s)]
		#DClabels = np.char.add(DCnames,DCtypes)
		#pattern3 = re.compile('image.pbcor.smoo(.*?)_Nyquist')
		#DCsmoo = [pattern3.search(s).group(1) for s in DCimage_final_all if pattern3.search(s)]
		#DClabels = np.char.add(DClabels,DCsmoo)

		# PB images (for weighting)
		#PBimage =  listofMSfiles("./products/",listexpr="tclean.pb",ends=True)
		#PBimage = "./products/" + str(PBimage[0])

		# Attempt to estimate SD image RMS (use first as reference) in first 1/5th of the BW
		#freqaxis = imhead(SDimage)['axisnames'] == "Frequency"	# Identify Freq axis
		#Ntot = imhead(SDimage)['shape'][freqaxis]			# Find Nchans
		#SDimageRMS = imstat(SDimage,chans="2~"+str(int(Ntot/10.)))["sigma"][0]	# stats 

		# Run IQA for smoothed images
		#for mytargetima,mytaregtname in zip(DCimage_final_all,DClabels):
		for mytargetima,mytargetPB in zip(DCimage_final_all,DCimage_final_allPBs):

			# Names
			myimagename=mytargetima
			myplotname=mytargetima.split("/")[-1]

			# Attempt to estimate Int image RMS in first 1/5th of the BW
			freqaxis = imhead(mytargetima)['axisnames'] == "Frequency"	# Identify Freq axis
			Ntot = imhead(mytargetima)['shape'][freqaxis]			# Find Nchans
			DCimageRMS = imstat(mytargetima,chans="2~"+str(int(Ntot/10.)))["sigma"][0]	# stats

			# get IQAs
			iqa.get_IQA(ref_image=SDimage,
				target_image=[mytargetima],
				pb_image=mytargetPB,
		#		masking_RMS=0.0)
				masking_RMS=DCimageRMS)
				
			iqa.Compare_Apar(ref_image=SDimage,
				target_image=[mytargetima],
				save=True,
				labelname=myimagename,
				plotname="./reports/Apar_" + str(myplotname))  

			iqa.Compare_Apar_cubes(ref_image=SDimage,
				target_image=[mytargetima],
				save=True,
				labelname=myimagename,
				plotname="./reports/AparSpectrogram_" +  str(myplotname))  

			iqa.Compare_Apar_signal(ref_image=SDimage,
				target_image=[mytargetima],
				save=True,
				labelname=myimagename,
				plotname="./reports/Apar_vs_signal_" +  str(myplotname)) 

	# Cleaning up
	os.system("rm -rf convo2ref* temp.mask* fake*")
	os.system("rm -rf ./products/*convo2ref* ./products/*thrsh")


	print(" ... do_getproducts() DONE")
	print("==================================================================")


	# Cleaning up
	if (debugging == False):
		print("==================================================================")
		print(" Cleaning all intermediate and tmp datasets...")

		os.system("rm -rf *_tmp *_raw")

		print(" ... Finalcleaning DONE")
		print("==================================================================")




# Copy data/images to /product folder
def do_getSDproducts(workdir="./FITS_tmp/",debugging=False):
	'''

	do_getTPproducts (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Function to get all final SD images into ./products folder

	Arguments
	----------
	  workdir : str()
		(Default = ./FITS_tmp/)
		Working folder where TP images/FITS are located
	  debugging : boolean
		Copy list of images to /product folder

	Outputs
	----------
	  Copies images to /product folder

	'''
	print("==================================================================")
	print(" Starting do_getSDproducts()...")

	# Copy SD image
	os.system("cp -r " + str(workdir) + "SDmap* ./products/.")

	print(" ... do_getSDproducts() DONE")
	print("==================================================================")


	# Cleaning up
	if (debugging == False):
		print("==================================================================")
		print(" Cleaning all intermediate and tmp datasets...")

		os.system("rm -rf *_tmp *_raw")

		print(" ... Finalcleaning DONE")
		print("==================================================================")



def FITSlinemosaic(myFITSima,myFITSpb):

	if (len(myFITSima)>1):
		num = 0
		imalist = []
		pblist = []
		
		# Import list of FITS images + pb into CASA images
		for ima,pb in zip(myFITSima,myFITSpb):
			num = num+1
			importfits(fitsimage=str(ima),imagename="tile_im"+str(num),overwrite=True)
			importfits(fitsimage=str(pb),imagename="tile_pb"+str(num),overwrite=True)
			imalist = np.append(imalist,"tile_im"+str(num))
			pblist = np.append(pblist,"tile_pb"+str(num))

		return imalist, pblist





	
