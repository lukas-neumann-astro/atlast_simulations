from scipy import fftpack
#from astropy.io import fits
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
import os
import copy
from numpy import inf
from matplotlib.colors import LogNorm 
from scipy.stats import kurtosis, skew
from scipy.optimize import curve_fit
from scipy import stats
import argparse
#from turbustat.statistics import PowerSpectrum
#import astropy.units as u
import matplotlib 
import matplotlib.pyplot as pyplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import glob

#import casatools# as cto
#import casatasks# as cta
import casatasks as cta
from casatasks import *

##-------------------------------------------------------
##-------------------------------------------------------

IQA_colours = ["red", "blue", "orange", "green" , "cyan", "pink", "brown","yellow","magenta","black", "grey", "purple", "darkorchid", "peru", "royalblue", "olive", "lightsalmon", "paleturquoise", "mocassin"]


## Taken from ALvaro's script
def CASA2fits(CASAfile):
    print(CASAfile)
    os.system("rm -rf "+  CASAfile+".fits")
    exportfits(imagename=CASAfile,fitsimage=CASAfile+".fits")
    
def get_simulation(myskymodel, myproject):

    ##-------------------------------------------------------
    # 12m ARRAY OBSERVATIONS

    #ALMA Full Ops Configuration: #1 the most copact, #28 the most extended
    #fullops=["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28"]

    print("Running the simulation for each Cycle5 configuration: the execution of the task can take some time")
    
    
    #simobserve
    simobserve(project=myproject, skymodel=myskymodel, direction = "J2000 10h00m00.0s -30d00m00.0s", obsmode = "int", antennalist = 'aca.cycle5.cfg', mapsize='10arcmin', incell='2.34arcsec', thermalnoise = '', setpointings=True, maptype='hex',
               pointingspacing='nyquist', integration='100s', totaltime='28800s')
    #simanalyze
    simanalyze(project=myproject, image = True, vis= myproject+".aca.cycle5.ms", imdirection = "J2000 10h00m00.0s -30d00m00.0s", niter = 5000, threshold = '0.0000000004Jy/beam', imsize=[256,256],
               analyze = True, showuv=False, showresidual=False, showconvolved=False, graphics = "", verbose = True, overwrite = True)

            
   # print("Running the simulation for each FullOps configuration: the execution of the task can take some time")
    #for i in fullops:
            #simobserve
        #    simobserve(project=myproject, skymodel=myskymodel, setpointings=True, direction="J2000 18h00m00.031s -22d59m59.6s", mapsize="0.76arcsec", obsmode="int", totaltime="1200s", antennalist="alma.out" + str(i) + ".cfg", thermalnoise = '',
          #             graphics='file')
            #simanalyze
           # simanalyze(project=myproject, image = True, modelimage=myproject+ ".alma.out" + str(i) + ".skymodel", vis=myproject+ ".alma.out" + str(i) + ".ms", imsize = [192, 192], niter = 10000, threshold = "1e-7Jy", weighting = "natural",
             #          analyze=True, showuv=False, showresidual=True, showconvolved=True, graphics = "file", verbose = True, overwrite = True)

##-------------------------------------------------------
##-------------------------------------------------------

## Convolve input image to the CASA simulated image resolution
## Taken from Alvaro's script (I've modified the output name)
def get_convo2target(convo_file,ref_image):
    # Get beam info from refence image
    hdr = imhead(ref_image,mode='summary')
    beam_major = hdr['restoringbeam']['major']
    beam_minor = hdr['restoringbeam']['minor']
    beam_PA = hdr['restoringbeam']['positionangle']
    ref_unit = hdr['unit']
    # Convolution into the same beam as the reference
    os.system("rm -rf tmp.tmp")
    imsmooth(imagename= convo_file,outfile= "tmp.tmp",kernel='gauss',major=beam_major,minor=beam_minor,pa=beam_PA,targetres=True)
    #imhead(convo_file, mode='put', hdkey='Bunit', hdvalue=ref_unit)
    os.system("rm -rf convo2ref")
    imregrid(imagename= "tmp.tmp",template= ref_image,output= convo_file + '_conv')
    imhead(convo_file + '_conv', mode='put', hdkey='Bunit', hdvalue=ref_unit)   # Lydia's modification to avoid lost bunits
    os.system("rm -rf tmp.tmp")

##-------------------------------------------------------
## Mask data (typically reference image)
## Taken from Alvaro's script (I've modified the output name)
def mask_image(myimage,threshold=0.0,relative=False):
    # Create a copy of your image
    #os.system('rm -rf masked.tmp')
    os.system('cp -r ' + str(myimage) + ' masked.tmp')
    # Create your mask
    ia.open('masked.tmp')
    if (relative == False):
        ia.calcmask(mask= 'masked.tmp >= '+str(threshold), name='mymask')
    if (relative == True):
        ima_sigma = imstat(myimage)["rms"][0]
        ia.calcmask(mask= 'masked.tmp >= '+str(threshold*ima_sigma), name='mymask')
    ia.close()
    os.system('mv masked.tmp ' + str(myimage) + '_masked')
    #makemask(mode='copy',inpimage='masked.tmp',inpmask=['masked.tmp:mymask'],output=str(myimage) + '_masked',overwrite=True)
    print(" New masked image : " + str(myimage) + '_masked')
    print("-----------------------------------------")
    print(" mask_image(): DONE ")
    print("=========================================")

##-------------------------------------------------------
## Drop Stokes axis
## Taken from Alvaro's script
def drop_axis(myimage):
    print("=================================================")
    print(" drop_axis(): drop additional axis (e.g. Stokes) ")
    print("=================================================")
    # reference: check axis 
    os.system("rm -rf " + myimage + "_subimage")
    imsubimage(imagename=myimage,outfile=myimage + "_subimage",dropdeg=True)
    print(" Reference image: " + str(myimage))
    print(" New image: " + str(myimage) +"_subimage")
    print("-----------------------------------------")
    print(" drop_axis(): DONE ")
    print("=========================================")

##-------------------------------------------------------
## Mask the target it similar to reference
## Taken from Alvaro's script and created a task
def mask_target(myimage, ref_image):
    #immath(imagename='convo2ref',mode='evalexpr',expr='IM0',outfile='convo2ref_masked',mask='mask('+str(ref_image.replace('-','\-').replace('_','\_'))+')')
    immath(imagename=myimage,mode='evalexpr',expr='IM0',outfile=myimage + "_masked",  mask='mask("'+str(ref_image)+'")')
    #os.system("mv convo2ref_masked " + target_image[j] + "_convo2ref")



# get values from FITS (ALL)
def get_ALLvalues(FITSfile,xmin,xmax,xstep):
    # FITS file
    image = fits.open(FITSfile)
    # Histogram
    bins_histo = np.arange(xmin,xmax,xstep)
    bins_mids = bins_histo[1:]-xstep/2.
    subimage = image[0].data.flatten()
    idxs = np.isfinite(subimage)
    hist, bin_edges = np.histogram(subimage[idxs],bins=bins_histo)
    ## values in log scale
    ##values_log = np.log10(hist.T)
    # return
    return 0. , bins_histo, bins_mids, hist.T


##-------------------------------------------------------
## Quality estimators
## taken from Alvaro's script

## Calculate Image Accuracy parameter (Apar)
def image_Apar(image,ref_image):

    # Resampling
    os.system('rm -rf tmp_resampled')
    imregrid(imagename= image,
         template= ref_image,
         #axes=[0, 1],
         output= 'tmp_resampled')
    os.system('rm -rf ' + image + '_Apar')
    # Q parameter
    immath(imagename=['tmp_resampled',ref_image], 
        outfile= image + '_Apar',
        expr='(IM0-IM1)/abs(IM1)')
    # Clean-up
    os.system('rm -rf tmp_resampled')

## Calculate image Fidelity
def image_Fidelity(image,ref_image):

    # Resampling
    os.system('rm -rf tmp_resampled')
    imregrid(imagename= image,
         template= ref_image,
         #axes=[0, 1],
         output= 'tmp_resampled')
    # Fidelity parameter
    os.system('rm -rf ' + image + '_Fidelity')
    immath(imagename=['tmp_resampled',ref_image], 
        outfile= image + '_Fidelity',
        expr='abs(IM1)/abs(IM1-IM0)')
    # Clean-up
    os.system('rm -rf tmp_resampled')


## Calculate image Difference
def image_Diff(image,ref_image):

    # Resampling
    os.system('rm -rf tmp_resampled')
    imregrid(imagename= image,
         template= ref_image,
         #axes=[0, 1],
         output= 'tmp_resampled')
    # Fidelity parameter
    os.system('rm -rf ' + image + '_Diff')
    immath(imagename=['tmp_resampled',ref_image], 
        outfile= image + '_Diff',
        expr='IM1-IM0')
    # Clean-up
    os.system('rm -rf tmp_resampled')

def noise_image(fitsfile,noise=0.1,noisefile="noise"):

    # Read Templeate
    hdu = fits.open(fitsfile)
    # Copy maks
    mask = np.isnan(hdu[0].data)
    # Create a noise dataset
    hdu[0].data = np.random.normal(0.,noise,(hdu[0].data.shape[0],hdu[0].data.shape[1])) # Noise image

def show_Apar_map(ref_image,target_image,#pathnametodrop,
                  channel=0, 
                  save=False, plotname='',
                  labelname='', titlename=''
                  ):

    # Figure
    fig = plt.figure(figsize=(10,15))
    #fig = plt.figure(figsize=(15,10))
    if titlename=='':   
        fig.suptitle('Accurray map', fontsize=16)
        #fig.suptitle('Accurray map for \n'+target_image.replace(pathnametodrop,''), fontsize=16)
    else:    
        fig.suptitle(titlename, fontsize=16)

    grid = plt.GridSpec(ncols=2,nrows=2, wspace=0.5, hspace=0.3)
    #grid = plt.GridSpec(ncols=2,nrows=2, wspace=0.3, hspace=0.3)

    # Panel #1: Reference
    ax1 = plt.subplot(grid[0, 0])
    image = fits.open(ref_image+".fits")
    # get min/max values
    #vmin , vmax = np.min(image[0].data[-np.isnan(image[0].data)]), np.max(image[0].data[-np.isnan(image[0].data)])
    vmin , vmax = np.min(image[0].data[~np.isnan(image[0].data)]), np.max(image[0].data[~np.isnan(image[0].data)])
    # Continuum or cube?
    Ndims = np.shape(np.shape(image[0].data))
    channel=int(channel)
    if (Ndims[0] == 2):
        # Continuum
        im = ax1.imshow(image[0].data,vmin=vmin,vmax=vmax,cmap='jet')
    else:
        # Cubes
        im = ax1.imshow(image[0].data[channel],vmin=vmin,vmax=vmax,cmap='jet')
    # Plot parameters, limits, axis, labels ...
    plt.gca().invert_yaxis()
    cbar = plt.colorbar(im, ax=ax1,orientation='vertical')
    cbar.ax.set_ylabel('Flux (image units)', fontsize=15)
    plt.text(0.1,0.1,"Reference", bbox={'facecolor': 'white', 'pad': 10},transform=ax1.transAxes)
    #plt.text(0.1,0.1,"Ref: " + str(ref_image.replace(pathnametodrop,'')), bbox={'facecolor': 'white', 'pad': 10},transform=ax1.transAxes)
    plt.xlabel("X (pixel units)",fontsize=15)
    plt.ylabel("Y (pixel units)",fontsize=15)
    plt.title(" Reference (Chan.# " + str(channel) + ")")

    # Panel #2: Target image at Reference resolution
    ax1 = plt.subplot(grid[0, 1])
    image = fits.open(target_image+".fits")
    # Continuum or cube?
    Ndims = np.shape(np.shape(image[0].data))
    if (Ndims[0] == 2):
        # Continuum
        im = ax1.imshow(image[0].data,vmin=vmin,vmax=vmax,cmap='jet')
    else:
        # Cubes
        im = ax1.imshow(image[0].data[channel],vmin=vmin,vmax=vmax,cmap='jet')
    # Plot parameters, limits, axis, labels ...
    plt.gca().invert_yaxis()
    cbar = plt.colorbar(im, ax=ax1,orientation='vertical')
    cbar.ax.set_ylabel('Flux (image units)', fontsize=15)
    plt.text(0.1,0.1,"Target", bbox={'facecolor': 'white', 'pad': 10},transform=ax1.transAxes)
    #plt.text(0.1,0.1,"Target: " + str(target_image.replace(pathnametodrop,'')), bbox={'facecolor': 'white', 'pad': 10},transform=ax1.transAxes)
    plt.xlabel("X (pixel units)",fontsize=15)
    plt.ylabel("Y (pixel units)",fontsize=15)
    plt.title(" Target at ref. resolution (Chan.# " + str(channel) + ")")

    # Panel #3: A-par map
    ax1 = plt.subplot(grid[1, 0])
    image = fits.open(target_image+"_Apar.fits")
    # Number of axis = Dimentions (Cont vs cubes)
    Ndims = np.shape(np.shape(image[0].data))
    contours = np.array([-1.0,-0.5,-0.25,0.25,0.5,1.0])
    contours = np.arange(-1.,1,0.1)
    # Continuum or cube?
    if (Ndims[0] == 2):
        # Continuum
        im = ax1.imshow(image[0].data,vmin=-1.1,vmax=1.1,cmap='bwr_r')
        cp = ax1.contour(image[0].data,levels=contours,colors="grey")
        ax1.clabel(cp,fontsize=10,colors="grey",fmt="%.2f",inline=1)
    else:
        # Cubes
        im = ax1.imshow(image[0].data[channel],vmin=-1.1,vmax=1.1,cmap='bwr_r')
        cp = ax1.contour(image[0].data[channel],levels=contours,colors="grey")
        ax1.clabel(cp,fontsize=10,colors="grey",fmt="%.2f")
    # Plot parameters, limits, axis, labels ...
    cbar = plt.colorbar(im, ax=ax1,orientation='vertical')
    cbar.ax.set_ylabel('Accuracy parameter', fontsize=15)
    plt.gca().invert_yaxis()
    #plt.show()
    plt.xlabel("X (pixel units)",fontsize=15)
    plt.ylabel("Y (pixel units)",fontsize=15)
    plt.title(" Accuracy map (Chan.# " + str(channel) + ")")

    # Panel #4: Histogram
    ax2 = plt.subplot(grid[1, 1])
    nchans, b, mids, h = get_ALLvalues(FITSfile=target_image+"_Apar.fits",xmin=-1.525,xmax=1.525,xstep=0.05)
    plt.plot(mids,h,label="ALL pixels",linewidth=3,c="red")
    # Mean value
    meanvalue = np.round(np.average(mids,weights=h),2)
    sigmavalue = np.round(np.sqrt(np.cov(mids, aweights=h)),2)
    plt.vlines(meanvalue,np.min(h[h>0]),np.max(h),linestyle="dotted",color="red",linewidth=3,label="A = "+str(meanvalue)+ " +/- " + str(sigmavalue),alpha=1.,zorder=-2)
    # Print results on screen
    print(" Accuracy = " + str(meanvalue) + " +/- " + str(sigmavalue))
    # Continuum or cube
    if (Ndims[0] > 2): # cubes only
        nchans_chan, b_chan, mids_chan, h_chan = get_CHANvalues(FITSfile=target_image+"_Apar.fits",xmin=-1.525,xmax=1.525,xstep=0.05,channel=channel)
        plt.plot(mids_chan,h_chan,label="Channel # " +str(channel),c="blue",linewidth=3,linestyle="dotted")
    # Plot limits, labels, axis...
    plt.vlines(0.,np.min(h[h>0]),np.max(h),linestyle="--",color="black",linewidth=3,label="Goal",alpha=1.,zorder=-2)
    plt.xlim(-1.1,1.1)
    plt.yscale('log')   # Make y axis in log scale
    plt.legend(loc="lower right")
    plt.xlabel("Accuracy parameter",fontsize=20)
    plt.ylabel(r'Number of pixels',fontsize=20)
    # Save plot?
    if save == True:
        #shortname=target_image.replace(pathnametodrop,'').replace('.image','')
        #plt.savefig('Accuracy_map_'+shortname+'.png')
        #print(' See results: Accuracy_map_'+shortname+'.png')
        if plotname == '':
            plotname="Accuracy_map_tmp"
        plt.savefig(plotname+'.png')
        print(" See results: "+plotname+".png")
        
        
        
        
        plt.close()
    # out
    print("---------------------------------------------")
    return True


def show_Fidelity_map(ref_image,target_image,
                       #pathnametodrop,
                       channel=0, save=False, plotname='',
                       labelname='', titlename=''
                       ):
    # Figure
    fig = plt.figure(figsize=(10,15))
    #fig = plt.figure(figsize=(15,10))
    if titlename=='':   
        fig.suptitle('Fidelity map', fontsize=16)
        #fig.suptitle('Fidelity map for \n'+target_image.replace(pathnametodrop,''), fontsize=16)
    else:    
        fig.suptitle(titlename, fontsize=16)
    grid = plt.GridSpec(ncols=2,nrows=2, wspace=0.5, hspace=0.3)
    #grid = plt.GridSpec(ncols=2,nrows=2, wspace=0.3, hspace=0.3)

    # Panel #1: Reference
    ax1 = plt.subplot(grid[0, 0])
    image = fits.open(ref_image+".fits")
    # get min/max values
    #vmin , vmax = np.min(image[0].data[-np.isnan(image[0].data)]), np.max(image[0].data[-np.isnan(image[0].data)])
    vmin , vmax = np.min(image[0].data[~np.isnan(image[0].data)]), np.max(image[0].data[~np.isnan(image[0].data)])
    # Continuum or cube?
    Ndims = np.shape(np.shape(image[0].data))
    channel=int(channel)
    if (Ndims[0] == 2):
        # Continuum
        im = ax1.imshow(image[0].data,vmin=vmin,vmax=vmax,cmap='jet')
    else:
        # Cubes
        im = ax1.imshow(image[0].data[channel],vmin=vmin,vmax=vmax,cmap='jet')
    # Plot parameters, limits, axis, labels ...
    plt.gca().invert_yaxis()
    cbar = plt.colorbar(im, ax=ax1,orientation='vertical')
    cbar.ax.set_ylabel('Flux (image units)', fontsize=15)
    plt.text(0.1,0.1,"Reference", bbox={'facecolor': 'white', 'pad': 10},transform=ax1.transAxes)
    #plt.text(0.1,0.1,"Ref: " + str(ref_image.replace(pathnametodrop,'')), bbox={'facecolor': 'white', 'pad': 10},transform=ax1.transAxes)
    plt.xlabel("X (pixel units)",fontsize=15)
    plt.ylabel("Y (pixel units)",fontsize=15)
    plt.title(" Reference (Chan.# " + str(channel) + ")")

    # Panel #2: Target image at Reference resolution
    ax1 = plt.subplot(grid[0, 1])
    image = fits.open(target_image+".fits")
    # Continuum or cube?
    Ndims = np.shape(np.shape(image[0].data))
    if (Ndims[0] == 2):
        # Continuum
        im = ax1.imshow(image[0].data,vmin=vmin,vmax=vmax,cmap='jet')
    else:
        # Cubes
        im = ax1.imshow(image[0].data[channel],vmin=vmin,vmax=vmax,cmap='jet')
    # Plot parameters, limits, axis, labels ...
    plt.gca().invert_yaxis()
    cbar = plt.colorbar(im, ax=ax1,orientation='vertical')
    cbar.ax.set_ylabel('Flux (image units)', fontsize=15)
    plt.text(0.1,0.1,"Target", bbox={'facecolor': 'white', 'pad': 10},transform=ax1.transAxes)
    #plt.text(0.1,0.1,"Target: " + str(target_image.replace(pathnametodrop,'')), bbox={'facecolor': 'white', 'pad': 10},transform=ax1.transAxes)
    plt.xlabel("X (pixel units)",fontsize=15)
    plt.ylabel("Y (pixel units)",fontsize=15)
    plt.title(" Target at ref. resolution (Chan.# " + str(channel) + ")")

    # Panel #3: Fidelity map
    ax1 = plt.subplot(grid[1, 0])
    image = fits.open(target_image+"_Fidelity.fits")
    # Number of axis = Dimentions (Cont vs cubes)
    Ndims = np.shape(np.shape(image[0].data))
    # Continuum or cube?
    if (Ndims[0] == 2):
        # Continuum
        im = ax1.imshow(image[0].data,vmin=1.,vmax=100,cmap='hot',norm=LogNorm())
    else:
        # Cubes
        im = ax1.imshow(image[0].data[channel],vmin=1.,vmax=100,cmap='hot',norm=LogNorm())
    # Plot parameters, limits, axis, labels ...
    cbar = plt.colorbar(im, ax=ax1,orientation='vertical')
    cbar.ax.set_ylabel('Fidelity', fontsize=15)
    plt.gca().invert_yaxis()
    plt.show()
    plt.xlabel("X (pixel units)",fontsize=15)
    plt.ylabel("Y (pixel units)",fontsize=15)
    plt.title(" Fidelity map (Chan.# " + str(channel) + ")")

    # Panel #4: Histogram
    ax2 = plt.subplot(grid[1, 1])
    nchans, b, mids, h = get_ALLvalues(FITSfile=target_image+"_Fidelity.fits",xmin=0,xmax=100,xstep=0.5)
    plt.plot(mids,h,label="ALL pixels",linewidth=3,c="red")
    # Continuum or cube
    if (Ndims[0] > 2): # cubes only
        nchans_chan, b_chan, mids_chan, h_chan = get_CHANvalues(FITSfile=target_image+"_Fidelity.fits",xmin=0,xmax=100,xstep=0.5,channel=channel)
        plt.plot(mids_chan,h_chan,label="Channel # " +str(channel),c="blue",linewidth=3,linestyle="dotted")
    # Plot limits, labels, axis...
    plt.xlim(1,100)
    plt.xscale('log')
    plt.yscale('log')   # Make y axis in log scale
    plt.legend()
    plt.xlabel("Fidelity",fontsize=20)
    plt.ylabel(r'Number of pixels',fontsize=20)
    # Save plot?
    if save == True:
        #shortname=target_image.replace(pathnametodrop,'').replace('.image','')
        #plt.savefig('Fidelity_map_'+shortname+'.png')
        #print(' See results: Fidelity_map_'+shortname+'.png')
        if plotname == '':
            plotname="Fidelity_map_tmp"
        plt.savefig(plotname+'.png')        
        plt.close()
    # out
    print("---------------------------------------------")
    return True


# Accuracy parameter comparisons
def Compare_Apar(ref_image = '',target_image=[''],
                  #pathnametodrop = '', 
                  save=False, plotname='', 
                  labelname=[''], titlename=''):
    # Reference image
    print("=============================================")
    print(" Accuracy parameter: comparisons")
    print(" Reference : "+str(ref_image))
    flux0 = np.round(imstat(ref_image)["flux"][0])
    print(" Total Flux = " + str(flux0) + " Jy")
    print("---------------------------------------------")
    # Number of plots
    Nplots = np.shape(target_image)[0]
    # Global comparisons 
    plt.figure(figsize=(8,11))
    grid = plt.GridSpec(ncols=1,nrows=5, wspace=0.3, hspace=0.3)
    ax1 = plt.subplot(grid[0:4, 0])
    # Loop over all images
    for m in np.arange(Nplots):
        # Get total flux in image
        flux = np.round(imstat(target_image[m]+".fits")["flux"][0])
        # Extract values from file
        nchans, b, mids, h = get_ALLvalues(FITSfile=target_image[m]+"_Apar.fits",xmin=-1.525,xmax=1.525,xstep=0.05)
        # Get mean and std
        meanvalue = np.round(np.average(mids,weights=h),2)
        sigmavalue = np.round(np.sqrt(np.cov(mids, aweights=h)),2)
        # Get skewness and kurtosis of the Apar image
        hdu = fits.open(target_image[m]+"_Apar.fits")
        Adist = hdu[0].data.flatten()
        Adist = Adist[(Adist <= 10.) & (Adist >= -10.)] # remove big outlayers
        skewness = np.round(skew(Adist),3)
        kurt = np.round(kurtosis(Adist),3)          
        ax1.plot(mids,h,label=target_image[m] + "; A = "+ str(meanvalue) + " +/- " + str(sigmavalue),linewidth=3,c=IQA_colours[m])
        # Print results on screen
        print(" Target image " + str(m+1) + " : " + str(target_image[m]))
        print(" Total Flux = " + str(flux) + " Jy ("+str(np.round(flux/flux0,2))+"\%)")
        print(" Accuracy:")
        print(" Mean +/- Std. = " + str(meanvalue) + " +/- " + str(sigmavalue))
        print(" Skewness, Kurtosis = " + str(skewness) + " , " + str(kurt) )
        print("................................................")
    # Add Goal line
    plt.vlines(0.,np.min(h[h>0]),np.max(h),linestyle="--",color="black",linewidth=3,label="Goal",alpha=1.,zorder=-2)
    # Plot limits, legend, labels...
    plt.xlim(-1.5,1.5)
    plt.yscale('log')   # Make y axis in log scale
    #plt.legend(loc='lower right')
    plt.legend(bbox_to_anchor=(0.5, -0.1),loc='upper center', borderaxespad=0.)
    plt.xlabel("Accuracy",fontsize=20)
    plt.ylabel(r'# pixels',fontsize=20)
    if titlename=='':
        plt.title("Accuracy Parameter: comparisons",fontsize=16)
    else: 
        plt.title(titlename,fontsize=16)

    # Save plot?
    if save == True:
        if plotname == '':
            plotname="AparALL_tmp"
        plt.savefig(plotname+'.png')
        print(" See results: "+plotname+".png")
        plt.close()
    # out
    print("---------------------------------------------")
    print(" Accuracy parameter comparisons... DONE")
    print("=============================================")
    return True
    # 
    hdu[0].data[mask] = np.nan
    # Write to file
    fits.writeto(noisefile+".fits",data=hdu[0].data,header=hdu[0].header,overwrite=True)


# Fidelity comparisons
def Compare_Fidelity(ref_image = '',target_image=[''],
                    save=False, plotname='', 
                    labelname=[''], titlename=''):
    """
    Compare all Fidelity images (continuum or mom0 maps) (A. Hacar, Univ. of Vienna)
    Arguments:
      ref_image - image used as reference
      target_image - list of images to be compared with reference
      save - save plot? (default = False)
    Requires:
      The script will look for target_image[i]_convo2ref_Fidelity.fits images produced by the get_IQA() script
    
    Results:
      1- Histogram including the Fidelity distributions for all input images
      2- Numerical results: Total flux in the image, mean + std + kurtosis + skewness of each histogram 
    Example:
      Compare_Fidelity(ref_image = 'TP_image',target_image=['Feather.image','TP2vis.image'])
    """
    # Reference image
    print("=============================================")
    print(" Fidelity comparisons:")
    print(" Reference : "+str(ref_image))
    print("---------------------------------------------")
    # Number of plots
    Nplots = np.shape(target_image)[0]
    # Global comparisons 
    plt.figure(figsize=(8,11))
    grid = plt.GridSpec(ncols=1,nrows=5, wspace=0.3, hspace=0.3)
    ax1 = plt.subplot(grid[0:4, 0])
    for m in np.arange(Nplots):
        print(" Target image " + str(m+1) + " : " + str(target_image[m]))
        nchans, b, mids, h = get_ALLvalues(FITSfile=target_image[m]+"_convo2ref_Fidelity.fits",xmin=0.,xmax=100.,xstep=0.5)
        # Calculate mean value
        hdu = fits.open(target_image[m]+"_convo2ref_Fidelity.fits")
        Fdist = hdu[0].data.flatten()
        Fdist = Fdist[(Fdist < 100.) & (Fdist > 0)]
        meanvalue = np.round(np.mean(Fdist),1)
        medianvalue = np.round(np.median(Fdist),1)
        q1value = np.round(np.percentile(Fdist, 25),1)  # Quartile 1st
        q3value = np.round(np.percentile(Fdist, 75),1)  # Quartile 3rd
        #meanvalue = np.round(np.average(mids,weights=h),1)
        if labelname[m]=='':
            plt.plot(mids,h,label=target_image[m] + "; <Fidelity> = "+ str(meanvalue),linewidth=3,c=IQA_colours[m])
        else:    
            plt.plot(mids,h,label=labelname[m] + "; <Fidelity> = "+ str(meanvalue),linewidth=3,c=IQA_colours[m])
        # Display on screen
        print(" Fidelity")
        print("  Mean = " + str(meanvalue))
        print("  [Q1,Median,Q3] = ["+str(q1value)+" , "+ str(medianvalue)+" , "+str(q3value)+"]")
    # plot lims, axis, labels, etc...
    plt.xlim(1,100.)
    plt.xscale('log')
    plt.yscale('log')   # Make y axis in log scale
    #plt.ylim(1,)
    #plt.legend(loc="lower left")
    plt.legend(bbox_to_anchor=(0.5, -0.1),loc='upper center', borderaxespad=0.)
    plt.xlabel("Fidelity",fontsize=20)
    plt.ylabel(r'# pixels',fontsize=20)
    if titlename=='':
        plt.title("Fidelity Comparisons",fontsize=16)
    else:
        plt.title(titlename,fontsize=16)
    if save == True:
        if plotname == '':
            plotname="FidelityALL_tmp"
        plt.savefig(plotname+'.png')
        print(" See results: "+plotname+".png")
        plt.close()
    print("---------------------------------------------")
    print(" Fidelity comparisons... DONE")
    print("=============================================")




def Compare_Apar_signal(ref_image = '',target_image=[''],
            #pathnametodrop = '', 
            save=False, noise=0.0, plotname='', 
            labelname=[''], titlename=''
            ):
    """
    Compare_Apar_signal (A. Hacar, Univ. of Vienna) 
    
    Compare all Apar images vs signal (continuum or mom0 maps).
    If No. of targets = 1, the mean and std of A-par will be calculated.
    This function can be applied in both cont/mom0 and cubes FITS files.
    
    Arguments:
      ref_image - image used as reference
      target_image - list of images to be compared with reference 
        (recommended to <= 4 targets)
      save - (optional) save plot? (default = False)
      noise - (optional) if noise > 0.0 the evolution of the noise level
        will be displayed  
    
    Requires:
      The script will look for target_image[i]_convo2ref_Apar.fits images produced by the get_IQA() script
    
    Results:
      Apar as function of reference & target signals
    
    Example 1: compare a list of targets
      Compare_Apar_signal(ref_image = 'TP_image',target_image=['Feather.image','TP2vis.image'])
    
    Example 2: investigate singel target (incl. A-par statistics)
      Compare_Apar_signal(ref_image = 'TP_image',target_image=['Feather.image'],noise=0.5)
    
    """
    # Reference image
    print("=============================================")
    print(" A-par vs Signal")
    print("---------------------------------------------")
    # Number of plots
    #Nplots = np.shape(target_image)[0]
    Nplots=1
    # These plots get too crowded with a high number of targets. If No. targets > 4 then exit this function
    #if (Nplots > 4):
     #   print(" Too many targets. Please use <= 4.")
      #  print(" No plot shown.")
       # print("---------------------------------------------")
        #print(" A-par vs Signal... DONE")
        #print("=============================================")
        #return None

    # Figure parameters 
    plt.figure(figsize=(8,14))
    grid = plt.GridSpec(ncols=1,nrows=7, wspace=0.3, hspace=0.3)
    
    # Plot #1: Reference vs A-par
    ax0 = plt.subplot(grid[0:2, 0])
    # Loop over all images
    xmax0 = 0.0; xmin0 = 1E6    # Dummy values
    for m in np.arange(Nplots):
        # Images
        im1 = fits.open(ref_image+".fits")
        im2 = fits.open(target_image +"_Apar"+".fits")

        # Define plot limits
        xmin = np.min(im1[0].data[np.isnan(im1[0].data)==False])
        xmax = np.max(im1[0].data[np.isnan(im1[0].data)==False])
        if (xmax > xmax0):
            xmax0=xmax+xmax/10. # Slightly larger
        if (xmin < xmin0):       
            xmin0=xmin
        if (xmin < 0.0):          #Lydia's modification to avoid negative values!
            xmin0=0.0001
        # Plot results
        ax0.scatter(im1[0].data,im2[0].data,c=IQA_colours[m],marker="o",rasterized=True,edgecolor='none')

    # Goal (A-par = 0)
    ax0.hlines(0.,xmin,xmax0,linestyle="--",color="black",linewidth=3,alpha=1.,zorder=2)

    # Calculate mean and sigma if there is only one target
    if (Nplots == 1):
        print("---------------------------------------------")
        print(" A-par values per bin: ")
        # Calculate bins & step in log-scale
        steplog=(np.log10(xmax0)-np.log10(xmin0))/10.
        xvalueslog=np.arange(np.log10(xmin0),np.log10(xmax0),steplog)
        # back to linear scale
        step=10.**steplog
        xvalues=10.**xvalueslog
        # Define stats vectors
        means=np.zeros(len(xvalues))    # Mean
        stds=np.zeros(len(xvalues)) # STD
        medians=np.zeros(len(xvalues))  # Median
        q1values=np.zeros(len(xvalues)) # Q1
        q3values=np.zeros(len(xvalues)) # Q2
  
        # helpers for debugging
        #print(xmin, xmax)        
        #print(xmin0, xmax0)
        #print(xvalues)

        count=0
        for j in xvalueslog:
            # Define bin ranges in log-space
            idx = (im1[0].data >= 10.**j) & (im1[0].data < 10.**(j+steplog)) & (np.isnan(im1[0].data)==False) & (np.isfinite(im1[0].data)==True)
            values = im2[0].data[idx]
            values = values[ (np.isnan(values)==False) & (np.isfinite(values)==True) ] # remove Nan & Inf.
            # Stats
            if (np.shape(values)[0] > 0):
                means[count] = np.mean(values)  # Mean
                stds[count] = np.std(values)    # STD
                medians[count] = np.median(values)  # Median
                q1values[count] = np.percentile(values, 10) # 10% Quartile
                q3values[count] = np.percentile(values, 90) # 90% Quartile
            # Show results on screen
            print("Bin "+str(count+1)+": Ref.Flux = " + str(np.round(10.**(j+steplog/2.),2)) + " ; A = " + str(np.round(means[count],2)) + " +/- " + str(np.round(stds[count],2)) + " ; [Q10,Q90] = ["+ str(np.round(q1values[count],3)) + " , " + str(np.round(q3values[count],3)) +"]")
            # Counter +1
            count+=1
            #
        # Display mean and STD
        ax0.errorbar(10.**(xvalueslog+steplog/2.),means, yerr=stds, fmt='o',c="blue",label=r"|y|$\pm 1 \sigma$ ",linewidth=2,markersize=10,zorder=2,capsize=5)
        #ax0.errorbar(10.**(xvalueslog+steplog/2.),medians, yerr=[q1values,q3values], fmt='o',c="cyan",label=r"[Q1,Median,Q3]",linewidth=2)
        ax0.vlines(10.**(xvalueslog+steplog/2.),q1values,q3values,color="cyan",label=r"[Q10,Q90]",linewidth=5,zorder=1)

    # Show noise effects?
    if (noise > 0):
        xvalues=np.arange(np.log10(xmin0),np.log10(xmax0),(np.log10(xmax0)-np.log10(xmin0))/20.)
        xvalues=10.**xvalues
        ax0.plot(xvalues,noise/np.abs(xvalues),c="blue",zorder=2,linewidth=4,linestyle="dotted")
        ax0.plot(xvalues,-noise/np.abs(xvalues),c="blue",zorder=2,linewidth=4,linestyle="dotted")

    # Plot limits, legend, labels...
    ax0.legend()
    ax0.set_yticks(np.arange(-2.,2.,0.25))
    ax0.set_xlim(xmin0,xmax0)
    ax0.set_ylim(-0.65, 0.65)
    # Adjust ylims if the results are really bad!
    if (np.mean(im2[0].data[np.isnan(im2[0].data)==False]) <= -0.5):
        ax0.set_ylim(-1.5, 0.5)
    ax0.set_xscale('log')
    ax0.set_ylabel(r" A-par",fontsize=20)
    if titlename=='':
        plt.title("Accuracy vs. Signal",fontsize=16)
    else: 
        plt.title(titlename,fontsize=16)


    # Plot #2: Reference vs Target
    ax1 = plt.subplot(grid[2:6, 0])
    # Loop over all images
    xmax0 = 0.0; xmin0 = 1E6    # Dummy values
    for m in np.arange(Nplots):
        # Images
        im1 = fits.open(ref_image+".fits")
        im2 = fits.open(target_image+".fits")
        # Define plot limits
        xmin = np.min(im1[0].data[np.isnan(im1[0].data)==False])
        xmax = np.max(im1[0].data[np.isnan(im1[0].data)==False])
        if (xmax > xmax0):
            xmax0=xmax+xmax/10. # Slightly larger
        if (xmin < xmin0):
            xmin0=xmin
        if (xmin < 0.0):          #Lydia's modification to avoid negative values!
            xmin0=0.0001

        # Plot results
        if labelname[m]=='':
            ax1.scatter(im1[0].data,im2[0].data,c=IQA_colours[m],marker="o",rasterized=True,label=target_image[m],edgecolor='none')
        else:    
            ax1.scatter(im1[0].data,im2[0].data,c=IQA_colours[m],marker="o",rasterized=True,label=labelname[m],edgecolor='none')

    # Show A values lines
    xvalues=np.arange(xmin0,xmax0,(xmax0-xmin0)/20.)
    plt.plot(xvalues,xvalues,c="k",zorder=2,linewidth=3,linestyle="--",label="Goal (linear correlation; A-par = 0.0)")
    plt.text((xmax0-xmin0)/3.,(xmax0-xmin0)/3.,"A=0",rotation=45,ha='center',va='center',rotation_mode="anchor",bbox=dict(boxstyle='square',facecolor='white', edgecolor='black'))
    # Note that the value of A=-1 needs values of Target=0, which is not allowed in ylog-plots
    for k in np.array([-0.75,-0.5,-0.25,0.25,0.5,0.75,1.0]):
        def Avalues(A,x):
            return A*x+x
        yvalues=Avalues(A=k,x=xvalues)
        ax1.plot(xvalues,yvalues,c="grey",zorder=2,linestyle="dashed",alpha=0.5)
        ax1.text((xmax0-xmin0)/3.,Avalues(A=k,x=(xmax0-xmin0)/3.),"A="+str(k),rotation=45,ha='center',va='center',rotation_mode="anchor",clip_on=True,size=10.,color="grey",zorder=2)

    # Show noise effects?
    if (noise > 0):
        xvalues=np.arange(np.log10(xmin0),np.log10(xmax0),(np.log10(xmax0)-np.log10(xmin0))/20.)
        xvalues=10.**xvalues
        plt.plot(xvalues,xvalues-noise,c="blue",zorder=2,linewidth=4,linestyle="dotted",label=r"White noise:  $\sigma = $"+str(np.round(noise,2))+" (image units)")
        plt.plot(xvalues,xvalues+noise,c="blue",zorder=2,linewidth=4,linestyle="dotted")

    # Plot limits, legend, labels...
    plt.xlim(xmin0,xmax0)
    plt.ylim(xmin0,xmax0)
    plt.xscale('log')
    plt.yscale('log')
    # Legend and labels
    plt.legend(bbox_to_anchor=(0.5, -0.15),loc='upper center', borderaxespad=0.)
    plt.ylabel(r" Target flux (image units)",fontsize=20)
    plt.xlabel(r" Reference flux (image units)",fontsize=20)


    # Save plot?
    if save == True:
        if plotname == '':
            plotname="Apar_signal_ALL_tmp"
        plt.savefig(plotname+'.png')
        print(" See results: "+plotname+".png")
        plt.close()
    # out
    print("---------------------------------------------")
    print(" A-par vs Signal... DONE")
    print("=============================================")
    return True


def Compare_Fidelity_signal(ref_image = '',target_image=[''],
             save=False,noise=0.0, plotname='', 
             labelname=[''], titlename=''
             ):
    """
    Compare all Fidelity images vs signal (continuum or mom0 maps) (A. Hacar, Univ. of Vienna)
    Arguments:
      ref_image - image used as reference
      target_image - list of images to be compared with reference (better only one)
      save - save plot? (default = False)
    Requires:
      The script will look for target_image[i]_convo2ref_Apar.fits images produced by the get_IQA() script
    
    Results:
      Fidelity vs signal in the original image
    Example:
      Compare_Fidelity_signal(ref_image = 'TP_image',target_image=['Feather.image','TP2vis.image'])
    """
    # Reference image
    print("=============================================")
    print(" Fidelity vs Signal")
    print("---------------------------------------------")
    # Number of plots
    Nplots = np.shape(target_image)[0]
    # Global comparisons 
    plt.figure(figsize=(8,11))
    grid = plt.GridSpec(ncols=1,nrows=5, wspace=0.3, hspace=0.3)
    ax1 = plt.subplot(grid[0:4, 0])
    # Loop over all images
    xmin0 = 100; xmax0 = 0.
    for m in np.arange(Nplots):
        # Images
        ##im1 = fits.open(target_image[m]+"_convo2ref.fits")    # signal
        im1 = fits.open(ref_image+".fits")
        im2 = fits.open(target_image + ".fits") # parameter
        xmin = np.min(im1[0].data[np.isnan(im1[0].data)==False])
        xmax = np.max(im1[0].data[np.isnan(im1[0].data)==False])
        ymin = np.min(im2[0].data[np.isnan(im2[0].data)==False])
        ymax = np.max(im2[0].data[np.isnan(im2[0].data)==False])
        if (xmax > xmax0):
            xmax0=xmax
        if (xmin < xmin0):
            xmin0=xmin
        if (xmin < 0.0):          #Lydia's modification to avoid negative values!
            xmin0=0.0001
        # Plot results
        if labelname[m]=='':            ax1.scatter(im1[0].data,im2[0].data,c=IQA_colours[m],marker="o",rasterized=True,label=target_image[m],edgecolor='none')
        else:    
            ax1.scatter(im1[0].data,im2[0].data,c=IQA_colours[m],marker="o",rasterized=True,label=labelname[m],edgecolor='none')
    # Plot limits, legend, labels...
    plt.xlim(xmin0,xmax0+xmax0/5)
    plt.ylim(xmin0,xmax0+xmax0/5)
    plt.xscale('log')
    plt.yscale('log')   
    # Get A values lines
    xvalues=np.arange(xmin0,xmax0,(xmax0-xmin0)/20.)
    plt.plot(xvalues,xvalues,linestyle="dashed",c="k",label="Goal (linear correlation)",zorder=2,linewidth=3)
    # Fidelities
    for k in np.array([2.,5.,10.]):
        def Fvalues(F,x):
            return F*x/(F+1.)
        yvalues=Fvalues(F=k,x=xvalues)
        plt.plot(xvalues,yvalues,c="grey",zorder=2,linestyle="dashed")
        plt.text((xmax0-xmin0)/2.,Fvalues(F=k,x=(xmax0-xmin0)/2.),"Fid.="+str(k),rotation=45,ha='center',va='center',rotation_mode="anchor",bbox=dict(boxstyle='square',facecolor='white', edgecolor='black'))
        def Fvalues(F,x):
            return F*x/(F-1.)
        yvalues=Fvalues(F=k,x=xvalues)
        plt.plot(xvalues,yvalues,c="grey",zorder=2,linestyle="dashed")
        plt.text((xmax0-xmin0)/2.,Fvalues(F=k,x=(xmax0-xmin0)/2.),"Fid.="+str(k),rotation=45,ha='center',va='center',rotation_mode="anchor",bbox=dict(boxstyle='square',facecolor='white', edgecolor='black'))
    #plt.legend(loc='lower right')
    plt.legend(bbox_to_anchor=(0.5, -0.1),loc='upper center', borderaxespad=0.)
    plt.ylabel(r" Target flux (image units)",fontsize=20)
    plt.xlabel(r'# Reference flux (image units)',fontsize=20)
    if titlename=='':
        plt.title("Fidelity vs. Signal",fontsize=16)
    else: 
        plt.title(titlename,fontsize=16)
    
    # Save plot?
    if save == True:
        if plotname == '':
            plotname="Fidelity_signal_ALL_tmp"
        plt.savefig(plotname+'.png')
        print(" See results: "+plotname+".png")
        plt.close()
    # out
    print("---------------------------------------------")
    print(" Fidelity vs Signal... DONE")
    print("=============================================")
    return True



##------------------------------------------------------------
##Power spectrum

## define linear fitting function
def linFunc(x, slope, b):
    return x*slope+b

## define function to compute MAD
MAD = lambda x: np.median(np.abs(x - np.median(x)))

## define function to return 1D SPS
def compute_1D_SPS(image):
    ## mask out image NaN's
    nan_inds = np.isnan(image)
    image[nan_inds] = 0.

    ## compute 2D power spectrum
    modulusImage = np.abs(np.fft.fftshift(np.fft.fft2(image)))**2

    ## obtain center pixel coordinates (where the DC power component is located)
    center = np.where(modulusImage == np.nanmax(modulusImage))
    center = [center[0][0], center[1][0]]


    ## define pixel grid
    y = np.arange(-center[0], modulusImage.shape[0] - center[0])
    x = np.arange(-center[1], modulusImage.shape[1] - center[1])
    yy, xx = np.meshgrid(y, x, indexing='ij')
    dists = np.sqrt(yy**2 + xx**2)
    
    ## define spatial frequency grid
    yfreqs = np.fft.fftshift(np.fft.fftfreq(modulusImage.shape[0]))
    xfreqs = np.fft.fftshift(np.fft.fftfreq(modulusImage.shape[1]))
    yy_freq, xx_freq = np.meshgrid(yfreqs, xfreqs, indexing='ij')
    freqs_dist = np.sqrt(yy_freq**2 + xx_freq**2)
    zero_freq_val = freqs_dist[np.nonzero(freqs_dist)].min() / 2.
    freqs_dist[freqs_dist == 0] = zero_freq_val
    
    ## define bin spacing
    max_bin = 0.5
    min_bin = 1.0 / np.min(modulusImage.shape)
    nbins = int(np.round(dists.max()) + 1)
    bins = np.linspace(min_bin, max_bin, nbins + 1)
    finite_mask = np.isfinite(modulusImage)
    ## compute the radial profile using median values within each annulus
    ps1D, bin_edges, cts = stats.binned_statistic(freqs_dist[finite_mask].ravel(), modulusImage[finite_mask].ravel(), bins=bins, statistic='median')
    ## compute MAD uncertainties
    ps1D_MADErrs, bin_edges, cts = stats.binned_statistic(freqs_dist[finite_mask].ravel(), modulusImage[finite_mask].ravel(), bins=bins, statistic=MAD)
    bin_cents = (bin_edges[1:] + bin_edges[:-1]) / 2.
    ## return profile & errors (in linear space) and bin centers for x-axis
    return ps1D, ps1D_MADErrs, bin_cents*1/u.pix


def genmultisps(fitsimages, save=False, plotname='', labelname='',
                titlename=''):
    """
    gensps
    Script to plot the power spectra of a user provided general 2D FITS file
    Arguments: 
      fitsimages - list of names of input images
      save      - save plot? (default False)
    """

    matplotlib.rc('font', family='sans-serif') 
    matplotlib.rc('font', serif='Helvetica Neue') 
    matplotlib.rc('text', usetex='false') 

    
    # Note: we use new IQA_colours instead

    if type(fitsimages) != list or len(fitsimages)==0:
        print("ERROR: fitsimages needs to be non-empty list")
        return False


    # initialise plotting
    fig, ax = pyplot.subplots(figsize = (8,8))

    ## set plotting parameters
    majorYLocFactor = 2
    minorYLocFactor = majorYLocFactor/2.
    majorXLocFactor = 0.5
    minorXLocFactor = majorXLocFactor/4.
    
    majorYLocator = MultipleLocator(majorYLocFactor)
    majorYFormatter = FormatStrFormatter('%d')
    minorYLocator = MultipleLocator(minorYLocFactor)
    majorXLocator = MultipleLocator(majorXLocFactor)
    majorXFormatter = FormatStrFormatter('%.1f')
    minorXLocator = MultipleLocator(minorXLocFactor)

    ## set tick mark parameters
    ax.yaxis.set_major_locator(majorYLocator)
    ax.yaxis.set_major_formatter(majorYFormatter)
    ax.yaxis.set_minor_locator(minorYLocator)
    ax.xaxis.set_major_locator(majorXLocator)
    ax.xaxis.set_major_formatter(majorXFormatter)
    ax.xaxis.set_minor_locator(minorXLocator)


    imnumber = 0
    for fitsimage in fitsimages:

        if type(fitsimage) != str or len(fitsimage)==0:
            print('ERROR: empty fitsimage name')
            return False

        fitsName = fitsimage

        ## open the header data unit
        hdu = fits.open(fitsName)

        ## open up 2D image
        image = hdu[0].data

        ## set useful variables
        pixUnits = False
        noBeamInfo = False
        fit = False

        ## check for degenerate stokes axes and remove
        if len(image.shape) == 4:
            image = image[0, 0, :, :]
        elif len(image.shape) == 3:
            image = image[0, :, :]
        else:
            print('No degenerate axes...')

        ## check for relevant header info
        try: 
            pixSize = hdu[0].header['CDELT2']*3600 *u.arcsec
        except (KeyError):
            pixUnits = True
            print('Missing pixel size in header, we will continue with units of pixels')

        try:
            bmaj = hdu[0].header['BMAJ']*3600 * u.arcsec
            bmin = hdu[0].header['BMIN']*3600. * u.arcsec 
        except (KeyError):
            noBeamInfo = True

            print('Missing beam information in header, we will apply no spatial frequency cut-off')

        if not noBeamInfo:
            print('pixSize bmaj bmin ', pixSize, bmaj, bmin)

            ## convert image to Jy/pixel 
            image*=(pixSize.value)**2/(1.1331*bmaj.value*bmin.value)

        ## compute power spectrum
        SPS_1D, MADErrs, spatialFreqs = compute_1D_SPS(image)
        

        ## compute power spectrum
        if pixUnits == True and noBeamInfo == False: ## use pixel units; however, set high spatial frequency cutoff from BMAJ
            angularScales = spatialFreqs**(-1)*u.pix

            ## determine low/high cutoffs indices (where angular scales < 0.5*map/3 & angular scales > BMAJ)
            lowCut_idx = np.where(angularScales < angularScales[0]/3)[0][0]
            bmaj_pix = bmaj/pixSize*u.pix/u.arcsec
            highCut_idx = np.where(angularScales > bmaj_pix)[0][-1]
        elif noBeamInfo == True and pixUnits == True: ## use pixel units with no cutoff from beam (i.e., simulated image)
            angularScales = spatialFreqs**(-1)

            ## determine low/high cutoffs indices (where angular scales < 0.5*map/3 & high cut is the pixel size
            lowCut_idx = np.where(angularScales < angularScales[0]/3)[0][0]
            highCut_idx = -1
        elif noBeamInfo == True and pixUnits == False: ## use angular units but no cutoff from beam (i.e., simulated image)
            angularScales = spatialFreqs**(-1) * pixSize * (1/u.pix)

            ## determine low/high cutoffs indices (where angular scales < 0.5*map/3 & high cut is the pixel size
            lowCut_idx = np.where(angularScales < angularScales[0]/3)[0][0]
            highCut_idx = -1
        else: # i.e.  noBeamInfo == False and pixUnits == False; use angular units and set high-freq cutoff from beam info
            angularScales = spatialFreqs**(-1) * pixSize * (1/u.pix) ## in arcsec
            ## determine low/high cutoffs indices (where angular scales < 0.5*map/3 & angular scales > BMAJ)
            lowCut_idx = np.where(angularScales < angularScales[0]/3)[0][0]
            highCut_idx = np.where(angularScales > bmaj)[0][-1]

            if lowCut_idx == highCut_idx: ## likely dealing with single-dish/tp image => DON'T FIT
                fit = False
        print('low/high spatial freq. cutoffs [arcsec]: ', angularScales[lowCut_idx].value, angularScales[highCut_idx].value)

        ## convert to log space
        logErrors= MADErrs/(SPS_1D*np.log(10))
        logPwr = np.log10(SPS_1D)
        logAngScales= np.log10(angularScales.value)

        if fit:
            ## fit from large-scale to brk    
            coeffs, matcov = curve_fit(linFunc, logAngScales[lowCut_idx:highCut_idx], logPwr[lowCut_idx:highCut_idx], [1,1])
            perr = np.sqrt(np.diag(matcov))

        ## finally, do the plotting

        ## plot full power spectra ##
        ######ax.plot(logAngScales, logPwr, label='Data '+str(imnumber), color=IQA_colours[imnumber])
        #if labelname[fitsimages.index(fitsimage)]=='':
        #    ax.plot(logAngScales, logPwr, label=str(fitsName), color=IQA_colours[imnumber])
        #else:    
        #    ax.plot(logAngScales, logPwr, label=labelname[fitsimages.index(fitsimage)], color=IQA_colours[imnumber])
        ax.plot(logAngScales, logPwr, label=str(fitsName), color=IQA_colours[imnumber])

        ## fill inbetween for errors ##
        #ax.fill_between(logAngScales, logPwr-logErrors, logPwr+logErrors, color=IQA_colours[1])

        if fit:
            ax.plot(logAngScales[lowCut_idx:highCut_idx], linFunc(logAngScales[lowCut_idx:highCut_idx], coeffs[0], coeffs[1]), color='black', linestyle='--', label = 'PL: %.1f+/-%.1f' % (coeffs[0], perr[0]))
        
        ## set x and y limits ##
        ax.set_xlim(np.max(logAngScales)+0.25, np.min(logAngScales) - 0.25)

        ## show fit limits ##
        ax.axvline(logAngScales[lowCut_idx], color = IQA_colours[imnumber], linestyle = '-.')
        ax.axvline(logAngScales[highCut_idx], color= IQA_colours[imnumber], linestyle = '--')

        imnumber += 1
    # end for


    if pixUnits == True:
        ax.set_xlabel(r'$log_{10}\left(\rm 1/pix\right)$')
    else:
        ax.set_xlabel(r'$log_{10}\left(\rm arcsec\right)$')

    ax.set_ylabel(r'$log_{10}\left(\rm Power\right)$')
    #pyplot.legend(loc='upper right', fontsize=6)
    pyplot.legend(loc='lower left')
    if titlename=='':
        plt.title("Power spectra",fontsize=16)
    else:
        plt.title(titlename,fontsize=16)

    pyplot.show()

    if save == True:
        #pyplot.savefig(fitsimages[0]+'_and_others.sps.pdf')
        if plotname == '':
            plotname="Power_spectra_ALL_tmp"
        plt.savefig(plotname+'.png')
        print(" See results: "+plotname+".png")
        plt.close()


    return True


def convert_JypB_JypP(sdimage):

    """
    convert image brightness unit from Jy/beam to Jy/pixel  (Moser-Fischer, L.)
    a helper for runWSM to prepare the startmodel format
    usemask - masking mode parameter as for tclean
    mask - file name of mask
    pbmask - PB mask cut-off level 
    niter - number of iterations spent on this mask
    """

    myimhead = cta.imhead(sdimage)


    print('Checking SD units...')
    if myimhead['unit']=='Jy/beam': 
        print('SD units {}. OK, will convert to Jy/pixel.'.format(myimhead['unit']))
        ##CHECK: header units
        SingleDishResolutionArcsec = myimhead['restoringbeam']['major']['value'] #in Arcsec
        CellSizeArcsec = abs(myimhead['incr'][0])*206265. #in Arcsec
        toJyPerPix = CellSizeArcsec**2/(1.1331*SingleDishResolutionArcsec**2)
        SDEfficiency = 1.0 #--> Scaling factor
        fluxExpression = "(IM0 * {0:f} / {1:f})".format(toJyPerPix,SDEfficiency)
        #scaled_name = sdimage.split('/')[-1]+'.Jyperpix'
        scaled_name = sdimage+'.Jyperpix'
        
        os.system('rm -rf '+scaled_name)
        cta.immath(imagename=sdimage,
                   outfile=scaled_name,
                   mode='evalexpr',
                   expr=fluxExpression)
        hdval = 'Jy/pixel'
        dummy = cta.imhead(imagename=scaled_name,
                           mode='put',
                           hdkey='BUNIT',
                           hdvalue=hdval)
        ### TO DO: MAY NEED TO REMOVE BLANK 
        ### and/or NEGATIVE PIXELS IN SD OBSERVATIONS
        return scaled_name
    elif myimhead['unit']=='Jy/pixel': 
        print('SD units {}. SKIP conversion. '.format(myimhead['unit']))
        return sdimage
    else: 
        print('SD units {}. NOT OK, needs conversion by user to Jy/beam or Jy/pixel. '.format(myimhead['unit']))
        return sys.exit()

def check_prep_tclean_param(  
                vis,     
                spw, 
                field, 
                specmode,                                 
                imsize, 
                cell, 
                phasecenter,         
                start, 
                width, 
                nchan, 
                restfreq,
                threshold, 
                niter,
                cycleniter,
                usemask,
                sidelobethreshold,
                noisethreshold,
                lownoisethreshold, 
                minbeamfrac,
                growiterations,
                negativethreshold,                
                mask, 
                pbmask,
                interactive,               
                multiscale, 
                maxscale,
                loadmask,
                fniteronusermask
                ):

    """
    check validity of parameters and set up tclean parameters in a uniform manner
    (Moser-Fischer, L.)
    a helper for runsdintimg and runtclean
    Currently, it provides 'cube' and 'mfs' as spectral modes - 'mtmfs'
    might be implemented later.
    
    steps:
    - check 
    
    vis  
    spw - 
    field, 
    specmode,                                 
    imsize, 
    cell, 
    phasecenter,         
    start, 
    width, 
    nchan, 
    restfreq,
    threshold, 
    niter,
    cycleniter,
    usemask,
    sidelobethreshold,
    noisethreshold,
    lownoisethreshold, 
    minbeamfrac,
    growiterations,
    negativethreshold,                
    mask, 
    pbmask,
    interactive,               
    multiscale, 
    maxscale,
    loadmask,
    fniteronusermask      
    """

    # valid specmode?
    if specmode not in ['mfs', 'cube']:
        print('specmode \"'+specmode+'\" is not supported.')
        return sys.exit()


    # valid threshold?
    if not type(threshold) == str or 'Jy' not in threshold and niter>1:
        if not interactive:
            print("You must provide a valid threshold, example '1mJy'")
            return sys.exit()
        else:
            print("You have not set a valid threshold. Please do so in the graphical user interface!")
            threshold = '1mJy'
    
    
    # valid image and cell size?
    if imsize==[] or cell=='':
        cta.casalog.post('You need to provide values for the parameters imsize and cell.', 'SEVERE', origin='runsdintimg')
        return sys.exit()    
    

    if loadmask==True and fniteronusermask>1.0 or fniteronusermask<0.0:
        print('fniteronusermask is out of range: ' +fniteronusermask+' Please choose a value between 0 and 1 (inclusively)')
        return sys.exit()    
    else:
        pass    
    

    #   # specmode, deconvolver and multiscale setup
    #   if multiscale:
    #       if specmode == 'mfs':
    #           mydeconvolver = 'mtmfs'   # needed bc it's the only mfs mode implemented into sdint
    #       elif specmode == 'cube':
    #           mydeconvolver = 'multiscale'
    #           #numchan = nchan           # not really needed here?
    #       mycell = myqa.convert(myqa.quantity(cell),'arcsec')['value']
    #       myscales = [0]
    #       for i in range(0, int(math.log(maxscale/mycell,3))):
    #           myscales.append(3**i*5)
    #   
    #       print("My scales (units of pixels): "+str(myscales))
    #   
    #   else:    
    #       myscales = [0]
    #       if specmode == 'mfs':
    #           mydeconvolver = 'mtmfs'   # needed bc the only mfs mode implemented into sdint
    #       elif specmode == 'cube':
    #           mydeconvolver = 'hogbom'
    #           #numchan = nchan           # not really needed here?





    
    # specmode, deconvolver and multiscale setup
    if multiscale:
        mydeconvolver = 'multiscale'
        if maxscale==-1:
            maxscale=derive_maxscale(vis, restfreq=restfreq)
        myqa = qatool()
        mycell = myqa.convert(myqa.quantity(cell),'arcsec')['value']
        myscales = [0]
        for i in range(0, int(math.log(maxscale/mycell,3))):
            myscales.append(3**i*5)

        print("My scales (units of pixels): "+str(myscales))

    else:    
        myscales = [0]
        mydeconvolver = 'hogbom'


    # weighting schemes
    if specmode == 'mfs':
        weightingscheme ='briggs'     # cont mode 
    elif specmode == 'cube':
        if get_casa_version() >= '6.2.0':
            weightingscheme ='briggsbwtaper'   # special briggs for cubes --- CURRENTLY SWITCHED OFF FOR SDINT 
        else:
            weightingscheme ='briggs'#bwtaper'   # special briggs for cubes --- WAIT FOR IMPLEMENTATION IN SDINT 

        



    # others
    npnt = 0    
    if phasecenter=='':
        phasecenter = npnt

    if restfreq=='':
        therf = []
    else:
        therf = [restfreq]


    clean_arg=dict(vis=vis,
                   field = field,
                   phasecenter=phasecenter,
                   imsize=imsize,
                   cell=cell,                                   
                   spw=spw,
                   specmode=specmode,
                   deconvolver=mydeconvolver,
                   scales=myscales,
                   nterms=1,                  # nterms=1 turns mtmfs into mfs, CASA 6.2 needs nterms=2 to run (bug?)
                   start=start,
                   width=width,
                   nchan = nchan, #numchan, 
                   restfreq=therf,
                   gridder='mosaic',          
                   weighting = weightingscheme,
                   robust = 0.5,
                   restoringbeam = 'common',   # SD-cube has only one beam - INT-cube needs it, too, else feather etc. fail
                   niter=niter,
                   cycleniter=cycleniter,
                   cyclefactor=2.0,
                   threshold=threshold,
                   interactive = interactive,
                   pbcor=True,               
                   # Masking Parameters below this line 
                   # --> Should be updated depending on dataset
                   usemask=usemask,
                   sidelobethreshold=sidelobethreshold,
                   noisethreshold=noisethreshold,
                   lownoisethreshold=lownoisethreshold, 
                   minbeamfrac=minbeamfrac,
                   growiterations=growiterations,
                   negativethreshold=negativethreshold,
                   mask=mask,
                   pbmask=pbmask,
                   verbose=True)

    return clean_arg

def runtclean(vis, 
              imname, 
              startmodel='',
              spw='', 
              field='', 
              specmode='mfs', 
              imsize=[], 
              cell='', 
              phasecenter='',
              start=0, 
              width=1, 
              nchan=-1, 
              restfreq='',
              threshold='', 
              niter=0, 
              cycleniter=-1,
              usemask='auto-multithresh' ,
              sidelobethreshold=2.0, 
              noisethreshold=4.25, 
              lownoisethreshold=1.5, 
              minbeamfrac=0.3, 
              growiterations=75, 
              negativethreshold=0.0,
              mask='', 
              pbmask=0.4, 
              interactive=True, 
              multiscale=False, 
              maxscale=0.,
              #restart=True,
              continueclean = False,
              loadmask=False,
              fniteronusermask=0.3
              ):

    """
    runtclean (A. Plunkett, NRAO, D. Petry, ESO)
    a wrapper around the CASA task "TCLEAN,"
    steps:
    - return mask setup
    - derive multiscale sizes for 'multiscale=True'
    - clean parameter definition as known from tclean (analoguously to 'runsdintimg')
    --- important fixed parameters you should be aware of: 
    restoring beam = common, gridder = 'mosaic',
    weighting='briggs', robust = 0.5
    - exportfits pcbor 
    vis - the MS containing the interferometric data
    imname - the root name of the output images
    startmodel - start model for cleaning
          default: '' i.e. no start model, example: 'TP.scaled.image'
    spw - the standard selection parameter spw of tclean
           default: '' i.e. all SPWs
    field - the standard selection parameter field of tclean
           default: '' i.e. all fields
    specmode - the standard tclean specmode parameter: supported are msf or cube
           default: msf 
    imsize - (optional) the standard tclean imsize parameter 
            should correspond to the imagesize for the most extended
            interferometer config.
           default: determine from the input MS with aU.pickCellSize
    cell - (optional) the standard tclean cell parameter
            should correspond to the cell size for the most extended
            interferometer config, i.e. smallest beam / 5.
           default: determine from the input MS with aU.pickCellSize
    phasecenter - the standard tclean phasecenter parameter
           e.g. 'J2000 12:00:00 -35.00.00.0000'
           default: '' - determine from the input MS with aU.pickCellSize
    start - the standard tclean start parameter
             default: 0
    width - the standard tclean width parameter
             default: 1
    nchan - the standard tclean nchan parameter
             default: -1
    restfreq - the restfrequency to write to the image for velocity calculations
             default: None, example: '115.271GHz'
    threshold - the threshold for cleaning
             default: None, example: '12mJy'
    niter - the standard tclean niter parameter
             default: 0, example: niter=1000000
    cycleniter -          
    usemask - the standard tclean mask parameter.  If usemask='auto-multithresh', can specify:
             sidelobethreshold, noisethreshold, lownoisethreshold, minbeamfrac, growiterations - 
             if usemask='user', must specify mask='maskname.mask' 
             if usemask='pb', can specify pbmask=0.4, or some level.
             default: 'auto-multithresh'
    interactive - the standard tclean interactive option
             default: True
    multiscale - if False (default) use hogbom cleaning, otherwise multiscale
    maxscale - for multiscale cleaning, use scales up to this value (arcsec)
             Recommended value: 10 arcsec
             default: 0.
    restart - True (default): Re-use existing images. False: Increment imagename
    continueclean - True: same as 'restart', False(default): Delete old version
    loadmask - run sdintimaging with user-specified mask for fniteronusermask*niter iterations 
              and continue with auto-masking (usemask='auto-multithresh') for the remaining 
              niter*(1-fniteronusermask) iterations 
             default: False
    fniteronusermask - adjusting the amount of iterations spend on a usermask for loadmask=True
              allowed values: 0.0 (none in theory, in fact: 1 iteration) - 1.0 (all)
             default: 0.3 
    Example: runtclean('gmc_120L.alma.all_int-weighted.ms', 
                'gmc_120L', phasecenter='J2000 12:00:00 -35.00.00.0000', 
                spw='0', field='0~68', imsize=[1120,1120], cell='0.21arcsec',
                threshold='12mJy',niter=100000,
                usemask='auto-multithresh')
    """


    #myvis = file_check_vis(vis)

    myvis = vis


    #   if type(vis) is str:
    #       myvis = file_check(vis)
    #   
    #           
    #   if isinstance(vis, list):
    #       #acceptinput=True
    #       for i in range(0,len(vis)):
    #           file_check(vis[i])          # if one of the files does not exist, script will exit here
    #           #if os.path.exists(vis[i]):
    #           #    #myvis = vis 
    #           #    pass 
    #           #else:
    #           #    acceptinput=False 
    #           #    print(vis[i]+' does not exist')
    #       #if acceptinput==False:
    #       #    return False             
    #       #else:
    #       #    myvis = vis 
    #       myvis = vis    
            
    #   if loadmask==True and fniteronusermask>1 or fniteronusermask<0:
    #       print('fniteronusermask is out of range: ' +fniteronusermask+' Please choose a value between 0 and 1 (inclusively)')
    #       return False        
    #   else:
    #       pass
            
        
    #print('')

    print('Start tclean ...')

    #mymaskname = ''
    #   if usemask == 'auto-multithresh':
    #       #print('Run with {0} mask'.format(usemask))
    #       mymask = usemask
    #       if os.path.exists(mask):
    #          mymaskname = mask
    #          print('Run with {0} mask and {1} '.format(mymask,mymaskname))
    #       else: print('Run with {0} mask'.format(mymask))        
    #   
    #   elif usemask =='pb':
    #       print('Run with {0} mask on PB level {1}'.format(usemask,pbmask))
    #   elif usemask == 'user':
    #       if os.path.exists(mask):
    #          print('Run with {0} mask {1}'.format(usemask,mask))
    #   
    #       else:
    #          print('### WARNING:   mask '+mask+' does not exist, or not specified')
    #          #return False
    #   else:
    #       print('check the mask options')
    #       return False

    #   # specmode and deconvolver
    #   if multiscale:
    #       mydeconvolver = 'multiscale'
    #       myqa = qatool()
    #       mycell = myqa.convert(myqa.quantity(cell),'arcsec')['value']
    #       myscales = [0]
    #       for i in range(0, int(math.log(maxscale/mycell,3))):
    #           myscales.append(3**i*5)
    #   
    #       print("My scales (units of pixels): "+str(myscales))
    #   
    #   else:    
    #       myscales = [0]
    #       mydeconvolver = 'hogbom'

    if niter==0:
        cta.casalog.post('You set niter to 0 (zero, the default). Only a dirty image will be created.', 'WARN', origin='runtclean')


    #   if specmode == 'mfs':
    #       weightingscheme ='briggs'   # cont mode
    #   elif specmode == 'cube':
    #       weightingscheme ='briggs' #bwtaper'   # special briggs for cubes --- WAIT FOR IMPLEMENTATION IN SDINT 
    #   

    #   tclean_arg=dict(vis = myvis,
    #                   imagename = imname, #+'.TCLEAN',
    #                   startmodel = startmodel,
    #                   field = field,
    #                   phasecenter = phasecenter,
    #                   imsize = imsize,
    #                   cell = cell,
    #                   spw = spw,
    #                   specmode = specmode,
    #                   deconvolver = mydeconvolver,   
    #                   scales = myscales,             
    #                   #nterms = 1,                    # needed by sdint for mtmfs
    #                   start = start, 
    #                   width = width, 
    #                   nchan = nchan, 
    #                   restfreq = restfreq,
    #                   gridder = 'mosaic',
    #                   weighting = weightingscheme,
    #                   robust = 0.5,
    #                   restoringbeam = 'common',   # SD-cube has only one beam - INT-cube needs it, too, else feather etc. fail
    #                   niter = niter,
    #                   cyclefactor=2.0,
    #                   threshold = threshold,
    #                   interactive = interactive,
    #                   pbcor = True,
    #                   # Masking Parameters below this line 
    #                   # --> Should be updated depending on dataset
    #                   usemask=usemask,
    #                   sidelobethreshold=sidelobethreshold,
    #                   noisethreshold=noisethreshold,
    #                   lownoisethreshold=lownoisethreshold, 
    #                   minbeamfrac=minbeamfrac,
    #                   growiterations=growiterations,
    #                   negativethreshold=negativethreshold,
    #                   mask=mask,
    #                   pbmask=pbmask,    # used by all usemasks! perhaps needed for fastnoise-calc!
    #                   verbose=True)#, 
    #                   #restart=restart)  # should switch off bc default



    tclean_arg = check_prep_tclean_param( 
                myvis,       
                spw, 
                field, 
                specmode,                                 
                imsize, 
                cell, 
                phasecenter,         
                start, 
                width, 
                nchan, 
                restfreq,
                threshold, 
                niter,
                cycleniter,
                usemask,
                sidelobethreshold,
                noisethreshold,
                lownoisethreshold, 
                minbeamfrac,
                growiterations,
                negativethreshold,                
                mask, 
                pbmask,
                interactive,               
                multiscale, 
                maxscale,
                loadmask,
                fniteronusermask
                )

    #tclean_arg['vis']        = myvis
    tclean_arg['imagename']  = imname #+'.TCLEAN',
    tclean_arg['startmodel'] = startmodel


    #if os.path.exists(imname+'.TCLEAN.image'):
    #    casalog.post('Image '+imname+'.TCLEAN already exists.  Running with restart='+str(restart), 'WARN')        
    ##os.system('rm -rf '+imname+'.TCLEAN.*')

    if continueclean == False:
        os.system('rm -rf '+imname+'.*') #+'.TCLEAN.*')   
        # if to be switched off add command to delete "*.pbcor.fits"
    

    if loadmask==True:

        tclean_arg['usemask']='user'
        if fniteronusermask==0 or fniteronusermask==0.0:
            tclean_arg['niter']=1
        else:   
            tclean_arg['niter']=int(niter*fniteronusermask)
        # load mask into tclean with limited iterations
        print('')
        print('### Load mask into tclean with iterations = fniteronusermask*niter = ', fniteronusermask, ' * ', niter)
        report_mask(tclean_arg['usemask'],tclean_arg['mask'],tclean_arg['pbmask'],tclean_arg['niter'])
        tcleansresults = cta.tclean(**tclean_arg)
        
        tclean_arg['usemask']=usemask
        tclean_arg['mask']=''
        # if startmodel used, it would have been loaded in tclean step before 
        # -> clear startmodel parameter for next tclean call, else crash!
        tclean_arg['startmodel']=''        
        
        tclean_arg['niter']=niter-tclean_arg['niter']
        
        if tclean_arg['niter']<=0:    #avoid negative niter values and pointless executions 
            pass
        else:               
            # clean and get tclean-feedback 
            print('')
            print('### Continue tclean with iterations = (1-fniteronusermask)*niter = ', (1.0-fniteronusermask), ' * ', niter)
            report_mask(tclean_arg['usemask'],tclean_arg['mask'],tclean_arg['pbmask'],tclean_arg['niter'])
            tcleansresults = cta.tclean(**tclean_arg)
        
        #### store feedback in a file 
        ###pydict_to_file2(tcleansresults, imname)
        ###
        ###os.system('cp -r summaryplot_1.png '+imname+'.png')   

    else: 
        # clean and get tclean-feedback 
        report_mask(tclean_arg['usemask'],tclean_arg['mask'],tclean_arg['pbmask'],tclean_arg['niter'])
        tcleansresults = cta.tclean(**tclean_arg)
        
    # store feedback in a file 
    pydict_to_file2(tcleansresults, imname)
    
    os.system('cp -r summaryplot_1.png '+imname+'.png') 



    # print('Exporting final pbcor image to FITS ...')
    # #exportfits(imname+'.TCLEAN.image.pbcor', imname+'.TCLEAN.pbcor.fits')
    # os.system('rm -rf '+imname+'.image.pbcor.fits') #+'.TCLEAN.*')   
    # os.system('rm -rf '+imname+'.pb.fits') #+'.TCLEAN.*')   
    # cta.exportfits(imname+'.image.pbcor', imname+'.image.pbcor.fits')
    # cta.exportfits(imname+'.pb', imname+'.pb.fits')

    export_fits(imname)


    return True

def report_mask(usemask, mask, pbmask, niter       
                ):

    """
    report selected mask used for tclean/sdint (Moser-Fischer, L.)
    a helper for runsdintimg and runtclean
    usemask - masking mode parameter as for tclean
    mask - file name of mask
    pbmask - PB mask cut-off level 
    niter - number of iterations spent on this mask
    """
     
    print('')
    if usemask == 'auto-multithresh':
        print('### Run with {0} mask for {1} iterations ###'.format(usemask,niter))        
    elif usemask =='pb':
        print('### Run with {0} mask on PB level {1} for {2} iterations ###'.format(usemask,pbmask,niter))
    elif usemask == 'user':
        if os.path.exists(mask):
           print('### Run with {0} mask {1} for {2} iterations ###'.format(usemask,mask,niter))
        else:
           print('### WARNING:   mask '+mask+' does not exist, or is not specified. ###')
           #return False
    else:
        print("### Invalid usemask '"+usemask+"'. Please, check the mask options. ###")
        return False
    print('---------------------------------------------------------')

def pydict_to_file2(pydict, filename):   
    """ 
    save dictionaries (L. Moser-Fischer)
    pydict - a python dictionary
    filename - will store dict under this string +'.txt' 
    """ 
    import pickle as pk 
    with open(filename+'.pickle', 'wb') as handle:
        pk.dump(pydict, handle, protocol=pk.HIGHEST_PROTOCOL)
    return True 
def pydict_to_file2(pydict, filename):   
    """ 
    save dictionaries (L. Moser-Fischer)
    pydict - a python dictionary
    filename - will store dict under this string +'.txt' 
    """ 
    import pickle as pk 
    with open(filename+'.pickle', 'wb') as handle:
        pk.dump(pydict, handle, protocol=pk.HIGHEST_PROTOCOL)
    return True 


def export_fits(imname, clean_origin=''):

    """
    standardized output from a combination method  (Moser-Fischer, L.)
    imname - file name base
    clean_origin - file name base of intermediate cleaning products that imname is based on (e.g. feather products is based on tclean products)
                   
    """

    print('')
    print('Exporting following final image products to FITS for '+imname+':')

    for suffix in ['.image.pbcor', '.pb']:
        if clean_origin!='' and suffix=='.pb':
            os.system('cp -r '+clean_origin+suffix+' '+imname+suffix)
        os.system('rm -rf '+imname+suffix+'.fits')   
        print('-', suffix)
        cta.exportfits(imname+suffix, imname+suffix+'.fits')
    print('Export done.')
    print('')
    
    return True


#Imaging: Using a wrapper of DC scripts (by Plunkett et al 2023)--Adapted from the EMERGE pipeline
def do_imaging_singlechan(myINTvis,mySD,myline,Niters,DCsteps=[1,2,3,5],DCdeconv="HB"):
	'''
do_imaging_singlechan (EMERGE)
	
	Author: A.Hacar + EMERGE team
	
	Wrapper of the DC scripts (see Plunkett+2023) for line imaging and data combination

	Arguments
	----------
	  myINTvis: str 
		Path to target visibilities.
	  mySD: str
		Path to SD data
	  myline : str or float
		Reference frequency for deconvolution. The sky frequency is selected after correction by
		the source velocity according to the MS file.
		If "str" the line ID will be identified in emerge.linecat (e.g. "N2H+10").
		If else "float" the line will be taken as frequency in Hz 
	  Niters : int
		Number of iterations for decolvolution
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
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('Running now DC script usign the parameters selected in set_parameters.py')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++') 

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
	with open("./sim_DCpars_template.py", "r") as file:
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
	exec(open("./my_DCpars.py").read(),globals())
	exec(open("./DataComb/DC_run.py").read(),globals())


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

