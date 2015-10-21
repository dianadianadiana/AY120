# -*- coding: utf-8 -*-

#==================================
#Import Packages
#==================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import astropy.io.fits as pf
import glob
#import os #it tries to be IDL's SPAWN and falls short but it still useful
import datetime as dt
#==================================

# How to use fits and astropy
# http://docs.astropy.org/en/stable/io/fits/index.html
#==================================
# Global Variables
#==================================
datadir='/Users/Diana/Desktop/Astro 120/Lab3/'
datestr='151013/' #probably want all directory strings storred in an array
datapath = datadir + datestr
fits_list=glob.glob(datadir+datestr+'*.fits') #build list of images
fig_path = '/Users/Diana/Desktop/Astro 120/Lab3/Figures/'
#==================================

#==========================
# making arrays for the diff types of objects
# each element of the array has the syntax of datapath + 'd100.fits'
# DOES INCLUDE the actual file path
#==========================
bias_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(100,115)]
flatB_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(115,118)]
flatV_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(118,121)]
flatR_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(121,124)]
flatI_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(124,127)]
WratR_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(131,142)]
WratI_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(142,152)]
DanaeR_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(152,163)]
DanaeI_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(163,172)]
#==========================
#===========================
def dead_pixels(img):
    avg = np.median(img)    
    img[:,255]=avg; img[:,256]=avg; img[:,1002]=avg # dead pixels -- set to zero
    img[:,1023]=avg; img[:,1024]=avg; img[:,1022]=avg 
    img[:,783]=avg; img[:,669]=avg
    return img
    
def loadfits(fil):
    """ Input: 
            fil: the filename, with the syntax of 'd100.fits'
        Output: 
            img: the matrix of the fits file
            hdr: the header of the fits file """
    hdulist = pf.open(fil)
    img = hdulist[0].data # the matrix data
    hdr= hdulist[0].header # the header
    avg = np.median(img)
    print "***", avg, '***'
    print hdr['object']
    print img[:,256]
    img = dead_pixels(img)
    return [img, hdr]
    
#==========================
# average the bias arrays
#===========================
bias_data_arr = [pf.getdata(bias_file_arr[i]) for i in range(len(bias_file_arr))]
def avg_bias(bias_file_arr):
    """ takes all the bias fits and makes them into one average bias fits array """
    #apparently this does the same thing, but I wanted to write it out myself
    #np.sum(bias_data_arr, axis = 0)/len(bias_data_arr)
    bias_arr = bias_data_arr[0]
    for i in range(1, len(bias_file_arr)):
        bias_arr += bias_data_arr[i]
    return bias_arr/(len(bias_file_arr))
bias = avg_bias(bias_file_arr)
bias = dead_pixels(bias)
#===========================

#===========================
# Flats
#===========================
flatB_img, flatB_hdr = loadfits(flatB_file_arr[0])
flatV_img, flatV_hdr = loadfits(flatV_file_arr[0])
flatR_img, flatR_hdr = loadfits(flatR_file_arr[0])
flatI_img, flatI_hdr = loadfits(flatI_file_arr[0])
#===========================

#===========================
    
def display_info(fil):
    info = loadfits(fil)
    img, hdr = info[0], info[1]
    print '****************'
    print "object:", hdr['object']
    print "exposure time:", hdr['exptime']
    print 'date begin:', hdr['DATE-BEG']
    print 'dec:', hdr['dec']
    print 'ra:', hdr['ra']
    print 'start time temperature (deg C):', hdr['tempdet']
    print 'filter:', hdr['filtnam']
    print 'cover:',hdr.get('cover')
    print '****************'

display_info(WratR_file_arr[0])
display_info(WratR_file_arr[1])
display_info(WratR_file_arr[2])
display_info(WratR_file_arr[3])

#An HDU (Header Data Unit) is the highest level component of the FITS file structure, 
#consisting of a header and (typically) a data array or table.
#hdulist[0] is the primary HDU
# header includes useful info like size of the image (NAXIS1 & NAXIS2), the data
# (DATE-OBS), the exposure time (EXPTIME), where the celestial coordiantes of where 
# telescopes was pointing (RA & DEC) and the optical filter in the light path (FILTNAM)
    
#shape, (1024, 1056) == 1056 elements in each array; 1024 arrays
#on the xaxis = 1056, yaxis = 1024
def bsub(bias, fil):
    info = loadfits(fil)
    img, hdr = info[0], info[1]
    #cov = hdr['cover']
    #average_bias= np.mean(bias[:,-cov:])
    #img-=average_bias
    img-=bias
    return img
def correct_fits(fil, bias, flat):
    img, hdr = loadfits(fil)
    return (img - bias)/(flat-bias) * (np.average(flat-bias))
def display_fits(fil, flat_img):
    """Display fits using matplotlib"""    
    img, hdr = loadfits(fil)
    img = correct_fits(fil, bias, flat_img)
    fig = plt.figure()
    plt.title(hdr['object'], fontsize=18)
    plt.imshow(img, origin='lower', interpolation='nearest', cmap='gray_r', 
               vmin=0.95*np.median(img), vmax=3*np.median(img))
    plt.colorbar()
    return fig
    
plt.show(display_fits(DanaeR_file_arr[5], flatR_img))
#display_fits(DanaeR_file_arr[5], flatR_img).savefig(fig_path+'hi1.png',dpi=300)
##============================
## raw-bias (average of flat-bias)
##----------
## flat - bias


#xb = bsub(bias, DanaeR_file_arr[5])
#print np.mean(xb)
#flatb = bsub(bias,flatR_file_arr[0])
#flatb = flatb/np.median(flatb)
#
#fig = plt.figure()
#ax = fig.add_subplot(121)
#cax = ax.imshow(flatb,cmap='gray_r',vmin=0.8,vmax=1.2)
#fig.colorbar(cax)
#ax1 = fig.add_subplot(122)
#cax1 = ax1.imshow(xb/flatb,cmap='gray_r',vmin=50,vmax=500)
#fig.colorbar(cax1)
##plt.show()









#============================
# Centroids
#============================
def max_peaks(arr, width, lower_limit):
    """ Returns an array of the indexes where the indexes indicate where the peaks are"""
    #return [i for i in range(1, len(arr)-1) if arr[i-1]<arr[i] and arr[i+1]<arr[i] and arr[i]>lower_limit]
    return [i for i in range(35, len(arr)-1-32) if all(arr[i] > arr[i-width:i]) and all(arr[i]>arr[i+1:i+width]) and arr[i]>lower_limit]
def limit_applier(index_arr, arr, lower_limit):
    """ Makes sure the considered indexes are above a certain value"""    
    return [i for i in index_arr if arr[i] > lower_limit]
def centroids1(fil, flat_img):
    info = loadfits(fil)
    img, hdr = info[0], info[1]
    #img = bsub(bias, fil)
    img = correct_fits(fil, bias, flat_img)
    avg = np.median(img)
    lim = 2.65*avg
    #[x][y]
    img_arr = []
    print "%%%%%"
    print 'avg', avg
    print 'limit', lim
    #print max_peaks(img[510], 2.5*avg)
    #print img[510][max_peaks(img[510], 2.5*avg)]
    
    for i in range(len(img)):
        print i, max_peaks(img[i], 10,lim)
        img_arr.append(max_peaks(img[i], 10, lim))
    print '****', img_arr[0]
    
    #for i in range(1, len(img_arr)-1):
    #    print i, img_arr[i]
    #    if len(img_arr[i-1]>0) or len(img_arr[i-1]<0):
    img_arr_tuple =[]
    for i in range(len(img_arr)):
        if len(img_arr)>0:
            for elem in img_arr[i]:
                img_arr_tuple.append([i,elem])
    print img_arr_tuple
    
    cluster_arr = []
    for i in range(1,len(img_arr_tuple)-1):
        
        
        
                
    
    #for i in range(len(img)):
    #    #print np.where(img[i][max_peaks(img[i])]>3*avg)[0]
    #    img_arr.append(np.where(img[i][max_peaks(img[i])]>3*avg)[0])

    #print img_arr
def centroids(index_arr, x_arr, y_arr, peak_width):
        n = peak_width/2
	centroids= []
	for peak_index in index_arr:
		x_range = x_arr[peak_index-n:peak_index+n]
		y_range = y_arr[peak_index-n:peak_index+n]
		centroid = np.sum(x_range*y_range)/np.sum(y_range) #<x>

		numerator = []
		for i in x_range:
		    numerator.append(y_arr[i]*(x_arr[i]-centroid)**2)
		error = np.sum(numerator) / (np.sum(y_range))**2 
		centroids.append([centroid, error])
		print error
		#centroids.append(centroid)
	centroids, centroid_errors = np.transpose(centroids)[0], np.transpose(centroids)[1]
	return centroids, centroid_errors
centroids1(DanaeR_file_arr[5], flatR_img)


    
