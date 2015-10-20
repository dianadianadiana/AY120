
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
fits_list=glob.glob(datadir+datestr+'*.fits') #build list of images
#==================================

#### average the bias arrays
bias_file_arr = ['d'+str(i) + '.fits' for i in np.arange(100,115)]
flatB_file_arr = ['d'+str(i) + '.fits' for i in np.arange(115,118)]
flatV_file_arr = ['d'+str(i) + '.fits' for i in np.arange(118,121)]
flatR_file_arr = ['d'+str(i) + '.fits' for i in np.arange(121,124)]
flatI_file_arr = ['d'+str(i) + '.fits' for i in np.arange(124,127)]
WratR_file_arr = ['d'+str(i) + '.fits' for i in np.arange(131,142)]
WratI_file_arr = ['d'+str(i) + '.fits' for i in np.arange(142,152)]
DanaeR_file_arr = ['d'+str(i) + '.fits' for i in np.arange(152,163)]
DanaeI_file_arr = ['d'+str(i) + '.fits' for i in np.arange(163,172)]



print bias_file_arr
print (pf.getdata(datadir+datestr+bias_file_arr[0])+pf.getdata(datadir+datestr+bias_file_arr[1]))/2
print pf.getdata(datadir+datestr+bias_file_arr[0])

def avg_bias(bias_file_arr):
    pf.getdata(fits1)

bias = fits_list[0]
fits1 = fits_list[25]
def display_info(fits_file):
    hdulist=pf.open(datadir+datestr+fits_file)
    img=hdulist[0].data
    hdr=hdulist[0].header
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

display_info(0)
display_info(12)
display_info(25)
display_info(55)


x = pf.getdata(fits1)
hdr = pf.getheader(fits1)

print x.mean()
print x.min()
print x.max()

#An HDU (Header Data Unit) is the highest level component of the FITS file structure, 
#consisting of a header and (typically) a data array or table.
#hdulist[0] is the primary HDU
# header includes useful info like size of the image (NAXIS1 & NAXIS2), the data
# (DATE-OBS), the exposure time (EXPTIME), where the celestial coordiantes of where 
# telescopes was pointing (RA & DEC) and the optical filter in the light path (FILTNAM)
hdulist = pf.open(fits1)
hdulist.info()
img = hdulist[0].data
print img
def display_fits(fits_index):
    """Display fits using matplotlib"""
    global datadir ; global datestr ; global fits_list
    
    #open random image
    hdulist=pf.open(fits_list[fits_index])
    img=hdulist[0].data
    hdr=hdulist[0].header
    
    print img.shape
    print np.argmax(img[0])
    print img[256]
    print img[0]
    print img[256][256]
    print img[0][256]
    print np.mean(img[0])

    fig = plt.figure()
    plt.title(hdulist[0].header['object'])
    plt.imshow(img[:], origin='lower', interpolation='nearest', cmap='afmhot', 
               vmin=0.95*np.median(img), vmax=1.4*np.median(img))
    plt.colorbar()
    return fig
plt.show(display_fits(25))
    
#shape, (1024, 1056) == 1056 elements in each array; 1024 arrays
#on the xaxis = 1056, yaxis = 1024
hdulist.close()

    #bias_arr = range(len(img))
    #for i in range(len(img[0])):
    #    temp_arr = np.array([])
    #    for j in range(len(img)):
    #        temp_arr = np.append(temp_arr, img[j][i])
    #    bias_arr[i] = np.mean(img[i])
    #print bias_arr[256]