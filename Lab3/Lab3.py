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
    print "***", np.median(img), '***'
    print hdr['object']
    print img[:,256]
    img = dead_pixels(img)
    img = img[:,:1024] #ignore the 32 columns
    return [img, hdr]
    
#==========================
# average the bias arrays
#===========================
def avg_bias(bias_file_arr):
    """ takes all the bias fits and makes them into one average bias fits array """
    #apparently this does the same thing, but I wanted to write it out myself
    #np.sum(bias_data_arr, axis = 0)/len(bias_data_arr)
    bias_data_arr = [pf.getdata(bias_file_arr[i]) for i in range(len(bias_file_arr))]
    bias_arr = bias_data_arr[0]
    for i in range(1, len(bias_file_arr)):
        bias_arr += bias_data_arr[i]
    return bias_arr/(len(bias_file_arr))
bias = avg_bias(bias_file_arr)
bias = dead_pixels(bias)
bias = bias[:,:1024]
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
               vmin=0.95*np.median(img), vmax=1.7*np.median(img))
    plt.colorbar()
    return fig
    
plt.show(display_fits(DanaeR_file_arr[5], flatR_img))
#display_fits(DanaeR_file_arr[5], flatR_img).savefig(fig_path+'DanaeR_5_corrected.png',dpi=300)
##============================


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
    return [i for i in range(36, len(arr)-1) if all(arr[i] > arr[i-width:i]) and all(arr[i]>arr[i+1:i+width]) and arr[i]>lower_limit] # 36 is arbitrary
    #assuming that there is nothing before 36
def limit_applier(index_arr, arr, lower_limit):
    """ Makes sure the considered indexes are above a certain value"""    
    return [i for i in index_arr if arr[i] > lower_limit]
#def peak_verifier(index_arr, arr, n, max = True):
#    """
#    """
#    k = 0
#    delete_arr = []
#    cluster_arr =[]
#    while k < len(index_arr) - 1:
#   
#        curr_index_x = index_arr[k][0]
#        curr_index_y = index_arr[k][1]
#        next_index_x = index_arr[k+1][0]
#        next_index_y = index_arr[k+1][1]
#        curr_len = len(index_arr[k])
#        next_len = len(index_arr[k+1])
#        curr_arr = index_arr[k]
#        next_arr = index_arr[k+1]
#
#        temp_arr=[]
#        #if np.abs(curr_index_x-next_index_x) <= 3) and np.abs(curr_index_x-next_index_x) <= 3)
#        if curr_len ==0:
#            continue
#        m = k+1
#        elif curr_len > 0 and next_len > 0:
#            for elem in curr_arr:
#                temp_arr.append([k,elem])
#                for m_elem in index_arr[m]:
#                    if np.abs(m_elem - elem) <= 4:
#                        temp.arr.append([m,m_elem])
#                m+=1
#                cluster_arr.appened(temp_arr)
#
#                while index_arr[m]
#                np.abs(index_arr[m] - elem) <= 4
#                temp_arr.append([m,elem])
#                
#        
#        if np.abs(curr_index-next_index) <= n:
#            curr_lower, curr_upper = curr_index - n/2., curr_index + n/2.
#            if curr_lower < 0:
#                curr_lower = 0
#            if max:
#                curr_max = np.amax(arr[curr_lower:curr_upper])
#            else:
#                curr_min = np.amin(arr[curr_lower:curr_upper])
#            if (max and arr[curr_index] == curr_max) or (not max and arr[curr_index] == curr_min):
#                delete_arr.append(k+1)
#            else:
#                delete_arr.append(k)  
#        k+=1
#    return np.delete(index_arr, delete_arr)
def centroids1(fil, flat_img):
    info = loadfits(fil)
    img, hdr = info[0], info[1]
    #img = bsub(bias, fil)
    img = correct_fits(fil, bias, flat_img)
    avg = np.median(img)
    lim = 2.35*avg
    #[x][y]
    img_arr = [[]] #ignore the 0th row as not having anything
    print "%%%%%"
    print 'avg', avg
    print 'limit', lim
    #print max_peaks(img[510], 2.5*avg)
    #print img[510][max_peaks(img[510], 2.5*avg)]
    
    for i_row in range(1, len(img)): #range starts at 1 bc of the assumption that img[0] has nothing
        print i_row, max_peaks(img[i_row], 10, lim)
        img_arr.append(max_peaks(img[i_row], 10, lim))
    print '****', img_arr[0]
    
    #for i in range(1, len(img_arr)-1):
    #    print i, img_arr[i]
    #    if len(img_arr[i-1]>0) or len(img_arr[i-1]<0):
    img_arr_pair =[]
    for i_row in range(len(img_arr)):
        if len(img_arr)>0:
            for elem in img_arr[i_row]:
                img_arr_pair.append([i_row,elem]) # [x_index, y_index]
    #print img_arr_pair # [ [x_1,y_1], [x_2,y_2] ... [x_n,y_n] ]
    def getKey(item):
        return item[1]
    img_arr_pair = sorted(img_arr_pair, key = getKey)
    cluster_arr = []
    k=0
    while k < (len(img_arr_pair)-2):
        curr_pair, next_pair, next_next_pair = img_arr_pair[k], img_arr_pair[k+1], img_arr_pair[k+2]
        temp_arr=[curr_pair]
        if curr_pair[0] == 541:
            print curr_pair, next_pair
        if not ((np.abs(curr_pair[1]-next_pair[1])<=10 and np.abs(curr_pair[0]-next_pair[0])<=10) or 
        (np.abs(curr_pair[0]-next_next_pair[0])<=10 and np.abs(curr_pair[1]-next_next_pair[1]))<=10):
            k+=1
        else:
            while (np.abs(curr_pair[1]-next_pair[1])<=10 and np.abs(curr_pair[0]-next_pair[0])<=10) or (np.abs(curr_pair[0]-next_next_pair[0])<=10 and np.abs(curr_pair[1]-next_next_pair[1])<=10):
                if (np.abs(curr_pair[1]-next_pair[1])<=10 and np.abs(curr_pair[0]-next_pair[0])<=10):
                    temp_arr.append(next_pair)
                elif (np.abs(curr_pair[0]-next_next_pair[0])<=10 and np.abs(curr_pair[1]-next_next_pair[1])<=10):
                    temp_arr.append(next_pair)
                k+=1
                curr_pair, next_pair, next_next_pair = img_arr_pair[k], img_arr_pair[k+1], img_arr_pair[k+2]
            
        cluster_arr.append(temp_arr)
        k+=1

    for elem in cluster_arr:
        print elem
    #for i in range(1,len(img_arr_tuple)-1):

        
        
  
    
    #for i in range(len(img)):
    #    #print np.where(img[i][max_peaks(img[i])]>3*avg)[0]
    #    img_arr.append(np.where(img[i][max_peaks(img[i])]>3*avg)[0])

    #print img_arr
#def centroids(index_arr, x_arr, y_arr, peak_width):
#        n = peak_width/2
#	centroids= []
#	for peak_index in index_arr:
#		x_range = x_arr[peak_index-n:peak_index+n]
#		y_range = y_arr[peak_index-n:peak_index+n]
#		centroid = np.sum(x_range*y_range)/np.sum(y_range) #<x>
#
#		numerator = []
#		for i in x_range:
#		    numerator.append(y_arr[i]*(x_arr[i]-centroid)**2)
#		error = np.sum(numerator) / (np.sum(y_range))**2 
#		centroids.append([centroid, error])
#		print error
#		#centroids.append(centroid)
#	centroids, centroid_errors = np.transpose(centroids)[0], np.transpose(centroids)[1]
#	return centroids, centroid_errors
centroids1(DanaeR_file_arr[5], flatR_img)

    


    
