# -*- coding: utf-8 -*-
#==================================
#Import Packages
#==================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as col
import astropy.io.fits as pf
import glob
import math
#import os #it tries to be IDL's SPAWN and falls short but it still useful
import datetime as dt
#==================================
# How to use fits and astropy
# http://docs.astropy.org/en/stable/io/fits/index.html
#==================================
# Global Variables
#==================================
datadir='/Users/Diana/Desktop/Astro 120/Lab3/'
date_folders = ['151013/', '151019/', '151020/']
fig_path = '/Users/Diana/Desktop/Astro 120/Lab3/Figures/'
#==================================

#==================================
# making arrays for the diff types of objects
# each element of the array has the syntax of datapath + 'd100.fits'
# DOES INCLUDE the actual file path
#==================================
def get_fitsfiles(date_folder, a, b):
    return [datadir + date_folder + 'd'+str(i) + '.fits' for i in np.arange(a,b)]
# from 10 13 15
bias_file_arr = get_fitsfiles(date_folders[0],100,115)
flatB_file_arr = get_fitsfiles(date_folders[0],115,118)
flatV_file_arr = get_fitsfiles(date_folders[0],118,121)
flatR_file_arr = get_fitsfiles(date_folders[0],121,124)
flatI_file_arr = get_fitsfiles(date_folders[0],124,127)
WratR_file_arr = get_fitsfiles(date_folders[0],131,142)
WratI_file_arr = get_fitsfiles(date_folders[0],142,152)
DanaeR_file_arr = get_fitsfiles(date_folders[0],152,163)
DanaeI_file_arr = get_fitsfiles(date_folders[0],163,173)
AlineR_file_arr = get_fitsfiles(date_folders[0],173,183)
AlineI_file_arr = get_fitsfiles(date_folders[0],183,193)
GiselaR_file_arr = get_fitsfiles(date_folders[0],193,203)
GiselaI_file_arr = get_fitsfiles(date_folders[0],203,213)
GalateaR_file_arr = get_fitsfiles(date_folders[0],213,223)
GalateaI_file_arr = get_fitsfiles(date_folders[0],223,233)
flat_arr = [flatB_file_arr,flatV_file_arr,flatR_file_arr,flatI_file_arr]

# from 10 19 15
bias_file_arr1 = get_fitsfiles(date_folders[1],100,115)
flatB_file_arr1 = get_fitsfiles(date_folders[1],115,118)
flatV_file_arr1 = get_fitsfiles(date_folders[1],118,121)
flatR_file_arr1 = get_fitsfiles(date_folders[1],121,124)
flatI_file_arr1 = get_fitsfiles(date_folders[1],125,128) # 124 was BAD
AlineR_file_arr1 = get_fitsfiles(date_folders[1],128,138)
AlineI_file_arr1 = get_fitsfiles(date_folders[1],138,148)
GiselaR_file_arr1BAD = get_fitsfiles(date_folders[1],148,156) # bad files
WratR_file_arr1 = get_fitsfiles(date_folders[1],156,167)
WratI_file_arr1 = get_fitsfiles(date_folders[1],167,177)
DanaeR_file_arr1 = get_fitsfiles(date_folders[1],177,188)
DanaeI_file_arr1 = get_fitsfiles(date_folders[1],188,198)
GiselaR_file_arr1 = get_fitsfiles(date_folders[1],198,208)
GiselaI_file_arr1 = get_fitsfiles(date_folders[1],208,218)
GalateaR_file_arr1 = get_fitsfiles(date_folders[1],218,228)
GalateaI_file_arr1 = get_fitsfiles(date_folders[1],228,238)
flat_arr1 = [flatB_file_arr1,flatV_file_arr1,flatR_file_arr1,flatI_file_arr1]

# from 10 20 15
bias_file_arr2 = get_fitsfiles(date_folders[2],100,115)
flatB_file_arr2 = get_fitsfiles(date_folders[2],115,118)
flatV_file_arr2 = get_fitsfiles(date_folders[2],118,121)
flatR_file_arr2 = get_fitsfiles(date_folders[2],122,125) #121 was BAD
flatI_file_arr2 = get_fitsfiles(date_folders[2],125,128) 
WratR_file_arr2 = get_fitsfiles(date_folders[2],128,138)
WratI_file_arr2 = get_fitsfiles(date_folders[2],138,148)
DanaeR_file_arr2 = get_fitsfiles(date_folders[2],148,158)
DanaeI_file_arr2 = get_fitsfiles(date_folders[2],158,168)
GalateaR_file_arr2 = get_fitsfiles(date_folders[2],168,178)
GalateaI_file_arr2 = get_fitsfiles(date_folders[2],178,188)
AlineR_file_arr2 = get_fitsfiles(date_folders[2],188,198)
AlineI_file_arr2 = get_fitsfiles(date_folders[2],198,208)
GiselaR_file_arr2 = get_fitsfiles(date_folders[2],208,218)
GiselaI_file_arr2 = get_fitsfiles(date_folders[2],218,228)
flat_arr2 = [flatB_file_arr2, flatV_file_arr2, flatR_file_arr2, flatI_file_arr2]

all_bias_arr = [bias_file_arr, bias_file_arr1, bias_file_arr2]
all_flat_arr = [flat_arr, flat_arr1, flat_arr2]
#==================================
#==================================
#####  ======= WAYS THERE ARE ERRORS IN OUR FLAT FIELD CORRECTION
#####  ======= TIPS/THOUGHTS TO KEEP IN MIND
# (*) the bias is 0s exposure time -- can be used in any filter
# (*) the flat is 5-10s exposure time and points at the sunset or sunrise -- when
#     using this, make sure you have corresponding filters
# (*) the last 32 cover pixels are there to take the bias if an actual bias img is
#     not available; otherwise if a bias img is provided, the last 32 cover pixels
#     are essentially not needed, so we can ignore them
# (*) shape, (1024, 1056) == 1056 elements in each array; 1024 arrays
#     on the xaxis = 1056, yaxis = 1024
# (1) in the dead_pixels fn, we are taking the avg of the [1024,1024] img and not
#     including the last 32 cover pixels. (I think it is safe to do so because the 
#     actual data lies within [1024,1024]. 
# (2) in the dead_pixels fn, we are setting arbitrary pixels to the average value;
#     the pixels were chosen based on manual selection, and i chose to set them to 
#     the avg and not 0 because i wanted a more "likely" value for the dead pixels
# (3) in the loadfits fn, we return the img matrix corrected for the dead pixels and 
#     has the dimensions of [1024,1024] (so we ignore/get rid of the cover pixels)
# (4) in the avg_fits_img fn, we are trying to minimalize errors. <ex> When we have 14 bias 
#     fits files, it is better to sum them all up and then take the average to get an
#     averaged bias fits. Likewise, we can use this with the flat images. NOTE: the flat
#     fits files have diff notes -- the diff notes may indicate that we can't just avg
#     the flats, so what we do instead is noted in the next step. EDIT: we do deal with the
#     flats differently by taking the gain of each one and then averaging those gains
# (5) in the find_gain fn, we are finding the gain of each flat and then taking the avg
#     gain. The reason we do this is due to the fact that the notes are diff for the flats.
#     If the flats had different brightness levels / saturations [I think that's what
#     those numbers, the 41k, etc. are?] it wouldn't skew the correction in any way.
#     gain = np.mean(flat-bias)/(flat-bias) == this is a matrix
# (6) the bias portion is just applying the functions to get the averaged bias
# (7) the flats portion can be a little confusing, i will explain by sections
#     (i)   i just loaded the img and hdr for each filter. only one file (instead of
#           multiple) is used. this is not corrected at all except for the dead pixels
#           and that i made them [1024,1024]
#     (ii)  i then decided to avg all the flats using the same fn i used on the bias,
#           they are not corrected in any way
#     (iii) i then corrected the avg flats in (ii)
#     (iv)  i then approached the flat by doing the gain methed (see (5) for more info)
# (8) Findings!
#     I then plotted the differences in correcting the img by using the bias and flats
#     I am very certain that finding the bias was a good way to go, I just wanted to see
#     and test out how the different approaches in handling the flats would yield diff
#     corrected displays. 
#     honestly, i did not find THAT much of a difference, if any when using flat (i) (iii)
#     or (iv). Maybe the error is not visible to us, but to the computer? Despite the fact
#     that the displays look relatively the same, we decided to go with flat option (iv)
#     because it really made little to no room for error for future calculations
# (9) the correct_fits fn corrects the raw img, and by raw, i mean the img loaded in after
#     the loadfits fn, so it is already dead pixel corrected and has the shape [1024,1024]
#     the correction is: (img - bias) * gain
#     and then we just display the beautiful graph!
# *edit: I decided to make a get_bias() function in order to get the correct bias since it
#       is surprisingly different everyday
#==================================
#==================================
def dead_pixels(img):
    """ There are many dead pixels in the CCD detector, so we just set them to avg
    of the img when it is [1024,1024]"""
    avg = np.median(img[:,:1024]) 
    img[:,255]=avg; img[:,256]=avg; img[:,1002]=avg # dead pixels -- set to avg
    img[:,1023]=avg; img[:,1024]=avg; img[:,1022]=avg 
    img[:,783]=avg; img[:,784]=avg; img[:,669]=avg
    return img
    
def loadfits(fil):
    """ we take care of the dead pixels and we return the img W/O the cover pixels
        Input: 
            fil: the filename, with the syntax of 'd100.fits'
        Output: 
            img: the matrix of the fits file [1024,1024]
            hdr: the header of the fits file """
    hdulist = pf.open(fil)
    img = hdulist[0].data # the matrix data
    hdr= hdulist[0].header # the header
    print "Loading file for", hdr['object']
    img = dead_pixels(img)
    img = img[:,:1024] #ignore the 32 columns
    return [img, hdr]
    
def avg_fits_img(files_arr):
    """ takes all the fits in the 'files_arr' and makes them into one average fits array (img) 
        Note: no corrections are applied, just the raw fits files are used """
    #apparently this does the same thing, but I wanted to write it out myself
    #np.sum(bias_data_arr, axis = 0)/len(bias_data_arr)
    data_arr = [pf.getdata(fil) for fil in files_arr]
    avg_arr = data_arr[0]
    for arr in data_arr[1:]:
        avg_arr += arr
    return avg_arr/(len(files_arr))
    
def find_gain(flat_files, bias):
    def gain(x,y):
        """ x: bias, y:flat -- both are the img matricies 
            gain = np.mean(flat-bias)/(flat-bias) """
        return np.mean(y-x)/(y-x) #### makes an error if x-y
    data_arr = [dead_pixels(pf.getdata(fil))[:,:1024] for fil in flat_files] 
    gain_arr = gain(bias, data_arr[0])
    for flat in data_arr[1:]:
        gain_arr += gain(bias, flat)
    return gain_arr/(len(flat_files))
def get_gain(hdr):
    """ returns the average gain based on the fits' filter """
    filt, day = hdr['filtnam'], hdr['DATE-BEG'][8:10]
    if filt == 'B': index = 0
    elif filt == 'V': index = 1
    elif filt == 'R': index = 2
    elif filt == 'I': index = 3
    if day == '14': return find_gain(all_flat_arr[0][index], get_bias(hdr))
    elif day == '20': return find_gain(all_flat_arr[1][index], get_bias(hdr))
    elif day == '21': return find_gain(all_flat_arr[2][index], get_bias(hdr))
#==================================
# Bias 
#==================================
def get_bias(hdr):
    """  returns the bias (corrected for pixel and 1024x1024) depending on the day """
    day = hdr['DATE-BEG'][8:10]
    if day == '14': index = 0 #for day 10 / 13 / 15
    elif day == '20': index = 1 #for day 10 / 19 / 15
    elif day == '21': index = 2 #for day 10 / 20 / 15
    bias = avg_fits_img(all_bias_arr[index]) # take all 14(?) bias files to create averaged one
    return dead_pixels(bias)[:,:1024] # correct for dead pixels and then make it [1024,1024]
#==================================

#==================================
# Flats
#==================================
## load the files normally (dead_pixels, and [:,:1024] are applied)
#flatB_img, flatB_hdr = loadfits(flatB_file_arr[0])
#flatV_img, flatV_hdr = loadfits(flatV_file_arr[0])
#flatR_img, flatR_hdr = loadfits(flatR_file_arr[0])
#flatI_img, flatI_hdr = loadfits(flatI_file_arr[0])
## average the fits -- does not correct them in anyway [1024,1056]; NO correction for dead pixels
#flatB_img_avg = avg_fits_img(flatB_file_arr)
#flatV_img_avg = avg_fits_img(flatV_file_arr)
#flatR_img_avg = avg_fits_img(flatR_file_arr)
#flatI_img_avg = avg_fits_img(flatI_file_arr)
## have to correct the avg fits [1024,1024]; YES correction for dead pixels
#flatB_img_avg = dead_pixels(flatB_img_avg)[:,:1024]
#flatV_img_avg = dead_pixels(flatV_img_avg)[:,:1024]
#flatR_img_avg = dead_pixels(flatR_img_avg)[:,:1024]
#flatI_img_avg = dead_pixels(flatI_img_avg)[:,:1024]
# gain 
#gainB = find_gain(flatB_file_arr, get_bias())
#gainV = find_gain(flatV_file_arr, get_bias())
#gainR = find_gain(flatR_file_arr, get_bias())
#gainI = find_gain(flatI_file_arr, get_bias())
#==================================

#==================================
def display_info(fil, display = False):
    """ fn that just outputs all the info on a file that may be useful """
    info = loadfits(fil)
    img, hdr = info[0], info[1]
    print "Getting info for", hdr['object']
    if display:
        print '****************'
        print "object:", hdr['object']
        print "exposure time:", hdr['exptime']
        print 'date begin:', hdr['DATE-BEG']
        print 'dec:', hdr['dec']
        print 'ra:', hdr['ra']
        print 'start time temperature (deg C):', hdr['tempdet']
        print 'filter:', hdr['filtnam']
        print 'cover:',hdr.get('cover')
        print 'epoch:', hdr['EQUINOXU']
        print '****************'
    return [hdr['object'], hdr['exptime'], hdr['DATE-BEG'], hdr['dec'], hdr['ra'], 
    hdr['tempdet'],  hdr['filtnam'], hdr.get('cover'), hdr['EQUINOXU']]
#info = display_info(AlineR_file_arr1[5])

#An HDU (Header Data Unit) is the highest level component of the FITS file structure, 
#consisting of a header and (typically) a data array or table.
#hdulist[0] is the primary HDU
# header includes useful info like size of the image (NAXIS1 & NAXIS2), the data
# (DATE-OBS), the exposure time (EXPTIME), where the celestial coordiantes of where 
# telescopes was pointing (RA & DEC) and the optical filter in the light path (FILTNAM)
#==================================
#==================================    
def correct_fits(fil):
    """Returns the corrected image"""
    img, hdr = loadfits(fil) # img is corrected for deadpixels and [1024,1024]
    return (img - get_bias(hdr)) * get_gain(hdr)
def normalize_img(img):
    return img / np.median(img)
    
def display_fits_raw_c(fil):
    """ Display raw and corrected FITS file """    
    hdr = loadfits(fil)[1]
    img = pf.open(fil)[0].data
    img = normalize_img(img)
    img_c = normalize_img(correct_fits(fil))
    
    fig = plt.figure(figsize = (12,5))
    fig.suptitle(hdr['object'] + ', Filter: ' + hdr['filtnam'], fontsize = 20)
    ax = fig.add_subplot(121)
    ax.set_title('Raw', fontsize = 18)
    ax.set_xlabel('Pixel', fontsize = 18);  ax.set_ylabel('Pixel', fontsize = 18)
    cax = ax.imshow(img,cmap='gray_r',origin='lower',vmin=0.8,vmax=1.8)
    fig.colorbar(cax, fraction=0.046, pad=0.04)
    ax1 = fig.add_subplot(122)
    ax1.set_title('Corrected', fontsize = 18)
    ax1.set_xlabel('Pixel', fontsize = 18)
    cax1 = ax1.imshow(img_c,cmap='gray_r',origin='lower',vmin=.8,vmax=1.8)
    fig.colorbar(cax1, fraction=0.046, pad=0.04)
    return fig


def display_all(fil):
    """ displays bias, flat, raw, AND corrected """
    hdr = loadfits(fil)[1]

    ### finding the flat
    def find_flat():
        filt, day = hdr['filtnam'], hdr['DATE-BEG'][8:10]
        print filt
        if filt == 'B': index = 0
        elif filt == 'V': index = 1
        elif filt == 'R': index = 2
        elif filt == 'I': index = 3
        if day == '14': flat_arr = all_flat_arr[0][index]
        elif day == '20': flat_arr = all_flat_arr[1][index]
        elif day == '21': flat_arr = all_flat_arr[2][index]
        return pf.open(flat_arr[0])[0].data
    flat_img = find_flat()
    flat_img = normalize_img(flat_img)
    ### finding the bias
    bias_img = get_bias(hdr)
    bias_img = normalize_img(bias_img)
    
    ### raw
    img_raw = pf.open(fil)[0].data
    img_raw = normalize_img(img_raw)
    
    ### correct
    img_correct = correct_fits(fil)
    img_correct = normalize_img(img_correct)
    
    fig = plt.figure(figsize = (10,10))
    fig.suptitle(hdr['object'], fontsize = 20)
    
    # bias
    ax = fig.add_subplot(221)
    ax.set_title('Bias', fontsize = 18)
    cax = ax.imshow(bias_img,cmap='gray_r',origin='lower',vmin=0.8,vmax=1.2)
    fig.colorbar(cax, fraction=0.046, pad=0.04)
    
    #flat
    ax1 = fig.add_subplot(222)
    ax1.set_title('Filter: ' + hdr['filtnam'], fontsize = 18)
    cax1 = ax1.imshow(flat_img,cmap='gray_r',origin='lower',vmin=.8,vmax=1.15)
    fig.colorbar(cax1, fraction=0.046, pad=0.04)
    
    # raw
    ax2 = fig.add_subplot(223)  
    ax2.set_title('Raw', fontsize = 18)
    cax2 = ax2.imshow(img_raw,cmap='gray_r',origin='lower',vmin=.8,vmax=1.8)
    fig.colorbar(cax2, fraction=0.046, pad=0.04)
    
    # corrected
    ax3 = fig.add_subplot(224)
    ax3.set_title('Corrected', fontsize = 18)
    cax3 = ax3.imshow(img_correct,cmap='gray_r',origin='lower',vmin=.8,vmax=2.5)
    fig.colorbar(cax3, fraction=0.046, pad=0.04)
    
    #one big plot colorbar
    #fig.subplots_adjust(right=0.85)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    #fig.colorbar(cax, cax=cbar_ax)

    plt.tight_layout()
    return fig
#plt.show(display_all(DanaeI_file_arr[0]))
#display_all(DanaeI_file_arr[0]).savefig(fig_path+'display_all_DanaeI0.png',dpi=200)

#==================================
#==================================

#==================================
#==================================
#=========== Centroids ============
#==================================
#==================================
def max_peaks(arr, width, lower_limit):
    """ Returns an array of the indexes where the indexes indicate where the peaks are"""
    #return [i for i in range(1, len(arr)-1) if arr[i-1]<arr[i] and arr[i+1]<arr[i] and arr[i]>lower_limit]
    return [i for i in range(1, len(arr)-1) if all(arr[i] > arr[i-width:i]) and all(arr[i]>arr[i+1:i+width]) and arr[i]>lower_limit] # 36 is arbitrary
    #assuming that there is nothing before 36

def find_peaks(fil):
    print 'Calculating the peaks for', pf.open(fil)[0].header['object']
    img = correct_fits(fil)
    avg = np.median(img)
    lim = 2.2*avg
    #arr[i][j] == [row][col] ~~~ [y][x] so when we plot them, we need to plot the y values, x values
    # find the initial pixel values that are above the limit (=constant*avg)
    # img_arr has the following syntax:
    # the index of img_arr corresponds to the row value of the pixel
    # img_arr[782] represents the 782th row of the array
    # when img_arr[782] = [653, 897], that means that pixels [782,653] and [782,897] are above the limit
    # if img_arr[456] = [], that means that in the 456th row, there are no pixels that are above the limit
    img_arr = [[]] #ignore the 0th row as not having anything
    for i_row in range(1, len(img)): #range starts at 1 bc of the assumption that img[0] has no peaks
        # apply max_peaks to each row (similar to Lab2)
        img_arr.append(max_peaks(img[i_row], width=10, lower_limit=lim))

    def cluster(img_arr):
        """ The return value cluster_arr is an array where it takes in img_arr and
        sees something like this (where the number correlates to the index of img_arr 
        and the array at that index)
        287 []; 288 []; 289 []; 290 [654, 786]; 291 [656]; 292 []; 293 []; 294 []
        and groups it like so [ [290, 654], [290,786], [291,656] ] -- this is one cluster
        cluster_arr is an array of multiple 'cluster's like above
        """
        cluster_arr, k = [], 0 # index k represents the current index in img_arr
        while k < len(img_arr)-1:
            curr = img_arr[k]
            if len(curr) >= 1:
                temp_arr = [[k,elem] for elem in curr]
                j=k+1
                while j < len(img_arr) and len(img_arr[j]) >= 1:
                    for elem in img_arr[j]:
                        temp_arr.append([j,elem])
                    j+=1
                k=j-1
                cluster_arr.append(temp_arr)
            k+=1
        return cluster_arr
        
    #for index, elem in enumerate(img_arr):
    #    print index, elem
        
    def boundary_limits(img_arr):  
    # to account for the top left corner and bottom right false positives  (first if)
    # to account for the right side where lots of false positives are recorded (second if)    
        index=0
        count=0
        while index < len(img_arr):
            if img_arr[index] and (index<=35 or index>=(1024-40)) and (img_arr[index][0]<=35):
                img_arr[index]=[]
            elif img_arr[index] and [i for i in img_arr[index] if i >=1010]:
                count+=1
            index+=1
        index=0
        if count > 50:
            while index < len(img_arr):
                if img_arr[index] and [i for i in img_arr[index] if i >=1010]:
                    for i in img_arr[index]:
                        if i>=1010: img_arr[index].remove(i)
                index+=1
        return img_arr
    img_arr = boundary_limits(img_arr)
    clusters = cluster(img_arr)

    #for cluster in clusters:
    #    print cluster
    
    def go_thru_clusters(clusters):
        new_cluster_arr =[]
        for cluster in clusters:
            #print '*** orig cluster', cluster
            def group(cluster):
                temp_arr = [cluster[0]]
                cluster.remove(cluster[0])
                curr_i = 0
                while cluster:
                    if np.abs(cluster[curr_i][1] - temp_arr[-1][1])<20 and np.abs(cluster[curr_i][0] - temp_arr[-1][0]) <=20:
                        # first cond: deals with the diff in cols
                        # second cond: deals with the diff in rows
                        temp_arr.append(cluster[curr_i])
                        cluster.remove(cluster[curr_i])
                        #print "TEMP", temp_arr
                        #print "CLUST", cluster
                        curr_i=0
                    else:
                        curr_i+=1
                        #print curr_i
                    if len(cluster) <= curr_i:
                        break
                #print '**temp ', temp_arr
                #print '***cluster', cluster
                new_cluster_arr.append(temp_arr)
            while cluster:
                group(cluster)
        return new_cluster_arr  
    
    # new cluster arr just groups the clusters way better and applies certain restrictions
    new_cluster_arr = go_thru_clusters(clusters)
    
    # === get rid of the clusters that are length one because those are just systematic errors
    # that were picked up/or pixels that were above the limit but have no correlating pixels
    # near them
    delete_arr=[]
    for index, elem in enumerate(new_cluster_arr):
        if len(elem) == 1:
            delete_arr.append(index)
    new_cluster_arr = np.delete(new_cluster_arr, delete_arr)

    peaks = [] # each element of peaks has the syntax of [[x,y],size]
    for cluster in new_cluster_arr:
        intensity_arr = [img[coords[0]][coords[1]] for coords in cluster]
        max_index = np.argmax(intensity_arr)
        peaks.append([cluster[max_index], len(cluster)])
    return peaks

def centroid(fil):
    """ 
    Purpose:
        To find the centroids of each peak, given the coordinates of each peak and 
        their corresponding size. To find the centroid, we simply compute the center
        of mass in both the x and y direction. To find x_cm and y_cm (center of mass),
        we find the x_cm of each row (or the y_cm of each col) and take the average.
        For example, if we have 3 rows and 3 columns, then we take the x_cm of the first,
        second, and third row individually and then take the average of the x_cm's (likewise 
        for y_cm)
    Input:
        fil: the uncorrected, raw file, that will be corrected
        bias: the bias 
    Output: 
        centroids_x: correlates w cols, all the x_cm's
        centroids_y: correlates w rows, all the y_cm's
        [centroids_x[i], centroids_y[i]] correlate to a centroid pixel 
    peaks has the form of [[x,y],size] 
    deciding to take into account +- size for centroiding
    """
    centroids =[]
    img = correct_fits(fil)
    peaks = find_peaks(fil) # find the peaks, elem of peaks is [[x,y],size]
    print 'Calculating the centroids for', pf.open(fil)[0].header['object']
    for pair, size in peaks:
        x_lower = 0 if pair[0]-size <= 0 else pair[0]-size
        x_upper = 1023 if pair[0]+size >=1024 else pair[0]+size
        y_lower = 0 if pair[1]-size <= 0 else pair[1]-size
        y_upper = 1023 if pair[1]+size >=1024 else  pair[1]+size
        x_vals = range(x_lower, x_upper) # just the ranging x values
        y_vals = range(y_lower, y_upper) # just the ranging y values
        x_arr =[]
        for y in y_vals:
            x_temp = []
            for x in x_vals:
                x_temp.append([x,y]) # keeps y constant for varying x; will be used for the x_cm
            x_arr.append(x_temp)
        
        y_arr=[]   
        for x in x_vals:
            y_temp = []
            for y in y_vals:
                y_temp.append([x,y]) # keeps y constant for varying y; will be used for the y_cm
            y_arr.append(y_temp)
        
        def find_cm(arr, vals):
            cm_arr=[]
            for elem in arr:
                range_img = [img[i[0],i[1]] for i in elem] # vals of img at constant y (or x), varying x (or y)
                centroid =  np.sum([a*b for a,b in zip(vals,range_img)])/np.sum(range_img) #<x> or <y>
                cm_arr.append(centroid)
            return np.sum(cm_arr)/len(cm_arr) # return the avg vals of the center of masses

        x_cm = find_cm(x_arr, x_vals)
        y_cm = find_cm(y_arr, y_vals)
        centroids.append([x_cm,y_cm])
    
    centroids = np.transpose(centroids)
    centroids_x, centroids_y = centroids[0], centroids[1]
    return centroids_x, centroids_y  # centroids_x ~ cols, centroids_y ~ rows, so when plotting in plt, 
    # need to plot plt.plot(centroids_y, centroids_x) - i think i just switched them along the way
        
def display_corrected_fits(fil, bool_centroids = False):
    """ Display fits with centroids (if bool_centroids is true)"""    
    hdr = pf.open(fil)[0].header
    img = normalize_img(correct_fits(fil))
    fig = plt.figure()
    plt.title(hdr['object']+ ', Filter: ' + hdr['filtnam'] + '\n# of centroids: ' + str(len(centroids_x)), fontsize=18)
    plt.imshow(img, origin='lower', interpolation='nearest', cmap='gray_r', 
               vmin=0.95*np.median(img), vmax=1.7*np.median(img))  
    if bool_centroids:
        centroids_x, centroids_y = centroid(fil) # centroids_x ~ cols, centroids_y ~ rows
        #plt.scatter(centroids_y, centroids_x, marker='x',s=100) #x's
        plt.scatter(centroids_y, centroids_x, s=150, facecolors='none', edgecolors='b') #circles
    plt.colorbar()
    return fig  


def display_corrected_fits2(fil, fil1, bool_centroids = False, bool_centroids1 = False):
    """Display 2 fits files with their centroids on the same figure"""    
    hdr, hdr1 = pf.open(fil)[0].header, pf.open(fil1)[0].header
    img, img1 = normalize_img(correct_fits(fil)), normalize_img(correct_fits(fil1)) # normalize and correct images
    
    fig = plt.figure(figsize = (12,5))    

    #first image
    ax = fig.add_subplot(121)
    ax.set_title(hdr['object'] + ' Filter: ' + hdr['filtnam'] + '\n# of centroids: ' + str(len(centroids_x)), fontsize=18)
    ax.set_xlabel('Pixel'); ax.set_ylabel('Pixel')
    cax = ax.imshow(img,cmap='gray_r',origin='lower',vmin=.8,vmax=1.8)
    fig.colorbar(cax, fraction=0.046, pad=0.04)
    if bool_centroids:
        centroids_x, centroids_y = centroid(fil) # centroids_x ~ cols, centroids_y ~ rows
        #ax.scatter(centroids_y, centroids_x, marker='x',s=100) #x's
        ax.scatter(centroids_y, centroids_x, s=170, facecolors='none', edgecolors='b') #circles
    
    #second image
    ax1 = fig.add_subplot(122)
    ax1.set_title(hdr1['object'] + ' Filter: ' + hdr1['filtnam'] + '\n# of centroids: ' + str(len(centroids_x1)), fontsize=18)
    ax1.set_xlabel('Pixel')
    cax1 = ax1.imshow(img1,cmap='gray_r',origin='lower',vmin=.8,vmax=1.8)
    fig.colorbar(cax1, fraction=0.046, pad=0.04)
    if bool_centroids1:
        centroids_x1, centroids_y1 = centroid(fil1) # centroids_x ~ cols, centroids_y ~ rows
        #ax1.scatter(centroids_y1, centroids_x1, marker='x',s=100) #x's
        ax1.scatter(centroids_y1, centroids_x1, s=170, facecolors='none', edgecolors='b') #circles
    
    return fig  
        
#plt.show(display_corrected_fits(GiselaI_file_arr[5]))
#plt.show(display_corrected_fits2(DanaeI_file_arr[5],DanaeI_file_arr1[5]))
#==================================
#==================================









#==================================
#==================================
# Star Matching
#==================================
#==================================
def usno(radeg,decdeg,fovam,epoch):
 # James R. Graham 2013/10/13 UC Berkeley
 # get USNO B-1 stars centered at radeg and decdeg (J2000.0) in degrees centered 
 # in a square field of view (arc min). Corrects for proper motion to current epoch.
    import urllib as url
    import string as str
    import pdb

    str1 = 'http://webviz.u-strasbg.fr/viz-bin/asu-tsv/?-source=USNO-B1'
    str2 = '&-c.ra={:4.6f}&-c.dec={:4.6f}&-c.bm={:4.7f}/{:4.7f}&-out.max=unlimited'.format(radeg,decdeg,fovam,fovam)
    string = str1+str2
    print 'Calling Vizier', string
    f = url.urlopen(string)
    # Read from the object, storing the page's contents in 's'.
    s = f.read()
    f.close()
 
    sl = s.splitlines()[45:-1] # get rid of header - updated Oct 2013
    name = np.array([])
    rad = np.array([]) # RA in degrees
    ded = np.array([]) # DEC in degrees
    rmag = np.array([]) # rmage
    for k in sl:
        kw = k.split('\t')
        ded0 = float(kw[2])
        pmrad = float(kw[6])/3600e3/np.cos(np.deg2rad(ded0)) # comvert from mas/yr to deg/year
        pmded = float(kw[7])/3600e3 
        name = np.append(name,kw[0])
        rad = np.append(rad,float(kw[1]) + pmrad*(epoch-2000.0) ) 
        ded = np.append(ded,float(kw[2]) + pmded*(epoch-2000.0) )
        if kw[12] != '     ': # case when no mag is reported
            rmag = np.append(rmag,float(kw[12])) 
        else:
            rmag = np.append(rmag,np.nan)
    return name,rad,ded,rmag
    
#==================================
#looking at star catalogs
#Read the position fromt the FITS file, and convert RA/DEC to degrees
#==================================
def USNOhelper(fil):
    info= display_info(fil)
    ra, dec, epoch = info[4], info[3], info[8]
    print 'ra', ra
    print 'dec', dec
    print 'epoch', epoch
    
    radeg = 15*(float(ra[0:2])+float(ra[3:5])/60.+float(ra[6:])/3600.)
    dsgn = np.sign(float(dec[0:2]))
    dedeg = float(dec[0:2])+dsgn*float(dec[3:5])/60. + dsgn*float(dec[6:])/3600.
    fovam = 6.3
    return radeg, dsgn, dedeg, fovam, epoch

def plotUSNO(fil):
    """ This function just plots the USNO catalog of stars, given a fits file"""
    radeg, dsgn, dedeg, fovam, epoch = USNOhelper(fil)
    name, rad, ded, rmag = usno(radeg,dedeg,fovam,epoch)
    #print 'name', name
    #print 'rad', rad
    #print 'ded', ded
    #print 'rmag', rmag
    hdr = loadfits(fil)[1]
    w = np.where(rmag <18.)[0]
    fig = plt.figure()
    plt.plot(rad[w],ded[w],'g.')
    plt.locator_params(axis='x',nbins=4); plt.locator_params(axis='y',nbins=4)
    plt.tick_params('x',pad=10)
    plt.xlabel('RA [Deg]'); plt.ylabel('Dec [Deg]')
    plt.ticklabel_format(useOffset=False)
    plt.axis('scaled')
    plt.gca().invert_xaxis()
    return fig

#plt.show(plotUSNO(DanaeI_file_arr[5]))

def plotUSNO2(fil,fil1):
    """ This function just plots the USNO catalog of stars, given a fits file"""
    radeg, dsgn, dedeg, fovam, epoch = USNOhelper(fil)
    name, rad, ded, rmag = usno(radeg,dedeg,fovam,epoch)
    radeg1, dsgn1, dedeg1, fovam1, epoch1 = USNOhelper(fil1)
    name1, rad1, ded1, rmag1 = usno(radeg1,dedeg1,fovam1,epoch1)

    w = np.where(rmag <18.)[0]
    fig = plt.figure(figsize=(12,5))
    ax = fig.add_subplot(121)
    ax.plot(rad[w],ded[w],'g.')
    ax.locator_params(axis='x',nbins=4); ax.locator_params(axis='y',nbins=4)
    ax.tick_params('x',pad=10)
    ax.set_xlabel('RA [Deg]'); ax.set_ylabel('Dec [Deg]')
    ax.ticklabel_format(useOffset=False)
    ax.axis('scaled')
    plt.gca().invert_xaxis()
    
    w1 = np.where(rmag1 <18.)[0]
    ax1 = fig.add_subplot(122)
    ax1.plot(rad1[w1],ded1[w1],'g.')
    ax1.locator_params(axis='x',nbins=4); ax1.locator_params(axis='y',nbins=4)
    ax1.tick_params('x',pad=10)
    ax1.set_xlabel('RA [Deg]'); ax1.set_ylabel('Dec [Deg]')
    ax1.ticklabel_format(useOffset=False)
    ax1.axis('scaled')
    plt.gca().invert_xaxis()
    return fig
plt.show(plotUSNO2(DanaeI_file_arr[5], DanaeI_file_arr1[5]))

###################### MELS CODE for mapcoord ########################
#catalog RA & dec
def mapcoord(fil):
    radeg, dsgn, dedeg, fovam,epoch = USNOhelper(fil)
    name, rad, ded, rmag = usno(radeg,dedeg,fovam,epoch)

    #finding ra0, dec0 - RA and DEC of center, which should be the object's coords.
    hdr = loadfits(fil)[1]
    #convert from time format to degrees
    ra0deg = 15*(float(hdr['RA'][0:2]) + float(hdr['RA'][3:5])/60. + float(hdr['RA'][6:])/3600.)
    findsign = np.sign(float(hdr['DEC'][0:2]))
    dec0deg = float(hdr['DEC'][0:2]) + findsign*float(hdr['DEC'][3:5])/60. + findsign*float(hdr['DEC'][6:])/3600.
    
    #convert from degrees to radians
    ra = np.deg2rad(rad); dec = np.deg2rad(ded)
    ra0 = np.deg2rad(ra0deg); dec0 = np.deg2rad(dec0deg)

    CCDxnum = np.cos(dec)*np.sin(ra - ra0)
    CCDxden = np.cos(dec0)*np.cos(dec)*np.cos(ra - ra0) + np.sin(dec)*np.sin(dec0)
    CCDynum = np.sin(dec0)*np.cos(dec)*np.cos(ra - ra0) - np.cos(dec0)*np.sin(dec)
    CCDyden = np.cos(dec0)*np.cos(dec)*np.cos(ra - ra0) + np.sin(dec)*np.sin(dec0)
    CCDx = -(CCDxnum/CCDxden)
    CCDy = -(CCDynum/CCDyden)

    f = 16840 #focal length in mm, from handout, what's our focal length?
    p = 0.015 #pixel size in mm, NOTE: DOES NOT ACCOUNT FOR X2 BINNING
    littlex0 = 512; littley0 = 512 # center on CCD (total 1024 pixels)
    
    littlex = f*(CCDx/(p*2))+littlex0
    littley = f*(CCDy/(p*2))+littley0
    return littlex, littley

def USNO_CCD(fil): 
    """ Plots the corrected image, with the centroids on top AND with the USNO catalog of stars
        Includes annotations
        Returns:
            the figure
            the x and y values from the USNO catalog
            the x and y values from the CCD/Centroids"""   
    radeg, dsgn, dedeg, fovam,epoch = USNOhelper(fil)
    name, rad, ded, rmag = usno(radeg,dedeg,fovam,epoch)
    hdr = loadfits(fil)[1]
    
    y_cenvals, x_cenvals = centroid(fil) # centroids_x (y_cenvals) ~ cols (correlate w y), centroids_y (x_cenvals)~ rows (correlate with x)
    
    littlex, littley = mapcoord(fil)
    mag_lim = 18
    w = np.where(rmag <mag_lim)[0]
    x_usnovals = littlex[w]
    y_usnovals = (1024-littley[w]) 
    
    img_c = correct_fits(fil)
    img_c = normalize_img(img_c)

    fig = plt.figure()
    plt.plot(x_usnovals, y_usnovals, 'ro', ms=6, label='USNO points')
    plt.plot(x_cenvals, y_cenvals, 'bo', ms=6, label = 'CCD/Centroid points')
    plt.imshow(img_c, origin='lower', interpolation='nearest', cmap='gray_r', 
                vmin=.8, vmax=1.8)
    plt.locator_params(axis='x',nbins=4); plt.locator_params(axis='y',nbins=4)
    plt.tick_params('x',pad=10) 
    plt.xlabel('Pixel', fontsize = 18); plt.ylabel('Pixel', fontsize = 18)
    plt.title(hdr['OBJECT'] + '\nUSNO catalog (mag. <'+str(mag_lim) + ') vs. CCD centroids', fontsize = 20)
    plt.ticklabel_format(useOffset=False)
    plt.axis('scaled')
    plt.ylim([0,1024]); plt.xlim([0,1024])
    plt.legend(loc = 'best', fontsize = 14, framealpha=0.5)
    for i, txt in enumerate(range(len(x_usnovals))):
        plt.annotate(txt, (x_usnovals[i],y_usnovals[i]), fontsize = 12, color = 'r')
    for i, txt in enumerate(range(len(y_cenvals))):
        plt.annotate(txt, (x_cenvals[i],y_cenvals[i]), fontsize = 12, color = 'b')
    return [fig, x_usnovals, y_usnovals, x_cenvals, y_cenvals]

#USNO_CCD_info = USNO_CCD(DanaeI_file_arr[-1])
#plt.show(USNO_CCD_info[0])

def USNO_CCD_display(fil):
    USNO_CCD_info = USNO_CCD(fil)
    x_usnovals = USNO_CCD_info[1]; y_usnovals = USNO_CCD_info[2] 
    x_cenvals = USNO_CCD_info[3];  y_cenvals = USNO_CCD_info[4]
    img_c = correct_fits(fil)
    img_c = normalize_img(img_c)
    hdr = loadfits(fil)[1]
    fig = plt.figure()
    plt.plot(x_usnovals, y_usnovals, 'ro', ms=6, label='USNO points')
    plt.plot(x_cenvals, y_cenvals, 'bo', ms=6, label = 'CCD/Centroid points')
    #plt.imshow(img_c, origin='lower', interpolation='nearest', cmap='gray_r', 
    #            vmin=.8, vmax=1.8)
    plt.locator_params(axis='x',nbins=4); plt.locator_params(axis='y',nbins=4)
    plt.tick_params('x',pad=10) 
    plt.xlabel('Pixel', fontsize = 18); plt.ylabel('Pixel', fontsize = 18)
    plt.title(hdr['OBJECT'] + '\nUSNO catalog vs. CCD centroids', fontsize = 20)
    plt.ticklabel_format(useOffset=False)
    plt.axis('scaled')
    plt.ylim([0,1024]); plt.xlim([0,1024])
    plt.legend(loc = 'best', fontsize = 14, framealpha=0.5)
    for i, txt in enumerate(range(len(x_usnovals))):
        plt.annotate(txt, (x_usnovals[i],y_usnovals[i]), fontsize = 12, color = 'r')
    for i, txt in enumerate(range(len(y_cenvals))):
        plt.annotate(txt, (x_cenvals[i],y_cenvals[i]), fontsize = 12, color = 'b')
    return fig


#plt.show(USNO_CCD_display(DanaeI_file_arr[5]))



























#==================================
#==================================
#MATCHING THE COORDINATES
#==================================
#==================================

#==================================
# Choosing which points to fit
#==================================
def point_chooser():
    """
    Purpose: To be able to visually match the CCD centroids with the USNO objects 
            (this is a semi-automatic/semi-manual process)
            Once you input all the information, the program automatically makes 
            it into a list to be used later on in determining which centroids and 
            USNO objects will be used for the LLS
    Input: the fn prompts the user to input the number for the CCD centroids and 
            the USNO object. 
    Output: 
            chosenpoints_arr - a list of pairs of indices, ex. [[0,5], [5,3]] CCD[0] 
                correlates w USNO[5] and CCD[5] w USNO[3]. the points/indices 
                correspond to the annotated indices from the plot in USNO_CCD(fil)[0] 
    """
    print '**********************************************'
    print 'Correlate points with eachother'
    print 'First number is the CCD/centroid number (blue)'
    print 'Second number is the USNO number (red)'
    print 'When you have no more correlations, input \'done\' in lowercase to finish.'
    chosenpoints_arr =[] 
    x = raw_input('Input CCD: ')
    temp_arr = []
    while x != 'done':
        #assume that only numbers or 'done' are inputted -- sort of crashes otherwise
        temp_arr.append(int(x))
        if len(temp_arr) == 2:
            chosenpoints_arr.append(temp_arr)
            temp_arr =[]
            x = raw_input('Input CCD: ')
        else:
            x = raw_input('Input USNO: ')
    return chosenpoints_arr

def fullfit(fil, chosenpoints_arr = []):
    """
    Purpose: To obtain a general Least Squares fit for USNO catalog and CCD centroids
    Input:
        fil - the filename
        chosenpoints_arr - the array of the points we end up considering for the star matching
                            if [], then the user will be asked to input this arr
    Output: [a11, a12, x0, a21, a22, y0, T_fp, T, chosen_ccdx, chosen_ccdy, origx, origy]
        a11, a12, x0, a21, a22, y0 - the plate constants
        T_fp, T - the T matrix, with and without fp
        chosen_ccdx, chosen_ccdy - the x and y pixels of the CCD centroids that we matched
        origx, origy - the x and y pixels of the USNO objects that we matched
    """
    print 'Calling fullfit on', pf.open(fil)[0].header['OBJECT']
    USNO_CCD_info = USNO_CCD(fil)
    plt.show(USNO_CCD_info[0]) # so we can match the points visually
    USNOx = USNO_CCD_info[1] # all USNO x values
    USNOy = USNO_CCD_info[2] #all USNO y values
    CCDx = USNO_CCD_info[3] # all CCD x values
    CCDy = USNO_CCD_info[4] # all CCD y values
    
    if not chosenpoints_arr:
        chosenpoints_arr = point_chooser() #choose the points we want to fit (visually/manually)
        print chosenpoints_arr # so we can copy and paste if need be
    chosen_centroids_indexes = np.transpose(chosenpoints_arr)[0] # to get the indexes for the CCD points
    chosen_usno_indexes = np.transpose(chosenpoints_arr)[1]      # to get the indexes for the USNO points

    chosen_ccdx = CCDx[chosen_centroids_indexes] # apply the indexes above 
    chosen_ccdy = CCDy[chosen_centroids_indexes] # to get the chosen indexes 
    chosen_usnox = USNOx[chosen_usno_indexes]
    chosen_usnoy = USNOy[chosen_usno_indexes]
    origx = np.array(chosen_usnox) - 512 #not quite sure why we have to subtract 512, but it works when we do
    origy = np.array(chosen_usnoy) - 512
    
    fp = 16840./(2*.015) #f/2p
    fovx = origx; fovy = origy
    ax = np.transpose(np.matrix(chosen_ccdx))
    Bx = np.transpose(np.matrix([fovx,fovy,np.ones(len(fovx))]))
    Btransposex = np.transpose(Bx)
    btbx = Btransposex*Bx
    pseudoinvx = np.linalg.inv(btbx)
    cx = pseudoinvx*Btransposex*ax

    ay = np.transpose(np.matrix(chosen_ccdy))
    By = np.transpose(np.matrix([fovx,fovy,np.ones(len(fovy))]))
    Btransposey= np.transpose(By)
    btby = Btransposey*By
    pseudoinvy = np.linalg.inv(btby)
    cy = pseudoinvy*Btransposey*ay
    
    chisqx = np.transpose(ax - Bx*cx) * (ax - Bx*cx)
    chisqy = np.transpose(ay - By*cy) * (ay - By*cy)

    cfitx = np.linalg.inv((np.transpose(Bx))*Bx)*(np.transpose(Bx)*ax)
    cfity = np.linalg.inv((np.transpose(By))*By)*(np.transpose(By)*ay)
    a11 = float(cx[0]); a12 = float(cx[1]); x0 = float(cx[2])
    a21 = float(cy[0]); a22 = float(cy[1]); y0 = float(cy[2])
    T_fp = [[fp*a11,fp*a12,x0],[fp*a21,fp*a22,y0],[0,0,1]]
    T = [[a11,a12,x0],[a21,a22,y0],[0,0,1]]
    print 'T', T
    print 'T_fp', T_fp
    print 'sq of det of T', np.sqrt(np.linalg.det(T))
    print 'sq of det of T_fp', np.sqrt(np.linalg.det(T_fp))
    return [a11, a12, x0, a21, a22, y0, T_fp, T, chosen_ccdx, chosen_ccdy, origx, origy]

def plotfinalfit_res(fil, chosenpoints_arr=[]):
    """
    Purpose: plot the results of using a general LLS, and plot the residuals of it
    Input:
        fil - name of file
        chosenpoints_arr - the array of the points we end up considering for the star matching
                            if [], then the user will be asked to input this arr
    Output:
        fig_fit - the fitted figure
        fig_residuals - the figure with all the residuals
        finalxarr/finalyarr - the fitting of the USNO values on our CCD
        resx/resy - the residuals for our fit, in x and y
        rmsx/rmsy - the rms for x and y
        fig3 - the combined fig_fit and fig_residuals
    """
    a11, a12, x0, a21, a22, y0, T_fp, T, chosen_ccdx, chosen_ccdy, origx, origy = fullfit(fil,chosenpoints_arr)
    #fp = 16840./(2*.015) #not using this, using 1, otherwise accounts for it twice
    finalxarr, finalyarr = [], [] # fitting the USNO values to fit with our CCD values
    k = 0
    while k < len(origx):
        finalxarr.append(1*(a11*origx[k] + a12*origy[k] + x0))
        finalyarr.append(1*(a21*origx[k] + a22*origy[k] + y0))
        k+=1
    # First figure is the fitted plot (pixel vs pixel)
    fig_fit = plt.figure() 
    plt.plot(finalxarr, finalyarr, 'ko', ms = 4, label = 'Fitted USNO values')
    plt.plot(chosen_ccdx, chosen_ccdy,'rx', ms = 10, label = 'CCD values')
    plt.xlim([0,1024]); plt.ylim([0,1024])
    plt.title('Final fit for USNO and CCD centroids', fontsize = 18)
    plt.legend()
    
    # residuals
    resx = chosen_ccdx - finalxarr
    resy = chosen_ccdy - finalyarr
    rms_x = np.sqrt(np.mean(resx**2))
    rms_y = np.sqrt(np.mean(resy**2))
    
    # Second figure is the residuals plot
    fig_residuals = plt.figure(figsize=(12,5))
    plt.plot(chosen_ccdx,resx,'r^', ms = 6, label = 'x')
    plt.plot(chosen_ccdy,resy,'bo', ms = 6, label = 'y')
    plt.title('Residuals for USNO catalog and CCD values, ' + str(len(resx)) + ' points', fontsize = 18)
    plt.ylabel('Pixel offset', fontsize = 16)
    plt.xlabel('Pixel, either x or y', fontsize = 16)
    plt.xlim([0,1024]); #plt.ylim([-2,2])
    plt.axhline(y=0, color='k', ls='--')
    plt.legend()
    
    # third plot is both the fitted figure and the residuals figure
    fig3 = plt.figure(figsize=(12,5))
    fig3.suptitle("USNO vs. CCD centroids", fontsize = 20)
    ax = fig3.add_subplot(121)
    ax.plot(finalxarr, finalyarr, 'ko', ms = 4, label = 'Fitted USNO values')
    ax.plot(chosen_ccdx, chosen_ccdy,'rx', ms = 10, label = 'CCD values')
    ax.set_xlim([0,1024]); ax.set_ylim([0,1024])
    ax.set_title('Final fit', fontsize = 18)
    ax.set_ylabel('Pixel', fontsize = 16); ax.set_xlabel('Pixel', fontsize = 16)
    ax.legend()
    
    ax1 = fig3.add_subplot(122)
    ax1.plot(chosen_ccdx,resx,'r^', ms = 6, label = 'x')
    ax1.plot(chosen_ccdy,resy,'bo', ms = 6, label = 'y')
    ax1.set_title('Residuals for ' + str(len(resx)) + ' points', fontsize = 18)
    ax1.set_ylabel('Pixel offset', fontsize = 16); ax1.set_xlabel('Pixel, either x or y', fontsize = 16)
    ax1.set_xlim([0,1024])
    ax1.axhline(y=0, color='k', ls='--')
    ax1.legend()
    fig3.tight_layout()
    
    return fig_fit, fig_residuals, finalxarr, finalyarr, resx, resy, rms_x, rms_y, fig3

def mapback(fil, index, chosenpoints_arr = []):
    """ 
    Purpose:
        Find RA and DEC of the asteroid
    Input:
        fil - the filename
        index - the index of where the asteroid is: the index is in the original centroids_x/centroids_y arr
        chosenpoints_arr - the array of the points we end up considering for the star matching
                            (only need it since the fullfit fn is called)
    Output:
        alpha - RA in degrees
        delta - DEC in degrees
    """    
    hdr = loadfits(fil)[1]
    #convert the hdr center from time format to degrees
    ra0deg = 15*(float(hdr['ra'][0:2]) + float(hdr['ra'][3:5])/60. + float(hdr['ra'][6:])/3600.)
    findsign = np.sign(float(hdr['dec'][0:2]))
    dec0deg = float(hdr['dec'][0:2]) + findsign*float(hdr['dec'][3:5])/60. + findsign*float(hdr['dec'][6:])/3600.
    #convert from degrees to radians
    ra0 = np.radians(ra0deg);     dec0 = np.radians(dec0deg)
    
    T_fp = fullfit(fil, chosenpoints_arr)[6] # find T_fp 
    invT_fp = np.linalg.inv(T_fp)
    x_arr, y_arr = USNO_CCD(fil)[3:5] # gets all the original x and y centroids of the CCD
    asterx = x_arr[index]     #CCD x
    astery = y_arr[index]     #CCD y
    transform = invT_fp*np.transpose(np.matrix([asterx, astery, 1]))
    bigx = transform[0]
    bigy = transform[1]
    
    tanalp = np.arctan((-bigx)/(np.cos(dec0) - bigy*np.sin(dec0)))
    alpha = tanalp + ra0
    alpha = np.degrees(alpha)

    sindelnum = (np.sin(dec0) + bigy*np.cos(dec0))
    sindelden = (np.sqrt(1+bigx**2+bigy**2))
    delta = np.arcsin(sindelnum/sindelden)
    delta = np.degrees(delta)
    return alpha.item(0), delta.item(0) #some reason they are in a double matrix? # in degrees

chosenpoints_arrDI = [[19, 28], [16, 29], [20, 26], [17, 25], [15, 27], [14, 15], 
                        [18, 30], [1, 20], [5, 18], [3, 17], [0, 13], [2, 11], 
                        [4, 1], [11, 0], [9, 6], [8, 4], [12, 7]]
chosenpoints_arrDI1 = [[18, 14], [17, 11], [14, 16], [9, 13], [4, 15], [2, 9], [12, 8], 
                        [24, 18], [21, 0], [16, 4], [15, 2], [10, 6], [8, 7], [3, 5]]
chosenpoints_arrDI2 = [[18, 9], [19, 0], [12, 12], [15, 4], [16, 2], [17, 1], [6, 24], 
                        [0, 25], [1, 18], [4, 17], [9, 16]]
chosenpoints_arrAI = [[14, 19], [15, 17], [2, 3], [10, 12], [8, 11], [9, 7], 
                        [5, 5], [3, 9], [0, 10], [4, 1], [16, 15]] # not [11,13]
chosenpoints_arrAI1 = [[0, 5], [1, 3], [2, 0]] #not really good, residuals at 6pixels
chosenpoints_arrAI2 = [[17,5], [9, 21], [6, 19], [0, 20], [2, 15], [16, 2], [21, 0], 
                        [20, 1], [19, 4], [18, 3],  [8, 17]]          
chosenpoints_arrGI = [[12, 2], [14, 4], [18, 1], [13, 0], [2, 9], [4, 8], [6, 6], 
                        [1, 7], [8, 3], [7, 5]]
chosenpoints_arrGI1 = [[22, 15], [23, 13], [21, 10], [20, 4], [19, 8], [16, 7],
                        [15, 5], [13, 9], [5, 6], [6, 3]]
chosenpoints_arrGI2 = [[9, 14], [8, 5], [7, 7], [2, 8], [5, 3], [3, 2], [0, 1],
                        [1, 0], [12, 12], [14, 11]]

#===============================================================================
#===============================================================================
#===============================================================================
#=============================== PARALLAX ======================================
#===============================================================================
#===============================================================================
#===============================================================================

#===================== constants =====================
m_sun = 1 # in solar masses
G = 4*np.pi**2 # AU^3 * yr^-2 * solarmass^-1
k = np.sqrt(G*m_sun) / 365. # AU^3/2 * days^-1
#===================== files =========================
firstfil = DanaeI_file_arr[5]
secondfil = DanaeI_file_arr1[5]
thirdfil = DanaeI_file_arr2[5]
#firstfil = AlineI_file_arr[5]
#secondfil = AlineI_file_arr1[5]
#thirdfil = AlineI_file_arr2[5]
#firstfil = GiselaI_file_arr[5]
#secondfil = GiselaI_file_arr1[5]
#thirdfil = GiselaI_file_arr2[5]
#===================== setting the chosenpoint arrays ==========================
chosenpoints_arr1 = chosenpoints_arrDI
chosenpoints_arr2 = chosenpoints_arrDI1
chosenpoints_arr3 = chosenpoints_arrDI2
#============== setting the indexes of where the asteroids are =================
index1, index2, index3 = 10, 13, 7 # DanaeI 5 the indexes of where the centroid is-input this index to the centroids_x, centroids_y
#index1, index2, index3 = 7,3,14 # AlineI 5 the indexes of where the centroid is-input this index to the centroids_x, centroids_y
#index1, index2, index3 = 11, 14, 6 #GiselaI 5
# =========== finding alpha and delta of the asteroids from the three days ============
#================ made diff variables for in RA/DEC and in radians =====================
'''
# the alpha and deltas for corresponding files
alpha1_deg, delta1_deg = mapback(firstfil, index1, chosenpoints_arr1)
alpha2_deg, delta2_deg = mapback(secondfil, index2, chosenpoints_arr2)
alpha3_deg, delta3_deg = mapback(thirdfil, index3, chosenpoints_arr3)
'''
#Danae I 5
alpha1_deg, delta1_deg =8.70125315478, 27.3414416801
alpha2_deg, delta2_deg =7.21409509162, 27.1301990039
alpha3_deg, delta3_deg  =7.00514806729, 27.082107484

print '$$$$$$$$$$$$$$$$$$$$'
print 'alpha and delta before radians'
print 'alpha1 and delta1', alpha1_deg,  delta1_deg
print 'alpha2 and delta2', alpha2_deg,  delta2_deg
print 'alpha3 and delta3', alpha3_deg,  delta3_deg

alpha1, delta1 = math.radians(alpha1_deg), math.radians(delta1_deg)#convert to radians
alpha2, delta2 = math.radians(alpha2_deg), math.radians(delta2_deg)
alpha3, delta3 = math.radians(alpha3_deg), math.radians(delta3_deg)
print '$$$$$$$$$$$$$$$$$$$$'
print 'alpha and delta after radians'
print 'alpha1', alpha1
print 'delta1', delta1
print 'alpha2', alpha2
print 'delta2', delta2
print 'alpha3', alpha3
print 'delta3', delta3
print '$$$$$$$$$$$$$$$$$$$$'
epsilon = np.radians(23.43929111) #based on being in the J2000 equinox
#============ Get all the info of each file, just to display it ================
# what display_info(fil) returns
#return [hdr['object'], hdr['exptime'], hdr['DATE-BEG'], hdr['dec'], hdr['ra'], 
#hdr['tempdet'],  hdr['filtnam'], hdr.get('cover'), hdr['EQUINOXU']]
one, two, three = [display_info(fil) for fil in [firstfil, secondfil, thirdfil]]
print '$$$$$$$$ INFO $$$$$$$$'
print one
print two
print three
print '$$$$$$$$$$$$$$$$$$$$$$'

#================ Find the dates of each file ==================
firstdate, seconddate, thirddate = [loadfits(fil)[1]['DATE'] for fil in [firstfil, secondfil, thirdfil]]
print '$$$$$$$$ DATE $$$$$$$$$'
print firstdate
print seconddate
print thirddate
print '$$$$$$$$$$$$$$$$$$$$$$$'
#====== Reformat the dates so we can mathematically manipulate them to =======
#============== them into days ~~ they have the format of ====================
#=============== [year, month, day, hour, minute, second] ====================
#================ example - [2015, 10, 14, 4, 13, 39.67] =====================
firstdate_new, seconddate_new, thirddate_new = [[int(arr[:4]), int(arr[5:7]), int(arr[8:10]), int(arr[11:13]), 
                    int(arr[14:16]), float(arr[17:])] for arr in [firstdate, seconddate, thirddate]]
print '$$$ Reformated Date $$$'
print firstdate_new 
print seconddate_new 
print thirddate_new
print '$$$$$$$$$$$$$$$$$$$$$$$'
#==== Set the base as the first day of October at time 00:00:00.0 =============
t_1, t_2, t_3 = [x[2] + x[3]/24. + x[4]/(60.*24) + x[5]/(3600*24.) for x in [firstdate_new, seconddate_new, thirddate_new]]
print '$$$$$$$$$ Time $$$$$$$$'
print t_1
print t_2
print t_3
#========= Finding the Taus (time differences) =========================
tau_1 = t_2 - t_1
tau_3 = t_3 - t_2
print 'tau1', tau_1
print 'tau3', tau_3
print '$$$$$$$$$$$$$$$$$$$$$$'

#============= Functions =====================
def get_xyzeq(alpha, delta):
    """ get x_eq, y_eq, and z_eq"""
    x_eq = np.cos(alpha)*np.cos(delta)
    y_eq = np.sin(alpha)*np.cos(delta)
    z_eq = np.sin(delta)
    return x_eq,y_eq,z_eq
def get_XYZeq(alpha0, delta0):
    """ get X_eq, Y_eq, and Z_eq (Earth to Sun)"""
    return get_xyzeq(alpha0, delta0)
    
def get_s(alpha,delta,ep=epsilon):
    """ returns the s vector """
    x_eq, y_eq, z_eq = get_xyzeq(alpha, delta)
    xyz_eq = np.transpose(np.matrix([x_eq,y_eq,z_eq]))

    cose = np.cos(ep)
    sine = np.sin(ep)
    T = np.matrix([[1.,0.,0.],[0.,cose,sine],[0.,-sine,cose]])

    xyz = T*xyz_eq
    x = float(xyz[0]); y = float(xyz[1]); z = float(xyz[2])
    s = np.array([x,y,z])
    return s

def s2dot(s1,s2,s3,tau_1,tau_3):
    s_dot = ((tau_3*(s2-s1))/(tau_1*(tau_1+tau_3))) + ((tau_1*(s3-s2))/(tau_3*(tau_1+tau_3)))
    return s_dot
def s2ddot(s1,s2,s3,tau_1,tau_3):
    s_ddot = ((2.*(s3-s2))/(tau_3*(tau_1 + tau_3))) - ((2.*(s2-s1)) / (tau_1*(tau_1 + tau_3)))
    return s_ddot
def findasteroid(R, s2, s_dot, s_ddot):
    eps = 1 
    R_len = np.sqrt(np.dot(R, R))
    rguess = 2
    r_cross_s= np.cross(R, s2)
    s_cross_s = np.cross(s_ddot, s2)
    numerator = np.dot(s_dot, r_cross_s)
    denom = np.dot(s_dot, s_cross_s)
    while eps > 10e-21:
        rho=k**2*((1/R_len**3 ) - (1./rguess**3))*(numerator/denom)
        rsq = rho**2+R_len**2 + 2*rho*np.dot(R,s2)
        eps = abs(np.sqrt(rsq)-rguess)
        rguess = np.sqrt(rsq)
    return rguess, rho

#========= Getting XYZ from the RA and DEC of the earth to sun =================
#====================== (didn't really work :/) ================================
#got all the ra0 and dec0 from HORIZONS http://ssd.jpl.nasa.gov/horizons.cgi#top 
#======== Observer Location = Sun -- 500@10; target body = Earth -- 399 ========
def ra_dec_to_deg(ra, dec):
    radeg = 15*(float(ra[0:2])+float(ra[3:5])/60.+float(ra[6:])/3600.)
    dsgn = np.sign(float(dec[0:2]))
    decdeg = float(dec[0:2])+dsgn*float(dec[3:5])/60. + dsgn*float(dec[6:])/3600.
    return radeg, decdeg
ra01, dec01 = '01:14:16.00', '07:51:31.2' # 2015-Oct-14 00:00
ra02, dec02 = '01:36:39.53', '10:03:45.5' # 2015-Oct-20 00:00 
ra03, dec03 = '01:40:25.50', '10:25:19.0' # 2015-Oct-21 00:00
alpha01, delta01 = ra_dec_to_deg(ra01, dec01)
alpha02, delta02 = ra_dec_to_deg(ra02, dec02)
alpha03, delta03 = ra_dec_to_deg(ra03, dec03)
alpha01, delta01 = math.radians(alpha01), math.radians(delta01)
alpha02, delta02 = math.radians(alpha02), math.radians(delta02)
alpha03, delta03 = math.radians(alpha03), math.radians(delta03)
#x_eq, y_eq, z_eq = get_s(alpha2, delta2)
#X_eq, Y_eq, Z_eq = get_s(alpha02, delta02)
#rho_sq = (x_eq-X_eq)**2. + (y_eq-Y_eq)**2. + (z_eq-Z_eq)**2.
#rho = np.sqrt(rho_sq)

#=========================== find the asteroid! ================================
s1 = get_s(alpha1, delta1, epsilon)
s2 = get_s(alpha2, delta2, epsilon)
s3 = get_s(alpha3, delta3, epsilon)
s2dot = s2dot(s1, s2, s3, tau_1, tau_3)
s2ddot = s2ddot(s1, s2, s3, tau_1, tau_3)
earth_sun_distance = [8.946655614799660E-01 ,4.376048778015129E-01,-2.270836795464193E-05] #10/20/15
r, rho = findasteroid(earth_sun_distance, s2, s2dot, s2ddot)
print '$$$$$$$$$$$'
print 's1', s1
print 's2', s2
print 's3', s3
print 's2dot', s2dot
print 's2ddot', s2ddot
print 'r:', r
print 'rho:', rho
print '$$$$$$$$$$$'

#===============================================================================
#================================ END OF PARALLAX ============================== 
#===============================================================================

#===============================================================================
#================== Fitting the star matching by first linear fit ==============
#=================== (will be able to see how there is rotation) ===============
#========================= used the LLS fn from Lab 2 ==========================
#===============================================================================

def linleastsquares(data, poly):
    """
    data: [x,y] what you are trying to fit
    poly: the number of terms you want i.e. poly = 3 --> ax**2 + bx + c
    """
    sum_x_arr = [np.sum(data[0]**i) for i in np.arange(1,(poly-1)*2+1)]
    sum_xy_arr = [np.sum(data[0]**i * data[1]) for i in np.arange(0, poly)]
    
    left = [[None for i in range(poly)] for j in range(poly)]
    right = [[i] for i in sum_xy_arr]
    for i in range(poly):
        for j in range(poly):
            if i == 0 and j == 0:
                left[i][j] = len(data[0])
            else:
                left[i][j] = sum_x_arr[i+j-1]
    inv_left = np.linalg.inv(left)
    final = np.dot(inv_left,right)
    return final

#========================= Plotting the first linear fit! ======================
#============= using the LLSgen, fit the x and y data separately ===============
chosenpoints_arr = chosenpoints_arrDI
fil = DanaeI_file_arr[5]
#NON-ROTATION
def plotmatch(fil, chosenpoints_arr):
    USNO_CCD_info = USNO_CCD(fil)

    USNOx = USNO_CCD_info[1] # all USNO x values
    USNOy = USNO_CCD_info[2] #all USNO y values
    CCDx = USNO_CCD_info[3] # all CCD x values
    CCDy = USNO_CCD_info[4] # all CCD y values
    
    chosen_centroids_indexes = np.transpose(chosenpoints_arr)[0] # to get the indexes for the CCD points
    chosen_usno_indexes = np.transpose(chosenpoints_arr)[1]      # to get the indexes for the USNO points
    
    chosen_ccdx = CCDx[chosen_centroids_indexes] # apply the indexes above 
    chosen_ccdy = CCDy[chosen_centroids_indexes] # to get the chosen indexes 
    chosen_usnox = USNOx[chosen_usno_indexes]
    chosen_usnoy = USNOy[chosen_usno_indexes]
    
    cfitx, slopex = linleastsquares([chosen_ccdx,chosen_usnox],2)
    cfity, slopey = linleastsquares([chosen_ccdy,chosen_usnoy],2)
    #newcent means the new fitted centroids from the CCD after matching
    newcentx = chosen_ccdx*slopex+cfitx
    newcenty = chosen_ccdy*slopey+cfity
    fig = plt.figure()
    plt.plot(newcentx, newcenty, 'b^', ms = 6, label = 'Fitted CCD centroids')
    plt.plot(chosen_usnox, chosen_usnoy, 'ro', ms = 6, label = 'Mapped USNO values')
    plt.title('Matched centroids after non-rotation LLS of x and y', fontsize = 20)
    plt.xlabel('Pixel', fontsize = 16); plt.ylabel('Pixel', fontsize = 16)
    plt.legend()
    plt.ylim([0,1024]); plt.xlim([0,1024])
    return fig
        
#plt.show(plotmatch(DanaeI_file_arr2[5], chosenpoints_arrDI2))

#===============================================================================
#===================== Plotting the trajectory of the asteroid =================
#===============================================================================
#Gisela I 5
#alpha1_deg, delta1_deg =2.36878624168, 7.28901377175
#alpha2_deg, delta2_deg =1.40495296043, 6.54222653246
#alpha3_deg, delta3_deg  =1.27175373155, 6.42141439152

#Danae I 5
alpha1_deg, delta1_deg =8.70125315478, 27.3414416801
alpha2_deg, delta2_deg =7.21409509162, 27.1301990039
alpha3_deg, delta3_deg  =7.00514806729, 27.082107484

def display_asteroid_traj(alpha1_deg, delta1_deg, alpha2_deg, delta2_deg, alpha3_deg, delta3_deg):
    fig = plt.figure(figsize=(12,5))
    plt.annotate('10/13', (alpha1_deg, delta1_deg), fontsize = 12)
    plt.annotate('10/19', (alpha2_deg, delta2_deg), fontsize = 12)
    plt.annotate('10/20', (alpha3_deg, delta3_deg), fontsize = 12)
    x = [alpha1_deg, alpha2_deg, alpha3_deg]
    y = [delta1_deg, delta2_deg, delta3_deg]
    plt.scatter(x, y)
    z = linleastsquares([x, y], 3) # fit the trajectory to a quadratic
    x1 = np.linspace(min(x),max(x),100)
    plt.plot(x1, z[0] +z[1]*x1 + z[2]*x1*x1)
    plt.gca().invert_xaxis()
    plt.title('Trajectory of Danae', fontsize =20)
    plt.xlabel('RA [deg]', fontsize =18)
    plt.ylabel('DEC [deg]', fontsize =18)
    return fig
#plt.show(display_asteroid_traj(alpha1_deg, delta1_deg, alpha2_deg, delta2_deg, alpha3_deg, delta3_deg))

#===============================================================================
#=============================== GET Functions =================================
#===============================================================================
def get_plate_constants(fil, chosenpoints_arr=[]):
    a11, a12, x0, a21, a22, y0, T_fp, T, chosen_ccdx, chosen_ccdy, origx, origy = fullfit(fil, chosenpoints_arr)
    print '####################'
    print 'a11', a11
    print 'a12', a12
    print 'a21', a21
    print 'a22', a22
    print 'x0', x0
    print 'y0', y0
    print 'sq of det of T', np.sqrt(np.linalg.det(T))
    print '####################'

def get_rms(fil, chosenpoints_arr = []):
    fig_fit, fig_residuals, finalxarr, finalyarr, resx, resy, rms_x, rms_y, fig3 = plotfinalfit_res(fil, chosenpoints_arr)
    print 'rms_x', rms_x
    print 'rms_y', rms_y
    
def plotUSNO_w_asteroid(fil, index, chosenpoints_arr):
    """ This function just plots the USNO catalog of stars, given a fits 
    file WITH the asteroid"""
    radeg, dsgn, dedeg, fovam, epoch = USNOhelper(fil)
    name, rad, ded, rmag = usno(radeg,dedeg,fovam,epoch)
    alpha, delta = conv_coords(fil, index, chosenpoints_arr)
    print 'HFASLFHA', alpha, delta
    hdr = loadfits(fil)[1]
    w = np.where(rmag <18.)[0]
    fig = plt.figure()
    plt.plot(rad[w],ded[w],'g.')
    plt.plot(alpha,delta,'ko')
    plt.locator_params(axis='x',nbins=4); plt.locator_params(axis='y',nbins=4)
    plt.tick_params('x',pad=10)
    plt.xlabel('RA [Deg]', fontsize = 18); plt.ylabel('Dec [Deg]',fontsize = 18)
    plt.ticklabel_format(useOffset=False)
    plt.axis('scaled')
    plt.gca().invert_xaxis()
    return fig
#plt.show(plotUSNO_w_asteroid(DanaeI_file_arr1[5] 13, chosenpoints_arrDI1))
