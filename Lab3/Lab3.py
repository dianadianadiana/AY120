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

#==================================
# making arrays for the diff types of objects
# each element of the array has the syntax of datapath + 'd100.fits'
# DOES INCLUDE the actual file path
#==================================
bias_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(100,115)]
flatB_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(115,118)]
flatV_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(118,121)]
flatR_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(121,124)]
flatI_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(124,127)]
WratR_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(131,142)]
WratI_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(142,152)]
DanaeR_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(152,163)]
DanaeI_file_arr = [datapath + 'd'+str(i) + '.fits' for i in np.arange(163,172)]
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
#     the flats, so what we do instead is noted in the next step.
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
    print "***", np.median(img), '***'
    print hdr['object']
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
        return np.mean(y-x)/(y-x)
    data_arr = [dead_pixels(pf.getdata(fil))[:,:1024] for fil in flat_files] 
    gain_arr = gain(bias, data_arr[0])
    for flat in data_arr[1:]:
        gain_arr += gain(bias, flat)
    return gain_arr/(len(flat_files))
def get_gain(hdr):
    """ returns the average gain based on the fits' filter """
    filt = hdr['filtnam']
    if filt == 'B':
        return find_gain(flatB_file_arr, bias)
    elif filt == 'V':
        return find_gain(flatV_file_arr, bias)
    elif filt == 'R':
        return find_gain(flatR_file_arr, bias)
    elif filt == 'I':
        return find_gain(flatI_file_arr, bias)
    
#==================================
#==================================
# Bias 
#==================================
bias = avg_fits_img(bias_file_arr) # take all 14(?) bias files to create averaged one
bias = dead_pixels(bias)[:,:1024]  # correct for dead pixels and then make it [1024,1024]
#==================================

#==================================
# Flats
#==================================
# load the files normally (dead_pixels, and [:,:1024] are applied)
flatB_img, flatB_hdr = loadfits(flatB_file_arr[0])
flatV_img, flatV_hdr = loadfits(flatV_file_arr[0])
flatR_img, flatR_hdr = loadfits(flatR_file_arr[0])
flatI_img, flatI_hdr = loadfits(flatI_file_arr[0])
# average the fits -- does not correct them in anyway [1024,1056]; NO correction for dead pixels
flatB_img_avg = avg_fits_img(flatB_file_arr)
flatV_img_avg = avg_fits_img(flatV_file_arr)
flatR_img_avg = avg_fits_img(flatR_file_arr)
flatI_img_avg = avg_fits_img(flatI_file_arr)
# have to correct the avg fits [1024,1024]; YES correction for dead pixels
flatB_img_avg = dead_pixels(flatB_img_avg)[:,:1024]
flatV_img_avg = dead_pixels(flatV_img_avg)[:,:1024]
flatR_img_avg = dead_pixels(flatR_img_avg)[:,:1024]
flatI_img_avg = dead_pixels(flatI_img_avg)[:,:1024]
# gain 
gainB = find_gain(flatB_file_arr, bias)
gainV = find_gain(flatV_file_arr, bias)
gainR = find_gain(flatR_file_arr, bias)
gainI = find_gain(flatI_file_arr, bias)
#==================================

#==================================
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
#==================================
#==================================    
def correct_fits(fil, bias,):
    img, hdr = loadfits(fil) # img is corrected for deadpixels and [1024,1024]
    return (img - bias) * get_gain(hdr)
def normalize_img(img):
    return img / np.median(img)
    
def display_fits(fil, bias):
    """Display fits using matplotlib"""    
    img, hdr = loadfits(fil)
    img = correct_fits(fil, bias)
    img = img/np.mean(img)
    fig = plt.figure()
    plt.title(hdr['object'], fontsize=18)
    plt.imshow(img, origin='lower', interpolation='nearest', cmap='gray_r', 
               vmin=0.95*np.median(img), vmax=1.7*np.median(img))
    plt.colorbar()
    return fig
def display_fits2(fil, bias):
    """Display fits using matplotlib"""    
    img, hdr = loadfits(fil)
    img_c = correct_fits(fil, bias)
    img = img/np.mean(img)
    img_c = img_c/np.mean(img_c)
    
    fig = plt.figure()
    ax = fig.add_subplot(121)
    cax = ax.imshow(img,cmap='gray_r',vmin=0.8,vmax=1.2)
    fig.colorbar(cax)
    ax1 = fig.add_subplot(122)
    cax1 = ax1.imshow(img_c,cmap='gray_r',vmin=.8,vmax=1.2)
    fig.colorbar(cax1)
    #plt.show()
    return fig

#plt.show(display_fits(DanaeR_file_arr[5], bias))
#plt.show(display_fits2(DanaeR_file_arr[5], bias))
#display_fits(DanaeR_file_arr[5], flatR_img).savefig(fig_path+'DanaeR_5_corrected.png',dpi=300)
#==================================
#==================================











###IGNORE

#==================================
# Centroids
#==================================
def max_peaks(arr, width, lower_limit):
    """ Returns an array of the indexes where the indexes indicate where the peaks are"""
    #return [i for i in range(1, len(arr)-1) if arr[i-1]<arr[i] and arr[i+1]<arr[i] and arr[i]>lower_limit]
    return [i for i in range(36, len(arr)-1) if all(arr[i] > arr[i-width:i]) and all(arr[i]>arr[i+1:i+width]) and arr[i]>lower_limit] # 36 is arbitrary
    #assuming that there is nothing before 36

def peaks(fil, bias):
    info = loadfits(fil)
    img, hdr = info[0], info[1]
    img = correct_fits(fil, bias)
    avg = np.median(img)
    lim = 2.2*avg
    #[x][y]
    img_arr = [[]] #ignore the 0th row as not having anything
    # find the initial pixel values that are above the limit (=constant*avg)
    # img_arr has the following syntax:
    # the index of img_arr corresponds to the x value of the pixel
    # img_arr[782] represents the 782th row of the array
    # when img_arr[782] = [653, 897], that means that pixels [782,653] and [782,897] are above the limit
    for i_row in range(1, len(img)): #range starts at 1 bc of the assumption that img[0] has nothing
        img_arr.append(max_peaks(img[i_row], 10, lim))

    

    def cluster(img_arr):
        cluster_arr=[]
        k=0
        while k <len(img_arr)-1:
            curr = img_arr[k]
            if len(curr) >= 1:
                temp_arr = [[k,elem] for elem in curr]
                j=k+1
                while len(img_arr[j]) >= 1:
                    for elem in img_arr[j]:
                        temp_arr.append([j,elem])
                    j+=1
                k=j-1
                cluster_arr.append(temp_arr)
            k+=1
        return cluster_arr
    # clusters is an array where it takes in img_arr and sees something like this
    # 289 []; 290 [654, 786]; 291 [656]; 292 [] 
    # and groups it like so [ [290, 654], [290,786], [291,656] ] -- this is one cluster
    # clusters is an array of multiple 'cluster' like above
    clusters = cluster(img_arr)
    
    def go_thru_clusters(clusters):
        new_cluster_arr =[]
        for cluster in clusters:
            #print '*** orig cluster', cluster
            def group(cluster):
                temp_arr = [cluster[0]]
                cluster.remove(cluster[0])
                curr_i = 0
                while cluster:
                    if np.abs(cluster[curr_i][1] - temp_arr[-1][1])<15 and np.abs(cluster[curr_i][0] - temp_arr[-1][0]) <=5:
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
    # that were picked up
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
    
    for winner in peaks:
        print winner
    return peaks




        
    
        
  
    
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

def centroid(fil, peaks):
    """ peaks has the form of [[x,y],size] 
    deciding to take into account +- size for centroiding"""
    centroids =[]
    img, hdr = loadfits(fil)
    img = correct_fits(fil, bias)
    for pair, size in peaks:
        x = range(pair[0]-size, pair[0]+size+1) # just the x values
        y =range(pair[1]-size, pair[1]+size+1) # just the y values
        x_range = [[x, pair[1]] for x in x] # y is constant val, x is changing
	y_range = [[pair[0], y] for y in y] # x is constant val, y is changing
	x_range_img = [img[i[0],i[1]] for i in x_range] # vals of img at constant y, varying x
	y_range_img = [img[i[0],i[1]] for i in y_range] # vals of img at constant x, varying y
	print pair
	print x
	print y

	x_centroid =  np.sum(x*x_range_img)/np.sum(x_range_img) #<x>
	y_centroid =  np.sum(y*y_range_img)/np.sum(y_range_img) #<y>
	
	centroids.append([x_centroid, y_centroid])
    centroids = np.transpose(centroids)
    centroids_x, centroids_y = centroids[0], centroids[1]
    return centroids_x, centroids_y
    
    
peaks = peaks(DanaeR_file_arr[5], bias)
print peaks

peaks1 = np.transpose(peaks)
print 'peaks1', peaks1
peaks2 = np.transpose(peaks1[0])
print 'peaks2', peaks2
print peaks2[0]
peak_x,peak_y=[],[]
for peak in peaks2:
    peak_x.append(peak[0])
    peak_y.append(peak[1])
print peak_x
print peak_y

centroids_x, centroids_y = centroid(DanaeR_file_arr[5],peaks)
#print centroids
#plt.scatter(centroids_x,centroids_y,marker='x',s=50)
plt.show()
def display_fits_w_centroids(fil, bias, centroids_x, centroids_y):
    """Display fits using matplotlib"""    
    img, hdr = loadfits(fil)
    img = correct_fits(fil, bias)
    #img = img/np.mean(img)
    fig = plt.figure()
    plt.title(hdr['object'], fontsize=18)
    plt.scatter(centroids_x,centroids_y,marker='x',s=50)
    plt.imshow(img, origin='lower', interpolation='nearest', cmap='gray_r', 
               vmin=0.95*np.median(img), vmax=1.7*np.median(img))
    plt.colorbar()
    return fig   

plt.show(display_fits_w_centroids(DanaeR_file_arr[5], bias, centroids_y, centroids_x))
	#plt.imshow(image,cmap='gray_r',vmin=50,vmax=500)
	


    
