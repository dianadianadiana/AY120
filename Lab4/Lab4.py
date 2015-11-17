#==================================
#Import Packages
#==================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as col
import astropy.io.fits as pf
import math

#==================================
# Global Variables
#==================================
datadir='/Users/Diana/Desktop/Astro 120/Lab4/Data/'
objects = ['baes_halogen_-0', 'baes_neon_m4_-0', 'baes_neon_2_-0', 'baes_redlaser_-0'] 
#halogen 01 - 20, neon_m4 01 -20, neon_2 01 - 20, redlaser 01 - 20
fig_path = '/Users/Diana/Desktop/Astro 120/Lab4/Figures/'
#================================== 
x = ["%.2d" % i for i in range(1,21)] # ['01','02',...,'19','20']
halogen_files = [datadir+objects[0]+i+'.fit' for i in x]
neon_files = [datadir+objects[1]+i+'.fit' for i in x]
neon_2_files = [datadir+objects[2]+i+'.fit' for i in x]
redlaser_files = [datadir+objects[3]+i+'.fit' for i in x]
#==================================
def avg_fits_img(files_arr): # taken from lab 3
    """ takes all the fits in the 'files_arr' and makes them into one average fits array (img) 
        Note: no corrections are applied, just the raw fits files are used """
    #apparently this does the same thing, but I wanted to write it out myself
    #np.sum(bias_data_arr, axis = 0)/len(bias_data_arr)
    data_arr = [pf.getdata(fil) for fil in files_arr]
    avg_arr = data_arr[0]
    for arr in data_arr[1:]:
        avg_arr += arr
    return avg_arr/(len(files_arr))
#==================================   
img_halogen = avg_fits_img(halogen_files)
img_neon = avg_fits_img(neon_files)
img_neon_2 = avg_fits_img(neon_2_files)
img_redlaser = avg_fits_img(redlaser_files)  
#==================================   
def normalize_img(img):
    return img / np.median(img)
    
    
def loadfits(fil):
    """ Input: 
            fil: the filename, with the syntax of 'd100.fits'
        Output: 
            img: the matrix of the fits file [1024,1024]
            hdr: the header of the fits file """
    hdulist = pf.open(fil)
    img = hdulist[0].data # the matrix data
    hdr= hdulist[0].header # the header
    print "Loading file for", hdr['object']
    return [img, hdr]

def display_fits(img, bool_centroids = False):
    """ Display fits with centroids (if bool_centroids is true)"""    
    #img, hdr = loadfits(fil)
    fig = plt.figure()
    img = normalize_img(img)
    #plt.title(hdr['object'], fontsize=18)
    plt.imshow(img, origin='lower', interpolation='nearest', cmap='gray_r', 
               vmin=0.95, vmax=1.2)  
    plt.axhline(y=500)
    if bool_centroids:
        centroids_x, centroids_y = centroid(fil) # centroids_x ~ cols, centroids_y ~ rows
        #plt.scatter(centroids_y, centroids_x, marker='x',s=100) #x's
        plt.scatter(centroids_y, centroids_x, s=150, facecolors='none', edgecolors='b') #circles
    plt.colorbar()
    return fig  
    
plt.show(display_fits(img_neon))

def find_bands(img):
    ''' 
    Purpose: This is so we can automatically determine where the bands are 
            (using the halogen image)
    Returns: an arr of pairs, where the first and second element of a pair
        are the starting and ending rows of a band
        example: [[268, 271], [298, 308], [334, 346], [373, 384], [414, 425],
        [457, 467], [503, 511], [550, 558], [601, 608], [652, 660], [707, 717],
        [765, 776], [826, 840], [891, 908], [960, 980]]

        ...img[268] and img[271] are the beginning and ending row of the band
        Careful: at lower indexes, the bands become really faint, so this may not
        work for really high orders of m
        Note: m = 1 corresponds to the last elem of the returned array, arr[-1]
              m = 2 -- arr[-2] and so on.
    '''
    avg = np.median(img)
    lim = 1.5 * avg #1.5 was arbitrary, worked pretty well 
    # produces an array where the each index corresponds to row, and the value
    # at each index is the number of pixels that are above the limit
    count_arr = [len(np.where(img[row]>lim)[0]) for row in range(len(img))]
    print count_arr
    row_arr = []
    for row, elem in enumerate(count_arr):
        if elem > 40: # if over 50 pixels reach the limit, then we consider them
            row_arr.append(row)
    print row_arr
    # below, the rows are grouped
    k = 0
    temp_arr = []
    group_arr = []
    while k < len(row_arr):
        curr = row_arr[k]
        if temp_arr == [] or (temp_arr != [] and np.abs(curr - temp_arr[0]) <= 20):
            temp_arr.append(curr)
            k += 1
        elif temp_arr != [] and np.abs(curr - temp_arr[0]) > 20:
            group_arr.append(temp_arr)
            temp_arr = []
        if k == len(row_arr) - 1: #idk how to simplify this out of the code
            group_arr.append(temp_arr)
    print group_arr
    # below, the group_arr is filtered so that we only have the beginning and ending
    # row of each group
    arr = []
    for elem in group_arr:
        arr.append([elem[0], elem[-1]])
    print arr
    return arr

def display_fits(img, m_arr = [], bool_centroids = False):
    """ Display the fits image with the band lines (if m_arr) 
        IGNORE the bool_centroids for now ~~ from lab 3"""    
    #img, hdr = loadfits(fil)
    fig = plt.figure()
    img = normalize_img(img)
    #plt.title(hdr['object'], fontsize=18)
    plt.imshow(img, origin='lower', interpolation='nearest', cmap='gray_r', 
               vmin=0.95, vmax=1.2)  
    if m_arr:
        for elem in m_arr:
            plt.axhline(y = .5*(elem[0] + elem[1])) #plot the average
    if bool_centroids: ### IGNORE (from lab3)
        centroids_x, centroids_y = centroid(fil) # centroids_x ~ cols, centroids_y ~ rows
        #plt.scatter(centroids_y, centroids_x, marker='x',s=100) #x's
        plt.scatter(centroids_y, centroids_x, s=150, facecolors='none', edgecolors='b') #circles
    plt.colorbar()
    return fig  
    
m_arr = find_bands(img_halogen)
plt.show(display_fits(img_neon,m_arr))

def avg_img(img, lower_bound, upper_bound): # taken and modified 'vg_fits_img' from lab 3
    """ This takes the img and the lower_bound and upper_bound for the rows
    It adds each row from lower_bound to upper_bound and takes the average
    Return: one array with the averaged values, 'band array' """
    #apparently this does the same thing, but I wanted to write it out myself
    #np.sum(bias_data_arr, axis = 0)/len(bias_data_arr)
    row_arr = [img[row, :] for row in np.arange(lower_bound, upper_bound + 1)]
    avg_arr = row_arr[0]
    for arr in row_arr[1:]:
        avg_arr += arr
    return avg_arr/(upper_bound - lower_bound + 1)

def display_band(img, order_num):
    """ Displays the band at the order_num and returns the figure and the band itself """
    lower_bound, upper_bound = m_arr[-order_num][0], m_arr[-order_num][1]
    current_band = avg_img(img, lower_bound, upper_bound)
    fig = plt.figure()
    plt.title('m = ' + str(order_num))
    plt.plot(current_band)
    return [fig, current_band]
    
img = img_neon
order_num = 4
plt.show(display_band(img, order_num)[0])
#img_curr = img[lower_bound:upper_bound, :]