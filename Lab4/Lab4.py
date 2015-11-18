#==================================
#Import Packages
#==================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as col
import astropy.io.fits as pf
import math
import Lab4_centroids as l4cen # import the centroid function
import LLS as LLS

#==================================
# Global Variables
#==================================
datadir='/Users/Diana/Desktop/Astro 120/Lab4/Data/'
objects = ['baes_halogen_-', 'baes_neon_m4_-', 'baes_neon_2_-', 'baes_redlaser_-'] 
sun_objects = ['sun-1-', 'sun-2-', 'sun-scan-1-', 'sun-scan-2-', 'sun-scan-3-', 
                'sun-3-', 'sun-4-']
#halogen 001 - 020, neon_m4 001 - 020, neon_2 001 - 020, redlaser 001 - 020
# sun_1 001 - 087, sun_2 001 - 100, sun_scan_1 001 - 073, sun_scan_1 001 - 076
# sun_scan_3 001 - 100; 11/16 sun_scan_1 001 - 069, 11/16 sun_scan_2 001 - 063,
# 11/16 sun_scan_3 001 - 076; 
fig_path = '/Users/Diana/Desktop/Astro 120/Lab4/Figures/'
#================================== 
x = ["%.3d" % i for i in range(1,21)] 
def get_num(lower, upper):
    # ['001','002',...,'019','020']
    return ["%.3d" % i for i in range(lower,upper)]
x = get_num(1,21)
halogen_files = [datadir+objects[0]+i+'.fit' for i in x]
neon_files = [datadir+objects[1]+i+'.fit' for i in x]
neon_2_files = [datadir+objects[2]+i+'.fit' for i in x]
redlaser_files = [datadir+objects[3]+i+'.fit' for i in x]

sun_1_files = [datadir+'sun_data/11-11/'+sun_objects[0]+i+'.fit' for i in get_num(1,88)]
sun_2_files = [datadir+'sun_data/11-11/'+sun_objects[1]+i+'.fit' for i in get_num(1,101)]
'''
sun_scan_1_files = [datadir+'sun_data/11-12/'+sun_objects[2]+i+'.fit' for i in get_num(1,74)]
sun_scan_2_files = [datadir+'sun_data/11-12/'+sun_objects[3]+i+'.fit' for i in get_num(1,77)]
sun_scan_3_files = [datadir+'sun_data/11-12/'+sun_objects[4]+i+'.fit' for i in get_num(1,101)]

sun_scan_1_files1 = [datadir+'sun_data/11-16/'+sun_objects[2]+i+'.fit' for i in get_num(1,70)]
sun_scan_2_files1 = [datadir+'sun_data/11-16/'+sun_objects[3]+i+'.fit' for i in get_num(1,64)]
sun_scan_3_files1 = [datadir+'sun_data/11-16/'+sun_objects[4]+i+'.fit' for i in get_num(1,77)]

sun_1_files1 = [datadir+'sun_data/11-17/'+sun_objects[0]+i+'.fit' for i in get_num(1,85)]
sun_3_files1 = [datadir+'sun_data/11-17/'+sun_objects[5]+i+'.fit' for i in get_num(1,101)]
sun_4_files1 = [datadir+'sun_data/11-17/'+sun_objects[6]+i+'.fit' for i in get_num(1,101)]
'''

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

img_sun_1 = avg_fits_img(sun_1_files) # 11/11
img_sun_2 = avg_fits_img(sun_2_files)
'''
img_sun_scan_1 = avg_fits_img(sun_scan_1_files) #11/12
img_sun_scan_2 = avg_fits_img(sun_scan_2_files)
img_sun_scan_3 = avg_fits_img(sun_scan_3_files)
img_sun_scan_1_1 = avg_fits_img(sun_scan_1_files1) #11/16
img_sun_scan_2_1 = avg_fits_img(sun_scan_2_files1)
img_sun_scan_3_1 = avg_fits_img(sun_scan_3_files1)
img_sun_1_1 = avg_fits_img(sun_1_files1) # 11/17
img_sun_3_1 = avg_fits_img(sun_3_files1)
img_sun_4_1 = avg_fits_img(sun_4_files1)
'''
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

m_arr = [[268, 271], [298, 308], [334, 346], [373, 384], [414, 425],
        [457, 467], [503, 511], [550, 558], [601, 608], [652, 660], 
        [707, 717], [765, 776], [826, 840], [891, 908], [960, 980]] #1.5lim, >40

#plt.show(display_fits(img_neon, m_arr = []))

def avg_img(img, lower_bound, upper_bound): # taken and modified 'avg_fits_img' from lab 3
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
#plt.show(display_band(img, order_num)[0])
#for order_num in range(1,9):
#    plt.show(display_band(img_halogen, order_num)[0])

#img_curr = img[lower_bound:upper_bound, :]

def display_all_bands(img):
    fig = plt.figure()
    plt.title('All bands, inverted xaxis')
    for m in range(1,9):
        lower_bound, upper_bound = m_arr[-m][0], m_arr[-m][1]
        current_band = avg_img(img, lower_bound, upper_bound)
        plt.plot(current_band, label = str(m))
    plt.legend()
    plt.gca().invert_xaxis()
    return fig

#plt.show(display_all_bands(img_neon))
#plt.show(display_all_bands(img_neon_2))

def display_fits(img, m_arr = [], bool_centroids = False):
    """ Display the fits image with the band lines (if m_arr) """    
    fig = plt.figure()
    img = normalize_img(img)
    #plt.title(hdr['object'], fontsize=18)
    plt.imshow(img, origin='lower', interpolation='nearest', cmap='gray', 
               vmin=0.8, vmax=1.8)  
    plt.gca().invert_yaxis()
    if m_arr:
        for elem in m_arr:
            plt.axhline(y = .5*(elem[0] + elem[1])) #plot the average
    if bool_centroids: ### IGNORE (from lab3)
        centroids_x, centroids_y = l4cen.centroid(img) # centroids_x ~ cols, centroids_y ~ rows
        #plt.scatter(centroids_y, centroids_x, marker='x',s=100) #x's
        plt.scatter(centroids_y, centroids_x, s=150, facecolors='none', edgecolors='b') #circles
        return [fig, centroids_y, centroids_x]
    #plt.colorbar()
     

    return fig

fig_centroids, x_cen_vals, y_cen_vals = display_fits(img_neon, [], True)  
print x_cen_vals
print y_cen_vals
centroid_pair_arr = np.transpose([x_cen_vals, y_cen_vals])
print centroid_pair_arr
print m_arr   
plt.show(fig_centroids)
plt.show(display_fits(img_redlaser))
'''
 tel = {'jack': 4098, 'sape': 4139}
>>> tel['guido'] = 4127
>>> tel
{'sape': 4139, 'guido': 4127, 'jack': 4098}
>>> tel['jack']
4098
>>> del tel['sape']
>>> tel['irv'] = 4127
>>> tel
{'guido': 4127, 'irv': 4127, 'jack': 4098}
>>> tel.keys()
['guido', 'irv', 'jack']
>>> 'guido' in tel
'''

k = 1
pair_dict = {} # a dict to hold all the pixel coordinates of the centroids at each order
while k <= len(m_arr):
    lower_bound, upper_bound = m_arr[-k][0], m_arr[-k][1]
    temp_arr = []
    for pair in centroid_pair_arr:
        if pair[1] >= lower_bound and pair[1] <= upper_bound:
            temp_arr.append(pair)
        elif pair[1] > upper_bound:
            break # to save time
    pair_dict[k] = temp_arr
    k+=1
print pair_dict
#==================================
# Fitting
#==================================
# focus on the 5th order
order5_pairs = pair_dict[5]
print order5_pairs
neonspec=  np.array([585.249,588.189,594.483,597.553,603.0,607.434,609.616,
614.306,616.359,621.72,626.649,630.479,633.443,638.299,640.225,650.653,
653.288,659.895,667.828,671.704])
neonspec5 = np.array([621.72,626.649,630.479,633.443,638.299,640.225,])
neonspec4 = np.array([638.299,640.225,650.653,653.288,659.895])
x_centroids5, y_centroids5 = np.transpose(order5_pairs)

order4_pairs = pair_dict[4]
x_centroids4, y_centroids4 = np.transpose(order4_pairs)
#x_centroids4 = np.append(x_centroids4, 356.314) # green laser
#neonspec4 = np.append(neonspec4, 650)
x_centroids4 = [36.623323379750964, 352.71187030023879, 475.51865542624859, 
943.4625849415421, 1027.5432586885437]
#x_centroids4 = np.append(x_centroids4, 356.314) # green laser


#==================================
centroids = x_centroids5
spec = neonspec5
centroids, spec = sorted(centroids), sorted(spec)
print centroids
print spec
fit = [centroids, spec]
linear = LLS.linleastsquares(fit, 2)
quad = LLS.linleastsquares(fit, 3)
#cubic = LLS.linleastsquares(fit, 4)

ideal_linear = linear[1]*centroids+linear[0]
ideal_quad = quad[2]*centroids*centroids + quad[1]*centroids + quad[0]
#ideal_cubic = cubic[3]*centroids*centroids*centroids + cubic[2]*centroids*centroids + cubic[1]*centroids + cubic[0]

def idealfigs(): #taken from Lab1
    fig = plt.figure(figsize = (12,5))
    plt.scatter(centroids,spec,s=70,c='r',label='experimental')
    
    plt.plot(centroids, ideal_linear,'--',label='linear best fit')
    plt.plot(centroids, ideal_quad,'k--',label='quadratic best fit')
    #plt.plot(centroids, ideal_cubic, 'o--', label = 'cubic best fit')
    plt.title('Least Squares Fitting', fontsize =20);    
    plt.xlabel('Pixels',fontsize=20);    plt.ylabel('Wavelength [nm]',fontsize=20)
    #for i in range(len(centroids_neon)):
    #    y, x0,x1 = ideal_quad[i], centroids_neon[i]-errors_neon[i]/2., centroids_neon[i]+errors_neon[i]/2.
    #    plt.plot((x0, x1), (y, y),c='k',ls ='-', linewidth =50)
    plt.xlim([np.amin(centroids), np.amax(centroids)])
    plt.legend(numpoints=1, loc = 4, fontsize = 18)
    plt.tick_params(labelsize=16)
    plt.tight_layout()
    return fig
    
plt.show(idealfigs())

def risiduals(): # taken from Lab2
    fig = plt.figure(figsize = (12,5))
    lin_r = ideal_linear - spec 
    quad_r = ideal_quad -  spec
    #cubic_r = ideal_cubic - spec
    plt.scatter(centroids, lin_r,90,'k','^',label = 'linear residuals')
    plt.scatter(centroids, quad_r,70,'r', 'o',label = 'quadratic residuals')
    #plt.scatter(centroids, cubic_r,70,'b', '2',label = 'cubic residuals')
    plt.title('Residuals', fontsize = 18)
    plt.axhline(y = 0, ls='--',c='k')
    plt.xlim([0,1024])
    plt.xlabel('Pixels',fontsize=16);    plt.ylabel('Difference of Fit and Actual Value',fontsize=16)
    plt.legend(numpoints=1, loc = 9, fontsize = 14)
    plt.tick_params(labelsize=14)
    return fig
    
plt.show(risiduals()) # off by like .5/1, want it to be like .1

#wavelength_arr = quad[2]*pixel_arr*pixel_arr + quad[1]*pixel_arr + quad[0]

