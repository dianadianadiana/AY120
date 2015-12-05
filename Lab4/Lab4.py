#==================================
#Import Packages
#==================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import astropy.io.fits as pf
import math
import Lab4_centroids as cen # import the centroid and such functions
import LLS as LLS
import heapq # for finding the max 3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETTING UP ALL THE VARIABLES
# (fig4 in https://drive.google.com/file/d/0B40Ynk22SiBpX2dvRFZzbEdZQmc/view)
# - In 'Global Variables', I set up where all our data is coming from
#   For calibration, we took 20 observations of: halogen lamp, neon, and redlaser
# - Then I 'obtain all the folder arrays' because it will be easier to call a folder
#   later on
# - Then I average out the images for the calibration data, I do it for the 
#   sun as well, but I don't think I should have but I left it there anyway
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
# Obtaining all the folder arrays
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

sun_1_files1 = [datadir+'sun_data/11-7/'+'sun01-'+i+'.fit' for i in get_num(1,191)]
'''
sun_scan_1_files = [datadir+'sun_data/11-12/'+sun_objects[2]+i+'.fit' for i in get_num(1,74)]
sun_scan_2_files = [datadir+'sun_data/11-12/'+sun_objects[3]+i+'.fit' for i in get_num(1,77)]
sun_scan_3_files = [datadir+'sun_data/11-12/'+sun_objects[4]+i+'.fit' for i in get_num(1,101)]

sun_scan_1_files1 = [datadir+'sun_data/11-16/'+sun_objects[2]+i+'.fit' for i in get_num(1,70)]
sun_scan_2_files1 = [datadir+'sun_data/11-16/'+sun_objects[3]+i+'.fit' for i in get_num(1,64)]
sun_scan_3_files1 = [datadir+'sun_data/11-16/'+sun_objects[4]+i+'.fit' for i in get_num(1,77)]

sun_1_files2 = [datadir+'sun_data/11-17/'+sun_objects[0]+i+'.fit' for i in get_num(1,85)]
sun_3_files2 = [datadir+'sun_data/11-17/'+sun_objects[5]+i+'.fit' for i in get_num(1,101)]
sun_4_files2 = [datadir+'sun_data/11-17/'+sun_objects[6]+i+'.fit' for i in get_num(1,101)]
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
# averaging the images 
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
img_sun_1_2 = avg_fits_img(sun_1_files1) # 11/17
img_sun_3_2 = avg_fits_img(sun_3_files1)
img_sun_4_2 = avg_fits_img(sun_4_files1)
'''
#==================================
def normalize(img):
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#==================================
#================================== 
# Calibration
#==================================
#==================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#(1) Finding the Bands using Halogen Lamp; 
#    - I wrote the 'find_bands' program to automatically find 
#    all the echelle bands: The way it works is that I input the Halogen averaged 
#    image and then found all the rows where there were over 40 (maybe I changed 
#    the number later) pixels that were over the limit (like 1.5*avg). Neighboring 
#    rows were put into the same band. This already showed us how the lower orders
#    (like m=4/m=34) have a much thicker band (~20 pixels) whereas higher orders
#    have a much thinner band (~10 pixels). Additionally, higher orders became much
#    closer to eachother than the lower orders
#    - Note: once I ran this function, I decided to copy n paste the returned array
#    so I don't need to keep on calling the function (aka wasting time)
#    - The thicker, spread apart bands (aka the lower order bands) are on the bottom
#    of the image and they are the longer wavelengths (redder) and the thinner,
#    closer together bands (aka the higher order bands) are on the top and they
#    are the shorter wavelengths (bluer) -- we figure this out after matching the 
#    pixels to corresponding wavelengths
#(2) Obtaining the flux of each echelle order using the bands from part 1
#    - the avg_img fn takes in an img and the lower and upper bound, sums up the 
#    flux of each row in that bound and then takes the average. the output is an
#    array that when plotted, it is a spectrum (note: this is a helper fn)
#    - the get_band_flux takes the above fn (avg_img) to calculate the band flux
#    - then we just plot the pixel vs intensity of the band
#    - the display_all_bands shows all the bands on top of eachother, which in 
#    reality should not be as interpreted as true but rather, it gives us a better 
#    idea of how the pixels (horizontally) correspond with wavelength -- I found
#    that longer wavelengths (redder) are on the left, and shorter wavelengths
#    (bluer) are on the right 
#        - this helps us when figuring out which wavelength (longer or shorter) do
#        we assign to each centroid
#        - the reason this method of displaying all the bands does not work, is
#        because of two reasons: (i) the end (right) of one band is the beginning
#        (left) of the band with an order higher; (ii) sounds great that they 
#        overlap, but they don't do see desirably, the distance between centroids
#        in the lower order one is less than the distance in the higher order one
#        (because of diffraction)
#    - I then included a fn that displays the actual img (aka the 2D array)
#(3a) Finding the centroids of each band and then find the max peaks and centroids 
#   - deciding whether to use the 2D centroid alogrithm from lab3 or the 1D 
#   centroid algorithm from lab2 (decided to do 1D), why -> because it was easy 
#   to obtain a 1D spectra of each order and calculate everything from there
#   - made variables for the neon spec (3 to 8) and then made variables for the 
#   centroids (3 to 7) and also found that the error in centroids was around 10e-5
#   - created a dictionary where you could just call a band order (like 2) and it 
#   would give you the 2 element array with the corresponding neonspec and centroids
#(3b) Figuring out which pixels correspond to which wavelengths using the red laser
#(had to match them up) did this for orders 3-8 (I think)
#    - we know that the red laser is at 653nm so we figure out that it is in the 
#    4th (34th order) and from there we can deduce the other corresponding 
#    wavelengths for the other pixels
#(4) Fitting the centroids to the known wavelengths using the LLS method from lab2, 
#and finding that the best fit is the quadratic *surprise surprise* And getting the
#residuals for each band
#    - we fit the centroids to the known wavelengths (similar to lab2) yeah
#(5) Finally applying the fit to the orders
#    - now i just have two functions that apply a quad fit and get the wavelength
#    calibration for a given order, respectfully. I created a wavelength dict so 
#    I can just get the calibrated wavelength when I input an order_num (3 to 7)
#    - we only need to take into account one order calibration, so i chose m=34 
#    because it had a strong Halpha line

#==================================
def find_bands(img):
    '''  Purpose: This is so we can automatically determine where the bands are 
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
              m = 2 -- arr[-2] and so on.'''
    avg = np.median(img)
    lim = 1.5 * avg #1.5 was arbitrary, worked pretty well 
    # produces an array where the each index corresponds to row, and the value
    # at each index is the number of pixels that are above the limit
    count_arr = [len(np.where(img[row]>lim)[0]) for row in range(len(img))]
    print count_arr
    row_arr = []
    for row, elem in enumerate(count_arr):
        if elem > 40: # if over 40 pixels reach the limit, then we consider them
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
    
#m_arr = find_bands(img_halogen)
#decided to copy n paste the m_arr just to save time
m_arr = [[268, 271], [298, 308], [334, 346], [373, 384], [414, 425],
        [457, 467], [503, 511], [550, 558], [601, 608], [652, 660], 
        [707, 717], [765, 776], [826, 840], [891, 908], [960, 980]] #1.5lim, >40
#==================================
def avg_img(img, lower_bound, upper_bound): # taken and modified 'avg_fits_img' from lab 3
    """ This takes the img and the lower_bound and upper_bound for the rows that we
            want to average
        It adds each row from lower_bound to upper_bound and takes the average
        Return: one array with the averaged values, 'band array', such that when you plot
            it, it is a spectrum (where x axis is in pixels)"""
    #apparently this does the same thing, but I wanted to write it out myself
    #np.sum(bias_data_arr, axis = 0)/len(bias_data_arr)
    row_arr = [img[row, :] for row in np.arange(lower_bound, upper_bound + 1)]
    avg_arr = row_arr[0]
    for arr in row_arr[1:]:
        avg_arr += arr
    return avg_arr/(upper_bound - lower_bound + 1)
def get_band_flux(img, order_num):
    lower_bound, upper_bound = m_arr[-order_num][0], m_arr[-order_num][1]
    current_band = avg_img(img, lower_bound, upper_bound)
    return current_band
def display_band(img, order_num, centroids = []):
    """ Displays the band at the order_num and returns the figure and the band itself """
    current_band = get_band_flux(img, order_num)
    fig = plt.figure()
    plt.title('Neon, m = ' + str(order_num + 30), fontsize =20) # 1 corresponds to 31
    plt.plot(current_band)
    plt.xlabel('Pixel [count]', fontsize = 18); plt.ylabel('Intensity [ADU]', fontsize =18)
    plt.xlim([0,1048])
    if len(centroids) != 0:
        for centroid in centroids:
            plt.axvline(x=centroid, c='r',ls ='--', linewidth =.75)#label = str(centroid))
    return fig
  
img = img_neon
order_num = 4

#plt.show(display_band(img, order_num))
#for order_num in range(1,9): #display all the bands
#    plt.show(display_band(img_halogen, order_num))

def display_all_bands(img):
    fig = plt.figure()
    plt.title('All bands, inverted xaxis')
    for m in range(1,9):
        lower_bound, upper_bound = m_arr[-m][0], m_arr[-m][1]
        current_band = avg_img(img, lower_bound, upper_bound)
        plt.plot(current_band, label = str(m + 30))
    plt.legend()
    plt.gca().invert_xaxis()
    return fig

#plt.show(display_all_bands(img_neon))
#plt.show(display_all_bands(img_neon_2))

def display_fits(img, m_arr = [], bool_centroids = False, band = []):
    """ Display the fits image with the band lines (if m_arr) """    
    fig = plt.figure()
    img = normalize(img)
    #plt.title(hdr['object'], fontsize=18)
    plt.imshow(img, origin='lower', interpolation='nearest', cmap='gray', 
               vmin=0.8, vmax=2)
    plt.xlabel('Pixel', fontsize = 18); plt.ylabel('Pixel', fontsize = 18) 
    plt.gca().invert_yaxis()
    if band:
        print '%%%%%%%%%%%%%%', band
        plt.axhline(y = band[0])
        plt.axhline(y = band[1])
    if m_arr:
        for elem in m_arr:
            plt.axhline(y = .5*(elem[0] + elem[1])) #plot the average
    if bool_centroids: ### IGNORE (from lab3)
        centroids_x, centroids_y = cen.centroid_2D(img) # centroids_x ~ cols, centroids_y ~ rows
        #plt.scatter(centroids_y, centroids_x, marker='x',s=100) #x's
        plt.scatter(centroids_y, centroids_x, s=150, facecolors='none', edgecolors='b') #circles
        return [fig, centroids_y, centroids_x]
    #plt.colorbar()
    return fig

fig_centroids, x_cen_vals, y_cen_vals = display_fits(img_neon, [], True)  
fig_wo_centroids = display_fits(img_neon, [], False)  
fig_w_bands =  display_fits(img_neon, m_arr, False)  
fig_w_band4 =  display_fits(img_neon, [], False, m_arr[-4])  
#plt.show(fig_w_bands)
centroid_pair_arr = np.transpose([x_cen_vals, y_cen_vals])
#print centroid_pair_arr
#print m_arr   
#plt.show(fig_centroids)
#plt.show(display_fits(img_redlaser + img_neon))
#plt.show(display_band(img_neon, 3))
#==================================

#==================================
#==================================
# Fitting
#==================================
#==================================
#==================================
def idealfigs_residuals(cen_spec): #taken from Lab1
    centroids, spec = cen_spec
    centroids, spec = sorted(centroids), sorted(spec)[::-1] # reverse the spec cause short lambda (left) to long lambda (right)
    print "IDEALFIGS_RESIDUALS"
    print 'centroids', centroids
    print 'neonspec', spec 
    fit = [centroids, spec]
    linear = LLS.linleastsquares(fit, 2)
    quad = LLS.linleastsquares(fit, 3)
    #cubic = LLS.linleastsquares(fit, 4)
    
    ideal_linear = linear[1]*centroids+linear[0]
    ideal_quad = quad[2]*centroids*centroids + quad[1]*centroids + quad[0]
    #ideal_cubic = cubic[3]*centroids*centroids*centroids + cubic[2]*centroids*centroids + cubic[1]*centroids + cubic[0]
   
    idealfig = plt.figure()
    plt.scatter(centroids,spec,s=70,c='r',label='experimental')
    
    #plt.plot(centroids, ideal_linear,'--',label='linear best fit')
    plt.plot(centroids, ideal_quad,'k--',label='quadratic best fit')
    #plt.plot(centroids, ideal_cubic, 'o--', label = 'cubic best fit')
    plt.title('Least Squares Fitting', fontsize =20);    
    plt.xlabel('Pixels',fontsize=20);    plt.ylabel('Wavelength [nm]',fontsize=20)
    #for i in range(len(centroids_neon)):
    #    y, x0,x1 = ideal_quad[i], centroids_neon[i]-errors_neon[i]/2., centroids_neon[i]+errors_neon[i]/2.
    #    plt.plot((x0, x1), (y, y),c='k',ls ='-', linewidth =50)
    plt.xlim([np.amin(centroids), np.amax(centroids)])
    plt.legend(numpoints=1, loc = 'best', fontsize = 18)
    plt.tick_params(labelsize=16)
    plt.tight_layout()
    
    residuals_fig = plt.figure()
    lin_r = ideal_linear - spec 
    quad_r = ideal_quad -  spec
    #cubic_r = ideal_cubic - spec
    plt.scatter(centroids, lin_r,90,'k','^',label = 'linear residuals')
    plt.scatter(centroids, quad_r,70,'r', 'o',label = 'quadratic residuals')
    #plt.scatter(centroids, cubic_r,70,'b', '2',label = 'cubic residuals')
    plt.title('Residuals', fontsize = 18)
    plt.axhline(y = 0, ls='--',c='k')
    plt.xlim([0,1024])
    plt.xlabel('Pixels',fontsize=20);    plt.ylabel('Difference of Fit and Actual Value',fontsize=20)
    plt.legend(numpoints=1, loc = 'best', fontsize = 14)
    plt.tick_params(labelsize=14)
    return idealfig, residuals_fig, quad
#==================================
def get_centroids_1D(img, order_num, width = 10, lower_limit = 300):
    intensity_arr = get_band_flux(img, order_num)
    peaks = cen.max_peaks(intensity_arr, width, lower_limit)
    return cen.centroid_1D(peaks, range(0,1048), intensity_arr, width)
#==================================
# Get all the neon values for each band
neonspec=  np.array([585.249,588.189,594.483,597.553,603.0,607.434,609.616,
614.306,616.359,621.72,626.649,630.479,633.443,638.299,640.225,650.653,
653.288,659.895,667.828,671.704]) 
neonspec7 = np.array([585.249,588.189,594.483,597.553])
neonspec6 = np.array([597.553,603.0,607.434,609.616,614.306,616.359,621.72])
neonspec5 = np.array([621.72,626.649,630.479,633.443,638.299,640.225,])
neonspec4 = np.array([638.29917,640.2248,650.65281,653.28822,659.89529])
neonspec3 = np.array([659.895, 667.828, 671.704])
#==================================
centroids3, centroids3_error = get_centroids_1D(img_neon, 3)
centroids4, centroids4_error = get_centroids_1D(img_neon, 4, width = 6)
centroids5, centroids5_error = get_centroids_1D(img_neon, 5)
centroids6, centroids6_error = get_centroids_1D(img_neon, 6)
centroids7, centroids7_error = get_centroids_1D(img_neon, 7)
centroids8, centroids8_error = get_centroids_1D(img_neon, 8) # don't need cause dont have neonspec8
centroid_redlaser, certoid_redlaser_error = get_centroids_1D(img_redlaser, 4)
#centroids4 = np.append(centroids4, centroid_redlaser) # red laser
#neonspec4 = np.append(neonspec4, 653)

# made a dict to easily call centroids and neonspec based on which order I want
cen_spec_dict = {3:[centroids3,neonspec3],4:[centroids4,neonspec4],
                5:[centroids5,neonspec5], 6:[centroids6,neonspec6],
                7:[centroids7,neonspec7]}
#==================================
print "========= centroids (using 1D) 8,7,6,5,4,3 ========="
print sorted(centroids8)
print sorted(centroids7)
print sorted(centroids6)
print sorted(centroids5)
print sorted(centroids4)
print sorted(centroids3)
print 'red laser', centroid_redlaser
print "=================================================="
#==================================
m = 4
#plt.show(fig_wo_centroids)
#plt.show(display_band(img_neon, m, cen_spec_dict[m][0]))
#ideal_fig, residuals_fig, quad = idealfigs_residuals(cen_spec_dict[m])
#plt.show(ideal_fig)
#plt.show(residuals_fig)
#==================================
# FINISH THE CALLIBRATION 
#==================================
def calibrate_quad(quad, x):
    """An easy way to just get automatically apply the quad fit"""
    return quad[2]*x*x + quad[1]*x + quad[0]
def get_wavelength(order_num):
    """Gets the calibration of pixel to wavelength for the xaxis for a given order"""
    centroids, spec = cen_spec_dict[order_num]
    centroids, spec = sorted(centroids), sorted(spec)[::-1] # reverse the spec cause short lambda (left) to long lambda (right)
    fit = [centroids, spec]
    quad = LLS.linleastsquares(fit, 3)
    wavelength = calibrate_quad(quad, np.arange(1048))
    return wavelength # size = 1048
#set up a dictionary that contains all the wavelength calibrations
cal_wavelengths_dict = {}
for i in range(3, 8): # orders 3 to 7
    cal_wavelengths_dict[i] = get_wavelength(i)
    
def plot_calibration(img, m):
    wavelength = cal_wavelengths_dict[m]
    fig = plt.figure()
    flux_arr = get_band_flux(img, m)
    plt.plot(wavelength, flux_arr)
    plt.title('Sun, m = ' + str(m + 30) + '\nWavelength Calibration', fontsize=20)
    plt.xlabel('Wavelength [nm]', fontsize = 18); plt.ylabel('Intensity [ADU]', fontsize = 18)
    plt.xlim([np.min(wavelength), np.max(wavelength)])
    return fig
#plt.show(plot_calibration(img_sun_1,4))
def plot_calibration_sun(fil, m): #for one file!
    img, hdr = loadfits(fil)
    wavelength = cal_wavelengths[m]

    exptime, date, filename, jd = hdr['EXPOSURE'], hdr['DATE-OBS'], hdr['OBJECT'], hdr['JD']
    date = date[5:7]+'/'+date[8:10]+'/'+date[:4] #format the date to look prettier
    fig = plt.figure()
    plt.title('Sun Spectrum on ' + date + '\nm = ' + str(m), fontsize = 20)
    flux_arr = get_band_flux(img, m)

    plt.plot(wavelength, flux_arr)
    plt.xlabel('Wavelength [nm]', fontsize = 18); plt.ylabel('Intensity [ADU]', fontsize = 18)
    plt.xlim([np.min(wavelength), np.max(wavelength)])
    return fig

#plt.show(plot_calibration_sun(sun_1_files[20],4))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
#==================================
#==================================
#SUNS FLUX OVER TIME
#==================================
#==================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#(1) Sum up the total flux of a file/image & find the time of each file observation
#    - to get the TIME of each observation, we can not just assume that the 
#    exposure time is the delta t between observations; so we must obtain the 
#    julian day [JD] of each file and convert that to seconds and then that is
#    our time array
#    - to get the FLUX of each of each observation, we look at the img of one file,
#    and add up the fluxes of each band in that image (I summed up over the first
#    15 bands). Might've been easier to just sum up the total flux of the whole
#    image but I was afraid of picking up false flux readings, and so I think
#    defining the total flux to be the fluxes from the bands is a safe route for
#    minimal errors in calculating total flux
#(2) Fixing the flux array and normalizing it
#    - So I just had the raw total flux of each file in the flux array at the start
#    - Then I found out which files are part of the transit of the sun; I took the 
#    non-transit fluxes, and took the average. I subtracted this average to account
#    for the scattering light the optical fiber cable picks up when it is not facing
#    the Sun. 
#    - Then I divided the flux array by its max value 
#    - So now the flux array is I/I_max (I/I_0)
#(3) Find delta t, t_0, and I_0 by 
#    - eqn 8 (https://drive.google.com/file/d/0B40Ynk22SiBpX2dvRFZzbEdZQmc/view)
#    - so this is a manual iterative process in finding the delta t (the time 
#    from the center of the transit to the edge), t_0 (the center of the transit
#    time), and I_0 (though I found that I_0 will more or less be around 1 because
#    of how I rescaled the flux array to be I/I_max)
#    - To start, I just found a ballpark/intial guess for delta t by simply 
#    taking the difference in time between the start of the transit and the end
#    and divided it by two. 
#    To find t_0 I simply took the difference in indexes of the transit start and 
#    end and divided it by two and found the time at that index.
#    To find I_0, I took the intensity at where the t_0 location index is.
#    - Now, I put in these parameters to the Intensity equation (aka the fit) and
#    then I plotted it against the actual data, visually I could see how it was 
#    a little off: so I plotted a residuals graph and calculated the rms. I was 
#    changing the parameters for the fit to see what parameters give me the smallest
#    rms error. However much the residuals did help, I felt like it was hard to 
#    find a fit that exactly matched the limb darkening portion of the transit, so
#    that's why I think the residuals/rms error were skewed
#================================== 
#~ Results for 11/11/15 Sun-1
#    spectrum files sun transit = [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
#    34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 
#    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65]
#    t0 loc 44 (location of where t_0 is)
#    Initial guesses:
#        delt 134.999969482
#        t_0 141.499984741
#        I_0 0.994419608416
#        rms:  0.0262863329707
#    After guessing and checking:
#        delt = 136.4
#        t_0 = 142.5
#        I_0 = 0.994419608416
#        rms: 0.0241672990943
#    After guessing and checking:
#        delt = 136.4
#        t_0 = 142.5
#        I_0 = 0.994419608416
#        rms: 0.0237264323039
#==================================

def get_flux_over_time(files):
    """ This takes in a directory (files) and finds the total flux of each file.
    To find the total flux of each file, I summed up the fluxes from each band 
    (in order to ensure no "false" flux readings were considered). 
    Returns: an array with the summed fluxes of each file in the 'files' dir"""
    flux_arr, time_arr = [], []
    hdr_0 = loadfits(files[0])[1]
    startjd = hdr_0['JD']*24*3600 # in seconds
    for k in range(len(files)):
        fil = files[k]
        img, hdr = loadfits(fil)
	jd = hdr['JD']*24*3600 # in seconds
	time_arr.append(jd-startjd)
        temp_flux = 0
        for num in range(len(m_arr)):
            band_flux = get_band_flux(img, num)
            temp_flux += np.sum(band_flux)
        flux_arr.append(temp_flux)
    return [time_arr, flux_arr]

def get_delt(time_arr, flux_arr):
    arr=[] #holds all the indexes where the sun is shining
    for index, value in enumerate(flux_arr):
        if value > (np.max(flux_arr) - np.min(flux_arr))*.1 + np.min(flux_arr):
            arr.append(index)
    delt_loc = arr[-1] - arr[0]
    t_0_loc = delt_loc/2 + arr[0]
    I_0 = flux_arr[t_0_loc]
    
    delt = time_arr[arr[-1]] - time_arr[arr[0]]
    t_0 = delt/2. + time_arr[arr[0]]
    print 't0 loc', t_0_loc
    print 'arr', arr
    print 'delt', delt
    print 't_0', t_0
    print 'I_0', I_0
    return arr, delt, t_0, I_0, t_0_loc

def intensity(t_arr, I_0, t_0, delt, flux_arr, arr):
    I_arr = []
    for t in t_arr:
        I = I_0 * (2./5 + 3./5 * np.sqrt(1 - (t-t_0)**2/delt**2))
        if math.isnan(I):
            I_arr.append(np.median(np.delete(flux_arr,arr)))
        else:
            I_arr.append(I)
    return I_arr
    
def sun_flux_over_time_fig(files):
    """ This funtion just plots the total relative flux over time """
    time_arr, flux_arr = get_flux_over_time(files)
    #flux_arr = flux_arr/ np.median(flux_arr)
    #flux_arr = flux_arr / np.max(flux_arr)
    hdr = loadfits(files[0])[1]
    exptime, date, filename, jd = hdr['EXPOSURE'], hdr['DATE-OBS'], hdr['OBJECT'], hdr['JD']
    date = date[5:7]+'/'+date[8:10]+'/'+date[:4] #format the date to look prettier
    transit_arr, delt, t_0, I_0, t_0_loc = get_delt(time_arr, flux_arr)
    flux_arr -= np.median(np.delete(flux_arr, transit_arr))
    flux_arr = flux_arr / np.max(flux_arr)
    I_0 = flux_arr[t_0_loc] #because I just rescaled the flux_arr
    
    # Figure of just the flux vs time
    delt_str = "{:.6f}".format(delt)
    fig = plt.figure(figsize = (12,5))
    plt.title('Flux versus time for ' + filename + ' on ' + date + '\nwith exposure time: ' 
                + str(exptime) + ' seconds', fontsize = 20)
    plt.xlabel('Time [s]', fontsize = 18) ; plt.ylabel('Total Relative Flux', fontsize = 18)
    plt.plot(time_arr, flux_arr , 'o', label = 'delta t: ' + delt_str + 's')
    plt.axvline(x=t_0, c='r',ls ='--', linewidth =.75, label = 't_0 = ' + str(t_0))
    for point in transit_arr:
         plt.axvline(x=time_arr[point], c='r',ls ='--', linewidth =.5)
    plt.tight_layout()
    plt.legend(loc = 'best')
    plt.xlim([np.min(time_arr), np.max(time_arr)]); 
    plt.ylim([-.1,1.1])
    
    
    t=time_arr
    I_0 = 1
    delt = 136.4
    t_0 = 142.5
    delt_str = "{:.1f}".format(delt)
    t_0_str = "{:.1f}".format(t_0)
    I_0_str = "{:.4f}".format(I_0)

    i_arr = intensity(t, I_0, t_0, delt/2, flux_arr, transit_arr) #set t_0 to be 0
    # Fitted figure
    fit_fig = plt.figure()
    plt.title('Flux versus time for ' + filename + ' on ' + date + '\nwith exposure time: ' 
                + str(exptime) + ' seconds', fontsize = 20)
    #t = np.arange(0 - delt , 0 + delt + 1 , .1) #go from -delt to delt with .1 increments for the fit
    #i_arr = intensity(t, I_0, 0, delt/2, flux_arr, transit_arr) #set t_0 to be 0

    plt.plot(t, i_arr, label = 'Fit\n' +
            r'$\Delta t$' + ': ' + delt_str + r'$\pm$' + '1.4s' + '\n' + 
            r'$t_0$' + ': ' + t_0_str + r'$\pm$' + '.8s' + '\n' + 
            r'$I_0$' + ': ' + I_0_str)
    plt.plot(np.array(time_arr), flux_arr, 'o', label = 'Data') #plot the original data but offset of t_0 to make t_0 be at 0
    jd_t0 = loadfits(files[t_0_loc])[1]['JD']
    plt.xlabel('Time - ' + str(jd_t0*24*3600)+' [JD in s]', fontsize=18); 
    plt.ylabel(r'$I / I_0$', fontsize = 18)
    plt.tight_layout()
    plt.legend(loc = 'best', framealpha=0.5, fancybox=True)
    #plt.xlim([np.min(np.array(time_arr)-t_0), np.max(np.array(time_arr)-t_0)]); 
    plt.ylim([-.1,1.1])
    
    residuals_fig = plt.figure()
    t_new = time_arr
    i_arr_new = intensity(t_new, I_0, t_0, delt/2, flux_arr, transit_arr) #set t_0 to be 0
    res = i_arr_new - flux_arr
    rms = np.sqrt(np.mean(res**2))
    plt.plot(t_new, res, 'o', label = str(rms))
    plt.title('delt: ' +str(delt) + ' t_0: ' + str(t_0) + ' I_0: ' +str(I_0))
    plt.legend()
    print 'RMS', rms
    
    return fig, fit_fig, residuals_fig
    
#plt.show(display_fits(img_sun_1)) #the fits img 
plt.show(sun_flux_over_time_fig(sun_1_files)) #the flux vs time
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#==================================
#==================================
# PLOTTING THE SUN'S DIFFERENCE IN SPECTRA
#==================================
#==================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FOR PLOTTING THE DOPPLER SHIFT
# (fig4 in https://drive.google.com/file/d/0B40Ynk22SiBpX2dvRFZzbEdZQmc/view)
# - I call in 'files' and 'ordernum' as parameters, where the files are the sun files
# - I then find the fluxes of each file in files and put that into flux_arr
#   - transit_arr is the array of indexes that are when the sun is passing over
#   To find out which two files out of files i will be using.. the first one
#       will be around t_0 (the midpoint of the transit) and the second one
#       will be the third from last one (arbitrary - i just didnt want one
#       that is really close to the edge because of limb darkening effects)
# - I then look for the noise in the fluxes by calculating the median of the non-transit
#   points; I subtract this median from the band fluxes for each index later on
# - I load in the img for each index, and then I get the flux of the band I input into 
#   the function and subtract the noise from the nontransit points and then normalize 
#   by dividing by the max of the band flux of each respective index
#   (I had to ignore 10 points because they were just not good)
# - Then I find the time change between the two indexes that I'm looking at by
#   obtaining their JD (julian day) and converting it to seconds
# - Then I need to calibrate the x axis, so I just call the corresponding quad
#   fit for the order_num I input -- and apply the fit and also ignore the annoying 
#   points
# - I just plot everything to finish it off!
#==================================

def get_current_band(img, order_num, noise, ignore_index = 10):
    # this is a helper function
    band = get_band_flux(img, order_num) - noise
    band /= np.max(band)
    return band[:-ignore_index]
def get_wavelength_fit(wavelength, current_band):
    fit = [wavelength, current_band]
    quad = LLS.linleastsquares(fit, 3)
    wavelength_fit = calibrate_quad(quad, wavelength)
    return wavelength_fit
def display_band_sun(files, order_num):
    # to figure out the indicies of the transit
    time_arr, flux_arr = get_flux_over_time(files)
    flux_arr = flux_arr / np.max(flux_arr)
    transit_arr = get_delt(time_arr, flux_arr)[0]
    index1, index2 = transit_arr[len(transit_arr)/2], transit_arr[-3]
    nontransit_noise = np.median(np.delete(flux_arr,transit_arr))
    #========================
    # get the flux arrays for the two different files
    img_1, hdr_1 = loadfits(files[index1])
    img_2, hdr_2 = loadfits(files[index2])
    exptime, filename = hdr_1['EXPOSURE'], hdr_1['OBJECT']
    date = hdr_1['DATE-OBS']
    date = date[5:7]+'/'+date[8:10]+'/'+date[:4] #format the date to look prettier
    ignore_index = 10 # for some reason there are 10 annoying points
    current_band_1 = get_current_band(img_1, order_num, nontransit_noise)
    current_band_2 = get_current_band(img_2, order_num, nontransit_noise)
    #========================
    # to find the time change
    jd_1, jd_2 = hdr_1['JD'], hdr_2['JD']
    time_change = (jd_2-jd_1)*24*3600
    time_change_str = "{:.1f}".format(time_change)
    #========================
    # Calibration for the x axis
    wavelength = cal_wavelengths_dict[order_num]
    wavelength = wavelength[:-ignore_index]
    #========================
    # Plot everything Yay!
    fig = plt.figure()
    plt.title('Change in time between the two events: ' + time_change_str + ' [s]'
                +'\nExposure time: ' + str(exptime) + ' [s]; The Sun (m='+str(order_num+30)+') on ' + date,
                fontsize=20)
    plt.plot(wavelength, current_band_1, label = 'spectrum1: ' + str(index1))
    plt.plot(wavelength, current_band_2, label = 'spectrum2: ' + str(index2))
    plt.xlabel('Wavelength [nm]', fontsize = 18); plt.ylabel('Intensity [ADU]', fontsize =18)
    plt.xlim([np.min(wavelength),np.max(wavelength)])
    plt.ylim([np.min([np.min(current_band_1), np.min(current_band_2)])-.05,1.03])
    plt.gca().invert_xaxis()
    plt.legend(loc='best')
    #========================
    return [fig, nontransit_noise, wavelength]

def plot_wavelength_fit(wavelength, files, index, order_num = 4):
    #========================
    # get the corrected band flux (using the noise, and ignoring the last 10)
    # get the poly fit of degree 2 of the band (wavelength1)
    img, hdr = loadfits(files[index])
    current_band = get_current_band(img, order_num, nontransit_noise)
    wavelength1 = get_wavelength_fit(wavelength, current_band)
    #========================
    # plot the polynomial fit of deg 2 on top of the spectrum
    fig_fit = plt.figure()
    ax1 = fig_fit.add_subplot(111)
    ax1.plot(wavelength, current_band, label = 'Data')
    ax1.plot(wavelength, wavelength1,'k', label = 'Fit')
    ax1.set_title('Spectrum ' + str(index), fontsize=20)
    ax1.set_xlabel('Wavelength [nm]', fontsize=18)
    ax1.set_ylabel('I/I_0', fontsize=18)
    ax1.set_xlim([np.min(wavelength), np.max(wavelength)])
    fig_fit.gca().invert_xaxis()
    ax1.legend(loc='best')
    #========================
    # plot the spectrum minus the poly fit
    fig_flat = plt.figure()
    ax2 = fig_flat.add_subplot(111)
    flattened_flux = current_band-wavelength1
    ax2.plot(wavelength, flattened_flux)
    ax2.set_title('Spectrum ' + str(index) + ' after poly fit subtracted', fontsize=20)
    ax2.set_xlabel('Wavelength [nm]', fontsize=18)
    ax2.set_ylabel('I/I_0 - Fit', fontsize=18)
    ax2.set_xlim([np.min(wavelength), np.max(wavelength)])
    fig_flat.gca().invert_xaxis()
    #========================
    return [fig_fit, fig_flat]
    
def get_flux_for_correlation(wavelength, files, index, order_num = 4):
    img, hdr = loadfits(files[index])
    current_band = get_current_band(img, order_num, nontransit_noise)
    wavelength1 = get_wavelength_fit(wavelength, current_band)  
    return current_band - wavelength1

transit_arr_11_11 = [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 
37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 
56, 57, 58, 59, 60, 61, 62, 63, 64, 65]
order_num = 4
fig_two_indexes, nontransit_noise, wavelength = display_band_sun(sun_1_files, order_num)
plt.show(fig_two_indexes)
index = 63

#for index in transit_arr_11_11:
#    fig_fit, fig_flat = plot_wavelength_fit(wavelength, sun_1_files, index)
#    plt.show(fig_fit)
#    plt.show(fig_flat)

def ham_it_up(arr):
    return np.hamming(len(arr)) * arr
    
def get_crossed(flux1, flux2):
    flux1 = ham_it_up(flux1)
    flux2 = ham_it_up(flux2)
    return np.correlate(flux1,flux2,'same')

def get_shift_axis(cross_corr):
    n_pix=len(cross_corr)
    if n_pix % 2: 
        print 'y'
        shift_axis=np.arange(-n_pix/2.+1,n_pix/2.+1)
    else: 
        shift_axis=np.arange(-n_pix/2.,n_pix/2.)
    return shift_axis
def get_shift(shift_axis, cross_corr, num = 3):
    #gets the indexes of the max 3 points
    max_three = heapq.nlargest(num, range(len(cross_corr)), cross_corr.__getitem__)
    max_three = sorted(max_three)
    #print max_three
    #print shift_axis[max_three]
    #print cross_corr[max_three]
    fit = [shift_axis[max_three], cross_corr[max_three]]
    quad = LLS.linleastsquares(fit, 3)
    x = np.linspace(-1,1, 100)
    max_three_fit = calibrate_quad(quad, x)
    shift_max_y = np.amax(max_three_fit)
    shift_max_x = x[np.argmax(max_three_fit)] 
    return [shift_max_x, shift_max_y, x, max_three_fit]
    
def corr_plot(flux1, flux2):
    #calculate the cross correlation (y-axis)
    cross_corr = get_crossed(flux1, flux2)
    #calculate shift axis (x-axis)
    shift_axis = get_shift_axis(cross_corr)
    
    fig = plt.figure()
    plt.subplot(2,1,1)
    plt.plot(one,'b',linewidth=2)
    plt.plot(two,'r',linewidth=2)
    plt.xlabel('X Axis', fontsize=16) ; plt.ylabel('Y Axis [Unit]', fontsize=16)
    plt.title('Cross-Correlation', fontsize=18)

    shift_max_x, shift_max_y, x, max_three_fit = get_shift(shift_axis, cross_corr)
    
    plt.subplot(2,1,2)
    plt.plot(shift_axis,cross_corr,'g',linewidth=1)
    plt.plot(shift_axis,cross_corr,'bo',linewidth=1)
    plt.axvline(x = shift_max_x)
    plt.plot(x, max_three_fit,'k', linewidth = 2, label = str(shift_max_x))
    plt.xlim([-2,2])
    plt.xlabel('Shift Axis', fontsize=16) ; plt.ylabel('Cross-Correlation [Unit$^2$]', fontsize=16)
    plt.legend(loc=0)
    plt.tight_layout()
    return fig, shift_max_x
folder = sun_1_files
one = get_flux_for_correlation(wavelength, folder, 60)
two = get_flux_for_correlation(wavelength, folder, 44)
#plt.show(corr_plot(one, two))

#t_0_flux = get_flux_for_correlation(wavelength, folder, 44)
#for index in range(len(folder)):
#    t_0_flux += get_flux_for_correlation(wavelength, folder, index)
#t_0_flux /= len(folder)
#shift_arr = []
#for index in transit_arr_11_11:
#    curr_flux = get_flux_for_correlation(wavelength, folder, index)
#    plt.show(corr_plot(curr_flux, t_0_flux))

def shift(t_0_index, transit_arr, wavelength, folder):
    #=========================
    # get the shifted array
    t_0_flux = get_flux_for_correlation(wavelength, folder, t_0_index)
    for index in range(len(folder)):
        t_0_flux += get_flux_for_correlation(wavelength, folder, index)
    t_0_flux /= len(folder)
    shift_arr = []
    for index in transit_arr_11_11:
        curr_flux = get_flux_for_correlation(wavelength, folder, index)
        cross_corr = get_crossed(curr_flux, t_0_flux)
        shift_axis = get_shift_axis(cross_corr)
        x_shift = get_shift(shift_axis, cross_corr)[0]
        shift_arr.append(x_shift)
    print shift_arr
    shift_fig = plt.figure()
    plt.plot(shift_arr, 'o')
    plt.ylabel('Pixel shift')
    
    fit_fig = plt.figure()
    x = np.arange(len(shift_arr))
    yint, m = LLS.linleastsquares([x, shift_arr],2)
    shift_fit = x*m + yint
    plt.plot(x, shift_arr, 'o')
    plt.plot(x, shift_fit,'r', label = str(m) + '   ' + str(yint))
    plt.title('Doppler Shift ',size=20)
    plt.xlabel('Time', size=18)
    plt.ylabel('Pixel Shift', size=18)
    plt.legend()
    
    residuals_fig = plt.figure()
    res = shift_fit - shift_arr
    plt.plot(x,res,'go')
    return shift_arr, shift_fig, fit_fig, residuals_fig
# 44th spectrum is the t_0, so go from transit_arr[0] to transit_arr[-1]
t_0_index = 44
shift_arr, shift_fig, fit_fig, residuals_fig = shift(t_0_index, transit_arr_11_11, wavelength, folder)
#plt.show(shift_fig)

shift = np.max(shift_arr)
print 'SHIFT', shift

def calculate_vrot(lamb_0, lamb_shift):
    c = 299792 #km/s
    return (lamb_shift-lamb_0)/lamb_0 * c

def get_speed(shift, order_num = 4):
    centroids, spec = cen_spec_dict[order_num]
    centroids, spec = sorted(centroids), sorted(spec)[::-1] # reverse the spec cause short lambda (left) to long lambda (right)
    fit = [centroids, spec]
    quad = LLS.linleastsquares(fit, 3) #get the a_0 and a_1 and a_2
    print quad
    vel_arr = [] #keeps all the velocities which will be averaged later
    for i in range(1048):
        lambda_new = calibrate_quad(quad, shift+i) # calculate the wavelength at pixel = (shift+i)
        lambda_0 = calibrate_quad(quad, i) # calculate the wavelength at pixel = i
        vel_arr.append(calculate_vrot(lambda_0, lambda_new))
    return np.sum(vel_arr)/len(vel_arr) # take the average (v will be in m/s)

def get_sun_radius(vel, period):
    period = period * 24 * 3600 # get period from days to seconds
    return (1./(2*np.pi)) * vel * period