import numpy as np
import matplotlib.pyplot as plt
import os 

#from matplotlib import rc
#rc('font',**{'family':'serif','serif':['Computer Modern']})
#rc('text', usetex=True)
#from pylab import rcParams
#rcParams['figure.figsize'] = 10, 5

#greenlaser, sun, incan, LED, mecury, neon, redlaser,
#path = "/Users/Diana/Desktop/Astro 120/Lab2/data/baes_neon/"
#files = os.listdir(path) # read in all 100 files
#filename = files[50] # choose one file
#data = np.genfromtxt(path+filename,skip_header=17, skip_footer=1)
#data=np.transpose(data)

def get_folder(data_path, name):
    return data_path + 'baes_' + name + '/'
def my_mean(arr):
    return np.sum(arr)/len(arr)
    
#greenlaser, sun, incan, LED, mecury, neon, redlaser,
data_path = "/Users/Diana/Desktop/Astro 120/Lab2/data/"
figure_path = "/Users/Diana/Desktop/Astro 120/Lab2/Figures/"

folder_greenlaser = get_folder(data_path,'greenlaser')#all the file names in this folder
folder_neon = get_folder(data_path,'neon')
folder_incan = get_folder(data_path,'incan')
folder_mercury = get_folder(data_path, 'mercury')

def average_intensity(folder):
    files = os.listdir(folder)
    intensity_all_arr = []
    for fil in files:
        if fil == '.DS_Store':
            continue
        data = np.transpose(np.genfromtxt(folder+fil,skip_header=17, skip_footer=1))
        intensity_all_arr.append(data[1])
    intensity_all_arr = np.transpose(intensity_all_arr)
    intensity_arr = []
    for elem in intensity_all_arr:
        intensity_arr.append(np.mean(elem))
    return np.round(intensity_arr)


pixel_arr = np.arange(2048)
intensity_neon = average_intensity(folder_neon)
intensity_mercury = average_intensity(folder_mercury)

intensity_arr = intensity_neon
name = "neon"

def min_peaks(arr):
    """ For Absorption -- Returns an array of the indexes where the indexes indicate where the minima are"""
    return [i for i in range(1, len(arr)-1) if arr[i-1]>arr[i] and arr[i+1]>arr[i]]
def max_peaks(arr):
    """ For Emission -- Returns an array of the indexes where the indexes indicate where the peaks are"""
    return [i for i in range(1, len(arr)-1) if arr[i-1]<arr[i] and arr[i+1]<arr[i]]

def limit_applier(index_arr, arr, lower_limit = 500):
    """ Makes sure the considered indexes are above a certain value"""    
    return [i for i in index_arr if arr[i] > lower_limit]
    
def peak_verifier(index_arr, arr, n, max = True):
    """
        The purpose of this function is to recognize if two points are considered to be
            on the same "peak", and if they are, then the point that is the true
            max(emission) or min (absorption) of the peak will be chosen, 
            whereas the other one will be discarded.
        Parameters:
            index_arr: this is the peak_indexes or harmonics_indexes
            arr: freq*constant
            n: the max number of data points between two points
            max - bool: True if looking for the max, False for the min
        Returns:
            the index_arr with the indexes of the true peak maxes or mins
    """
    k = 0
    delete_arr = []
    while k < len(index_arr) - 1:
   
        curr_index = index_arr[k]
        next_index = index_arr[k+1]

        if np.abs(curr_index-next_index) <= n:
            curr_lower, curr_upper = curr_index - n/2., curr_index + n/2.
            if curr_lower < 0:
                curr_lower = 0
            if max:
                curr_max = np.amax(arr[curr_lower:curr_upper])
            else:
                curr_min = np.amin(arr[curr_lower:curr_upper])
            if (max and arr[curr_index] == curr_max) or (not max and arr[curr_index] == curr_min):
                delete_arr.append(k+1)
            else:
                delete_arr.append(k)  
        k+=1
    return np.delete(index_arr, delete_arr)

def peaks(intensity_arr, n = 10, lower_limit = 200, max = True):
    if max:
        peaks = max_peaks(intensity_arr)
    else:
        peaks = min_peaks(intensity_arr)
    peaks = limit_applier(peaks, intensity_arr, lower_limit)
    peaks = peak_verifier(peaks, intensity_arr, n, max)
    return peaks

peaks_neon = peaks(intensity_neon, n = 10, lower_limit = 200, max = True)
peaks_mercury = peaks(intensity_mercury, n = 10, lower_limit = 200, max = True)
peaks_mercury = np.delete(peaks_mercury,[3,4,7,8,11,14, 15,16,17, 18,19]) #deleting extra ones
#peaks_mercury = np.append(peaks_mercury, 1412)
#center of mass of peak=centroid
#http://ugastro.berkeley.edu/infrared09/PDF-2009/centroid-error.pdf


def centroids(index_arr, x_arr, y_arr, peak_width):
        n = peak_width/2
	centroids= []
	print n
	for peak_index in index_arr:
		x_range = x_arr[peak_index-n:peak_index+n]
		y_range = y_arr[peak_index-n:peak_index+n]
		centroid = np.sum(x_range*y_range)/np.sum(y_range) #<x>

		#numerator = []
		#for i in range(len(x_range)):
		#    numerator.append(y_arr[i]*(x_arr[i]-centroid)**2)
		#error = np.sqrt( np.sum(numerator) / (np.sum(y_range))**2 )
		#centroids.append([centroid, error])
		centroids.append(centroid)
	#centroids, centroid_errors = np.transpose(centroids)[0], np.transpose(centroids)[1]
	#return centroids, centroid_errors
	return centroids
    
centroids_neon = centroids(peaks_neon, pixel_arr, intensity_neon, peak_width = 10)
centroids_mercury= centroids(peaks_mercury, pixel_arr, intensity_mercury, peak_width = 10)
dark_correction_neon = my_mean(intensity_neon[:1100]) -50
dark_correction_mercury = my_mean(intensity_mercury[:200])-50
intensity_neon -= dark_correction_neon
intensity_mercury -= dark_correction_mercury
def spectra_fig(name):
    fig = plt.figure()
    plt.plot(pixel_arr, intensity_arr, linewidth =.75)
    plt.xlabel('Pixels [count]') ; plt.ylabel('Intensity [ADU]')
    plt.scatter(centroids_neon, intensity_neon[peaks_neon], color='k', s= 5)
    plt.title('Spectra of: ' + name.upper() + ' -- ' +str(len(centroids)) + ' peaks')
    for centroid in centroids_neon:
        plt.axvline(x=centroid, c='r',ls ='--', label = str(centroid), linewidth =.75)
    plt.xlim([0,2048]) ; plt.ylim(bottom=0)
    plt.legend(loc=2, fontsize =10)
    return fig
    
#plt.show(spectra_fig(name))
#spectra_fig(name).savefig(figure_path + name + "_spectra.png",dpi=300)

def spectra_neon_mercury():
    fig = plt.figure(figsize = (12,5))
    ax = fig.add_subplot(121)
    ax1 = fig.add_subplot(122)
    
    ax.plot(pixel_arr, intensity_neon, 'k',linewidth =.75, label='Number of Peaks: ' + str(len(centroids_neon)))
    ax1.plot(pixel_arr, intensity_mercury,'k', linewidth =.75, label='Number of Peaks: ' + str(len(centroids_mercury)))
    ax.set_ylabel('Intensity [ADU]', fontsize = 16)
    ax.set_title('Spectra of Neon (Ne I)', fontsize = 18)
    ax1.set_title('Spectra of Mercury (Hg I)', fontsize = 18)
    for centroid in centroids_neon:
        ax.axvline(x=centroid, c='k',ls ='--', linewidth =.75)
    for centroid in centroids_mercury:
        ax1.axvline(x=centroid, c='k',ls ='--', linewidth =.75)
    ax.set_xlim([1250,2000])
    ax1.set_xlim([0,1700])
    for ax in [ax,ax1]:
        ax.legend(loc=9, fontsize =13)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('Pixels [count]', fontsize = 16)
    return fig
#plt.show(spectra_neon_mercury())
#spectra_neon_mercury().savefig(figure_path + 'spectra_neon_mercury.png',dpi = 300)

neonspec=  np.array([585.249,588.189,594.483,597.553,603.0,607.434,609.616,
614.306,616.359,621.72,626.649,630.479,633.443,638.299,640.225,650.653,653.288,659.895,667.828,671.704])
mercuryspec = np.array([405.4,436.6,487.7,542.4,546.5,577.7,580.2,584.0,587.6,593.4,599.7,611.6,
625.7,631.1,650.8,662.6,687.7,693.7,707,712.3,760.0,811.0])
mercuryspec = mercuryspec[:9]

spec = np.array([585.249,588.189,594.483,597.553,603.0,607.434,609.616,
614.306,616.359,621.72,626.649,630.479,633.443,638.299,640.225,650.653,653.288,659.895,667.828,671.704,
405.4,436.6,487.7,542.4,546.5,577.7,580.2,584.0,587.6])#,593.4,599.7,611.6])
centroids = centroids_neon+ centroids_mercury
spec = sorted(spec)
centroids= sorted(centroids)

    
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

fit = [centroids, spec]
linear = linleastsquares(fit, 2)
quad = linleastsquares(fit, 3)

ideal_linear = linear[1]*centroids+linear[0]
ideal_quad = quad[2]*centroids*centroids + quad[1]*centroids + quad[0]


def idealfigs():
    fig = plt.figure(figsize = (12,5))
    plt.scatter(centroids_neon,neonspec,s=70,c='r',label='experimental, neon')
    plt.scatter(centroids_mercury,mercuryspec,s=70,c='b',label='experimental, mercury')
    
    plt.plot(centroids, ideal_linear,'--',label='linear best fit')
    plt.plot(centroids, ideal_quad,'k--',label='quadratic best fit')
    plt.title('Least Squares Fitting', fontsize =20);    
    plt.xlabel('Pixels',fontsize=20);    plt.ylabel('Wavelength [nm]',fontsize=20)
    #for i in range(len(centroids)):
    #    y, x0,x1 = ideal_quad[i], centroids[i]-centroid_errors[i]/2., centroids[i]+centroid_errors[i]/2.
    #    plt.plot((x0, x1), (y, y),c='k',ls ='-', linewidth =1)
    plt.xlim([np.amin(centroids), np.amax(centroids)])
    plt.legend(numpoints=1, loc = 4, fontsize = 18)
    plt.tick_params(labelsize=16)
    plt.tight_layout()
    return fig
    
#plt.show(idealfigs())
#idealfigs().savefig(figure_path + "least_squares_fitting.png",dpi=300)

def risiduals():
    fig = plt.figure(figsize = (12,5))
    lin_r = ideal_linear - spec 
    quad_r = ideal_quad-  spec
    plt.scatter(centroids, lin_r,90,'k','^',label = 'linear residuals')
    plt.scatter(centroids, quad_r,70,'r', 'o',label = 'quadratic residuals')
    plt.title('Residuals', fontsize = 18)
    plt.axhline(y = 0, ls='--',c='k')
    plt.xlim([0,2048])
    plt.xlabel('Pixels',fontsize=16);    plt.ylabel('Difference of Fit and Actual Value',fontsize=16)
    plt.legend(numpoints=1, loc = 9, fontsize = 14)
    plt.tick_params(labelsize=14)

    return fig
    
#plt.show(risiduals())
#risiduals().savefig(figure_path + "residuals.png",dpi=300)


