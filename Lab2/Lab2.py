import numpy as np
import matplotlib.pyplot as plt
import os 

#greenlaser, sun, incan, LED, mecury, neon, redlaser,
path = "/Users/Diana/Desktop/Astro 120/Lab2/data/baes_neon/"
files = os.listdir(path) # read in all 100 files
filename = files[50] # choose one file
data = np.genfromtxt(path+filename,skip_header=17, skip_footer=1)
data=np.transpose(data)

pixel_arr = data[0]
intensity_arr = data[1]
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
    
# use min_peaks for absorption and max_peaks for emission
#peaks = max_peaks(intensity_arr)
#peaks = limit_applier(peaks, intensity_arr, lower_limit = 250)
#peaks = peak_verifier(peaks, intensity_arr, 20, max = True)

def peaks(intensity_arr, n = 10, lower_limit = 500, max = True):
    if max:
        peaks = max_peaks(intensity_arr)
    else:
        peaks = min_peaks(intensity_arr)
    peaks = limit_applier(peaks, intensity_arr, lower_limit)
    peaks = peak_verifier(peaks, intensity_arr, n, max)
    return peaks

peaks = peaks(intensity_arr, n = 10, lower_limit = 200, max = True)

#center of mass of peak=centroid
def centroids(index_arr, x_arr, y_arr, peak_width):
        n = peak_width/2.
	centroids= []
	for peak_index in index_arr:
		x_range = x_arr[peak_index-n:peak_index+n]
		y_range = y_arr[peak_index-n:peak_index+n]
		centroid = np.sum(x_range*y_range)/np.sum(y_range)
		centroids.append(centroid)
	return centroids
    
print len(peaks)
centroids = centroids(peaks, pixel_arr, intensity_arr, 10)

def spectra_fig(name):
    fig = plt.figure()
    plt.plot(pixel_arr, intensity_arr, linewidth =.75)
    plt.xlabel('Pixels') ; plt.ylabel('Intensity')
    plt.scatter(centroids, intensity_arr[centroids], color='k', s= 5)
    plt.title('Spectra of: ' + name.upper())
    for centroid in centroids:
        plt.axvline(x=centroid, c='r',ls ='--', label = str(centroid), linewidth =.75)
    plt.grid(True)
    plt.xlim([0,2048]) ; plt.ylim(bottom=0)
    plt.legend(loc=2, fontsize =10)
    return fig
    
plt.show(spectra_fig(name))
figure_path = "/Users/Diana/Desktop/Astro 120/Lab2/Figures/"
#spectra_fig(name).savefig(figure_path + name + "_spectra.png",dpi=300)

neonspec=  np.array([585.249,588.189,594.483,597.553,603.0,607.434,609.616,
614.306,616.359,621.72,626.649,630.479,633.443,638.299,640.225,650.653,653.288,659.895,667.828,671.704])

fig = plt.figure()
plt.plot(neonspec, centroids, 'o')
#plt.plot(neonspec, centroids)

plt.show()

def linleastsquares(data): #data : x,y  data= np.array
    n = len(data[0])
    sum_x = np.sum(data[0])
    sum_y = np.sum(data[1])
    sum_xsq = np.sum(data[0] ** 2)
    sum_xy = np.sum(data[0]*data[1])
    
    left = np.array([[n,sum_x],[sum_x,sum_xsq]])
    right = np.array([[sum_y],[sum_xy]])

    inv_left = np.linalg.inv(left)
    final = np.dot(inv_left,right)
    return [final[0], final[1]]
    
def linleastsquares(data, poly):
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

def linleastsquares_quad(data):
    n = len(data[0])
    
    sum_x = np.sum(data[0])
    sum_x2 = np.sum(data[0] ** 2)
    sum_x3 = np.sum(data[0] ** 3)
    sum_x4 = np.sum(data[0] ** 4)
    sum_y = np.sum(data[1])
    sum_xy = np.sum(data[0]*data[1])
    sum_x2y = np.sum(data[0] ** 2 * data[1])
    
    left = np.array([[n,sum_x,sum_x2],[sum_x,sum_x2,sum_x3],[sum_x2,sum_x3,sum_x4]])
    right = np.array([[sum_y],[sum_xy],[sum_x2y]])

    inv_left = np.linalg.inv(left)
    final = np.dot(inv_left, right)
    return final
    #return [final[0], final[1], final[2]]
    
def linleastsquares_4(data):
    n = len(data[0])
    
    sum_x = np.sum(data[0])
    sum_x2 = np.sum(data[0] ** 2)
    sum_x3 = np.sum(data[0] ** 3)
    sum_x4 = np.sum(data[0] ** 4)
    sum_x5 = np.sum(data[0] ** 5)
    sum_x6 = np.sum(data[0] ** 6)
    sum_y = np.sum(data[1])
    sum_xy = np.sum(data[0]*data[1])
    sum_x2y = np.sum(data[0] ** 2 * data[1])
    sum_x3y = np.sum(data[0] ** 3 * data[1])
    
    left = np.array([[n,sum_x,sum_x2,sum_x3],[sum_x,sum_x2,sum_x3,sum_x4],[sum_x2,sum_x3,sum_x4,sum_x5],[sum_x3,sum_x4,sum_x5,sum_x6]])
    right = np.array([[sum_y],[sum_xy],[sum_x2y],[sum_x3y]])
	
    inv_left = np.linalg.inv(left)
    final = np.dot(inv_left, right)
    return [final[0], final[1], final[2], final[3]]
    
x = neonspec
y = centroids
data = [x,y]
linear = linleastsquares(data, 2)
ideal_linear = linear[1]*x+linear[0]

quad = linleastsquares(data, 3)
ideal_quad = quad[2]*x**2 + quad[1]*x + quad[0]

lls4 = linleastsquares(data, 4)
ideal_4 = lls4[3]*x**3 + lls4[2]*x**2 + lls4[1]*x + lls4[0]

def idealfigs():
    fig = plt.figure()
    plt.plot(neonspec, centroids,'o',label='experimental')
    plt.plot(neonspec, ideal_linear,'--',label='linear best fit')
    plt.plot(neonspec, ideal_quad,'--',label='quadratic best fit')
    plt.plot(neonspec, ideal_4,'--',label='4 best fit')

    plt.xlabel('Wavelength[nm]')
    plt.ylabel('Pixels')
    plt.legend()
    return fig
plt.show(idealfigs())