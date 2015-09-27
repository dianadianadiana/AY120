import numpy as np
import matplotlib.pyplot as plt
import os 

#greenlaser, sun, incan, LED, mecury, neon, redlaser,
path = "/Users/Diana/Desktop/Astro 120/Lab2/data/baes_sun/"
files = os.listdir(path) # read in all 100 files
filename = files[50] # choose one file
data = np.loadtxt(path+filename,skiprows=17,usecols=(0,1) )
data=np.transpose(data)

pixel_arr = data[0]
intensity_arr = data[1]
    
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
#print peaks

def peaks(intensity_arr, n = 20, lower_limit = 500, max = True):
    if max:
        peaks = max_peaks(intensity_arr)
    else:
        peaks = min_peaks(intensity_arr)
    peaks = limit_applier(peaks, intensity_arr, lower_limit)
    peaks = peak_verifier(peaks, intensity_arr, n, max)
    return peaks

peaks = peaks(intensity_arr, n = 130, lower_limit = 500, max = False)

fig = plt.figure()
plt.plot(pixel_arr, intensity_arr, linewidth =.5)
plt.xlabel('Pixels') ; plt.ylabel('Intensity')

for peak in peaks:
    plt.axvline(pixel_arr[peak], c='r',ls ='--', label = str(pixel_arr[peak]))
    
plt.legend()
plt.show()

#figure_path = "/Users/Diana/Desktop/"
#fig.savefig(figure_path + "test.png",dpi=200)
