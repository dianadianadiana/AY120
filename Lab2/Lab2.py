import numpy as np
import matplotlib.pyplot as plt
import os 

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

def limit_applier(index_arr, arr, lower_limit = 250):
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
            next_lower, next_upper = next_index - n/2., next_index + n/2.

            if curr_lower < 0:
                curr_lower = 0
            if next_lower < 0:
                next_lower = 0
            if max:
                #curr_maxpow and next_maxpow should be the same
                curr_max = np.amax(arr[curr_lower:curr_upper])
                next_max = np.amax(arr[next_lower:next_upper])
            else:
                curr_min = np.amin(arr[curr_lower:curr_upper])
                next_min = np.amin(arr[next_lower:next_upper])

            if (max and arr[curr_index] == curr_min) or (not max and arr[curr_index] == curr_min):
                delete_arr.append(k+1)
            else:
                delete_arr.append(k)  
        k+=1
    return np.delete(index_arr, delete_arr)
    
# use min_peaks for absorption and max_peaks for emission
peaks = min_peaks(intensity_arr)
peaks = limit_applier(peaks, intensity_arr, lower_limit = 250)
peaks = peak_verifier(peaks, intensity_arr, 100, max = False)
print peaks

fig = plt.figure()
plt.plot(pixel_arr, intensity_arr)
for peak in peaks:
    plt.axvline(pixel_arr[peak], c='r',ls ='--')
figure_path = "/Users/Diana/Desktop/"
#fig.savefig(figure_path + "test.png",dpi=200)
plt.show()
