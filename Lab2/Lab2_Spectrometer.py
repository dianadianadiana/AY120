import numpy as np
import matplotlib.pyplot as plt
import os 
def my_mean(arr):
    return np.sum(arr)/len(arr)
def get_folder(data_path, name):
    return data_path + 'baes_' + name + '/'
    
data_path = "/Users/Diana/Desktop/Astro 120/Lab2/data/"
figure_path = "/Users/Diana/Desktop/Astro 120/Lab2/Figures/"

folder_incan_dark_50ms = get_folder(data_path,'incan_dark_1000_50ms')
folder_incan_dark_100ms = get_folder(data_path,'incan_dark_1000_100ms')
folder_incan_100ms_80v = get_folder(data_path,'incan_1000_100ms_80v')
folder_incan_100ms_85v = get_folder(data_path,'incan_1000_100ms_85v')
folder_incan_100ms_90v = get_folder(data_path,'incan_1000_100ms_90v')
folder_incan_50ms_90v = get_folder(data_path,'incan_1000_50ms_90v')

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
############
#intensity_arr_100_dark = average_intensity(folder_incan_dark_100ms)
#intensity_arr_100_80v = average_intensity(folder_incan_100ms_80v)
#intensity_arr_100_85v = average_intensity(folder_incan_100ms_85v)
#intensity_arr_100_90v = average_intensity(folder_incan_100ms_90v)

def fig_9():
    fig = plt.figure(figsize = (10,5))
    plt.plot(pixel_arr, intensity_arr_100_80v - intensity_arr_100_dark, linewidth =1, label = '80V')
    plt.plot(pixel_arr, intensity_arr_100_85v - intensity_arr_100_dark, linewidth =1, label = '85V')
    plt.plot(pixel_arr, intensity_arr_100_90v - intensity_arr_100_dark, linewidth =1, label = '90V')
    plt.xlabel('Pixels [count]',fontsize=16) ; plt.ylabel('Intensity/Signal [ADU]', fontsize=16)
    plt.title('Spectra of Incandescent Lamp, 100ms', fontsize = 18)
    plt.xlim([0,2048]) ; plt.ylim(bottom=0)
    plt.legend(loc=2, fontsize =16)
    plt.tick_params(labelsize=14)
    plt.tight_layout()
    return fig
    
#plt.show(fig_9())
#fig_9().savefig(figure_path + "incan_80v85v90v_fig9.png",dpi=300)

def fig_10(folder, dark_folder, index):
    #the x axis is in a way the time where each event is 100ms
    intensity_arr = np.array([])
    dark_arr = np.array([])
    files = os.listdir(folder)
    files1 = os.listdir(dark_folder)
    for fil in files:
        if fil == '.DS_Store':
            continue
        data = np.transpose(np.genfromtxt(folder+fil,skip_header=17, skip_footer=1))
        intensity_arr = np.append(intensity_arr, data[1][index])
    for fil in files1:
        if fil == '.DS_Store':
            continue
        data = np.transpose(np.genfromtxt(dark_folder+fil,skip_header=17, skip_footer=1))
        dark_arr = np.append(dark_arr, data[1][index])
    dark_correction = np.sum(dark_arr)/len(dark_arr)
    fig = plt.figure(figsize = (10,5))
    plt.plot(pixel_arr[:len(intensity_arr)], intensity_arr - dark_correction,'k', linewidth =1)
    plt.xlabel('Sample', fontsize =16); plt.ylabel('Signal [ADU]', fontsize =16)
    mean = my_mean(intensity_arr)
    mean = np.mean(intensity_arr-dark_correction)
    plt.xlim([0,1000]) ; plt.ylim([mean - 70, mean+ 70])
    plt.tick_params(labelsize=14)
    plt.title('Measured Signals for Index ' + str(index) + ' (Time Series)\nmean ' + str(mean), fontsize=18)
    plt.tight_layout()
    return fig
    
#plt.show(fig_10(folder_incan_50ms_90v, folder_incan_dark_50ms, 1000))
#fig_10(folder_incan_50ms_90v, folder_incan_dark_50ms, 1000).savefig(figure_path + "incan_1000_50ms_90v_fig10.png",dpi=300)
#plt.show(fig_10(folder_incan_100ms_90v, folder_incan_dark_100ms, 1000))
#plt.show(fig_10(folder_incan_100ms_80v, folder_incan_dark_100ms, 1000))
#plt.show(fig_10(folder_incan_100ms_85v, folder_incan_dark_100ms, 1000))
#fig_10(folder_incan_100ms_90v, folder_incan_dark_100ms, 1000).savefig(figure_path + "incan_1000_100ms_90v_fig10.png",dpi=300)
#fig_10(folder_incan_100ms_80v, folder_incan_dark_100ms, 1000).savefig(figure_path + "incan_1000_100ms_80v_fig10.png",dpi=300)
#fig_10(folder_incan_100ms_85v, folder_incan_dark_100ms, 1000).savefig(figure_path + "incan_1000_100ms_85v_fig10.png",dpi=300)

#sampling_rate = 100 #in milliseconds
#n_freq = .5 * sampling_rate
#print np.float((np.amax(pixel_arr)-np.amin(pixel_arr)))/(len(pixel_arr)-1)

        
def my_std(arr):
   """Calculates the STD of the arr and returns mu (the mean)"""
   mu = my_mean(arr)
   sq_sum = np.sum([(elem-mu)**2 for elem in arr]) #sums (x_i-mu)^2
   return mu, np.sqrt(sq_sum/(len(arr)-1))  
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
def fig_12(folder, dark_folder):
    intensity_arr = np.array([])
    dark_arr = np.array([])
    files = os.listdir(folder)
    files1 = os.listdir(dark_folder)
    mean_arr = np.array([])
    
    filenum =1000 #1000 files
    indexnum = 2048 #2048 pixels
    intensity_arr = np.zeros(shape=(filenum,indexnum))
    for index, fil in enumerate(files[:filenum]):
        if fil == ".DS_Store":
            continue
        data = np.transpose(np.genfromtxt(folder+fil,skip_header=17, skip_footer=1))
        intensity_arr[index] = data[1][:indexnum]
    intensity_arr = np.transpose(intensity_arr)
    #print intensity_arr
    
    dark_arr = np.zeros(shape=(filenum,indexnum))
    for index, fil in enumerate(files1[:filenum]):
        if fil == ".DS_Store":
            continue
        data = np.transpose(np.genfromtxt(dark_folder+fil,skip_header=17, skip_footer=1))
        dark_arr[index] = data[1][:indexnum]
    dark_arr = np.transpose(dark_arr)
    #print dark_arr

    
    dark_corr_arr = np.array([])
    for elem in dark_arr:
        dark_corr_arr = np.append(dark_corr_arr, my_mean(elem))
    #print dark_corr_arr
   
    mean_arr = np.array([])
    variance_arr = np.array([])
    for index, elem in enumerate(intensity_arr):
        mean, std = my_std(elem - dark_corr_arr[index])
        mean_arr, variance_arr = np.append(mean_arr, mean), np.append(variance_arr, std**2)
    #print mean_arr
    #print variance_arr
    
    linear = linleastsquares([mean_arr, variance_arr], 2)

    ideal_linear = linear[1]*mean_arr+linear[0]
    #this graph shows the correlation of how when the brightness increases, the variance increases too
    fig = plt.figure(figsize = (10,5))
    plt.plot(mean_arr, variance_arr, 'o', linewidth =1)
    plt.plot(mean_arr, ideal_linear,'r-',label='gain = ' + str(linear[1]) + '\nb = ' + str(linear[0]), linewidth = 1.5)
    plt.xlabel("Mean, " + r'$\left(\overline{ADU - ADU_0}\right)$', fontsize =16); plt.ylabel("Variance, " + r'$\sigma^2$', fontsize =16)
    plt.tick_params(labelsize=14)
    plt.title('Correlation between Brightness and Variance', fontsize=18)
    plt.tight_layout()
    plt.legend(loc = 2)
    return fig
    
    
plt.show(fig_12(folder_incan_100ms_90v, folder_incan_dark_100ms))

fig_12(folder_incan_100ms_90v, folder_incan_dark_100ms).savefig("/Users/Diana/Desktop/astro temp/test.png", dpi =300)