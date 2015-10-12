import numpy as np
import matplotlib.pyplot as plt
import os 
def my_mean(arr):
    return np.sum(arr)/len(arr)
def my_std(arr):
   """Calculates the STD of the arr and returns mu (the mean)"""
   mu = my_mean(arr)
   sq_sum = np.sum([(elem-mu)**2 for elem in arr]) #sums (x_i-mu)^2
   return mu, np.sqrt(sq_sum/(len(arr)-1)) 
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

def get_intensity_index(folder, index):
    """returns an intensity array of all the calculated intensities at a given index"""
    intensity_arr = np.array([])
    files = os.listdir(folder)
    for fil in files:
        if fil == '.DS_Store':
            continue
        data = np.transpose(np.genfromtxt(folder+fil,skip_header=17, skip_footer=1))
        intensity_arr = np.append(intensity_arr, data[1][index])
    return intensity_arr
    
def fig_10(folder, dark_folder, index):
    #the x axis is in a way the time where each event is 100ms
    intensity_arr = get_intensity_index(folder, index)
    dark_arr = get_intensity_index(dark_folder, index)
    dark_correction = my_mean(dark_arr)
    
    fig = plt.figure(figsize = (10,5))
    plt.plot(pixel_arr[:len(intensity_arr)], intensity_arr - dark_correction,'k', linewidth =1)
    plt.xlabel('Sample', fontsize =16); plt.ylabel('Signal [ADU]', fontsize =16)
    mean, std = my_std(intensity_arr-dark_correction)
    limit = 70
    plt.xlim([0,1000]) ; plt.ylim([mean - limit, mean+ limit])
    plt.tick_params(labelsize=14)
    plt.title('Measured Signals for Index ' + str(index) + ' (Time Series)\nmean ' + str(mean), fontsize=18)
    plt.tight_layout()
    return fig
    
def fig_10_two(folder, dark_folder,index, index1):
    intensity_arr = get_intensity_index(folder, index)
    intensity_arr1 = get_intensity_index(folder, index1)
    dark_arr = get_intensity_index(dark_folder, index)
    dark_arr1 = get_intensity_index(dark_folder, index1)

    dark_correction = my_mean(dark_arr)
    dark_correction1 = my_mean(dark_arr1)
    
    fig = plt.figure(figsize = (12,5))
    ax = fig.add_subplot(121)
    ax1 = fig.add_subplot(122)
    mean, std = my_std(intensity_arr-dark_correction)
    mean1, std1 = my_std(intensity_arr1-dark_correction1)
    limit = 80
    ax.set_ylim([mean - limit, mean + limit]); ax1.set_ylim([mean1 - limit, mean1 + limit]) 
    plt.suptitle('Measured Intensity Signals for a Given Index (Time Series)', fontsize=18)
    ax.plot(pixel_arr[:len(intensity_arr)], intensity_arr - dark_correction,'k', linewidth =1, 
    label='Pixel: '+str(index)+'\nMean = '+str(mean)+'\nVariance = '+"{:.2f}".format(std**2))
    "{:.3f}".format(std**2)
    ax1.plot(pixel_arr[:len(intensity_arr1)], intensity_arr1 - dark_correction1,'k', linewidth =1,
    label='Pixel: '+str(index1)+'\nMean = '+str(mean1)+'\nVariance = '+"{:.2f}".format(std1**2))
    ax.set_ylabel('Signal [ADU]', fontsize =16)
   
    for ax in [ax, ax1]:
        ax.set_xlabel('Sample', fontsize =16)
        ax.legend(loc='upper right', fancybox = True, fontsize = 14)
        ax.tick_params(labelsize=14)
    #plt.tight_layout()
    return fig
    
#plt.show(fig_10_two(folder_incan_100ms_90v, folder_incan_dark_100ms,500,1000))
#fig_10_two(folder_incan_100ms_90v, folder_incan_dark_100ms,500,1000).savefig(figure_path + "fig10_incan_1000_100ms_90v_500vs1000.png",dpi=300)

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
        
 
def linleastsquares(data, poly):
    """data: [x,y] what you are trying to fit
       poly: the number of terms you want i.e. poly = 3 --> ax**2 + bx + c"""
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
    intensity_arr, dark_arr = np.array([]), np.array([])
    files, dark_files = os.listdir(folder), os.listdir(dark_folder)
    
    filenum =1000 #1000 files
    indexnum = 2048 #2048 pixels

    #this array will eventually contain an indexnum amount of rows;
    #where in each row, the varying intensities at a certain index are stored
    #rows correspond to pixel count
    intensity_arr = np.zeros(shape=(filenum,indexnum)) #an array that has filenum rows and indexnum cols
    for index, fil in enumerate(files[:filenum]):
        if fil == ".DS_Store": #sometimes this appears in the folder and causes problems
            continue
        data = np.transpose(np.genfromtxt(folder+fil,skip_header=17, skip_footer=1))
        intensity_arr[index] = data[1][:indexnum]
    intensity_arr = np.transpose(intensity_arr)
    
    #likewise as intensity_arr, but accounting for the dark subtraction
    dark_arr = np.zeros(shape=(filenum,indexnum)) #an array that has filenum rows and indexnum cols
    for index, fil in enumerate(dark_files[:filenum]):
        if fil == ".DS_Store":
            continue
        data = np.transpose(np.genfromtxt(dark_folder+fil,skip_header=17, skip_footer=1))
        dark_arr[index] = data[1][:indexnum]
    dark_arr = np.transpose(dark_arr)
    
    #dark_corr_arr takes the average of the dark_arr at each pixel
    dark_corr_arr = np.array([])
    for elem in dark_arr:
        dark_corr_arr = np.append(dark_corr_arr, my_mean(elem))
   
    mean_arr = np.array([])
    variance_arr = np.array([])
    #get the mean (of the intensity_arr - dark subtraction) and the variance to plot
    for index, elem in enumerate(intensity_arr):
        mean, std = my_std(elem - dark_corr_arr[index])
        mean_arr, variance_arr = np.append(mean_arr, mean), np.append(variance_arr, std**2)
    
    linear = linleastsquares([mean_arr, variance_arr], 2)
    gain, ADU_0 = linear[1][0], linear[0][0]
    ideal_linear = linear[1]*mean_arr+linear[0]
    #this graph shows the correlation of how when the brightness increases, the variance increases too
    fig = plt.figure(figsize = (10,5))
    plt.plot(mean_arr, variance_arr, 'o', linewidth =1)
    plt.plot(mean_arr, ideal_linear,'r-',label='gain = ' + str(gain) + '\nb = ' + str(ADU_0), linewidth = 1.5)
    plt.xlabel("Mean, " + r'$\left(\overline{ADU - ADU_0}\right)$', fontsize =16); plt.ylabel("Variance, " + r'$\sigma^2$', fontsize =16)
    plt.tick_params(labelsize=14)
    plt.xlim([np.amin(mean_arr), np.amax(mean_arr)])
    plt.title('Correlation between Brightness and Variance', fontsize=18)
    #plt.yscale('log')
    plt.tight_layout()
    plt.legend(loc = 2, fontsize=16)
    return fig
      
#plt.show(fig_12(folder_incan_100ms_90v, folder_incan_dark_100ms))
#fig_12(folder_incan_100ms_90v, folder_incan_dark_100ms).savefig(figure_path + "fig12_100ms_90v_log.png", dpi =300)
fig_12(folder_incan_100ms_80v, folder_incan_dark_100ms).savefig(figure_path + "fig12_100ms_80v.png", dpi =300)

