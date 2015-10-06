import numpy as np
import matplotlib.pyplot as plt
import os 

def get_folder(data_path, name):
    return data_path + 'baes_' + name + '/'
    
data_path = "/Users/Diana/Desktop/Astro 120/Lab2/data/"
folder_incan_50ms = get_folder(data_path,'incan_1000_50ms')#all the file names in this folder
folder_incan_dark_50ms = get_folder(data_path,'dark_incan_50ms')
folder_incan_100ms = get_folder(data_path,'incan_1000_100ms')
folder_incan_dark_100ms = get_folder(data_path,'dark_incan_100ms')

folder_incan_dark_100ms 
folder_incan_100ms_80v
folder_incan_100ms_85v
folder_incan_100ms_90v
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
intensity_arr_50 = average_intensity(folder_incan_50ms)
intensity_arr_50_dark = average_intensity(folder_incan_dark_50ms)
intensity_arr_100 = average_intensity(folder_incan_100ms)
intensity_arr_100_dark = average_intensity(folder_incan_dark_100ms)

def fig_9():
    fig = plt.figure(figsize = (10,5))
    plt.plot(pixel_arr, intensity_arr_50 - intensity_arr_50_dark, linewidth =1, label = '50ms')
    plt.plot(pixel_arr, intensity_arr_100 - intensity_arr_100_dark, linewidth =1, label = '100ms')
    plt.xlabel('Pixels [count]',fontsize=16) ; plt.ylabel('Intensity/Signal [ADU]', fontsize=16)
    plt.title('Spectra of Incandescent Lamp, 50ms vs. 100ms', fontsize = 18)
    plt.xlim([0,2048]) ; plt.ylim(bottom=0)
    plt.legend(loc=2, fontsize =16)
    plt.tick_params(labelsize=14)
    plt.tight_layout()
    return fig
    
plt.show(fig_9())
figure_path = "/Users/Diana/Desktop/Astro 120/Lab2/Figures/"
#fig_9().savefig(figure_path + "incan_50vs100ms_fig9.png",dpi=300)

def fig_10(folder, dark_folder, index):
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
    plt.xlim([0,1000]) #; plt.ylim(bottom=0)
    plt.tick_params(labelsize=14)
    plt.title('Measured Signals for Index ' + str(index), fontsize=18)
    plt.tight_layout()
    return fig
#plt.show(fig_10(folder_incan_50ms, folder_incan_dark_50ms, 1000))
fig_10(folder_incan_50ms, folder_incan_dark_50ms, 1000).savefig(figure_path + "incan_50ms_1000_fig10.png",dpi=300)
#plt.show(fig_10(folder_incan_100ms, folder_incan_dark_100ms, 1000))
fig_10(folder_incan_100ms, folder_incan_dark_100ms, 1000).savefig(figure_path + "incan_100ms_1000_fig10.png",dpi=300)

    