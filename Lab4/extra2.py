import numpy as np
import matplotlib.pyplot as plt

dates = ['1_11_12', '2_11_12', '1_11_16', '2_11_16']
rot_v_dict = {'1_11_12': 1.99472712, '2_11_12': 2.22551396, 
'1_11_16':1.80421032, '2_11_16': 2.39975034}
sol_radius_dict = {'1_11_12': 719749.20893956, '2_11_12': 803023.07993681, 
'1_11_16':651005.8146031, '2_11_16': 865892.07838828}
earth_sun_dist_dict = {'1_11_12': 1.54703354e+08, '2_11_12': 1.72602293e+08, 
'1_11_16':1.39927605e+08, '2_11_16': 1.86115395e+08}

def my_std(arr):
    """Calculates the STD of the arr"""
    mu = np.mean(arr)
    sq_sum = np.sum([(elem-mu)**2 for elem in arr]) #sums (x_i-mu)^2
    return np.sqrt(sq_sum/(len(arr)-1))
def my_sdom(arr):
    return my_std(arr)/(np.sqrt(len(arr)))
    
def plot_rot_v(dates, rot_v_dict, sol_radius_dict, earth_sun_dist_dict):
    rot_v_arr = []
    sol_radius_arr = []
    earth_sun_arr = []
    for elem in dates:
        rot_v_arr.append(rot_v_dict[elem])
        sol_radius_arr.append(sol_radius_dict[elem])
        earth_sun_arr.append(earth_sun_dist_dict[elem])
    print rot_v_arr
    print sol_radius_arr 
    print earth_sun_arr
    print "STD"
    print my_std(rot_v_arr)
    print my_std(sol_radius_arr)
    print my_std(earth_sun_arr)
    print "SDOM"
    print my_sdom(rot_v_arr)
    print my_sdom(sol_radius_arr)
    print my_sdom(earth_sun_arr)
    
    fig = plt.figure()
    plt.plot(rot_v_arr, 'o')
    plt.xlim([-.5, 3.5])
    plt.ylim([1.5,2.5])
    plt.annotate(dates[0], (0, rot_v_arr[0]), fontsize = 12)
    plt.annotate(dates[1], (1, rot_v_arr[1]), fontsize = 12)
    plt.annotate(dates[2], (2, rot_v_arr[2]), fontsize = 12)
    plt.annotate(dates[3], (3, rot_v_arr[3]), fontsize = 12)
    return fig
    
plt.show(plot_rot_v(dates, rot_v_dict, sol_radius_dict, earth_sun_dist_dict))

v_rot = 2.10604


T_rot = (26.24*24*3600)
radius = (1/(2*np.pi))*v_rot*T_rot
print 'solar radius = ', radius
angsun = np.radians(1919.3/3600.)
my_d = (2*radius)/np.sin(angsun) # km
print 'earth sun dist= ', my_d
