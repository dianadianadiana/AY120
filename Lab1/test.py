# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os

path = "/Users/Diana/Desktop/Astro 120/Lab1/data/"
file_arr = os.listdir(path) #import all the files into an array
print file_arr
filename = file_arr[3] #"baes20_150902_1956_43.csv"

##################### First Question ########################
x = np.loadtxt(path+filename, delimiter=',', dtype = np.int32)
t = x[:,1]

fig=plt.figure()
plt.plot(t,'o')
plt.xlabel('Event Number')
plt.ylabel('Clock Tick')
plt.show(fig)
#plt.savefig('eventnumber10000.png')
#plt.close(fig)

# Compute the intervals
dt = t[1:] - t[0:-1]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(dt,',')
ax1.set_xlabel('Event number')
ax1.set_ylabel('Interval [clock ticks]')
#fig.savefig('/Users/Diana/Desktop/Astro 120/Lab1/fig_int32.png', dpi = 50)
fig.show()

##################### Second Question ######################
# Compute the means in chunks
marr = np.array([])
i = np.arange(dt.size)
nstep = 1000
for j in i[0::nstep]:
    print j, j+nstep-1,np.mean(dt[j:j+nstep])
    m = np.mean(dt[j:j+nstep])
    marr = np.append(marr, m)

fig1 = plt.figure()
ax1_1 = fig1.add_subplot(111)
ax1_1.plot(marr, 'o')
ax1_1.set_xlabel("Start Index (10^3)")
fig1.show()

#################### Figure 5 -- ########################
mean_arr = np.array([])
datastep = 100
count = 100
while count < dt.size:
    mean = np.mean(dt[0:count+datastep])
    print count, count+datastep, mean
    mean_arr = np.append(mean_arr, mean)
    count += datastep
x_arr = np.arange(1,len(mean_arr)+1)*100
print mean_arr, x_arr

fig2 = plt.figure()
ax1_2 = fig2.add_subplot(111)
ax1_2.plot(x_arr, mean_arr, 'o')
ax1_2.set_xlabel('Number of intervals averaged')
ax1_2.set_ylabel('Mean Interval [clock ticks]')
#fig2.savefig('/Users/Diana/Desktop/Astro 120/Lab1/fig5.png', dpi = 50)
fig2.show()

################### Standard deviation part #############
mu = np.sum(marr)/np.float(marr.size)
std = np.sqrt(np.sum((marr-mu)**2.)/(np.float(marr.size)-1.)) #standard deviation
print mu, std, marr[0]

########################## s/sqrtN #########################
count = 10
rang = np.arange(10,500,10)
stddev_arr= np.array([])
for j in rang:
	ith_mean_arr=np.array([])
	for i in np.arange(j,10000,j):
		ith_mean=np.mean(dt[0:i])
		ith_mean_arr = np.append(ith_mean_arr,ith_mean)
	mu_i = np.sum(ith_mean_arr)/np.float(ith_mean_arr.size) 
	stddev_i = np.sqrt(np.sum((ith_mean_arr-mu_i)**2.)/(np.float(ith_mean_arr.size)-1.))
	stddev_arr = np.append(stddev_arr,stddev_i)

fig3 = plt.figure()
ax1_3 = fig3.add_subplot(111)
ax1_3.plot(np.arange(10,500,10),stddev_arr,'ko')
ax1_3.set_xlabel('Number of Events Averaged')
ax1_3.set_ylabel('Standard Deviation of the Mean')
fig3.show()

fig4 = plt.figure()
ax1_4 = fig4.add_subplot(111)
m, b = np.polyfit(1./np.sqrt(np.arange(10,500,10)),stddev_arr,1)
ax1_4.plot(1./np.sqrt(np.arange(10,500,10)),m/np.sqrt(np.arange(10,500,10)) +b )
ax1_4.plot(1./np.sqrt(np.arange(10,500,10)),stddev_arr,'ro')
ax1_4.set_xlabel('1/sqrt(N)')
ax1_4.set_ylabel('Standard deviation of the mean [ticks]')
fig4.show()

######################### Histogram #########################

N = 500
# define the lower and upper bin edges and bin width
bw = (dt.max()-dt.min())/(N-1.)
bin1 = dt.min() + bw*np.arange(N)

# define the array to hold the occurrence count
bincount = np.array([])
# loop through the bins
for bin in bin1:
    count = np.where((dt >= bin) & (dt < bin+bw))[0].size
    bincount = np.append(bincount, count)
    
# compute bin centers for plotting
binc = bin1 + .5*bw
fig = plt.figure()
plt.plot(binc, bincount, drawstyle='steps-mid')
plt.xlabel('Interval [ticks]')
plt.ylabel('Frequency')
plt.show(fig)
#plt.savefig(figpath + "")

def poisson(deltat, tau):
    return (i/tau)*exp(-deltat/tau)
    


########################### OTher #######################
#for filename in file_arr:
#    x = np.loadtxt(path+filename, delimiter=',', dtype = np.int32)
#    t = x[:,1]
#    marr = np.array([])
#    i = np.arange(dt.size)
#    nstep = 1000
#    for j in i[0::nstep]:
#        m = np.mean(dt[j:j+nstep])
#        marr = np.append(marr, m)
#    mu = np.sum(marr)/np.float(marr.size)
#    std = np.sqrt(np.sum((marr-mu)**2.)/(np.float(marr.size)-1.)) #standard deviation
#
#    print mu, std

# Compute the intervals
#x1 = np.loadtxt(path+filename, delimiter=',')
#t1 = x1[:,1]
#dt1 = t1[1:] - t1[0:-1]
#
#fig1 = plt.figure()
#ax2 = fig1.add_subplot(111)
#ax2.plot(dt1,',')
#ax2.set_xlabel('Event number')
#ax2.set_ylabel('Interval [clock ticks]')
#fig1.savefig('/Users/Diana/Desktop/Astro 120/Lab1/fig_no_int32.png', dpi = 50)