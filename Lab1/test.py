import numpy as np
import matplotlib.pyplot as plt

path = "/Users/Diana/Desktop/Astro 120/Lab1/"
filename = "baes20_150902_1956_43.csv"

x = np.loadtxt(path+filename, delimiter=',', dtype = np.int32)
t = x[:,1]

# Compute the intervals
dt = t[1:] - t[0:-1]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(dt,',')
ax1.set_xlabel('Event number')
ax1.set_ylabel('Interval [clock ticks]')
plt.show()

plt.savefig(path + 'figtest.png', dpi = 50)

## Compute the means in chunks
#marr = np.array([])
#i = np.arange(dt.size)
#nstep = 1000
#for j in i[0::nstep]:
#    m = np.mean(dt[j:j+nstep])
#    marr = np.append(marr, m)
#
#plt.plot(marr, 'o')
#plt.xlabel("Start Index (10^3)")
#plt.show()













#path = "/Users/Diana/Desktop/Astro 120/Lab1/"
#filename = "baes20_150902_1956_43.csv"
#
#x = np.loadtxt(path+filename, delimiter=',', dtype = np.int32)
#t = x[:,1]
#
## Compute the intervals
#dt = t[1:] - t[0:-1]
#
## Compute the means in chunks
#marr = np.array([])
#i = np.arange(dt.size)
#nstep = 1000
#for j in i[0::nstep]:
#    m = np.mean(dt[j:j+nstep])
#    marr = np.append(marr, m)
#
#plt.plot(marr, 'o')
#plt.xlabel("Start Index (10^3)")
#plt.show()




#mu = np.sum(marr)/np.float(marr.size)
#std = np.sqrt(np.sum((marr-mu)**2.)/(np.float(marr.size)-1.))
#
#print mu, std
#
#
#N = 500
## define the lower and upper bin edges and bin width
#bw = (dt.max()-dt.min())/(N-1.)
#bin1 = dt.min() + bw*np.arange(N)
#
## define the array to hold the occurrence count
#bincount = np.array([])
## loop through the bins
#for bin in bin1:
#    count = np.where((dt >= bin) & (dt < bin+bw))[0].size
#    bincount = np.append(bincount, count)
#    
## compute bin centers for plotting
#binc = bin1 + .5*bw
#plt.figure()
#plt.plot(binc, bincount, drawstyle='steps-mid')
#plt.show()