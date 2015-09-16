# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os

path = "/Users/Diana/Desktop/Astro 120/Lab1/data/"
file_arr = os.listdir(path) #import all the files into an array
print file_arr
filename = file_arr[3] #"baes20_150902_1956_43.csv"
print "chosen file", filename
num_events = 20000

##################### First Question ########################
x = np.loadtxt(path+filename, delimiter=',', dtype = np.int32)
t = x[:,1]

#### Figure 1: Times Series #####
def fig1():
    fig=plt.figure()
    plt.plot(t,'o')
    plt.xlabel('Event Number')
    plt.ylabel('Clock Tick')
    return fig

# Compute the intervals
dt = t[1:] - t[0:-1]

###### Figure 2: Interval subsequent events (in clock ticks) as a function of event number ####
def fig2():
    fig=plt.figure()
    plt.plot(dt)
    plt.xlabel('Event Number')
    plt.ylabel('Interval [clock ticks]')
    return fig
   
###### Figure 3: Interval subsequent events (in clock ticks) as a function of event number - single dot####
def fig3():
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(dt,',')
    ax1.set_xlabel('Event number')
    ax1.set_ylabel('Interval [clock ticks]')
    return fig    
       
##################### Second Question ######################
# Compute the means in chunks
def mean_interval(nstep):
    marr = np.array([])
    i = np.arange(dt.size)
    for j in i[0::nstep]:
        #print j, j+nstep-1,np.mean(dt[j:j+nstep])
        m = np.mean(dt[j:j+nstep])
        marr = np.append(marr, m)
    x_arr = np.arange(len(marr))*nstep
    return x_arr, marr
    
def mean_interval_fig(nstep):
    x_arr, marr = mean_interval(nstep)    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(x_arr, marr, 'o')
    ax1.set_title('Chunks of ' + str(nstep) + ' events; mean = ' + str(np.mean(marr)) + ' ticks')
    ax1.set_xlabel("Start Index")
    ax1.set_ylabel('Mean Interval [clock ticks]')
    return fig

######### Figure 4: The mean interval fromm Figure 3 computed for ten chunks of 1000 events.#####
def fig4():
    return mean_interval_fig(1000)    

###### Figure 6: Interval subsequent events (in clock ticks) as a function of event number ####
def fig6():
    return mean_interval_fig(100)
    
#################### Figure 5 ########################
def mean_datastep(datastep, dt):
    mean_arr = np.array([])
    count = 0
    while count < dt.size:
        mean = np.mean(dt[0:count+datastep])
        #print count, count+datastep, mean
        mean_arr = np.append(mean_arr, mean)
        count += datastep
    x_arr = np.arange(1,len(mean_arr)+1)*datastep
    return x_arr, mean_arr
    
def mean_datastep_fig(datastep, dt):
    x_arr, mean_arr = mean_datastep(datastep, dt)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(x_arr, mean_arr, 'o')
    ax1.set_title('Mean computed for Increasing Intervals of ' + str(datastep))
    ax1.set_xlabel('Number of intervals averaged')
    ax1.set_ylabel('Mean Interval [clock ticks]')
    return fig

#datastep = 100 
def fig5():
    return mean_datastep_fig(100)


################### Standard deviation part #############
marr = mean_interval(1000)[1]
mu = np.sum(marr)/np.float(marr.size)
std = np.sqrt(np.sum((marr-mu)**2.)/(np.float(marr.size)-1.)) #standard deviation
print mu, std, marr[0]

########################## s/sqrtN #########################
#
#count = 10
#rang = np.arange(10,500,10)
#stddev_arr= np.array([])
#for j in rang:
#	ith_mean_arr=np.array([])
#	for i in np.arange(j,10000,j):
#		ith_mean=np.mean(dt[0:i])
#		ith_mean_arr = np.append(ith_mean_arr,ith_mean)
#	mu_i = np.sum(ith_mean_arr)/np.float(ith_mean_arr.size) 
#	stddev_i = np.sqrt(np.sum((ith_mean_arr-mu_i)**2.)/(np.float(ith_mean_arr.size)-1.))
#	stddev_arr = np.append(stddev_arr,stddev_i)

####### Figure 7: Variation of the standard deviation of the mean with the size of the data chunk averaged
increment = 10
std_means = np.array([])
while increment < 500:
	ith_mean_arr = np.array([])
	dt_arr = np.arange(dt.size)
	nstep = increment
	for j in dt_arr[::nstep]:
		if dt[j:j+nstep].size != 0:
			m =np.mean(dt[j:j+nstep])
		ith_mean_arr = np.append(ith_mean_arr,m)
	std_m = np.std(ith_mean_arr)
	std_means = np.append(std_means,std_m)
	increment+=10
fig=plt.figure()
plt.plot(np.arange(10,500,10),std_means,'ko')
plt.xlabel('Number of Events Averaged')
plt.ylabel('Standard Deviation of the Mean')
plt.close(fig)
#plt.show()
#plt.savefig('/Users/Diana/Desktop/Astro 120/Lab1/fig7_20000.png')


############# Figure 8: Standard deviation of the mean vs. 1/sqrt(N) showing linear behavior. The green
############# line is the theoretical expectation with SDOM = s/sqrt(N), where s is the sample standard deviation
def fig8(dt):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    x_arr = 1./np.sqrt(np.arange(10,500,10))
    s = np.std(dt)
    ax1.plot(x_arr,np.float(s)*x_arr)
    ax1.plot(x_arr,std_means,'ro')
    ax1.set_xlabel('1/sqrt(N)')
    # r'$\frac{1}{\sqrt{N}}$'
    ax1.set_ylabel('Standard deviation of the mean [ticks]')
    ax1.set_xlim([np.amin(x_arr8), np.amax(x_arr8)])
    #fig4.show()
    #fig4.savefig('/Users/Diana/Desktop/Astro 120/Lab1/fig8_20000.png')


######################### Histogram #########################
def histogram(N, dt):
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
    return bw, binc, bincount
def histogram_fig(N,dt):
    bw, binc, bincount = histogram(N,dt)
    fig = plt.figure()
    plt.plot(binc, bincount, drawstyle='steps-mid')
    plt.title("Histogram with "+str(N)+" bins")
    plt.xlabel('Interval [ticks]')
    plt.ylabel('Frequency')
    return fig
    #plt.savefig('/Users/Diana/Desktop/astro temp/historgram_' + str(N)+'.png')
#plt.savefig(figpath + "")

#plt.show(histogram_fig(500,dt))

def poisson(bw, n, dt, tau):
    return bw*n*(1/tau)*np.exp(-dt/tau)  
def fig10(N, dt): #the zoomed in histogram
    return histogram_fig(N, dt[np.where(dt<4000)[0]])

def fig11a(N, dt):
    dt_no_afterpulse = dt[np.where(dt > 3000)[0]]
    tau = np.mean(dt_no_afterpulse)
    bw, binc, bincount = histogram(N, dt_no_afterpulse)
    fig = plt.figure()
    plt.plot(binc, bincount, drawstyle='steps-mid')
    y = poisson(bw, len(dt_no_afterpulse), dt_no_afterpulse,tau)
    plt.plot(dt_no_afterpulse, y,',',label='Poisson fit',)
    plt.title('Histogram with no afterpulses, N = '+str(N))
    plt.xlabel('Interval [ticks]')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    return fig

def fig11b(N, dt):
    dt_no_afterpulse = dt[np.where(dt > 3000)[0]]
    tau = np.mean(dt_no_afterpulse)
    bw, binc, bincount = histogram(N, dt_no_afterpulse)
    fig = plt.figure()
    plt.yscale('log')
    plt.plot(binc, bincount, drawstyle='steps-mid')
    y = poisson(bw, len(dt_no_afterpulse), dt_no_afterpulse,tau)
    plt.plot(dt_no_afterpulse, y,',', label='Poisson fit')
    plt.title('Histogram with no afterpulses (log scale), N = '+str(N))
    plt.xlabel('Interval [ticks]')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    return fig
    
plt.show(fig10(100, dt))
plt.show(fig11a(100, dt))
plt.show(fig11b(100,dt))


fig_path = '/Users/Diana/Desktop/Astro 120/Lab1/Figures/'
num_events_str = '_' + str(num_events)
fig_arrx = []#[fig1(),fig2()]
fig_arry = []#['fig1','fig2']
fig_arr = [fig_arrx, fig_arry]
for figure, figure_name in np.transpose(fig_arr):
    print figure_name
    plt.show(figure)
    figure.savefig(fig_path+figure_name+num_events_str+'.png')

########################### Other #######################
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