# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os


path = "/Users/Diana/Desktop/Astro 120/Lab1/data/"
file_arr = os.listdir(path) #import all the files into an array
print file_arr
filename = "baes20_150902_1956_43.csv"
print "chosen file", filename
num_events = 20000

def format_scinotation(z, decimal):
    zstr = str(z)
    e = 0
    while z>9:
        e +=1
        z/=10
    return str(np.int(z%(10**e)))+'.'+zstr[1:decimal]+'x10^'+str(e)

##################### First Question ########################
x = np.loadtxt(path+filename, delimiter=',', dtype = np.int32)
t = x[:,1]

#### Figure 1: Times Series #####
def fig1(t):
    fig=plt.figure()
    plt.plot(t,',')#,'marker',12)
    #plt.axhline(y=2**31,color='r', linewidth = .5)
    plt.title("Raw Event Time (Scattered)")
    plt.xlabel('Event Number')
    plt.ylabel('Clock Tick')
    return fig
    
def fig1a(t):
    fig=plt.figure()
    plt.plot(t)
    plt.title("Raw Event Time")
    plt.xlabel('Event Number')
    plt.ylabel('Clock Tick')
    return fig
    
def fig1_1a(t):
    fig=plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(t,'k', label = 'Points Joined')
    ax.set_title("Raw Event Time")
    #ax.set_xlabel('Event Number')
    ax.set_ylabel('Clock Tick')
    ax.legend(loc='upper right', frameon = True, fancybox = True, fontsize = 10)

    
    ax1 = fig.add_subplot(212)
    ax1.plot(t,'k,', label = 'Scattered')
    ax1.set_xlabel('Event Number')
    ax1.set_ylabel('Clock Tick')
    ax1.legend(loc='upper right', frameon = True, fancybox = True,fontsize = 10)
    return fig

# Compute the intervals
dt = t[1:] - t[0:-1]

###### Figure 2: Interval subsequent events (in clock ticks) as a function of event number ####
def fig2(dt):
    fig=plt.figure()
    plt.plot(dt, 'k')
    plt.title('Interval between Subsequent Events')
    plt.xlabel('Event Number')
    plt.ylabel('Interval [clock ticks]')
    return fig
   
###### Figure 3: Interval subsequent events (in clock ticks) as a function of event number - single dot####
def fig3(dt):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(dt,'k,')
    ax1.set_title('Interval between Subsequent Event')
    ax1.set_xlabel('Event number')
    ax1.set_ylabel('Interval [clock ticks]')
    return fig    

def fig2_3(dt):
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(dt, 'k',linewidth=.5, label = 'Points Joined')
    ax.set_title('Interval between Subsequent Events')
    #ax.set_xlabel('Event number', size = 12)
    ax.set_ylabel('Interval [clock ticks]', size = 12)
    ax.legend(loc='upper right', fancybox = True)
    
    ax1 = fig.add_subplot(212)
    ax1.plot(dt,'k,', label='Scattered')
    ax1.set_xlabel('Event number')
    ax1.set_ylabel('Interval [clock ticks]')
    ax1.legend(loc='upper right', fancybox = True)
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

def mean_interval_fig_lab():
    x_arr1000, marr1000 = mean_interval(1000)  
    x_arr100, marr100 = mean_interval(100)
    x_arr10, marr10 = mean_interval(10)
    
    #x_arr = [x_arr1000,x_arr100,x_arr10]
    #m_arr = [marr1000,marr100,marr10]
    #means = [np.mean(marr) for marr in m_arr]

    fig = plt.figure()
    ax = fig.add_subplot(311)
    ax.plot(x_arr1000, marr1000, 'ko', label = 'Chunks of 1000 events\nmean = ' + format_scinotation(np.mean(marr1000),5))
    ax.set_title('Mean Interval for Different Sized Chunks')
    ax.axhline(y = np.mean(marr1000), ls ='--', color = 'r')

    ax1 = fig.add_subplot(312)
    ax1.plot(x_arr100, marr100, 'ko', label = 'Chunks of 100 events\nmean = ' + format_scinotation(np.mean(marr100),5))
    ax1.axhline(y = np.mean(marr100), ls ='--', color = 'r')
    ax1.set_ylabel('Mean Interval [clock ticks]', fontsize = 12)
   
    ax2 = fig.add_subplot(313)
    ax2.plot(x_arr10, marr10, 'ko', label = 'Chunks of 10 events\nmean = ' + format_scinotation(np.mean(marr10),5))
    ax2.axhline(y = np.mean(marr10), ls ='--', color = 'r')
    ax2.set_xlabel("Start Index")
    
    ax_arr = [ax, ax1,ax2]
    for ax in ax_arr:
        ax.legend(loc='upper right', fancybox = True, fontsize = 8)
        ax.grid(True)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.tick_params(axis='y', labelsize=8)
        ax.tick_params(axis='x', labelsize=8)

    return fig
    
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
    return mean_datastep_fig(100, dt)
    
def mean_datastep_fig_lab(dt):
    x_arr1000, marr1000 = mean_datastep(1000, dt)  
    x_arr100, marr100 = mean_datastep(100, dt)

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(x_arr1000, marr1000, 'ko', label = 'Intervals of 1000 events\nmean = ' + format_scinotation(np.mean(marr1000),5))
    ax.set_title('Mean Interval for Different Sized Chunks')
    ax.axhline(y = np.mean(marr1000), ls ='--', color = 'r')
    ax.set_ylabel('Mean Interval [clock ticks]')

    ax1 = fig.add_subplot(212)
    ax1.plot(x_arr100, marr100, 'ko', label = 'Intervals of 100 events\nmean = ' + format_scinotation(np.mean(marr100),5))
    ax1.axhline(y = np.mean(marr100), ls ='--', color = 'r')
    ax1.set_ylabel('Mean Interval [clock ticks]')
    ax1.set_xlabel('Index Count')

    ax_arr = [ax, ax1]
    for ax in ax_arr:
        ax.legend(loc='lower right', fancybox = True, fontsize = 8)
        ax.grid(True)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.tick_params(axis='y', labelsize=8)
        ax.tick_params(axis='x', labelsize=8)
        #ax.xticklabel = np.arange(0,20000,2500)

    return fig

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
    ax1.set_xlim([np.amin(x_arr), np.amax(x_arr)])
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
    
################################################################################
###############################----Figures-----#################################
################################################################################
#plt.show(fig10(100, dt))
#plt.show(fig11a(100, dt))
#plt.show(fig11b(100,dt))

#plt.show(fig1(t))
#plt.show(fig1a(t))
#plt.show(fig2(dt))
#plt.show(fig3(dt))

plt.show(fig2_3(dt))
plt.show(fig1_1a(t))
plt.show(mean_interval_fig_lab())
plt.show(mean_datastep_fig_lab(dt))


fig_path = '/Users/Diana/Desktop/Astro 120/Lab1/Figures/'
num_events_str = '_' + str(num_events)
#fig1_1a(t).savefig(fig_path+'fig1_2graphs'+num_events_str+'.png', dpi = 150)
#fig2_3(dt).savefig(fig_path+'fig2&fig3'+num_events_str+'.png', dpi = 150)
mean_interval_fig_lab().savefig(fig_path+'mean_interval'+num_events_str+'.png', dpi = 150)
mean_datastep_fig_lab(dt).savefig(fig_path+'mean_datastep'+num_events_str+'.png', dpi = 150)


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