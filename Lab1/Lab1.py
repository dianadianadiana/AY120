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
data1= np.transpose(np.loadtxt(path+'baes_vary1_150914_1448_43.csv' ,delimiter=',',dtype='int32'))[1]
data2=np.transpose(np.loadtxt(path+'baes_vary2_150914_1451_43.csv' ,delimiter=',',dtype='int32'))[1]
data3=np.transpose(np.loadtxt(path+'baes_vary3_150914_1454_43.csv' ,delimiter=',',dtype='int32'))[1]
data4=np.transpose(np.loadtxt(path+'baes_vary4_150914_1457_43.csv' ,delimiter=',',dtype='int32'))[1]
data5=np.transpose(np.loadtxt(path+'baes_vary5_150914_1457_43.csv' ,delimiter=',',dtype='int32'))[1]
data6=np.transpose(np.loadtxt(path+'baes_vary6_150914_1459_43.csv' ,delimiter=',',dtype='int32'))[1]
datas = [data1,data2,data3,data4,data5,data6]

fig_path = '/Users/Diana/Desktop/Astro 120/Lab1/Figures/'
num_events_str = '_' + str(num_events)

def format_scinotation(z, decimal=3):
    zstr = str(z)
    e = 0
    while z>9:
        e +=1
        z/=10
    return str(np.int(z%(10**e)))+'.'+zstr[1:decimal]+'x'+r'$10^{e}$'.format(e=str(e))

def my_std(arr):
    """Calculates the STD of the arr"""
    mu = np.mean(arr)
    sq_sum = np.sum([(elem-mu)**2 for elem in arr]) #sums (x_i-mu)^2
    return np.sqrt(sq_sum/(len(arr)-1)) 
def my_sdom(arr):
    return my_std(arr)/(np.sqrt(len(arr)))

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
    """ Returns the figure 1 (Raw Event Time) in two graphs, top one showing the 
    lines connected, and bottom one showing the lines not connected """
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(t,'k', label = 'Points Joined')
    ax.set_title("Raw Event Time")
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
    m_arr = [marr1000,marr100,marr10]
    #means = [np.mean(marr) for marr in m_arr]

    fig = plt.figure()
    ax = fig.add_subplot(311)
    ax.plot(x_arr1000, marr1000, 'ko', label = 'Chunks of 1000 events, '+r'$\sigma$'+' = '+format_scinotation(my_std(marr1000))+'\nmean = '+str(np.mean(marr1000))+r'$\pm$'+"{:.3f}".format(my_sdom(marr1000)))
    #\nmean = ' + format_scinotation(np.mean(marr1000),5)+'
    ax.set_title('Mean Interval for Different Sized Chunks')
    ax.axhline(y = np.mean(marr1000), ls ='--', color = 'r')

    ax1 = fig.add_subplot(312)
    ax1.plot(x_arr100, marr100, 'ko', label = 'Chunks of 100 events, '+r'$\sigma$'+' = '+format_scinotation(my_std(marr100))+'\nmean = '+str(np.mean(marr100))+r'$\pm$'+"{:.3f}".format(my_sdom(marr100)))
    ax1.axhline(y = np.mean(marr100), ls ='--', color = 'r')
    ax1.set_ylabel('Mean Interval [clock ticks]', fontsize = 12)
   
    ax2 = fig.add_subplot(313)
    ax2.plot(x_arr10, marr10, 'ko', label = 'Chunks of 10 events, '+r'$\sigma$'+' = '+format_scinotation(my_std(marr10))+'\nmean = '+str(np.mean(marr10))+r'$\pm$'+"{:.3f}".format(my_sdom(marr10)))
    ax2.axhline(y = np.mean(marr10), ls ='--', color = 'r')
    ax2.set_xlabel("Start Index")
    
    ax_arr = [ax, ax1,ax2]
    for ax in ax_arr:
        ax.legend(loc='upper right', fancybox = True, fontsize = 12)
        ax.grid(True)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.tick_params(axis='y', labelsize=8)
        ax.tick_params(axis='x', labelsize=8)
    print [my_std(arr) for arr in m_arr]
    return fig
#plt.show(mean_interval_fig_lab())
#mean_interval_fig_lab().savefig(fig_path+'fig4&6_MeanOfDiffChunks'+num_events_str+'.png', dpi = 150)

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

def fig5():
    return mean_datastep_fig(100, dt)
    
def mean_datastep_fig_lab(dt):
    x_arr1000, marr1000 = mean_datastep(1000, dt)  
    x_arr100, marr100 = mean_datastep(100, dt)

    fig = plt.figure()
    ax = fig.add_subplot(211)
    #"{:.3f}".format(3.1415926))
    ax.plot(x_arr1000, marr1000, 'ko', label = 'Increasing by 1000 events\nmean = '+str(np.mean(marr1000))+r'$\pm$'+"{:.3f}".format(my_sdom(marr1000)))
    ax.set_title('Mean for Progressively Increasing Fractions of Data')
    ax.axhline(y = np.mean(marr1000), ls ='--', color = 'r')
    ax.set_ylabel('Mean Interval [clock ticks]')

    ax1 = fig.add_subplot(212)
    ax1.plot(x_arr100, marr100, 'ko', label = 'Increasing by 100 events\nmean = ' + str(np.mean(marr100))+r'$\pm$'+"{:.3f}".format(my_sdom(marr100)))
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
#plt.show(mean_datastep_fig_lab(dt))
#mean_datastep_fig_lab(dt).savefig(fig_path+'fig5_MeanOfDiffFractions'+num_events_str+'.png', dpi = 150)

################### Standard deviation part #############
marr = mean_interval(1000)[1]
mu = np.sum(marr)/np.float(marr.size)
std = np.sqrt(np.sum((marr-mu)**2.)/(np.float(marr.size)-1.)) #standard deviation
#print mu, std, marr[0]

########################## s/sqrtN  SDOM #########################
####### Figure 7: Variation of the standard deviation of the mean with the size of the data chunk averaged
lower, upper, nstep = 10, 500, 10
event_arr = np.arange(lower,upper+nstep,nstep)
sdom_arr = np.array([])
for nstep in event_arr:
    ith_mean_arr = np.array([])
    for i in np.arange(0,len(dt),nstep):
        if dt[i:i+nstep].size !=0: #ex [0:10][10:20]etc or [0:20][20:40]etc
            ith_mean = np.mean(dt[i:i+nstep]) #calculate the mean of that interval
            ith_mean_arr = np.append(ith_mean_arr, ith_mean) #add that mean to the array of means for that interval
    sdom_arr = np.append(sdom_arr,my_std(ith_mean_arr)) #add the std for the array of means
			
def fig7():
    fig=plt.figure()
    plt.plot(event_arr,sdom_arr,'ko')
    plt.title('SDOM for Different Sized Chunks')
    plt.xlabel('Number of Events Averaged')
    plt.ylabel('Standard Deviation of the Mean')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid(True)
    return fig
    
#plt.show(fig7())
#fig7().savefig(fig_path+'fig7_SDOMforDiffSizedChunks'+num_events_str+'.png', dpi = 150)

############# Figure 8: Standard deviation of the mean vs. 1/sqrt(N) showing linear behavior. The green
############# line is the theoretical expectation with SDOM = s/sqrt(N), where s is the sample standard deviation
def fig8():
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    s = my_std(dt)
    oneoverN = 1./np.sqrt(event_arr)
    ax1.plot(oneoverN,np.float(s)*oneoverN,'k',label='theorectical')
    ax1.plot(oneoverN,sdom_arr,'ko', label='experimental')
    ax1.set_title('Linear Relation\nbetween Standard Deviation of the Means and '+ r'$\left(\frac{1}{\sqrt{N}}\right)$')
    ax1.set_xlabel(r'$\left(\frac{1}{\sqrt{N}}\right)$') #1/sqrt(N)
    ax1.set_ylabel('Standard deviation of the mean [ticks]')
    ax1.set_xlim([np.amin(oneoverN), np.amax(oneoverN)])
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax1.legend(loc='lower right', fancybox = True, fontsize = 12)
    return fig
#plt.show(fig8())
#fig8().savefig(fig_path+'fig8_LinRelationSDOM'+num_events_str+'.png', dpi = 150)


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
def histogram_fig(N,dt, vertical_line=False):
    bw, binc, bincount = histogram(N,dt)
    fig = plt.figure()
    plt.plot(binc, bincount, 'k', drawstyle='steps-mid')
    plt.title("Histogram with "+str(N)+" bins")
    plt.xlabel('Interval [ticks]')
    plt.ylabel('Frequency')
    if vertical_line:
        plt.axvline(x=binc[np.argmax(bincount)],ls='--',color='r',label='Peak at '+str(np.round(binc[np.argmax(bincount)]))+' ticks')
        plt.legend()
    return fig


def poisson(bw, n, dt, tau):
    return bw*n*(1/tau)*np.exp(-dt/tau)  
def fig10(N, dt, vertical_line=False): #the zoomed in histogram
    return histogram_fig(N, dt[np.where(dt<4000)[0]], vertical_line)

#plt.show(histogram_fig(150,dt))
#histogram_fig(150,dt).savefig(fig_path+'fig9_HistogramWAfterpulse'+num_events_str+'.png', dpi = 150)
#plt.show(fig10(200,dt,True))
#fig10(200,dt,True).savefig(fig_path+'fig10_HistogramWAfterpulseZoomedIn'+num_events_str+'.png', dpi = 150)

def fig9_10(N, N1, dt):
    
    bw, binc, bincount = histogram(N,dt)
    bw1, binc1, bincount1 = histogram(N1, dt[np.where(dt<4000)[0]])
    fig = plt.figure()
    #fig.suptitle("Histogram")
    ax = fig.add_subplot(211)
    ax.plot(binc, bincount, 'k', drawstyle='steps-mid', label = str(N)+" bins")
    ax.set_title('Histogram')
    ax.set_ylabel('Frequency')
    ax.set_xlim([0,5.0e7])
    ax.legend(fancybox = True, fontsize = 10)
    ax.grid(True)
    
    ax1 = fig.add_subplot(212)
    ax1.plot(binc1, bincount1, 'k', drawstyle='steps-mid', label = str(N1)+" bins")
    ax1.set_xlabel('Interval [ticks]')
    ax1.set_ylabel('Frequency')
    #if vertical_line:
    ax1.axvline(x=binc1[np.argmax(bincount1)],ls='--',color='r',label='Peak at '+str(np.round(binc1[np.argmax(bincount1)]))+' ticks')
    ax1.legend(fancybox = True, fontsize = 10)
    ax1.grid(True)
    return fig
    
#plt.show(fig9_10(200,200,dt))
#fig9_10(200,200,dt).savefig(fig_path+'fig9_10_HistogramWAfterpulse'+num_events_str+'.png', dpi = 150)

def fig11a(N, dt):
    dt_no_afterpulse = dt[np.where(dt > 3000)[0]]
    tau = np.mean(dt_no_afterpulse)
    bw, binc, bincount = histogram(N, dt_no_afterpulse)
    fig = plt.figure()
    plt.plot(binc, bincount, c='k',drawstyle='steps-mid')
    y = poisson(bw, len(dt_no_afterpulse), dt_no_afterpulse,tau)
    plt.plot(dt_no_afterpulse, y,',',c='grey',label='Poisson fit',)
    plt.title('Histogram with no afterpulses, N = '+str(N))
    plt.xlabel('Interval [ticks]')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    plt.xlim([np.amin(dt_no_afterpulse), np.amax(dt_no_afterpulse)])
    plt.ylim(bottom=0)
    return fig
    
plt.show(fig11a(50,dt))

def fig11b(N, dt):
    dt_no_afterpulse = dt[np.where(dt > 3000)[0]]
    tau = np.mean(dt_no_afterpulse)
    bw, binc, bincount = histogram(N, dt_no_afterpulse)
    fig = plt.figure()
    plt.yscale('log')
    plt.plot(binc, bincount, c='k',drawstyle='steps-mid')
    y = poisson(bw, len(dt_no_afterpulse), dt_no_afterpulse,tau)
    plt.plot(dt_no_afterpulse, y,',', c='grey',label='Poisson fit')
    plt.title('Histogram with no afterpulses (log scale), N = '+str(N))
    plt.xlabel('Interval [ticks]')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    plt.xlim([np.amin(dt_no_afterpulse), np.amax(dt_no_afterpulse)])
    plt.ylim(bottom=0)
    return fig
plt.show(fig11b(50,dt))

def fig11(N, dt):
    dt_no_afterpulse = dt[np.where(dt > 3000)[0]]
    tau = np.mean(dt_no_afterpulse)
    bw, binc, bincount = histogram(N, dt_no_afterpulse)
    y = poisson(bw, len(dt_no_afterpulse), dt_no_afterpulse,tau)

    fig = plt.figure()
    fig.suptitle('Histogram with no afterpulses, N = '+str(N))

    ax = fig.add_subplot(121)
    ax1 = fig.add_subplot(122)
    for ax in [ax,ax1]:
        ax.plot(binc, bincount, c='k',drawstyle='steps-mid')
        ax.plot(dt_no_afterpulse, y,',',c='grey',label='Poisson fit',)
        ax.set_xlabel('Interval [ticks]')
        ax.set_ylabel('Frequency')
        ax.legend(loc='upper right')
    ax.set_xlim([np.amin(dt_no_afterpulse), np.amax(dt_no_afterpulse)])
    ax.set_title('Normal Histogram')
    ax1.set_title('Log Scale')
    ax.set_ylim(bottom=0)
    ax1.set_yscale('log')
    return fig

plt.show(fig11(50, dt))
#data1= np.transpose(np.loadtxt(path+'baes_vary1_150914_1448_43.csv' ,delimiter=',',dtype='int32'))[1]
#data2=np.transpose(np.loadtxt(path+'baes_vary2_150914_1451_43.csv' ,delimiter=',',dtype='int32'))[1]
#data3=np.transpose(np.loadtxt(path+'baes_vary3_150914_1454_43.csv' ,delimiter=',',dtype='int32'))[1]
#data4=np.transpose(np.loadtxt(path+'baes_vary4_150914_1457_43.csv' ,delimiter=',',dtype='int32'))[1]
#data5=np.transpose(np.loadtxt(path+'baes_vary5_150914_1457_43.csv' ,delimiter=',',dtype='int32'))[1]
#data6=np.transpose(np.loadtxt(path+'baes_vary6_150914_1459_43.csv' ,delimiter=',',dtype='int32'))[1]
#datas = [data1,data2,data3,data4,data5,data6]

###fig12
dt_arr = [t[1:] - t[0:-1] for t in datas]
dt_arr = [dt[np.where(dt>3000)[0]] for dt in dt_arr]

def fig12():
    means_arr = [np.mean(dt) for dt in dt_arr]
    std_arr = [my_std(dt) for dt in dt_arr]
    fig = plt.figure()
    plt.plot(means_arr, std_arr, 'ko',label='experimental')
    plt.plot([0,1.0e7],[0,1.0e7], c='grey', label='theoretical')
    fig.suptitle('hi')
    plt.title('Changing the LED Brightness')
    plt.xlabel('Interval sampe mean [ticks]')
    plt.ylabel('Interval Standard Deviation [ticks]')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.legend(loc='lower right', fancybox = True, fontsize = 12)
    return fig
plt.show(fig12())
    
cumsumdt_arr = [np.cumsum(dt) for dt in dt_arr]

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

#plt.show(fig2_3(dt))
#plt.show(fig1_1a(t))
#plt.show(mean_interval_fig_lab())
#plt.show(mean_datastep_fig_lab(dt))

#9/19/15



fig_path = '/Users/Diana/Desktop/Astro 120/Lab1/Figures/'
num_events_str = '_' + str(num_events)
#fig1_1a(t).savefig(fig_path+'fig1_2graphs'+num_events_str+'.png', dpi = 150)
#fig2_3(dt).savefig(fig_path+'fig2&fig3'+num_events_str+'.png', dpi = 150)
#mean_interval_fig_lab().savefig(fig_path+'mean_interval'+num_events_str+'.png', dpi = 150)
#mean_datastep_fig_lab(dt).savefig(fig_path+'mean_datastep'+num_events_str+'.png', dpi = 150)


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