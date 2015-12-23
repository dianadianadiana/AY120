import numpy as np


slope =  doppm = 0.00340463
delt = 136.4
calibm = -0.0218834

c = 6.60651408e+02
b = -2.03936794e-02
a = -1.32851462e-06

linear_m = -0.0218834 
linear_yint = 660.8859559

xi11 = np.radians(3.37) #from JPL 11
eta11 = np.radians(226.33)


def transform(xi, eta, doppm, delt, calibm):    
    cosfactor = np.cos(eta)*np.cos(xi)
    c = 299792.458 #km/s
    lambd = 524*linear_m + linear_yint
    
    vrot = (doppm*c*calibm*delt)/(lambd*cosfactor)
    print 'vrot: ', vrot
    Trot = (24.47*24*3600)
      
    radius = (1/(2*np.pi))*vrot*Trot
    # 696,300 should be in km
    
    print 'sun radius: ', radius
    #AU = 1.496e8 #km
    angsun = (1919.3/3600)     #angular diameter of sun in degrees (/3600)
    angsun = np.radians(angsun)
    my_AU = (2*radius)/np.sin(angsun)
    print 'AU: ', my_AU
    return my_AU
    
transform(xi11,eta11, .002, delt, calibm)

def transform1(zi, eta, doppm, time, calibm):
	corrfactor=np.cos(eta)*np.cos(zi)
	c = 299792.548
	lambd = 524*linear_m + linear_yint
	vrot = (doppm*c*calibm*time)/(lambd*corrfactor)
	print 'vrot: ',vrot
	T_rot = (26.24*24*3600)
	radius = (1/(2*np.pi))*vrot*T_rot
	print 'solar radius = ', radius
	angsun = np.radians(1919.3/3600.)
	AU = (2*radius)/np.sin(angsun)
	print 'earth sun dist= ', AU
	return vrot,radius, AU
zi = np.radians(3.16)
eta = np.radians(202.16)
ans = transform1(xi11, eta11, .0012, delt/2, calibm)


