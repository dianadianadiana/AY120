import numpy as np


sun_period = 24.47 # days sidereal
def get_sun_radius(vel, period):
    period = period * 24 * 3600 # get period from days to seconds
    return (1./(2*np.pi)) * vel * period
c = 299792 #km/s

slope = 0.00340463
delt = 134.999969482


c = 6.60651408e+02
b = -2.03936794e-02
a = -1.32851462e-06
xi11 = np.radians(3.37) #from JPL 11
eta11 = np.radians(226.33)
xi12 = np.radians(3.25)
eta12 = np.radians(213.15)
xi16 = np.radians(2.79)
eta16 = np.radians(160.42)
def transform(xi, eta, doppm, delt):
    wlsol = doppm**2*a + doppm*b + c
    print wlsol
    cosfactor = np.cos(eta)*np.cos(xi)
    print 'delt', delt
    Trot = (24.47*24*3600)
    print Trot
    
    vrot = (wlsol*delt)/cosfactor
    print 'vrot', vrot
    radius = (1/(2*np.pi))*vrot*Trot

    print 'sun radius: ', radius
    return radius
    
#radius111 = transform(xi11, eta11, doppm111, time111[-1], time111[0])
#radius112 = transform(xi12, eta12, doppm112, time112[-1], time112[0])

#radius = transform(xi11, eta11, slope, delt/2)

def get_vrot(beta, delt, nm_pix, xi, eta):
    lamb = 524**2*a+ 524*b + c # the center wavelength
    print lamb
    cosfactor = np.cos(eta)*np.cos(xi)
    return (beta*c*nm_pix*delt)/(lamb *cosfactor)
    
nm_pix = -2.18834018e-02
 
v_rot = get_vrot(slope, delt, nm_pix, xi11, eta11)

print v_rot
radius = get_sun_radius(v_rot, sun_period)
print radius


ang_dia_sun = 1919.3 #arc seconds
theta_sun = ang_dia_sun * (1./3600)
theta_sun = np.radians(theta_sun)
theta_sun = .533 # degrees angular diameter
def get_earthsun_dist(theta, radius):
    return 2*radius/np.sin(theta)
m_d_11_11 = get_earthsun_dist(theta_sun, radius)
print m_d_11_11
def dist(vals):
    x,y,z = vals
    return np.sqrt(x**2 + y**2 + z**2)
    
coords11_11 = [-6.619813967311905E-01, -7.364541638390490E-01, -1.468250727594430E-06]
coords11_12 = [-6.488009866335709E-01, -7.477827710351609E-01, -5.680459615135181E-07]

AU = 149597870.700 # km

real_d_11_11 = dist(coords11_11) * AU #km
my_d_11_11 = 124974985.27485408 #km

def accurate_error(calc, real):
    return (calc-real)/real * 100 # percent
    
print accurate_error(m_d_11_11, real_d_11_11)
