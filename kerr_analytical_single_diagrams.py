from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

###################################################################################################

#photon orbits
##########################################################
#quantities to change
mass = 1.0
angmom=0.5
E = 1.0/3.0
L =  1.0
Q = 1.157
filename = 'photon.txt'
######################################################

r_event_max = mass + np.sqrt(mass**2.0 - angmom**2.0)
sigma = lambda r,th,a: r**2.0 + (np.cos(th)*a)**2.0
delta = lambda M,r,a: r**2.0 + a**2.0 - 2*M*r

g_tt = lambda M,r,th,a: -1*(1- (2*M*r)/(sigma(r,th,a)))
g_tphi = lambda M,r,th,a: -(2*M*r*a*np.sin(th)**2.0)/(sigma(r,th,a))
g_rr = lambda M,r,th,a : sigma(r,th,a)/delta(M,r,a)
g_thth = lambda r,th,a: sigma(r,th,a)
g_phiphi = lambda M,r,th,a : (a**2.0 + r**2.0 + (2*M*r*(a*np.sin(th))**2.0)/(sigma(r,th,a)))*np.sin(th)**2.0


data = np.loadtxt(filename, delimiter=',')

V_eff = lambda E,L,M,r,th,a : (g_rr(M,r,th,a)**(-1.0))*( (g_phiphi(M,r,th,a)*E**2.0 + g_tt(M,r,th,a)*L**2.0 +2*E*L*g_tphi(M,r,th,a))/(-g_tphi(M,r,th,a)**2.0 + g_tt(M,r,th,a)*g_phiphi(M,r,th,a)))

test_r, test_th = np.meshgrid(np.linspace(r_event_max+0.1, 10, 100), np.linspace(0.1, np.pi-0.1, 100))
test_V_eff = V_eff(E,L,mass,test_r,test_th,angmom)

fig = plt.figure()
fig.suptitle(r'Single Photon Orbit and $V_{eff}$', fontsize=20)
x = data[:,0]
y = data[:,1]
z = data[:,2]
l = len(data)
spacetime=plt.subplot(1,2,1, projection='3d')
spacetime.scatter(0,0,0,c='black', s =500)
spacetime.scatter(x[l-1], y[l-1], z[l-1], c ='yellow', s =50)
spacetime.plot(x[:l],y[:l], z[:l], c ='red', linestyle='dashed', linewidth=3)
spacetime.set_axis_off()
spacetime.view_init(17, 45)
spacetime.set_xlim([np.min(x),np.max(x)])
spacetime.set_ylim([np.min(y),np.max(y)])
spacetime.set_zlim([np.min(z),np.max(z)])
spacetime.text2D(0.05, 0.90, "a:" + str("%.3f"%(angmom)), size=12, c='blue',transform=spacetime.transAxes)
spacetime.text2D(0.05, 0.825, "M:" + str("%.3f"%(mass)), size=12, c='blue',transform=spacetime.transAxes)
spacetime.text2D(0.05, 0.300, "Q:" + str("%.3f"%(Q)), size=12, c='blue',transform=spacetime.transAxes)
spacetime.text2D(0.05, 0.225, "E:" + str("%.3f"%(E)), size=12, c='blue',transform=spacetime.transAxes)
spacetime.text2D(0.05, 0.150, "L:" + str("%.3f"%(L)), size=12, c='blue',transform=spacetime.transAxes)

pot = plt.subplot(1,2,2)


pot.pcolormesh(test_r,test_th,test_V_eff, vmax=0.00, vmin=None)
pot.axvline(x=2.7, linestyle='dashed', c='black')
#fig.colorbar(pot)
pot.set_xlabel('r')
pot.set_ylabel(r'$\theta$')


plt.show()
###########################################################################################################################




#particle orbits
##########################################################
#quantities to change
mass = 1.0
angmom = 0.5
E = 1.0
L =  -2.0
Q = 4.0
filename = 'particle.txt'
######################################################

data = np.loadtxt(filename, delimiter=',')
x = data[:,0]
y = data[:,1]
z = data[:,2]
l = len(data)


fig = plt.figure()
spacetime=plt.subplot(1,1,1, projection='3d')
spacetime.set_title('Single Particle Orbit', fontsize=20)
spacetime.scatter(0,0,0,c='black', s =500)
spacetime.scatter(x[l-1], y[l-1], z[l-1], c ='yellow', s =50)
spacetime.plot(x[:l],y[:l], z[:l], c ='red', linestyle='dashed', linewidth=3)
spacetime.set_axis_off()
spacetime.view_init(17, 45)
spacetime.set_xlim([np.min(x),np.max(x)])
spacetime.set_ylim([np.min(y),np.max(y)])
spacetime.set_zlim([np.min(z),np.max(z)])
spacetime.text2D(0.05, 0.90, "a:" + str("%.3f"%(angmom)), size=12, c='blue',transform=spacetime.transAxes)
spacetime.text2D(0.05, 0.825, "M:" + str("%.3f"%(mass)), size=12, c='blue',transform=spacetime.transAxes)
spacetime.text2D(0.05, 0.300, "Q:" + str("%.3f"%(Q)), size=12, c='blue',transform=spacetime.transAxes)
spacetime.text2D(0.05, 0.225, "E:" + str("%.3f"%(E)), size=12, c='blue',transform=spacetime.transAxes)
spacetime.text2D(0.05, 0.150, "L:" + str("%.3f"%(L)), size=12, c='blue',transform=spacetime.transAxes)

plt.show()

##################################################################################################################
