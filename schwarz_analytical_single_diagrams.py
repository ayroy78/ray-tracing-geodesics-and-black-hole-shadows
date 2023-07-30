from matplotlib import pyplot as plt
import numpy as np

#massive particle ray tracing #############################################################################################

filename = 'particle.txt'


ray_data = np.loadtxt(filename, delimiter=',')
plt.title('Single Particle Orbit')
plt.plot(ray_data[:,0], ray_data[:,1], c='red', linestyle='dashed', label='particle')
phi = np.linspace(0,2*np.pi,100)
r = 2.0
x = r*np.cos(phi)
y =r*np.sin(phi)
plt.plot(x,y, c='black', label='event horizon')
r=3.0
x = r*np.cos(phi)
y =r*np.sin(phi)
plt.plot(x,y, c='blue', label='photonsphere')
plt.xlim([-np.max(ray_data[:,0])-2,np.max(ray_data[:,0])+2])
plt.ylim([-np.max(ray_data[:,1])-2,np.max(ray_data[:,1])+2])
plt.scatter(ray_data[-1,0], ray_data[-1,1], c='red')
plt.legend(loc=3)
plt.show()

########################################################################################################################

#photon ray tracing#######################################################################################################

filename = 'photon.txt'

ray_data = np.loadtxt(filename, delimiter=',')
plt.title('Single Photon Orbit')
plt.plot(ray_data[:,0], ray_data[:,1], linestyle='dashed', label='photon')
plt.scatter(ray_data[-1,0], ray_data[-1,1])
phi = np.linspace(0,2*np.pi,100)
r = 2.0
x = r*np.cos(phi)
y =r*np.sin(phi)
plt.plot(x,y, c='black', label='event horizon')
plt.xlim([-np.max(ray_data[:,0])-1,np.max(ray_data[:,0])+1])
plt.ylim([-np.max(ray_data[:,1])-1,np.max(ray_data[:,1])+1])
plt.legend()
plt.show()


##########################################################################################################################




