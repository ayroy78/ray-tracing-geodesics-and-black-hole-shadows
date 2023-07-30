from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

filename = 'schwarz_shadow_data_spline.txt'
filename2 = 'schwarz_shadow_data.txt'

##########################################################################################################
#comparison

data = pd.read_csv(filename, delimiter=',')
data2 = pd.read_csv(filename2, delimiter=',')
spline_colours = data.iloc[:,2]
analytical_colours = data2.iloc[:,2]

colour_diff = spline_colours - analytical_colours
print("The number of deviations was equal to: "+ str(np.count_nonzero(colour_diff != 0)))

#r
print(np.min(data.iloc[:,3]), np.max(data.iloc[:,3]), np.mean(data.iloc[:,3]))
#theta
print(np.min(data.iloc[:,4]), np.max(data.iloc[:,4]), np.mean(data.iloc[:,4]))
#hamiltonian
print(np.min(data.iloc[:,5]), np.max(data.iloc[:,5]))
#carter constant
print(np.min(data.iloc[:,6]), np.max(data.iloc[:,6]))
################################################################################################################
#shadow
fig, ax = plt.subplots(figsize=(5, 5), gridspec_kw={'left': 0.148})

# Colormap definition
colors = ['black', 'red', 'green', 'yellow', 'blue']
#colors = ['red', 'green', 'yellow', 'blue']
cmap_new = plt.cm.colors.ListedColormap(colors)
#data['flipped'] = -(data.iloc[:,3])
scatter = data.plot.scatter(x=0, y=1, c=2, cmap=cmap_new, s=7, figsize=(5,5), ax =ax)
# Get the colorbar object
colorbar = scatter.collections[0].colorbar

# Set the tick locations and labels
#colorbar.set_ticks(range(5))
#colorbar.set_ticklabels(['0', '1', '2', '3', '4'])
colorbar.remove()
plt.xlabel(r'$x=-rsin\beta$')
plt.ylabel(r'$y=rsin \alpha$')

plt.show()

################################################################################################################
#hamiltonian
fig, ax = plt.subplots(figsize=(5, 5), gridspec_kw={'left': 0.148})

data.plot.scatter(x=0, y=1, c=np.log10(data.iloc[:,5]), cmap ='magma', s=7, figsize=(5,5), ax=ax)

plt.xlabel(r'$x=-rsin\beta$')
plt.ylabel(r'$y=rsin \alpha$')

plt.show()
################################################################################################################
#carter constant
fig, ax = plt.subplots(figsize=(5, 5), gridspec_kw={'left': 0.148})

data.plot.scatter(x=0, y=1, c=np.log10(data.iloc[:,6]), cmap ='magma', s=7, figsize=(5,5), ax=ax)

plt.xlabel(r'$x=-rsin\beta$')
plt.ylabel(r'$y=rsin \alpha$')

plt.show()
