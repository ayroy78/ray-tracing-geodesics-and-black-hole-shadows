from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

filename = 'gc=0, q=0.txt'
filename2 = 'gc=0, q=0_analytical.txt'

##########################################################################################################
#comparison

data = pd.read_csv(filename, delimiter=',')
data2 = pd.read_csv(filename2, delimiter=',')

print("The horizon pixel count for lucas' data was equal to: "+ str(np.count_nonzero(data.iloc[:,2] == 0)))
print("The horizon pixel count for the analytical shadow was equal to: "+ str(np.count_nonzero(data2.iloc[:,2] == 0)))
########################################################################################################
#manipulation

# Create a new DataFrame by appending the original DataFrame to itself
data = data.append(data, ignore_index=True)

# Use NumPy to flip the y-values (2nd column) in the new DataFrame
data.iloc[31250:, 1] = -data.iloc[31250:, 1]

# Use NumPy to modify the colors (3rd column) in the new DataFrame
data.iloc[31250:, 2] = np.where(data.iloc[31250:, 2] == 1, 3, np.where(data.iloc[31250:, 2] == 2, 4, np.where(data.iloc[31250:, 2] == 3, 1, np.where(data.iloc[31250:, 2] == 4, 2,0))))   
    
#########################################################################################################
#error analysis
#r
print(np.min(data.iloc[:,3]), np.max(data.iloc[:,3]), np.mean(data.iloc[:,3]))
#th
print(np.min(data.iloc[:,6]), np.max(data.iloc[:,6]))
#initial H
print(np.min(data.iloc[:,4]), np.max(data.iloc[:,4]))
#final H
print(np.min(data.iloc[:,5]), np.max(data.iloc[:,5]))

#Number with H >1
print(np.count_nonzero(data.iloc[:,5]>=10**0))
########################################################################################################

# Colormap definition
colors = ['black', 'red', 'green', 'yellow', 'blue']
#colors = ['red', 'green', 'yellow', 'blue']
cmap_new = plt.cm.colors.ListedColormap(colors)
#data['flipped'] = -(data.iloc[:,3])
scatter = data.plot.scatter(x=0, y=1, c=2, cmap=cmap_new, s=7, figsize=(5,5))
# Get the colorbar object
colorbar = scatter.collections[0].colorbar


# Set the tick locations and labels
colorbar.set_ticks(range(5))
colorbar.set_ticklabels(['0', '1', '2', '3', '4'])
colorbar.remove()
plt.xlabel(r'$x=-rsin\beta$')
plt.ylabel(r'$y=rsin \alpha$')
plt.show()

###########################################################################################################
#Hamiltonian
data.plot.scatter(x=0,y=1, c=np.log10(data.iloc[:,5]), cmap='magma', s=7, figsize=(5,5))
plt.xlabel(r'$x=-rsin\beta$')
plt.ylabel(r'$y=rsin \alpha$')
plt.show()

plt.hist(np.log10(data.iloc[:,5]))
plt.xlabel(r'$log(\overline{H})$')
plt.ylabel('Counts')
plt.show()


