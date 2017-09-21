import numpy as np
import scipy.spatial as spatial
from scipy.spatial.distance import euclidean
import pandas as pd
import matplotlib.pyplot as plt

#universal variables
data = 'DM_2p8M.dat'
output_coords = 'coords_and_flags9.csv'
output_halos = 'halos9.csv'
cube_length = 141.3
particle_mass = 1.4*10**10

#get data
x,y,z = np.loadtxt(data,unpack=True)
xyz=list(zip(x,y,z))
numData = len(xyz)

#create dataframe to hold final values
#and to hold initial points and eventually flags
haloDF = pd.DataFrame(columns=['Mass', 'Center of Mass', 'RMS Radius', 'Potential Issue'], index=np.arange(1,len(xyz)))

d = {'x':x, 'y':y, 'z':z}
haloCoords = pd.DataFrame(d)

#r_link is the mean interparticle separation times .2
#which is the number density
r_link = .2*(cube_length/np.power(numData,1/3.))

#create flag array, tree, and halo count for flag
flag = np.linspace(0,0,len(x))
potential_particles = 0
tree = spatial.cKDTree(xyz)
halo_count = 0

for i in range(0,numData):
    if(i % 1000 == 0):
        print(i)
    #make sure the particle hasn't been used
    if (flag[i] == 0):
        point_list = [i]
        coords = [xyz[i]]
        halo_count = halo_count + 1
        flag[i] = halo_count

        #iterate through all points in point list
        for p in point_list:
            neighbors = tree.query_ball_point(xyz[p],r_link)
            #make sure there are some neighbors
            if neighbors:
                for j in neighbors:
                    if x[j] - r_link < 0 or x[j] + r_link > cube_length \
                            or y[j] - r_link < 0 or y[j] + r_link > cube_length \
                            or z[j] - r_link < 0 or z[j] + r_link > cube_length:
                        haloDF['Potential Issue'][halo_count] = 'yes'
                        potential_particles+=1
                    #make sure the particle hasn't been used already
                    if (flag[j] == 0):
                        point_list.append(j)
                        coords.append(xyz[j])
                        flag[j] = halo_count

        total_points = len(coords)
        center_of_mass = tuple(map(lambda y: sum(y) / float(len(y)), zip(*coords)))

        #get distance from center for rms
        dist_from_center = []
        for p in coords:
            dist_from_center.append(euclidean(p,center_of_mass))

        rms_radius = np.sqrt(np.mean(np.square(dist_from_center)))

        #store data for each halo
        haloDF['Mass'][halo_count] = 1.4*10**10*total_points
        haloDF['Center of Mass'][halo_count] = center_of_mass
        haloDF['RMS Radius'][halo_count] = rms_radius

haloDF = haloDF[haloDF.Mass > 1.4*10**11]
haloDF = haloDF.reset_index()


haloDF.to_csv(output_halos)
haloCoords['flag'] = flag
haloCoords.to_csv(output_coords)

print('done')
print(halo_count)
print(potential_particles)

f, plts = plt.subplots(2,1,figsize=(15,15))

#get number of halos with mass above a certain threshold

x = np.linspace(12,15,30)

count_mass = np.linspace(0,0,30)

i=0

for mass in x:
    for row in haloDF['Mass']:
        if np.log10(row) >= mass:
            count_mass[i]+=1
    i+=1

print(count_mass)

num_dens_halos = [np.log10(x/(np.power(cube_length, 3)*.1)) for x in count_mass]
log_rms = np.log10(haloDF['RMS Radius'].astype(float))
log_mass = np.log10(haloDF['Mass'].astype(float))

plts[0].plot(x,num_dens_halos)
plts[0].set_ylabel('log(dN/Vdm)')
plts[0].set_xlabel('log(M)')
plts[0].set_title('Number Density of Halos vs. Mass')

plts[1].scatter(log_rms,log_mass)
plts[1].set_ylabel('log(rms radius)')
plts[1].set_xlabel('log(M)')
plts[1].set_title('RMS Radius vs. Mass')


plt.show()

