import numpy as np
import scipy.spatial as spatial
import pandas as pd
import matplotlib as plt


x, y = np.mgrid[0:4, 0:4]
points = list(zip(x.ravel(), y.ravel()))
tree = spatial.KDTree(points)
print(tree.query_ball_point([2, 0], 1))

x = [1,2,3,4]
y = [1,1,1,1]
z = [1,1,2,3]
points2 = list(zip(x,y,z))
print(points2)
tree2 = spatial.KDTree(points2)
print(tree2.query_ball_point([0,2,0],2))