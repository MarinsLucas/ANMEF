import matplotlib.pyplot as plt
import numpy as np

x = [0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1]
solexata = [1,0.92388,0.707107,0.382683,6.12323e-17,-0.382683,-0.707107,-0.92388,-1]


#alpha = -1
aproximada0 = [9.87711,10.9877,8.40258,9.19188,5.64883,5.99664,2.03509,1.88847,-1.88847,-2.50721,-5.52452,-6.52118,-8.31952,-9.54237,-9.84795,-11.1108] #beta = 1 
aproximada = [1.03325,1.07581,0.947376,0.831211,0.717268,0.460068,0.377962,0.0188847,-0.0188847,-0.425174,-0.412857,-0.804504,-0.743975,-1.06136,-0.96183,-1.15662] #beta = 10 
aproximada2=[0.994781,0.930528,0.918985,0.712273,0.703283,0.385582,0.380512,0.000188847,-0.000188847,-0.385233,-0.380861,-0.712006,-0.70355,-0.930383,-0.91913,-1.00712] #beta = 100
plt.plot(x, solexata, label="Solução Exata")
for i in range(len(x)-1):
    plt.plot(x[i:i+2], aproximada[i*2:(2*i)+2], '-o', color='b')
    plt.plot(x[i:i+2], aproximada2[i*2:(2*i)+2], '-o', color='r')
    plt.plot(x[i:i+2], aproximada0[i*2:(2*i)+2], '-o', color='g')



plt.show()

#alpha = 0
aproximada0 = [0.38315,1.49377,0.353984,1.14329,0.270928,0.618742,0.146625,9.90035e-17,2.34612e-17,-0.618742,-0.146625,-1.14329,-0.270928,-1.49377,-0.353984,-1.61685] #beta = 1 
aproximada = [0.938315,0.980869,0.86689,0.750725,0.663489,0.406289,0.359078,6.50095e-17,5.74552e-17,-0.406289,-0.359078,-0.750725,-0.663489,-0.980869,-0.86689,-1.06169] #beta = 10 
aproximada2=[0.993831,0.929578,0.918181,0.711469,0.702745,0.385044,0.380323,6.16101e-17,6.08546e-17,-0.385044,-0.380323,-0.711469,-0.702745,-0.929578,-0.918181,-1.00617] #beta = 100
plt.plot(x, solexata, label="Solução Exata")
for i in range(len(x)-1):
    plt.plot(x[i:i+2], aproximada[i*2:(2*i)+2], '-o', color='b')
    plt.plot(x[i:i+2], aproximada2[i*2:(2*i)+2], '-o', color='r')
    plt.plot(x[i:i+2], aproximada0[i*2:(2*i)+2], '-o', color='g')

plt.show()

#alpha = 1
aproximada0 = [-9.11082,-8.00019,-7.69461,-6.90531,-5.10697,-4.75916,-1.74184,-1.88847,1.88847,1.26973,5.23127,4.23461,7.77767,6.55482,9.13998,7.87711] #beta = 1 
aproximada = [0.843375,0.885929,0.786404,0.670239,0.60971,0.35251,0.340193,-0.0188847,0.0188847,-0.387405,-0.305299,-0.696946,-0.583003,-0.900383,-0.77195,-0.966745] #beta = 10 
aproximada2=[0.992882,0.928629,0.917376,0.710664,0.702207,0.384506,0.380134,-0.000188847,0.000188847,-0.384855,-0.379785,-0.710931,-0.70194,-0.928774,-0.917231,-1.00522] #beta = 100
plt.plot(x, solexata, label="Solução Exata")
for i in range(len(x)-1):
    plt.plot(x[i:i+2], aproximada[i*2:(2*i)+2], '-o', color='b')
    plt.plot(x[i:i+2], aproximada2[i*2:(2*i)+2], '-o', color='r')
    plt.plot(x[i:i+2], aproximada0[i*2:(2*i)+2], '-o', color='g')

plt.show()