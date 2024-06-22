import matplotlib.pyplot as plt
import numpy as np

def ler_arquivo(nome_arquivo):
    try:
        with open(nome_arquivo, 'r') as arquivo:
            # Lendo o tamanho do primeiro vetor
            tamanho_vetor1 = int(arquivo.readline().strip())

            # Lendo o primeiro vetor (xs)
            vetor1 = list(map(float, arquivo.readline().strip().split(',')))

            # Lendo o segundo veto (solução exata)
            vetor2 = list(map(float, arquivo.readline().strip().split(',')))

            #(u)
            vetor3 = list(map(float, arquivo.readline().strip().split(',')))
            return tamanho_vetor1, vetor1, vetor2, vetor3
                   
    except FileNotFoundError:
        print("Arquivo não encontrado:", nome_arquivo)
        return None, None, None, None
    except Exception as e:
        print("Ocorreu um erro ao ler o arquivo:", e)
        return None, None, None,


""" x = [0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1]
solexata = [1,0.92388,0.707107,0.382683,6.12323e-17,-0.382683,-0.707107,-0.92388,-1]

#Questão C e D da lista 3

#alpha = -1
aproximada0 = [9.87711,10.9877,8.40258,9.19188,5.64883,5.99664,2.03509,1.88847,-1.88847,-2.50721,-5.52452,-6.52118,-8.31952,-9.54237,-9.84795,-11.1108] #beta = 1 
aproximada = [1.03325,1.07581,0.947376,0.831211,0.717268,0.460068,0.377962,0.0188847,-0.0188847,-0.425174,-0.412857,-0.804504,-0.743975,-1.06136,-0.96183,-1.15662] #beta = 10 
aproximada2=[0.994781,0.930528,0.918985,0.712273,0.703283,0.385582,0.380512,0.000188847,-0.000188847,-0.385233,-0.380861,-0.712006,-0.70355,-0.930383,-0.91913,-1.00712] #beta = 100
plt.plot(x, solexata, label="Solução Exata")

for i in range(len(x)-1):
    plot, = plt.plot(x[i:i+2], aproximada0[i*2:(2*i)+2], '-o', color='b')
    plot1, =plt.plot(x[i:i+2], aproximada[i*2:(2*i)+2], '-o', color='r')
    plot2, = plt.plot(x[i:i+2], aproximada2[i*2:(2*i)+2], '-o', color='g') 

plot.set_label("B=1")
plot1.set_label("B=10")
plot2.set_label("B=100")
plt.legend()
plt.title("a=-1")
plt.show() 

#alpha = 0
aproximada0 = [0.38315,1.49377,0.353984,1.14329,0.270928,0.618742,0.146625,9.90035e-17,2.34612e-17,-0.618742,-0.146625,-1.14329,-0.270928,-1.49377,-0.353984,-1.61685] #beta = 1 
aproximada = [0.938315,0.980869,0.86689,0.750725,0.663489,0.406289,0.359078,6.50095e-17,5.74552e-17,-0.406289,-0.359078,-0.750725,-0.663489,-0.980869,-0.86689,-1.06169] #beta = 10 
aproximada2=[0.993831,0.929578,0.918181,0.711469,0.702745,0.385044,0.380323,6.16101e-17,6.08546e-17,-0.385044,-0.380323,-0.711469,-0.702745,-0.929578,-0.918181,-1.00617] #beta = 100
plt.plot(x, solexata, label="Solução Exata")
for i in range(len(x)-1):
    plot, = plt.plot(x[i:i+2], aproximada0[i*2:(2*i)+2], '-o', color='b')
    plot1, =plt.plot(x[i:i+2], aproximada[i*2:(2*i)+2], '-o', color='r')
    plot2, = plt.plot(x[i:i+2], aproximada2[i*2:(2*i)+2], '-o', color='g')  

plot.set_label("B=1")
plot1.set_label("B=100")
plot2.set_label("B=100") 
plt.legend()
plt.title("a=0")
plt.show() 

#alpha = 1
aproximada0 = [-9.11082,-8.00019,-7.69461,-6.90531,-5.10697,-4.75916,-1.74184,-1.88847,1.88847,1.26973,5.23127,4.23461,7.77767,6.55482,9.13998,7.87711] #beta = 1 
aproximada = [0.843375,0.885929,0.786404,0.670239,0.60971,0.35251,0.340193,-0.0188847,0.0188847,-0.387405,-0.305299,-0.696946,-0.583003,-0.900383,-0.77195,-0.966745] #beta = 10 
aproximada2=[0.992882,0.928629,0.917376,0.710664,0.702207,0.384506,0.380134,-0.000188847,0.000188847,-0.384855,-0.379785,-0.710931,-0.70194,-0.928774,-0.917231,-1.00522] #beta = 100
plt.plot(x, solexata, label="Solução Exata")
for i in range(len(x)-1):
    plot, = plt.plot(x[i:i+2], aproximada0[i*2:(2*i)+2], '-o', color='b')
    plot1, =plt.plot(x[i:i+2], aproximada[i*2:(2*i)+2], '-o', color='r')
    plot2, = plt.plot(x[i:i+2], aproximada2[i*2:(2*i)+2], '-o', color='g')  
plot.set_label("B=1")
plot1.set_label("B=100")
plot2.set_label("B=100") 
plt.legend()
plt.title("a=1")
plt.show()



#Questão E e F
tam, xs, exata, aprox = ler_arquivo("output.txt")
plt.plot(xs, exata, label="Solução Exata")
plt.legend()
plt.plot(xs, aprox, '-o', label="Aproximada")

plt.show()  """

#Erros questão F lista 3
nels = [8,16,32,64,128,256,512]

#!A = -1
nint2 = [ 0.0085603,  0.00214637, 0.000536986 , 0.000134271, 3.35693e-05, 8.39244e-06, 2.09814e-06]
nint3 = [ 0.000201812, 2.527e-05, 3.16012e-06, 3.95057e-07, 4.93835e-08, 6.17318e-09, 7.83276e-10]
nint4 = [ 4.41474e-06,2.76298e-07,1.72745e-08, 1.07976e-09, 7.04876e-11, 8.14829e-11, 3.2539e-10]
nint5 = [ 8.16613e-08, 2.55482e-09, 7.98983e-11, 1.01174e-11, 3.92363e-11, 1.56918e-10, 6.27715e-10]

print("\na = -1")
print("k = 1: " + str((np.log(nint2[1]) - np.log(nint2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 2: " + str((np.log(nint3[1]) - np.log(nint3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 3: " + str((np.log(nint4[1]) - np.log(nint4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 4: " + str((np.log(nint5[1]) - np.log(nint5[0]))/(np.log(nels[1])- np.log(nels[0]))))


plt.title("a=-1")
plt.plot(np.log(nels), np.log(nint2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(nint3),'-o', label="k=2")
plt.plot(np.log(nels[:-2]), np.log(nint4[:-2]),'-o', label="k=3")
plt.plot(np.log(nels[:-3]), np.log(nint5[:-3]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()

#!A = 0
nint2 = [0.0110966,0.00378296,0.00159389,0.000755097,0.000372131,0.000185382,9.26055e-05]
nint3 = [0.000222713,2.78845e-05,3.48699e-06,4.35918e-07,5.44911e-08,6.81142e-09,8.51452e-10]
nint4 = [4.81569e-06,3.52578e-07,3.17802e-08,3.48375e-09,4.20338e-10,1.53007e-10,5.75738e-10]
nint5 = [8.5846e-08,2.68561e-09,8.40551e-11,1.71474e-11,6.77903e-11,2.71105e-10,1.08416e-09]

print("\na = 0")
print("k = 1: " + str((np.log(nint2[1]) - np.log(nint2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 2: " + str((np.log(nint3[1]) - np.log(nint3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 3: " + str((np.log(nint4[1]) - np.log(nint4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 4: " + str((np.log(nint5[1]) - np.log(nint5[0]))/(np.log(nels[1])- np.log(nels[0]))))


plt.title("a=0")
plt.plot(np.log(nels), np.log(nint2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(nint3),'-o', label="k=2")
plt.plot(np.log(nels[:-2]), np.log(nint4[:-2]),'-o', label="k=3")
plt.plot(np.log(nels[:-3]), np.log(nint5[:-3]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()

#!A=1
nint2 = [0.0154369,0.00642262,0.00302667,0.00148932,0.000741629,0.000370435,0.00018517]
nint3 = [0.000298184,4.77118e-05,9.53091e-06,2.20748e-06,5.40354e-07,1.34336e-07,3.34446e-08]
nint4 = [5.72515e-06,5.14477e-07,5.62907e-08,6.7625e-09,8.36545e-10,1.15894e-10,2.04317e-10]
nint5 = [9.99792e-08,3.56363e-09,1.54265e-10,7.9343e-12,8.11686e-13,5.16625e-12,2.08165e-11]

print("\na = 1")
print("k = 1: " + str((np.log(nint2[1]) - np.log(nint2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 2: " + str((np.log(nint3[1]) - np.log(nint3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 3: " + str((np.log(nint4[1]) - np.log(nint4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 4: " + str((np.log(nint5[1]) - np.log(nint5[0]))/(np.log(nels[1])- np.log(nels[0]))))


plt.title("a=1")
plt.plot(np.log(nels), np.log(nint2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(nint3),'-o', label="k=2")
plt.plot(np.log(nels[:-1]), np.log(nint4[:-1]),'-o', label="k=3")
plt.plot(np.log(nels[:-2]), np.log(nint5[:-2]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()


#!Galerkin Clássico
nint5 = [8.24125e-08,2.57833e-09,8.05958e-11,2.51988e-12,2.8097e-13,4.78381e-13,3.11229e-12]
nint4 = [4.4682e-06,2.79648e-07,1.7484e-08,1.09284e-09,6.8316e-11,6.58986e-12,2.15421e-11]
nint3 = [0.000205541,2.57381e-05,3.21869e-06,4.0238e-07,5.0299e-08,6.28743e-09,7.86338e-10]
nint2 = [0.00904308,0.00226901,0.000567769,0.000141975,3.54957e-05,8.87404e-06,2.21852e-06]


print("\nGalerkin Clássico")
print("k = 1: " + str((np.log(nint2[1]) - np.log(nint2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 2: " + str((np.log(nint3[1]) - np.log(nint3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 3: " + str((np.log(nint4[1]) - np.log(nint4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("k = 4: " + str((np.log(nint5[1]) - np.log(nint5[0]))/(np.log(nels[1])- np.log(nels[0]))))


plt.title("Galerkin Clássico")
plt.plot(np.log(nels), np.log(nint2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(nint3),'-o', label="k=2")
plt.plot(np.log(nels[:-1]), np.log(nint4[:-1]),'-o', label="k=3")
plt.plot(np.log(nels[:-2]), np.log(nint5[:-2]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()