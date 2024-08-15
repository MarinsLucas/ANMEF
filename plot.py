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

            #u
            vetor3 = list(map(float, arquivo.readline().strip().split(',')))

            return tamanho_vetor1, vetor1, vetor2, vetor3

           
    except FileNotFoundError:
        print("Arquivo não encontrado:", nome_arquivo)
        return None, None, None, None
    except Exception as e:
        print("Ocorreu um erro ao ler o arquivo:", e)
        return None, None, None, None

import subprocess

# Nome do arquivo a ser lido
nome_do_arquivo = "output.txt"

# Chamando a função para ler o arquivo
tamanho_vetor1, vetor1, vetor2, vetor3 = ler_arquivo(nome_do_arquivo)

plt.title("Comparação entre solução exata e aproximação")
plt.plot(vetor1, vetor3,'-o', label="aproximada")
plt.plot(vetor1, vetor2,'-o', label="exata", linestyle='--')

plt.grid()
plt.legend()
plt.show() 

#SIPG
nels = [8,16,32,64,128,256,512,1024]
""" erros2 = [0.00899566,0.00226545,0.000567528,0.000141959,3.54947e-05,8.87398e-06,2.21852e-06,5.54641e-07]
erros3 = [0.000186388,2.28792e-05,2.83317e-06,3.52474e-07,4.3955e-08,5.48795e-09,6.90571e-10,3.39351e-10]
erros4 = [4.46067e-06,2.79508e-07,1.74816e-08,1.09281e-09,6.86668e-11,2.89563e-11,1.34383e-10,4.49918e-10]
erros5 = [7.89157e-08,2.44767e-09,7.6221e-11,1.24112e-11,4.74672e-11,1.79354e-10,7.74392e-10,0]
print("u = 1: " + str((np.log(erros2[1]) - np.log(erros2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 2: " + str((np.log(erros3[1]) - np.log(erros3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 3: " + str((np.log(erros4[1]) - np.log(erros4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 4: " + str((np.log(erros5[1]) - np.log(erros5[0]))/(np.log(nels[1])- np.log(nels[0]))))

plt.title("Gráfico de convergência SIPG")
plt.plot(np.log(nels), np.log(erros2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(erros3),'-o', label="k=2")
plt.plot(np.log(nels[:-3]), np.log(erros4[:-3]),'-o', label="k=3")
plt.plot(np.log(nels[:-4]), np.log(erros5[:-4]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show() 


#Zero zero DG

erros2 =[0.0103201,0.00235094,0.000572926,0.000142298,3.55159e-05,8.8753e-06,2.2186e-06]
erros3 =[0.000819767,0.000129918,1.78763e-05,2.32834e-06,2.96552e-07,3.74011e-08,4.69549e-09]
erros4 = [4.64899e-06,2.82704e-07,1.75335e-08,1.09363e-09,6.83166e-11,4.28688e-12,1.94952e-12]
erros5 = [1.68843e-07,5.90148e-09,1.93395e-10,6.18461e-12,8.87581e-13,3.47602e-12,1.68548e-11]

print("u = 1: " + str((np.log(erros2[1]) - np.log(erros2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 2: " + str((np.log(erros3[1]) - np.log(erros3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 3: " + str((np.log(erros4[1]) - np.log(erros4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 4: " + str((np.log(erros5[1]) - np.log(erros5[0]))/(np.log(nels[1])- np.log(nels[0]))))

plt.title("Gráfico de convergência 00-DG")
plt.plot(np.log(nels), np.log(erros2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(erros3),'-o', label="k=2")
plt.plot(np.log(nels[:-3]), np.log(erros4[:-3]),'-o', label="k=3")
plt.plot(np.log(nels[:-4]), np.log(erros5[:-4]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()

#NIPG
erros2 =[0.0210547,0.006062,0.00161376,0.000415228,0.000105206,2.64676e-05,6.63678e-06,1.66159e-06]
erros3 =[0.00931808,0.00333415,0.000980207,0.000264308,6.85249e-05,1.74392e-05,4.39838e-06, 1.10443e-06]
erros4 = [1.71446e-05,1.1695e-06,7.64893e-08,4.8925e-09,3.09598e-10,2.06963e-11, 8.19958e-12]
erros5 = [2.20165e-06,1.72443e-07,1.19059e-08,7.78511e-10,4.6467e-11,9.11939e-12, 5.65376e-11]


print("u = 1: " + str((np.log(erros2[1]) - np.log(erros2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 2: " + str((np.log(erros3[1]) - np.log(erros3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 3: " + str((np.log(erros4[1]) - np.log(erros4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 4: " + str((np.log(erros5[1]) - np.log(erros5[0]))/(np.log(nels[1])- np.log(nels[0]))))

plt.title("Gráfico de convergência NIPG")
plt.plot(np.log(nels), np.log(erros2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(erros3),'-o', label="k=2")
plt.plot(np.log(nels[:-3]), np.log(erros4[:-2]),'-o', label="k=3")
plt.plot(np.log(nels[:-4]), np.log(erros5[:-3]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()"""


#NIPG0
erros2 = [0.027006,0.00679947,0.00170284,0.000425894,0.000106485,2.6622e-05,6.65555e-06,1.66389e-06]
erros3 = [0.0101949,0.00352115,0.00100949,0.000268355,6.90551e-05,1.7507e-05,4.40695e-06,1.1055e-06]
erros4 = [1.7518e-05,1.18158e-06,7.68761e-08,4.90468e-09,3.09767e-10,1.98743e-11,4.42772e-12,1.12924e-11]
erros5 = [2.24265e-06,1.74162e-07,1.19668e-08,7.80538e-10,4.68424e-11,1.12165e-11,5.81505e-11,2.45245e-10]

print("u = 1: " + str((np.log(erros2[1]) - np.log(erros2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 2: " + str((np.log(erros3[1]) - np.log(erros3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 3: " + str((np.log(erros4[1]) - np.log(erros4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 4: " + str((np.log(erros5[1]) - np.log(erros5[0]))/(np.log(nels[1])- np.log(nels[0]))))

plt.title("Gráfico de convergência NIPG 0")
plt.plot(np.log(nels), np.log(erros2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(erros3),'-o', label="k=2")
plt.plot(np.log(nels[:-3]), np.log(erros4[:-3]),'-o', label="k=3")
plt.plot(np.log(nels[:-4]), np.log(erros5[:-4]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()