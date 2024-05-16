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

            #(p)
            vetor3 = list(map(float, arquivo.readline().strip().split(',')))
            
            #(U)
            vetor4 = list(map(float, arquivo.readline().strip().split(',')))
            return tamanho_vetor1, vetor1, vetor2, vetor3, vetor4

    except FileNotFoundError:
        print("Arquivo não encontrado:", nome_arquivo)
        return None, None, None
    except Exception as e:
        print("Ocorreu um erro ao ler o arquivo:", e)
        return None, None, None

import subprocess

#subprocess.run("./out.exe")

# Nome do arquivo a ser lido
nome_do_arquivo = "output.txt"

# Chamando a função para ler o arquivo
tamanho_vetor1, vetor1, vetor2, vetor3, vetor4 = ler_arquivo(nome_do_arquivo)

# Verificando se a leitura foi bem sucedida
if tamanho_vetor1 is not None:
    print("Tamanho do vetor 1:", tamanho_vetor1)
    print("Vetor 1:", vetor1)
    print("Vetor 2:", vetor2)

plt.plot(vetor1, vetor2, label="Exata u")
plt.plot(vetor1, np.sin([x * np.pi for x in vetor1]), label = "Exata p")
plt.plot(vetor1, vetor3, label="u")
plt.plot(vetor1, vetor4, label="p")
plt.grid()
plt.legend()
plt.show()


elementos = [4, 16, 64, 256, 1024]
#Erru u:
#erros5 = [0.0011364, 2.37805e-05, 3.95218e-07, 6.26002e-09, 8.38519e-11]
#erros4 = [4.26415e-05, 1.6477e-07, 6.42919e-10, 4.93914e-12, 8.38519e-11]
#erros3 = [0.0011364,2.37805e-05, 3.95218e-07, 6.26002e-09, 1.04265e-10]
#erros2 = [0.0177512, 0.00113316, 7.09813e-05, 4.437e-06, 2.77316e-07]

erros4 = [0.00119622, 1.71928e-05, 2.66771e-07, 4.16664e-09, 4.31678e-10]
erros3 = [0.0215468, 0.00267676, 0.000192219, 1.23584e-05, 7.77729e-07]
erros2 = [0.077882, 0.00388342, 0.000227991, 1.40171e-05, 8.72417e-07]

print((np.log(erros4[1]) - np.log(erros4[0]))/(np.log(16)- np.log(4)))
print((np.log(erros3[1]) - np.log(erros3[0]))/(np.log(16)- np.log(4)))
print((np.log(erros2[1]) - np.log(erros2[0]))/(np.log(16)- np.log(4)))


plt.plot(np.log(elementos), np.log(erros4), label="k=4")
plt.plot(np.log(elementos), np.log(erros3), label="k=3")
plt.plot(np.log(elementos), np.log(erros2), label ="k=2")
plt.legend()
plt.show()
