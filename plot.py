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


elementos = [8, 16, 32, 64]

#nint = 2
erros2p = [0.00450433, 0.00113316, 0.000283793, 7.09813e-05]
erros2u = [0.0167836, 0.00388342, 0.000931678, 0.000227991]

#nint = 3
erros3p = [0.000172815, 2.37805e-05, 3.1012e-06, 3.95218e-07]
erros3u = [0.00850185, 0.00267676, 0.00073731, 0.000192219]
#nint = 4
erros4p = [2.64426e-06, 1.6477e-07, 1.02891e-08, 6.42919e-10]
erros4u = [0.000140355, 1.71928e-05, 2.13724e-06, 2.66771e-07]

#nint = 5
erros5p = [7.24192e-08, 2.42906e-09, 7.83637e-11, 1.46856e-11]
erros5u = [8.06921e-06, 5.8378e-07, 3.88357e-08, 2.50884e-09]

print("p 2" + str((np.log(erros5p[1]) - np.log(erros5p[0]))/(np.log(elementos[1])- np.log(elementos[0]))))
print("u 2" + str((np.log(erros5u[1]) - np.log(erros5u[0]))/(np.log(elementos[1])- np.log(elementos[0]))))

plt.title("Erros da solução de P")
plt.plot(np.log(elementos), np.log(erros2p),'-o', label="k=1")
plt.plot(np.log(elementos), np.log(erros3p), '-o',label="k=2")
plt.plot(np.log(elementos), np.log(erros4p), '-o',label="k=3")
plt.plot(np.log(elementos), np.log(erros5p),'-o', label="k=4")
plt.legend()
plt.show()

plt.title("Erros da solução de U")
plt.plot(np.log(elementos), np.log(erros2u),'-o', label="k=1")
plt.plot(np.log(elementos), np.log(erros3u), '-o',label="k=2")
plt.plot(np.log(elementos), np.log(erros4u), '-o',label="k=3")
plt.plot(np.log(elementos), np.log(erros5u),'-o', label="k=4")
plt.legend()
plt.show()
