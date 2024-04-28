import matplotlib.pyplot as plt
import numpy as np

def ler_arquivo(nome_arquivo):
    try:
        with open(nome_arquivo, 'r') as arquivo:
            # Lendo o tamanho do primeiro vetor
            tamanho_vetor1 = int(arquivo.readline().strip())

            # Lendo o primeiro vetor
            vetor1 = list(map(float, arquivo.readline().strip().split(',')))

            # Lendo o segundo vetor
            vetor2 = list(map(float, arquivo.readline().strip().split(',')))

            vetor3 = list(map(float, arquivo.readline().strip().split(',')))

            return tamanho_vetor1, vetor1, vetor2, vetor3

    except FileNotFoundError:
        print("Arquivo não encontrado:", nome_arquivo)
        return None, None, None
    except Exception as e:
        print("Ocorreu um erro ao ler o arquivo:", e)
        return None, None, None

import subprocess

subprocess.run("./out.exe")

# Nome do arquivo a ser lido
nome_do_arquivo = "output.txt"

# Chamando a função para ler o arquivo
tamanho_vetor1, vetor1, vetor2, vetor3 = ler_arquivo(nome_do_arquivo)

# Verificando se a leitura foi bem sucedida
if tamanho_vetor1 is not None:
    print("Tamanho do vetor 1:", tamanho_vetor1)
    print("Vetor 1:", vetor1)
    print("Vetor 2:", vetor2)

plt.plot(vetor1, vetor2, label="Exata")
plt.plot(vetor1, vetor3, label="Aproximada")
plt.legend()
plt.show()


elementos = [4, 16, 64, 256, 1024]
erros4 = [0.000436801, 1.7329e-06, 6.77568e-09, 5.85645e-10, 4.66978e-11  ]
erros3 = [0.00667044, 0.00010631, 1.66317e-06, 2.5989e-08, 4.06081e-10]
erros2 = [0.0963726, 0.00624492, 0.000391206, 2.44539e-05, 1.52838e-06]

print((np.log(erros4[1]) - np.log(erros4[0]))/(np.log(16)- np.log(4)))
print((np.log(erros3[1]) - np.log(erros3[0]))/(np.log(16)- np.log(4)))
print((np.log(erros2[1]) - np.log(erros2[0]))/(np.log(16)- np.log(4)))


plt.plot(np.log(elementos), np.log(erros4), label="k=4")
plt.plot(np.log(elementos), np.log(erros3), label="k=3")
plt.plot(np.log(elementos), np.log(erros2), label ="k=2")
plt.legend()
plt.show()
