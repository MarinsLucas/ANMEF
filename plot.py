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

            #u
            vetor4 = list(map(float, arquivo.readline().strip().split(',')))
            return tamanho_vetor1, vetor1, vetor2, vetor3, vetor4

           
    except FileNotFoundError:
        print("Arquivo não encontrado:", nome_arquivo)
        return None, None, None, None, None
    except Exception as e:
        print("Ocorreu um erro ao ler o arquivo:", e)
        return None, None, None, None, None

import subprocess

# Nome do arquivo a ser lido
nome_do_arquivo = "output.txt"

# Chamando a função para ler o arquivo
tamanho_vetor1, vetor1, vetor2, vetor3, vetor4 = ler_arquivo(nome_do_arquivo)

# Verificando se a leitura foi bem sucedida
if tamanho_vetor1 is not None:
    print("Tamanho do vetor 1:", tamanho_vetor1)
    print("Vetor 1:", vetor1)
    print("Vetor 2:", vetor2)

plt.title("Comparação entre solução exata e aproximação")
plt.plot(vetor1, vetor3,'-o', label="aproximada")
plt.plot(vetor1, vetor2,'-o', label="exata", linestyle='--')
plt.plot(vetor1, vetor4, '-o', label="p", linestyle='--')

plt.grid()
plt.legend()
plt.show()


def solexata(vetor1):
    vetor1 = np.pi*vetor1
    return np.cos(vetor1)

plt.plot(np.linspace(0, 1, 13), [  1,   0.965926,   0.866025,   0.707107,        0.5,   0.258819,-3.8382e-15,  -0.258819,       -0.5,  -0.707107,  -0.866025,  -0.965926,         -1], '-o')
plt.plot(np.linspace(0, 1, 13), solexata(np.linspace(0, 1, 13)))
plt.show()


