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

#Erros questão F lista 3
def solexatap(v):
    x = np.array(v)
    return np.cos(x * np.pi)

nels = [64,128,32,64,128,256,512]
nome_do_arquivo = "output.txt"
tamanho_vetor1, vetor1, vetor2, vetor3, vetor4 = ler_arquivo(nome_do_arquivo)

plt.plot(np.linspace(0, 1, 100),np.pi*np.sin(np.pi*np.linspace(0, 1, 100)), label="Solução Exata U")
plt.plot(np.linspace(0, 1, 100),solexatap(np.linspace(0, 1, 100)), label="Solução Exata P")
for i in range(0,len(vetor1)-1, 2):
    plot, = plt.plot(vetor1[i:i+2], vetor3[i:i+2], '-o', color='b')
    plot1, =plt.plot(vetor1[i:i+2], vetor4[i:i+2], '-o', color='r')
plot.set_label("Aproximada U")
plot1.set_label("Aproximada P")  
plt.legend()
plt.title("nel=8")
plt.show()

nels = [8,16,32,64,128,256,512,1024]
#u
#nint = 2
errosu2 = [0.263398, 0.131911, 0.0659819, 0.0329943, 0.0164975, 0.00824882, 0.00412442, 0.00206221]
errosu3 = [0.0133564, 0.00334394, 0.000836286, 0.00020909, 5.22738e-05, 1.30685e-05, 3.26713e-06, 8.16784e-07]
errosu4 = [0.00044341, 5.5495e-05, 6.93902e-06, 8.67445e-07, 1.08433e-07, 1.35542e-08, 1.69589e-09, 4.62811e-10]
errosu5 = [ 1.09714e-05, 6.86448e-07, 4.29145e-08, 2.68234e-09, 1.67727e-10, 2.29348e-11, 8.7639e-11, 2.33266e-10]

#p
#nint = 2
errosp2 = [0.0279973, 0.0134039, 0.00662498, 0.00330278, 0.00165017, 0.000824935, 0.000412448, 0.000206222]
errosp3 = [0.00135291, 0.00033548, 8.36966e-05, 2.09133e-05,5.22764e-06,  1.30687e-06, 3.26714e-07, 8.16785e-08]
errosp4 = [4.45877e-05, 5.55723e-06, 6.94144e-07, 8.6752e-08, 1.08435e-08, 1.35546e-09, 1.71028e-10, 1.30128e-10]
errosp5 = [1.10054e-06, 6.8698e-08, 4.29228e-09, 2.68247e-10, 1.6831e-11, 6.19553e-12, 2.75745e-11,7.42807e-11]

print("u = 1: " + str((np.log(errosu2[1]) - np.log(errosu2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 2: " + str((np.log(errosu3[1]) - np.log(errosu3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 3: " + str((np.log(errosu4[1]) - np.log(errosu4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 4: " + str((np.log(errosu5[1]) - np.log(errosu5[0]))/(np.log(nels[1])- np.log(nels[0]))))


plt.title("Gráfico de convergência de U")
plt.plot(np.log(nels), np.log(errosu2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(errosu3),'-o', label="k=2")
plt.plot(np.log(nels[:-2]), np.log(errosu4[:-2]),'-o', label="k=3")
plt.plot(np.log(nels[:-3]), np.log(errosu5[:-3]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()


print("p = 1: " + str((np.log(errosp2[1]) - np.log(errosp2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("p = 2: " + str((np.log(errosp3[1]) - np.log(errosp3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("p = 3: " + str((np.log(errosp4[1]) - np.log(errosp4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("p = 4: " + str((np.log(errosp5[1]) - np.log(errosp5[0]))/(np.log(nels[1])- np.log(nels[0]))))

plt.title("Gráfico de convergência de P")
plt.plot(np.log(nels), np.log(errosp2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(errosp3),'-o', label="k=2")
plt.plot(np.log(nels[:-2]), np.log(errosp4[:-2]),'-o', label="k=3")
plt.plot(np.log(nels[:-3]), np.log(errosp5[:-3]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()