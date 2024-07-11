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
#Questão F
""" #u
#nint = 2
errosu2 = [0.269971, 0.12825, 0.0632512, 0.0315154, 0.0157439, 0.0078702, 0.00393489, 0.001967424]
errosu3 = [0.0136897,0.00325112,0.000801676,0.000199718,4.98857e-05,1.24687e-05,3.117e-06,7.79239e-07]
errosu4 = [0.000454474,5.39547e-05,6.65185e-06,8.28563e-07,1.03479e-07,1.2932e-08,1.61648e-09,2.21117e-10]
errosu5 = [1.12452e-05,6.67395e-07,4.11385e-08,2.56211e-09,1.5999e-10,1.00052e-11,8.09577e-13,1.33423e-12]

#p
#nint = 2
errosp2 = [0.0351204,0.00834264,0.00205729,0.000512534,0.000128021,3.19983e-05,7.99914e-06,1.99976e-06]
errosp3 = [0.0017254,0.000204881,2.52603e-05,3.1465e-06,3.92966e-07,4.911e-08,6.13841e-09,7.6729e-10]
errosp4 = [5.70118e-05,3.38419e-06,2.08611e-07,1.29925e-08,8.11312e-10,5.07049e-11,4.77384e-12,1.41159e-11]
errosp5 = [1.40843e-06,4.17949e-08,1.28813e-09,4.01123e-11,1.25354e-12,1.32545e-13,1.82901e-13,4.84166e-13]

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
plt.show() """


#Questão G
""" 
#erros de U
errosu2 = [0.0284883, 0.00713322, 0.00178401,0.000446046, 0.000111514, 2.78787e-05, 6.96969e-06, 1.74242e-06]
errosu3 = [0.00291447, 0.000715269,  0.000177969, 4.44391e-05, 1.11064e-05, 2.7764e-06, 6.94094e-07, 1.73548e-07]
errosu4 = [0.000137, 1.70801e-05, 2.1336e-06, 2.66655e-07, 3.33306e-08, 4.1664e-09, 5.21042e-10, 6.64808e-11]
errosu5 = [4.1449e-06, 2.58985e-07, 1.61854e-08, 1.01139e-09, 6.27256e-11, 2.12921e-11, 8.57913e-11, 3.4245e-10]

#erros de P
errosp2 = [0.00458334, 0.00113838, 0.000284127, 7.10024e-05, 1.77488e-05, 4.43708e-06, 1.10926e-06, 2.77315e-07]
errosp3 = [0.000112135, 1.40396e-05, 1.75567e-06, 2.19481e-07, 2.74358e-08, 3.4295e-09, 4.287e-10, 5.46456e-11]
errosp4 = [2.6277e-06, 1.64488e-07, 1.02845e-08, 6.42846e-10,4.0179e-11, 2.51363e-12, 4.68165e-13, 1.75709e-12]
errosp5 = [5.15014e-08, 1.61141e-09, 5.0372e-11, 1.58848e-12,8.44115e-13, 3.35101e-12, 1.33904e-11, 5.34526e-11]


print("u = 1: " + str((np.log(errosu2[1]) - np.log(errosu2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 2: " + str((np.log(errosu3[1]) - np.log(errosu3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 3: " + str((np.log(errosu4[1]) - np.log(errosu4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 4: " + str((np.log(errosu5[1]) - np.log(errosu5[0]))/(np.log(nels[1])- np.log(nels[0]))))

print("p = 1: " + str((np.log(errosp2[1]) - np.log(errosp2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("p = 2: " + str((np.log(errosp3[1]) - np.log(errosp3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("p = 3: " + str((np.log(errosp4[1]) - np.log(errosp4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("p = 4: " + str((np.log(errosp5[1]) - np.log(errosp5[0]))/(np.log(nels[1])- np.log(nels[0]))))


plt.title("Gráfico de convergência de U")
plt.plot(np.log(nels), np.log(errosu2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(errosu3),'-o', label="k=2")
plt.plot(np.log(nels[:-1]), np.log(errosu4[:-1]),'-o', label="k=3")
plt.plot(np.log(nels[:-3]), np.log(errosu5[:-3]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()

plt.title("Gráfico de convergência de P")
plt.plot(np.log(nels), np.log(errosp2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(errosp3),'-o', label="k=2")
plt.plot(np.log(nels[:-2]), np.log(errosp4[:-2]),'-o', label="k=3")
plt.plot(np.log(nels[:-4]), np.log(errosp5[:-4]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()  """


#Questao H
#erros de U
errosu2 = [0.0284883,0.00713322,0.00178401,0.000446046,0.000111514,2.78787e-05,6.96963e-06,1.74242e-06]
errosu3 = [0.000645901,8.08642e-05,1.0112e-05,1.26412e-06,1.58019e-07,1.97525e-08,2.47695e-09,1.05536e-09]
errosu4 = [1.40318e-05,8.78453e-07,5.49262e-08,3.43325e-09,2.17196e-10,1.00263e-10,9.78533e-10,4.0443e-09]
errosu5 = [2.5887e-07,8.09977e-09,2.5329e-10,1.0811e-11,1.19501e-10,1.59567e-10,1.01108e-09,6.21091e-10]

#erros de P
errosp2 = [0.00458334,0.00113838,0.000284127,7.10024e-05,1.77488e-05,4.43708e-06,1.10927e-06,2.77328e-07]
errosp3 = [0.000205142,2.57256e-05,3.2183e-06,4.02368e-07,5.02986e-08,6.2874e-09,7.8609e-10,1.64885e-10]
errosp4 = [4.4703e-06,2.79681e-07,1.74845e-08,1.09285e-09,6.87715e-11,1.3227e-11,4.74665e-11,3.41797e-10]
errosp5 = [8.2415e-08,2.57835e-09,8.05975e-11,7.01079e-12,2.06815e-11,1.7182e-11,2.40453e-10,5.10061e-10]


print("u = 1: " + str((np.log(errosu2[1]) - np.log(errosu2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 2: " + str((np.log(errosu3[1]) - np.log(errosu3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 3: " + str((np.log(errosu4[1]) - np.log(errosu4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("u = 4: " + str((np.log(errosu5[1]) - np.log(errosu5[0]))/(np.log(nels[1])- np.log(nels[0]))))

print("p = 1: " + str((np.log(errosp2[1]) - np.log(errosp2[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("p = 2: " + str((np.log(errosp3[1]) - np.log(errosp3[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("p = 3: " + str((np.log(errosp4[1]) - np.log(errosp4[0]))/(np.log(nels[1])- np.log(nels[0]))))
print("p = 4: " + str((np.log(errosp5[1]) - np.log(errosp5[0]))/(np.log(nels[1])- np.log(nels[0]))))


plt.title("Gráfico de convergência de U")
plt.plot(np.log(nels), np.log(errosu2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(errosu3),'-o', label="k=2")
plt.plot(np.log(nels[:-3]), np.log(errosu4[:-3]),'-o', label="k=3")
plt.plot(np.log(nels[:-4]), np.log(errosu5[:-4]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show()

plt.title("Gráfico de convergência de P")
plt.plot(np.log(nels), np.log(errosp2),'-o', label="k=1")
plt.plot(np.log(nels), np.log(errosp3),'-o', label="k=2")
plt.plot(np.log(nels[:-2]), np.log(errosp4[:-2]),'-o', label="k=3")
plt.plot(np.log(nels[:-4]), np.log(errosp5[:-4]),'-o', label="k=4")
plt.legend()
plt.grid(True)
plt.xlabel("log(nel)")
plt.ylabel("log(erro)")
plt.show() 