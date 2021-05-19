from cmath import *
from math import *
from numpy import linalg
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg

def distancia_entre_pontos(x1, x2, y1, y2):
    """
    função que calcula a distância entre dois pontos
    recebe: coordenadas de dois pontos [inteiro]
    retorna: distância entre os pontos [inteiro]
    """
    
    return sqrt((x2 - x1)**2 + (y2 - y1)**2)

def matriz_conectividade(numero_do_membro, incidencia, numero_de_elementos=3):
    """
    função que calcula a matriz de conectividade de um elemento específico
    recebe: número do elemento [inteiro] e a matriz de incidência lida do excel [matriz]
    retorna: conectividade [lista]
    """
    
    conectividade = numero_de_elementos*[0]
    
    # O numero do membro-1 é a linha a matriz que eu tenho que utilizar
    no_1 = int(incidencia[numero_do_membro-1, 0])
    no_2 = int(incidencia[numero_do_membro-1, 1])
    
    conectividade[no_1-1] = -1
    conectividade[no_2-1] = 1

    return conectividade

def conec_global_T(nos, incidencia):
    """
    função responsável por devolver a matriz de conectividade global
    recebe: número de nós [inteiro], matriz de incidencia [matriz]
    retorna: matriz onde cada linha é uma conectividade de um membro [matriz]
    """
  
    C = []
    for i in range(nos):
        C.append(matriz_conectividade(i+1, incidencia)) #repare que i começa em 0
    return np.array(C).T


def calculate_Se(n_membro, m_incidencia, m_nos, m_membros):
    """ 
    Função responsável por calcular o valor de Se para um dado elemento
    recebe: número do membro [inteiro], matriz de incidência [matriz], matriz dos nós [matriz], matriz dos membros [matriz]
    retorna: valor de Se para o elemento escolhido [matriz/inteiro]
    """
    
    E = m_incidencia[n_membro-1, 2] #elemento linear elástico
    A = m_incidencia[n_membro-1, 3] #área da senção transversal
    
    no_1 = int(m_incidencia[n_membro-1, 0])
    no_2 = int(m_incidencia[n_membro-1, 1])
    
    # Pegar a coordenada do no_1
    x_no1 = m_nos[0, no_1-1]
    y_no1 = m_nos[1, no_1-1]
    
    # Pegar a coordenada do no_2
    x_no2 = m_nos[0, no_2-1]
    y_no2 = m_nos[1, no_2-1]
    
    l = distancia_entre_pontos(x_no2, x_no1, y_no2, y_no1)
    
    cordeenadas_membro = m_membros[:, n_membro-1] #pega as coordenadas do membro
    coornadas_matriz = np.array([cordeenadas_membro]) #transforma em matriz
    coornadas_matriz_T = coornadas_matriz.T #calcula a transposta
    me = sum(i**2 for i in cordeenadas_membro)
    segunda_parte = (np.dot(coornadas_matriz_T, coornadas_matriz))/me
    return ((E*A)/l)*segunda_parte

def  calculate_K(n_membro, m_incidencia, m_nos, m_membros):
    """
    função responsável por calcular a matriz K para um elemento
    recebe: número do membro [inteiro], matriz incidência, matriz de nós e matriz de membros
    retorna: matriz K para o dado elemento
    """
    
    MC = np.array([matriz_conectividade(n_membro, m_incidencia)]) #matriz conectividade
    MC_T = np.transpose(MC) #matriz conectividade transposta
    Se = calculate_Se(n_membro, m_incidencia, m_nos, m_membros)
    dot = MC * MC_T
    return np.kron(dot, Se) 

def matriz_global(elementos, m_incidencia, m_nos, m_membros):
    """
    função responsável por realizar a somatória das matrizes K de cada elementos
    recebe: número de membros [inteiro], matriz incidência, matriz de nós e matriz de membros
    retorna: somatória das matrizes K de cada elemento [matriz]
    """
    
    get_shape = calculate_K(1, m_incidencia, m_nos, m_membros) 
    
    x = get_shape.shape[0] #linhas
    y = get_shape.shape[1] #colunas
    
    kg = np.zeros((x, y)) #sempre vai ser num_membros*2
    for i in range(elementos):
        kg += calculate_K(i, m_incidencia, m_nos, m_membros)        
    return kg

def MR_para_solucao(matriz):
    """
    função responsável por excluir as colunas/linhas desejadas de uma matriz
    recebe: matriz
    retorna: matriz com as linhas e colunas desejadas apagadas
    """
    
    matriz_certa = np.zeros((6, 6))
    matriz_certa = np.delete(matriz, (0, 2, 3), axis=0) #apaga 1º, 3º e 4º linha
    matriz_certa = np.delete(matriz_certa, (0, 2, 3), axis=1) #apaga 1º, 3º e 4º coluna
    return matriz_certa

def calcula_deslocamentos(matriz_rigidez, matriz_força):
    L,U = scipy.linalg.lu(matriz_rigidez, permute_l=True)
    y = scipy.linalg.solve(L, matriz_força)
    x = scipy.linalg.solve_triangular(U, y)
    return x

def calculate_force(matriz_k, u, linha_number):
    """
    função responsável por calcular a força para um dado elemento
    recebe: matriz K do elemento, matriz u (completo), número da linha do elemento
    """
    
    linha = matriz_k[linha_number,:] #pega a linha desejada na matriz
    return np.dot(linha, u) #retorna multiplicação de linha X coluna

def tensao_e_deformacao(n_elemento, n_de_membros, matriz_u, m_incidencia, m_nos):
    """
    função responsável por calcular a tensão e deformação para cada membro
    recebe: número do membro desejado [inteiro], número total de membros [inteiro], matriz u calculada (completa), matriz de incidencia, matriz de nós.
    retorna: tensão [inteiro] e deformação [inteiro] calculadas para o membro desejado
    """
    
    if n_elemento == n_de_membros:
        matriz_aux = np.array((
            [matriz_u[(n_elemento*2) - 2]], 
            [matriz_u[(n_elemento*2) - 1]], 
            [matriz_u[0]], 
            [matriz_u[1]]))
    else:
         matriz_aux = np.array((
            [matriz_u[(n_elemento*2) - 2]], 
            [matriz_u[(n_elemento*2) - 1]], 
            [matriz_u[(n_elemento*2)]], 
            [matriz_u[(n_elemento*2)+1]]))
    
    E =  m_incidencia[n_elemento-1, 2]  
    
    no_1 = int(m_incidencia[n_elemento-1, 0])
    no_2 = int(m_incidencia[n_elemento-1, 1])
    
    # Pegar a coordenada do no_1
    x_no1 = m_nos[0, no_1-1]
    y_no1 = m_nos[1, no_1-1]
    
    # Pegar a coordenada do no_2
    x_no2 = m_nos[0, no_2-1]
    y_no2 = m_nos[1, no_2-1]
    
    l = sqrt((x_no2-x_no1)**2+(y_no2-y_no1)**2)
    
    sen = (y_no2-y_no1)/l #calcula seno do elemento
    cos = (x_no2-x_no1)/l #calcula coss do elemento
    
    c = np.array(([-cos, -sen, cos, sen]))
    
    tensao = (E/l) * np.dot(c, matriz_aux)
    deformacao = (1/l) * np.dot(c, matriz_aux)
    
    return tensao[0], deformacao[0]

def solucao_gauss(k, F, ite, tol):
    """
    função responsável por calcular a solução de Gauss para um sistema de equações
    recebe: matriz k, matriz de forças, número de iterações [inteiro], tolerância [float]
    retorna: solução do sistema de equações através da teoria de Gauss [matriz]
    """
    
    matriz_x = np.zeros((F.shape[0], 1)) #cria uma matriz nx1
    
    for iteracao in range(ite):
        for indice in range(matriz_x.shape[0]):
            b = F[indice]
            ax = sum(a*x for a,x in zip(k[indice, :], matriz_x[:,0])) - k[indice, indice]*matriz_x[indice,0]            
            matriz_x[indice] = (b - ax)/k[indice, indice]
            
    return matriz_x

def solucao_jacobi(k, F, ite, tol):
    """
    função responsável por calcular a solução de Gauss para um sistema de equações
    recebe: matriz k, matriz de forças, número de iterações [inteiro], tolerância [float]
    retorna: solução do sistema de equações através da teoria de Gauss [matriz]
    """
    
    matriz_x = np.zeros((F.shape[0], 1)) #cria uma matriz nx1
    matriz_x_auxiliar = np.zeros((F.shape[0], 1))
    
    for iteracao in range(ite):
        for indice in range(matriz_x.shape[0]):
            b = F[indice]
            ax = sum(a*x for a,x in zip(k[indice, :], matriz_x[:,0])) - k[indice, indice]*matriz_x[indice,0]            
            matriz_x_auxiliar[indice] = (b - ax)/k[indice, indice]
            
        matriz_x = matriz_x_auxiliar
        #print(matriz_x)
        #print("--------------") #Com esse print é possivel perceber que a de Gauss converte antes
            
    return matriz_x