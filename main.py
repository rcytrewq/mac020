#from mpl_toolkits import mplot3d 
#from mpl_toolkits import axes3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
from tqdm import tqdm

##################################################################
##############                                   #################
############## MAC020-TRABALHO MULTIDISCIPLINAR  #################
##############                                   #################
############## Transporte em Meios Porosos:      #################
############## Escoamento de Fluidos  Imiscíveis #################
##############                                   #################
############## João Victor Lopes Borges          #################
############## Yuri Ramos Correa                 #################
##############                                   #################
##################################################################

##################################################################
##################################################################
##################################################################
##################################################################



##################################################################
###############                                 ##################
############### Definições de variáveis globais ##################
###############                                 ##################
##################################################################

theta = 10.#permeabilidade 
phi = 0.63636364 #porosidade
mi_w = 1. #viscosidade da agua
mi_o = 5. #viscosidade do oleo
rho_w = 1. #densidade da agua
rho_o = 0.8 #densidade do oleo
swi=0.2 #saturacao de agua irreduzivel
sor=0.2 #saturacao de oleo residual

i=0 #contador
j=0 #contador
cont = 0 #contador
J = 2.18
dominioX = [0.0, 1.]
dominioY = [0.0, 1.]

nelX = 21 #tamanho da malha em x
dx = (dominioX[1]-dominioX[0])/(nelX-1) #delta x

dy = dx #delta y
h=dx
dt = h**2 #delta t

numit = (250) #numero de iteracoes
sig = 0 #contorno de fluxo = 0

nelY = int (np.ceil((dominioY[1]-dominioY[0])/dy)+1) #tamanho da malha em y

F = np.zeros((nelX*nelY)) #vetor fonte
S = np.zeros((nelX, nelY)) #matriz de s em cada ponto da malha

for i in range(nelX):
	for j in range (nelY):
		S[i][j] = 0 #condicao inicial nula

A = np.zeros((nelX*nelY, nelX*nelY)) #matriz A de Au=F
K = np.zeros((nelX, nelY)) # armazena o valor de K(kro/mi_o + krw/mi_w)
mag = np.zeros((nelX, nelY)) #matriz de magnitude da velocidade
S_new = np.zeros((nelX, nelY)) #matriz de s no passo de tempo seguinte
vely = np.zeros((nelX, nelY)) #matriz de velcidade em y
velx = np.zeros((nelX, nelY)) #matriz de velocidade em x

##################################################################
##################################################################
##################################################################
##################################################################


##################################################################
###############                                 ##################
###############       Definições de funcoes     ##################
###############                                 ##################
##################################################################


def permW(s): #permeabilidade relativa da agua
	return 0.4*((s-swi)/(1-sor-swi))**2
	
def permO(s): #permeabilidade relativa do oleo
	return 0.8*((1-sor-s)/(1-sor-swi))**2
		
def efe(s): #funcao f(s)
	lambda_w = permW(s)/mi_w
	lambda_o = permO(s)/mi_o
	lambda_tot = lambda_w+lambda_o
	
	return (lambda_w/lambda_tot)


def kaI(i,j): #retorna k_(i+1/2). Para k_(i-1/2), usar i-1 como parametro
	if i<0: return K[i+1][j]/2
	elif i+1>= nelX: return K[i][j]/2
	else: return (K[i][j]+K[i+1][j])/2
	
def kaJ(i,j): #retorna k_(j+1/2). Para k_(j-1/2), usar j-1 como parametro
	if j<0: return K[i][j+1]/2
	elif j+1>= nelX: return K[i][j]/2
	else: return (K[i][j]+K[i][j+1])/2



def func (i,j): #funcao que retorna o termo fonte
	if (i==np.round(nelX/2)-1 and j==0):
		return J
	elif (i==np.round(nelX/2)-1 and j ==nelY-1):
		return -J*efe(S[i][j])
	else: return 0
	
def func2 (i,j):
	if (i==np.round(nelX/2)-1 and j==0):
		return J*efe(S[i][j])
	elif (i==np.round(nelX/2)-1 and j ==nelY-1):
		return efe(S[i][j])*func(i,j)
	else: return 0


##################################################################
##################################################################
##################################################################
##################################################################	

##################################################################
###############                                 ##################
###############       Inicio da aproximacao     ##################
###############                                 ##################
##################################################################


itera = np.arange(0, numit, 1) #inicio da iteracao no tempo




for it in tqdm (itera): #barra de progresso
	
		
	for i in range (nelX):
		for j in range (nelY):
			K[i][j] = -theta*((permO(S[i][j])/mi_o) + (permW(S[i][j])/mi_w)) #atualiza a matriz K
			
		

		
	#contorno de fluxo nulo, nao precisa incluir funcao

	for i in range (nelX):
		for j in range (nelY):
			F[(nelY)*j + i] = func(i,j) #vetor fonte
			
	

	
	#######################################################
	###############    Criacao da matriz A    #############
	#######################################################



	for j in range (nelX):
		for i in range (nelX):
			if(j==0):
				if(i==0):
					A[i][j] = -(kaI(i,j) + kaI(i-1,j) + kaJ(i,j) + kaJ(i,j-1)) #Diagonal principal
					A[i][j+1] = kaI(i,j) #Diagonal superior
					A[i][nelX+j] = kaJ(i,j) #Diagonal direita
				elif(i==nelX-1):
					A[i][i] = -(kaI(i,j) + kaI(i-1,j) + kaJ(i,j) + kaJ(i,j-1)) #Diagonal principal
					A[i][nelX+i] = kaJ(i,j) #Diagonal direita
					A[i][i-1] = kaI(i-1,j) #Diagonal inferior
				else:
					A[i][i] = -(kaI(i,j) + kaI(i-1,j) + kaJ(i,j) + kaJ(i,j-1)) #Diagonal principal
					A[i][i+1] = kaI(i,j) #Diagonal superior
					A[i][i-1] = kaI(i-1,j) #Diagonal inferior
					A[i][nelX+i] = kaJ(i,j) #Diagonal direita

			elif(j==nelX-1):
				if(i==0):
					A[nelX*j][nelX*(j)] = -(kaI(i,j) + kaI(i-1,j) + kaJ(i,j) + kaJ(i,j-1)) #Diagonal principal
					A[nelX*j][nelX*(j)+1] = kaI(i,j) #Diagonal superior
					A[nelX*j][nelX*(j-1)] = kaJ(i,j-1)	#Diagonal esquerda
				elif(i==nelX-1):
					A[nelX*nelX-1][nelX*nelX-1] = -(kaI(i,j) + kaI(i-1,j) + kaJ(i,j) + kaJ(i,j-1)) #Diagonal principal
					A[nelX*nelX-1][nelX*nelX-2] = kaI(i-1,j) #Diagonal inferior
					A[nelX*nelX-1][nelX*j-1] = kaJ(i,j-1) #Diagonal esquerda
				else:
					A[nelX*j+i][nelX*j+i] = -(kaI(i,j) + kaI(i-1,j) + kaJ(i,j) + kaJ(i,j-1)) #Diagonal principal
					A[nelX*j+i][nelX*j+i+1] = kaI(i,j) #Diagonal superior
					A[nelX*j+i][nelX*j+i-1] = kaI(i-1,j) #Diagonal inferior
					A[nelX*j+i][nelX*(j-1)+i] = kaJ(i,j-1) #Diagonal esquerda
					
			else:
				if(i==0):
					A[nelX*j][nelX*j] = -(kaI(i,j) + kaI(i-1,j) + kaJ(i,j) + kaJ(i,j-1)) #Diagonal principal
					A[nelX*j][nelX*j+1] = kaI(i,j) #Diagonal superior
					A[nelX*j][nelX*(j+1)+i] = kaJ(i,j) #Diagonal direita
					A[nelX*j][nelX*(j-1)+i] = kaJ(i,j-1) #Diagonal esquerda
				elif(i==nelX-1):
					A[nelX*j+i][nelX*j+i] = -(kaI(i,j) + kaI(i-1,j) + kaJ(i,j) + kaJ(i,j-1)) #Diagonal principal
					A[nelX*j+i][nelX*j+i-1] = kaI(i-1,j) #Diagonal inferior
					A[nelX*j+i][nelX*(j+1)+i] = kaJ(i,j) #Diagonal direita
					A[nelX*j+i][nelX*(j-1)+i] = kaJ(i,j-1) #Diagonal esquerda
				else:
					A[nelX*j+i][nelX*j+i] = -(kaI(i,j) + kaI(i-1,j) + kaJ(i,j) + kaJ(i,j-1)) #Diagonal principal
					A[nelX*j+i][nelX*j+i+1] = kaI(i,j) #Diagonal superior
					A[nelX*j+i][nelX*j+i-1] = kaI(i-1,j) #Diagonal inferior
					A[nelX*j+i][nelX*(j+1)+i] = kaJ(i,j) #Diagonal direita
					A[nelX*j+i][nelX*(j-1)+i] = kaJ(i,j-1) #Diagonal esquerda
					
	##############################################################
	##############################################################
	##############################################################
	##############################################################	

	##############################################################
	##############                                 ###############
	##############        Solucao do sistema       ###############
	##############                                 ###############
	##############################################################

	F_ = h**2*F

	P = np.linalg.solve((A*h**2), F_) #Resolve o sistema
	P = np.reshape(P,(nelX, nelY)) #Ajusta a solucao em 2D
	P = P.T #Ajusta a solucao em 2d (.T = matriz transposta)

	
	##############################################################
	##############################################################
	##############################################################
	##############################################################	

		
	#######################################################
	###############    Calculo de vx e vy     #############
	#######################################################


	for i in range(0, nelX):
		for j in range (0, nelY):
			t =(-theta*((permO(S[i][j])/mi_o)+(permW(S[i][j])/mi_w)))
			if (i!=0 and j!=0 and i!=nelX-1 and j!=nelY-1):	velx[i][j] = t*((P[i+1][j]-P[i-1][j])/(2*dx))
			

	for i in range(0, nelX):
		for j in range (0, nelY):
			t =(-theta*((permO(S[i][j])/mi_o)+(permW(S[i][j])/mi_w)))
			if (i!=0 and j!=0 and i!=nelX-1 and j!=nelY-1):	vely[i][j] = t*((P[i][j+1]-P[i][j-1])/(2*dx))
			
	for i in range(nelX):
		velx[i][0] = velx[i][1]
		vely[i][0] = vely[i][1]
		velx[i][nelX-1] = velx[i][nelX-2]
		vely[i][nelX-1] = vely[i][nelX-2]
		
	for j in range(nelY):
		velx[0][j] = velx[1][j]
		vely[0][j] = vely[1][j]
		velx[nelX-1][j] = velx[nelX-2][j]
		vely[nelX-1][j] = vely[nelX-2][j]


	#######################################################
	###########  calculo da magnitude da vel  #############
	#######################################################

	for i in range (nelX):
		for j in range (nelY):
			mag[i][j] = np.sqrt((velx[i][j]**2)+(vely[i][j]**2))

	##############################################################
	##############################################################
	##############################################################
	##############################################################	

	
	#######################################################
	#####   calculo de S no proximo passo de tempo   ######
	#######################################################


	for i in range (1,nelX-1):
		for j in range (1,nelY-1):
			if (velx[i][j]>=0): 
				if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i-1][j]))*velx[i][j] + (efe(S[i][j])-efe(S[i][j-1]))*vely[i][j]   )/h )) + S[i][j]
				else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i-1][j]))*velx[i][j] + (efe(S[i][j+1])-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
			else:
				if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i+1][j])-efe(S[i][j]))*velx[i][j] + (efe(S[i][j])-efe(S[i][j-1]))*vely[i][j]   )/h )) + S[i][j]
				else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i+1][j])-efe(S[i][j]))*velx[i][j] + (efe(S[i][j+1])-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
				
	for i in range(1,nelX-1):
	
	#CONTORNO SUPERIOR
		if (velx[0][i]>=0): 
			if(vely[0][i]>=0):
				S_new[0][i] = ((dt/phi) * (func2(0,i) - ( (efe(S[0][i])-efe(S[1][i]-2*h*func2(0,i)))*velx[0][i] + (efe(S[0][i])-efe(S[0][i-1]))*vely[0][i]   )/h )) + S[0][i]
	
			else: S_new[0][i] = ((dt/phi) * (func2(0,i) - ( (efe(S[0][i])-efe(S[1][i]-2*h*func2(0,1)))*velx[0][i] + (efe(S[0][i+1])-efe(S[0][i]))*vely[0][i]   )/h )) + S[0][i]
				
		else:
			if(vely[0][i]>=0): S_new[0][i] = ((dt/phi) * (func2(0,i) - ( (efe(S[1][i])-efe(S[0][i]))*velx[0][i] + (efe(S[0][i])-efe(S[0][i-1]))*vely[0][i]   )/h )) + S[0][i]
			else: S_new[0][i] = ((dt/phi) * (func2(0,i) - ( (efe(S[1][i])-efe(S[0][i]))*velx[0][i] + (efe(S[0][i+1])-efe(S[0][i]))*vely[0][i]   )/h )) + S[0][i]
	
	
	#CONTORNO ESQUERDO
	
		if (velx[i][0]>=0):
			if(vely[i][0]>=0): S_new[i][0] = ((dt/phi) * (func2(i,0) - ( (efe(S[i][0])-efe(S[i-1][0]))*velx[i][0] + (efe(S[i][0])-efe(S[i][1]+2*h*func2(i,0)))*vely[i][0]   )/h )) + S[i][0]
			else: S_new[i][0] = ((dt/phi) * (func2(i,0) - ( (efe(S[i][0])-efe(S[i-1][0]))*velx[i][0] + (efe(S[i][1])-efe(S[i][0]))*vely[i][0]   )/h )) + S[i][0]
		
		else:
			if(vely[i][0]>=0): S_new[i][0] = ((dt/phi) * (func2(i,0) - ( (efe(S[i+1][0])-efe(S[i][0]))*velx[i][0] + (efe(S[i][0])-efe(S[i][1]+2*j*func2(i,0)))*vely[i][0]   )/h )) + S[i][0]
			else: S_new[i][0] = ((dt/phi) * (func2(i,0) - ( (efe(S[i+1][0])-efe(S[i][0]))*velx[i][0] + (efe(S[i][1])-efe(S[i][0]))*vely[i][0]   )/h )) + S[i][0]
		
	#CONTORNO INFERIOR		
		
		if (velx[nelX-1][i]>=0):
			if(vely[nelX-1][i]>=0): S_new[nelX-1][i] = ((dt/phi) * (func2(nelX-1,i) - ( (efe(S[nelX-1][i])-efe(S[nelX-2][i]))*velx[nelX-1][i] + (efe(S[nelX-1][i])-efe(S[nelX-1][i-1]))*vely[i][j]   )/h )) + S[nelX-1][j]
			else: S_new[nelX-1][i] = ((dt/phi) * (func2(nelX-1,i) - ( (efe(S[nelX-1][i])-efe(S[nelX-2][i]))*velx[nelX-1][i] + (efe(S[nelX-1][i+1])-efe(S[nelX-1][i]))*vely[nelX-1][i]   )/h )) + S[nelX-1][i]
			
		else:
			if(vely[nelX-1][i]>=0): S_new[nelX-1][i] = ((dt/phi) * (func2(nelX-1,i) - ( (efe(S[nelX-2][j]-2*h*func2(nelX-1,i))-efe(S[nelX-1][i]))*velx[nelX-1][i] + (efe(S[nelX-1][i])-efe(S[nelX-1][i-1]))*vely[i][j]   )/h )) + S[nelX-1][i]
			else: S_new[nelX-1][i] = ((dt/phi) * (func2(nelX-1,i) - ( (efe(S[nelX-2][i]-2*h*func2(nelX-1,i))-efe(S[nelX-1][i]))*velx[nelX-1][i] + (efe(S[nelX-1][i+1])-efe(S[nelX-1][i]))*vely[nelX-1][i]   )/h )) + S[nelX-1][i]
		
	#CONTORNO DIREITO
	
		if(velx[i][nelX-1]>=0):
			if(vely[i][nelX-1]>=0): S_new[i][nelX-1] = ((dt/phi) * (func2(i,nelX-1) - ( (efe(S[i][j])-efe(S[i-1][j]))*velx[i][j] + (efe(S[i][j])-efe(S[i][j-1]))*vely[i][j]   )/h )) + S[i][j]
			else: S_new[i][nelX-1] = ((dt/phi) * (func2(i,nelX-1) - ( (efe(S[i][nelX-1])-efe(S[i-1][nelX-1]))*velx[i][nelX-1] + (efe(S[i][nelX-2]- 2*h*func2(i,nelX-1))-efe(S[i][nelX-1]))*vely[i][nelX-1]   )/h )) + S[i][nelX-1]
		else:
			if(vely[i][nelX-1]>=0): S_new[i][nelX-1] = ((dt/phi) * (func2(i,nelX-1) - ( (efe(S[i+1][nelX-1])-efe(S[i][nelX-1]))*velx[i][nelX-1] + (efe(S[i][nelX-1])-efe(S[i][nelX-2]))*vely[i][nelX-1]   )/h )) + S[i][nelX-1]
			else: S_new[i][nelX-1] = ((dt/phi) * (func2(i,nelX-1) - ( (efe(S[i+1][nelX-1])-efe(S[i][nelX-1]))*velx[i][nelX-1] + (efe(S[i][nelX-2]-2*h*func2(i, nelX-1))-efe(S[i][nelX-1]))*vely[i][nelX-1]   )/h )) + S[i][nelX-1]
			
	#PONTO (0,0)
	
	i=0
	j=0
	
	if (velx[i][j]>=0): 
		if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i+1][j])+2*h*func2(i,j))*velx[i][j] + (efe(S[i][j])-efe(S[i][j+1]+2*h*func2(i,j)))*vely[i][j]   )/h )) + S[i][j]
		else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i+1][j]+2*h*func2(i,j)))*velx[i][j] + (efe(S[i][j+1])-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
	else:
		if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i+1][j])-efe(S[i][j]))*velx[i][j] + (efe(S[i][j])-efe(S[i][j+1]+2*h*func2(i,j)))*vely[i][j]   )/h )) + S[i][j]
		else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i+1][j])-efe(S[i][j]))*velx[i][j] + (efe(S[i][j+1])-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
	
	#PONTO (nelX-1,0)
	
	i=nelX-1
	j=0
	
	if (velx[i][j]>=0): 
		if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i-1][j]))*velx[i][j] + (efe(S[i][j])-efe(S[i][j+1]+2*h*func2(i,j)))*vely[i][j]   )/h )) + S[i][j]
		else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i-1][j]))*velx[i][j] + (efe(S[i][j+1])-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
	else:
		if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i-1][j]-2*h*func2(i,j))-efe(S[i][j]))*velx[i][j] + (efe(S[i][j])-efe(S[i][j+1]+2*h*func2(i,j)))*vely[i][j]   )/h )) + S[i][j]
		else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i-1][j]-2*h*func2(i,j))-efe(S[i][j]))*velx[i][j] + (efe(S[i][j+1])-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
	
	#PONTO (0, nelX-1)
	i=0
	j=nelX-1
	
	
	if (velx[i][j]>=0): 
		if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i+1][j]+2*h*func2(i,j)))*velx[i][j] + (efe(S[i][j])-efe(S[i][j-1]))*vely[i][j]   )/h )) + S[i][j]
		else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i+1][j]+2*h*func2(i,j)))*velx[i][j] + (efe(S[i][j-1]-2*h*func2(i,j))-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
	else:
		if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i+1][j])-efe(S[i][j]))*velx[i][j] + (efe(S[i][j])-efe(S[i][j-1]))*vely[i][j]   )/h )) + S[i][j]
		else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i+1][j])-efe(S[i][j]))*velx[i][j] + (efe(S[i][j-1]-2*h*func2(i,j))-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
		
		
	#PONTO (nelX-1, nelX-1)
	
	i=nelX-1
	j=nelX-1
	
	if (velx[i][j]>=0): 
		if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i-1][j]))*velx[i][j] + (efe(S[i][j])-efe(S[i][j-1]))*vely[i][j]   )/h )) + S[i][j]
		else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i][j])-efe(S[i-1][j]))*velx[i][j] + (efe(S[i][j-1]-2*h*func2(i,j))-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
	else:
		if(vely[i][j]>=0): S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i-1][j]-2*h*func2(i,j))-efe(S[i][j]))*velx[i][j] + (efe(S[i][j])-efe(S[i][j-1]))*vely[i][j]   )/h )) + S[i][j]
		else: S_new[i][j] = ((dt/phi) * (func2(i,j) - ( (efe(S[i-1][j]-2*h*func2(i,j))-efe(S[i][j]))*velx[i][j] + (efe(S[i][j-1]-2*h*func2(i,j))-efe(S[i][j]))*vely[i][j]   )/h )) + S[i][j]
	
	##############################################################
	##############################################################
	##############################################################
	##############################################################
	##############################################################
	##############################################################
	##############################################################
	##############################################################	

		
	#######################################################
	##########    Atualizacao das variaveis    ############
	##########    Sequencia do Loop            ############
	#######################################################
	
	
	
	for i in range (nelX):
		for j in range (nelY):
			S[i][j] = abs(S_new[i][j])
			S_new[i][j] = 0
			K[i][j] = 0
			F[nelX*j+i] = 0
	
	#######################################################
	##########    Salvar graficos da solucao   ############
	#######################################################
	
	if (it%10==0):
		cont += 1
		fig = plt.figure()
		ax = fig.gca(projection='3d')

		# Make data.
		X_ = np.arange(dominioX[0], dominioX[1], dx)
		Y_ = np.arange(dominioX[0], dominioY[1], dy)
		X_, Y_ = np.meshgrid(X_, Y_)

		colors=np.arange(-10,15,1)
		plt.rcParams['figure.figsize'] = [8, 6]
		fig = plt.figure()
		
		plt.imshow(S[0:nelX-2][0:nelX-2], cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

		#plt.contour(P, levels=15,linewidths=0.5)
		#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
		plt.clim(0,0.8)
		plt.colorbar()
		plt.title("Tempo = "+str(dt*it))

		plt.xlabel('Y')
		plt.ylabel('X')
		#plt.savefig("Imgs/"+str(nelX)+"X"+str(nelX)+"_"+str(theta)+"_"+str(swi)+"_"+str(sor)+"/"+str(cont)+".png")

		plt.close()
		plt.close()
			
	##############################################################
	##############################################################
	##############################################################
	##############################################################	

	##############################################################
	##############                                 ###############
	##############        Final da aproximacao     ###############
	##############                                 ###############
	##############################################################



##############################################################
##############################################################
##############################################################
##############################################################	
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################	
##############################################################
##############################################################

##############################################################
##############                                 ###############
##############   Graficos ao fim da execucao   ###############
##############                                 ###############
##############################################################


##############################################################
##############                                 ###############
##############         Grafico da pressao      ###############
##############                                 ###############
##############################################################

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X_ = np.arange(dominioX[0], dominioX[1]+dx, dx)
Y_ = np.arange(dominioX[0], dominioY[1]+dy, dy)
X_, Y_ = np.meshgrid(X_, Y_)

colors=np.arange(-10,15,1)
plt.rcParams['figure.figsize'] = [8, 6]
fig = plt.figure()

#ax = fig.gca(projection='3d')

#ax.plot_surface(X_, Y_, P, cmap=cm.coolwarm,linewidth=0, antialiased=False)
#ax.plot_surface(x, y, solucao, cmap=plt.cm.gray ,linewidth=0, antialiased=False, vmin = 0, vmax = amplitude)
#ax.set_title("t="+"{:.2f}".format(tempo)+"s")
#ax.set_zlim([0,amplitude])

#plt.contourf(P, levels=15,cmap='jet', linewidths=0.5)
plt.imshow(P, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

#plt.contour(P, levels=15,linewidths=0.5)
#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
plt.colorbar()
plt.title("TITLE")

plt.xlabel('Y')
plt.ylabel('X')
plt.show()

plt.close()

##############################################################
##############################################################
##############################################################
##############################################################	

##############################################################
##############                                 ###############
##############          Grafico de Vx          ###############
##############                                 ###############
##############################################################

plt.imshow(velx, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

#plt.contour(P, levels=15,linewidths=0.5)
#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
plt.colorbar()
plt.title("TITLE")

plt.xlabel('Y')
plt.ylabel('X')
plt.show()

plt.close()

##############################################################
##############################################################
##############################################################
##############################################################	

##############################################################
##############                                 ###############
##############           Grafico de Vy         ###############
##############                                 ###############
##############################################################


plt.imshow(vely, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

#plt.contour(P, levels=15,linewidths=0.5)
#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
plt.colorbar()
plt.title("TITLE")

plt.xlabel('Y')
plt.ylabel('X')
plt.show()

plt.close()


##############################################################
##############################################################
##############################################################
##############################################################	

##############################################################
##############                                 ###############
##############   Grafico da magnitude da vel   ###############
##############                                 ###############
##############################################################

plt.imshow(mag, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

#plt.contour(P, levels=15,linewidths=0.5)
#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
plt.colorbar()
plt.title("TITLE")

plt.xlabel('Y')
plt.ylabel('X')
plt.show()

plt.close()

##############################################################
##############################################################
##############################################################
##############################################################	

##############################################################
##############                                 ###############
##############       Grafico da saturacao      ###############
##############                                 ###############
##############################################################

plt.imshow(S, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

#plt.contour(P, levels=15,linewidths=0.5)
#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
plt.clim(0,0.8)
plt.colorbar()
plt.title("TITLE")

plt.xlabel('Y')
plt.ylabel('X')
plt.show()

plt.close()


