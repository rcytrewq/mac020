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
i=0
j=0

J = 2.18
dominioX = [0.0, 1.]
dominioY = [0.0, 1.]
#tempo = [0., .]




nelX = 21
dx = (dominioX[1]-dominioX[0])/(nelX-1)
print(dx)

dy = dx
h=dx
dt = h**2

#T_ = np.arange(tempo[0], tempo[1]+dt, dt)

numit = (500)
sig = 0 #contorno de fluxo = 0

nelY = int (np.ceil((dominioY[1]-dominioY[0])/dy)+1)

F = np.zeros((nelX*nelY))
S = np.zeros((nelX, nelY))
for i in range(nelX):
	for j in range (nelY):
		S[i][j] = 0

A = np.zeros((nelX*nelY, nelX*nelY))
K = np.zeros((nelX, nelY))
mag = np.zeros((nelX, nelY))
S_new = np.zeros((nelX, nelY))
vely = np.zeros((nelX, nelY))
velx = np.zeros((nelX, nelY))

def efe(s):
	lambda_w = permW(s)/mi_w
	lambda_o = permO(s)/mi_o
	lambda_tot = lambda_w+lambda_o
	
	return (lambda_w/lambda_tot)


def permW(s):
	return 0.4*((s-swi)/(1-sor-swi))**2
	
def permO(s):
	return 0.8*((1-sor-s)/(1-sor-swi))**2	

def kaI(i,j): #retorna k_(i+1/2). Para k_(i-1/2), usar i-1 como parametro
	if i<0: return K[i+1][j]/2
	elif i+1>= nelX: return K[i][j]/2
	else: return (K[i][j]+K[i+1][j])/2
	
def kaJ(i,j):
	if j<0: return K[i][j+1]/2
	elif j+1>= nelX: return K[i][j]/2
	else: return (K[i][j]+K[i][j+1])/2



def func (i,j):
	if (i==np.round(nelX/2)-1 and j==0):
		#print("AAAAAAAAAAAAAAAAAAAAA")
		return J#2*h*J/K[i][j]
	elif (i==np.round(nelX/2)-1 and j ==nelY-1):
		#print("BBBBBBBBBBBBBBBBB")
		return -J#2*h*J/K[i][j]
	else: return 0
	
def func2 (i,j):
	if (i==np.round(nelX/2)-1 and j==0):
		#print("AAAAAAAAAAAAAAAAAAAAA")
		return J
	elif (i==np.round(nelX/2)-1 and j ==nelY-1):
		#print("BBBBBBBBBBBBBBBBB")
		return 0#efe(S[i][j])*func(i,j)
	else: return 0


	




itera = np.arange(0, numit, 1)




for it in tqdm (itera):

	
	
	for i in range (nelX):
		for j in range (nelY):
			K[i][j] = -theta*((permO(S[i][j])/mi_o) + (permW(S[i][j])/mi_w))
			
		

	'''
	THETA = np.zeros((nelX, nelY))
	for i in range (nelX):
		for j in range (nelY):
			THETA[i][j] = -K*((permO(S[i][j])/mi_o)+(permW(S[i][j])/mi_w))

	#s=0.9 #saturacao inicial
	'''

	#theta = -K*((permO(s)/mi_o)+(permW(s)/mi_w))

	#print(THETA)


		
	#contorno de fluxo nulo, nao precisa incluir funcao

	for i in range (nelX):
		for j in range (nelY):
			F[(nelY)*j + i] = func(i,j)
			
	

	'''
	B = np.zeros((nelX, nelY))
	for i in range (nelX):
		for j in range (nelY):
			if (i==j): B[i][j] = 4
			elif(abs(i-j)==1):
				if((i==0 and j==1) or (i==nelX-1 and j==nelY-2)):
					B[i][j] = -2
				else:
					B[i][j] = -1
			else: B[i][j]=0
		
				
	I = np.zeros((nelX, nelY))
	for i in range (nelX):
		I[i][i] = -1

	O = np.zeros((nelX, nelY))
	'''

	##################################################################
	###############                                 ##################
	###############       Criacao da matriz A       ##################
	###############                                 ##################
	##################################################################




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
					

	'''
	mat_list = [] 
	for i in range(0,nelX): # generate row
		tmp = []
		for j in range(0,nelY): # loop through A^j*B
			if (j == i): tmp.append(B)
			elif (abs(j-i)== 1):
				if ((i==0 and j==1) or (i==nelX-1 and j==nelY-2)):
					tmp.append(2*I)
				else:
					tmp.append(I)
			else: tmp.append(O)
		print(tmp)
		mat_list.append(tmp)
		

	A = np.block(mat_list)
	print(len(A))
	print(len(F))
	A_=np.zeros((nelX*nelY, nelY*nelX))

	for q in range(nelX*nelY):
		for i in range (nelX):
			for j in range (nelY):
				A_[q][nelX*i+j] = -THETA[j][i]

	A_ = A*A_*(h**2)
	'''


	F_ = h**2*F

	P = np.linalg.solve((A*h**2), F_)
	P = np.reshape(P,(nelX, nelY))
	P = P.T




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




	for i in range (nelX):
		for j in range (nelY):
			mag[i][j] = np.sqrt((velx[i][j]**2)+(vely[i][j]**2))




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
	
	
	'''
		else:
			if(vely[i][0]>=0):
			
			else:		
				
				
		S_new[0][i] = S_new[1][i]
		S_new[nelX-1][i] = S_new[i][nelX-2]


	'''
	for i in range (nelX):
		for j in range (nelY):
			S[i][j] = abs(S_new[i][j])
			S_new[i][j] = 0
			K[i][j] = 0
			F[nelX*j+i] = 0
			
	
	fig = plt.figure()
	ax = fig.gca(projection='3d')

	# Make data.
	X_ = np.arange(dominioX[0], dominioX[1]+dx, dx)
	Y_ = np.arange(dominioX[0], dominioY[1]+dy, dy)
	X_, Y_ = np.meshgrid(X_, Y_)

	colors=np.arange(-10,15,1)
	plt.rcParams['figure.figsize'] = [8, 6]
	fig = plt.figure()
	
	plt.imshow(S, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

	#plt.contour(P, levels=15,linewidths=0.5)
	#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
	plt.clim(0,0.8)
	plt.colorbar()
	#plt.title("TITLE")

	plt.xlabel('Y')
	plt.ylabel('X')
	plt.savefig("Ok/nel21_tf60_"+str(it)+".png")

	plt.close()
	plt.close()
			
	'''
	fig = plt.figure()
	ax = fig.gca(projection='3d')

	# Make data.
	X_ = np.arange(dominioX[0], dominioX[1]+dx, dx)
	Y_ = np.arange(dominioX[0], dominioY[1]+dy, dy)
	X_, Y_ = np.meshgrid(X_, Y_)
	
	fig = plt.figure()
	colors=np.arange(-10,15,1)
	plt.rcParams['figure.figsize'] = [8, 6]
	
	plt.imshow(S, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

	#plt.contour(P, levels=15,linewidths=0.5)
	#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
	plt.clim(0,0.4)
	plt.colorbar()
	plt.title("TITLE")

	plt.xlabel('Y')
	plt.ylabel('X')
	plt.savefig("figura"+str(it)+".png")

	plt.close()
	
	
	'''
	
	
	
'''
	print(S)
	print(S_new)

'''


'''



velx = np.zeros((nelX, nelY))
for i in range(0, nelX):
	for j in range (0, nelY):
		t =(-K*((permO(P[i][j])/mi_o)+(permW(P[i][j])/mi_w)))
		if (i!=0 and j!=0 and i!=nelX-1 and j!=nelY-1):	velx[i][j] = t*((P[i+1][j]-P[i-1][j])/(2*dx))
		
vely = np.zeros((nelX, nelY))
for i in range(0, nelX):
	for j in range (0, nelY):
		t =(-K*((permO(P[i][j])/mi_o)+(permW(P[i][j])/mi_w)))
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


vel=np.zeros((nelX, nelY))
for i in range (nelX):
	for j in range (nelY):
		vel[i][j] = np.sqrt(velx[i][j]**2 + vely[i][j]**2)



print(P)
'''

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

plt.imshow(velx, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

#plt.contour(P, levels=15,linewidths=0.5)
#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
plt.colorbar()
plt.title("TITLE")

plt.xlabel('Y')
plt.ylabel('X')
plt.show()

plt.close()

plt.imshow(vely, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

#plt.contour(P, levels=15,linewidths=0.5)
#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
plt.colorbar()
plt.title("TITLE")

plt.xlabel('Y')
plt.ylabel('X')
plt.show()

plt.close()


plt.imshow(mag, cmap='jet',interpolation='nearest', extent=[0,X_[-1][-1], 0, Y_[-1][-1]])

#plt.contour(P, levels=15,linewidths=0.5)
#surf = ax.plot_surface(X_, Y_, P, cmap='jet')
plt.colorbar()
plt.title("TITLE")

plt.xlabel('Y')
plt.ylabel('X')
plt.show()

plt.close()

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

