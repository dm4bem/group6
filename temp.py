#Bibliothèques:
import numpy as np

#Variables :
nq, n0 = 40, 29 # number of flow-rates branches and of temperaure nodes
ho, hi = 25, 8 # en W/m2/K
lambda_air = 0.026 # en W/m/K
lambda_wind = 1.015 # en W/m/K
lambda_ext = 0.22 # mur exterieur ; en W/m/K
lambda_isol = 0.023 #isolation ; en W/m/K
lambda_int = 0.33 #mur interieur ; en W/m/K
lambda_door = 0.17 #mur interieur ; en W/m/K
H = 3 # hauteur du mur; en m 
w_mur, w_isol, w_wind, w_door = 0.1, 0.5, 0.025, 0.05 #epaisseur des couches de mur, d'isolant, et des fenètres et de la porte; en m
Li, Lw1, Lw2, Ld = [0.2, 0.2, 0.1, 1]
L1, L2, L3, L4, L5, L6 = [3, 3, 3, 3, 3, 3]  #largeur des différents murs, fenètres, et porte ; en m
rhoa = 1.2 #kg/m3, masse volumique de l'air 
ca = 1000 #J/kg/K, chaleur spécifique de l'air
ACH = 2 #2/hour, air change per hour in room 2, ventilation
To = 275 # outside; en K 
Ti = 292 # inside; en K
Th = 295 # heating; en K
E = 200 # en W/m2
alpham = 0.55 #coefficient d'absorbtion du mur pour du béton lisse
alphaw = 0.02 #coefficient d'absorbtion d'une fenetre
C_layer_out = 18216000
C_layer_in = 239580
Glass = 1089000
Air = 32400

#Matrix A

A=np.zeros([nq, n0])
for i in range(4):
    A[i,i]=1
    A[i+1,i]=-1
A[5,4], A[6,4]=1, -1
for i in range(7,11):
    A[i,i-2]=1
    A[i+1,i-2]=-1
A[12,9], A[13,9]= 1, -1
for i in range(12,15):
    A[i,i-3]=1
    A[i+1,i-3]=-1
for i in range(17,20):
    A[i,i-5]=1
    A[i+1,i-5]=-1    
for i in range(23,27):
    A[i,i-8]=1
    A[i+1,i-8]=-1     
for i in range(28,32):
    A[i,i-9]=-1
    A[i+1,i-9]=1
A[33,23], A[34,23]=-1, 1
for i in range(35,38):
    A[i,i-11]=-1
    A[i+1,i-11]=1
for i in [4, 6, 11, 15, 16, 20, 21]:
    A[i, 27]=1
A[22, 27], A[23, 27]= -1, -1
for i in [22, 27, 28, 33, 35, 39]:
    A[i, 28]=1

#Matrix G

G=np.zeros(40)

G[0: 5] = L4*H*np.array([ho, lambda_ext/w_mur, lambda_isol/w_isol, lambda_int/w_mur, hi]) #wall 4
G[5: 7] = Lw1*H*np.array([ho, hi+lambda_wind/w_wind]) #window 1
G[7: 12] = L3*H*np.array([ho, lambda_ext/w_mur, lambda_isol/w_isol, lambda_int/w_mur, hi]) #wall 3
G[12: 16] = L2*H*np.array([hi+lambda_air/w_mur, lambda_isol/w_isol, lambda_int/w_mur, hi]) #wall 2
G[16] = Ld*H*(2*hi+lambda_door/w_door) #main door
G[17: 21] = L3*H*np.array([ho+lambda_ext/w_mur, lambda_isol/w_isol, lambda_int/w_mur, hi]) #wall 1
G[21] = rhoa*ca*ACH/3600*L6*(Lw2-2*(2*w_mur+w_isol))*H #ventilation
G[22] = Ld*H*(2*hi+lambda_door/w_door) #door
G[23: 28] = Li*H*np.array([hi, lambda_air/w_mur, lambda_isol/w_isol, lambda_air/w_mur, hi]) #wall inside
G[28: 33] = L5*H*np.array([hi, lambda_ext/w_mur, lambda_isol/w_isol, lambda_int/w_mur, ho]) #wall 4
G[33: 35] = Lw2*H*np.array([hi+lambda_wind/w_wind, ho]) #window 2
G[35: 39] = L6*H*np.array([hi, lambda_int/w_mur, lambda_isol/w_isol, hi+lambda_air/w_mur]) #wall 2
G[39] = 10**5 #HVAC

G=np.diag(G)

#Matrix C

C=np.zeros(29)

C=np.array([0, C_layer_out, C_layer_in, 0, 0, 0, Glass, C_layer_out, 0, C_layer_out, C_layer_in, 0, C_layer_out, C_layer_in, 0, 0, C_layer_in, C_layer_in, 0, C_layer_in, C_layer_out, 0, Glass, 0, C_layer_in, C_layer_out, Air, Air])

#Matrix b 

b=np.zeros(40)
b[5], b[7], b[21], b[32], b[34] = To, To, To, To, To
b[16], b[17], b[38] = Ti, Ti, Ti
b[39] = Th

#Matrix f 

f=np.zeros(29)
f[0]=H*L4*alpham*E
f[5]=H*L3*alpham*E
f[22]=H*L5*alpham*E

f[3]=H*L4*E
f[8]=H*L3*E
f[11]=H*L2*E
f[14]=H*L1*E
f[15]=H*Li*E
f[18]=H*Li*E
f[19]=H*L5*E
f[24]=H*L6*E
f[27]=L2*(L1+Ld)*E
f[28]=L2*L6*E

f[4]=H*Lw1*alphaw*E
f[23]=H*Lw2*alphaw*E


#Matrix y

y=np.zeros(29)
y[27]=y[28]=1


# Temperature vector in steady-state
θ = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)

# Flow-rate vector in steady-state
q = G @ (-A @ θ + b)
