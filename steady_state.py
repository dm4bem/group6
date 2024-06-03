#Bibliothèques:
import numpy as np

#Variables :
nq, nθ = 24, 17 # number of flow-rates branches and of temperaure nodes
ho, hi = 25, 8 # en W/m2/K
lambda_wind = 1.015 #conductivité thermique de la fenêtre en W/m/K
lambda_wall = 0.22 #conductivité thermique du mur; en W/m/K
lambda_isol = 0.023 #conductivité thermique de l'isolant ; en W/m/K
lambda_door = 0.17 #conductivité thermique de la porte ; en W/m/K
H = 3 # hauteur sous plafond; en m 
w_wall, w_isol, w_wind, w_door = 0.2, 0.1, 0.025, 0.05 #epaisseurs des murs, isolant, fenêtre et portes; en m
Lwind, L_door, L1, L2, L3, L4, L5, L6 = 4, 1, 4, 3, 2.5, 3, 2, 6 #largeur des murs, fenêtre et portes; en m
Swind, S_door, S1, S2, S3, S4, S5, S6 = H*Lwind, H*L_door, H*L1, H*L2, H*L3, H*L4, H*L5, H*L6  #surfaces des murs, fenêtre et portes; en m

To = -5 # température extérieure; en °C 
Ti = 18 # température intérieure; en °C
E = 200 # en W/m2
alphawall_o = 0.55 #coefficient d'absorbtion du mur pour du béton lisse
alphawall_i = 1 #coefficient d'absorbtion du mur intérieur

P_int=1000+200 #puissance produite dans la maison par deux personnes et un radiateur en W

#Matrix A
A=np.zeros([nq, nθ])
A[0][0]=A[1][1]=A[2][6]=A[3][2]=A[4][3]=A[5][6]=A[6][6]=A[7][4]=A[8][5]=A[9][6]=A[10][7]=A[11][7]=A[12][8]=A[13][9]=A[14][10]=A[15][11]=A[17][12]=A[19][13]=A[20][14]=A[21][15]=A[22][16]=1
A[1][0]=A[2][1]=A[4][2]=A[5][3]=A[8][4]=A[9][5]=A[10][6]=A[11][6]=A[12][7]=A[13][8]=A[14][9]=A[15][10]=A[16][11]=A[17][7]=A[18][12]=A[19][7]=A[20][13]=A[21][14]=A[22][15]=A[23][16]=-1


#Matrix G
G=np.zeros(nq)
#wall2
G[0]=(hi+lambda_wall/(2*w_wall))*S2
G[1]=(lambda_wall/(2*w_wall)+lambda_isol/(w_isol*2))*S2
G[2]=(lambda_isol/(w_isol*2)+hi)*S2
#wall1
G[3]=(hi+lambda_wall/(2*w_wall))*S1
G[4]=(lambda_wall/(2*w_wall)+lambda_isol/(w_isol*2))*S1
G[5]=(lambda_isol/(w_isol*2)+hi)*S1
#wall5
G[7]=(hi+lambda_wall/(2*w_wall))*S5
G[8]=(lambda_wall/(2*w_wall)+lambda_isol/(w_isol*2))*S5
G[9]=(lambda_isol/(w_isol*2)+hi)*S5
#wall3
G[16]=ho*S3
G[15]=(lambda_wall/(w_wall*2))*S3
G[14]=(lambda_wall/(2*w_wall)+lambda_isol/(w_isol*2))*S3
G[13]=(lambda_isol/(w_isol*2))*S3
G[12]=hi*S3
#wall6
G[23]=ho*S6
G[22]=(lambda_wall/(w_wall*2))*S6
G[21]=(lambda_wall/(2*w_wall)+lambda_isol/(w_isol*2))*S6
G[20]=(lambda_isol/(w_isol*2))*S6
G[19]=hi*S6
#wall4
G[11]=(2*hi+lambda_wall/w_wall)*S4
#door1
G[6]=(2*hi+lambda_door/w_door)*S_door
#door2
G[10]=(2*hi+lambda_door/w_door)*S_door
#window
G[17]=(hi+lambda_wind/(w_wind*2))*Swind
G[18]=(ho+lambda_wind/(w_wind*2))*Swind
G=np.diag(G)


#Matrix b 
b=np.zeros(nq)
b[0]=b[3]=b[6]=b[7]=Ti
b[16]=b[18]=b[23]=-To

#Matrix f 
f=np.zeros(nθ)
f[11]=E*alphawall_o*S3
f[8]=E*alphawall_i*S3
f[16]=E*alphawall_o*S6
f[13]=E*alphawall_i*S6
f[7]=P_int


# Temperature vector in steady-state
θ = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)

# Flow-rate vector in steady-state
q = G @ (-A @ θ + b)


print(f"The temperature in the room 1 is : θ6 = {θ[6]:.2f} °C, and the temperature in the room 2 is :θ7 = {θ[7]:.2f} °C")