import numpy as np
#Variables of the model:
nq, nθ = 24, 17 # number of flow-rates branches and of temperaure nodes
H = 3 # height; in m
w_wall, w_insul, w_wind, w_door = 0.2, 0.1, 0.025, 0.05 # width of walls, insulation layers, windows and doors; in m
L_wind, L_door, L1, L2, L3, L4, L5, L6 = 4, 1, 4, 3, 2.5, 3, 2, 6 # length of walls, windows and doors; in m
S_wind, S_door, S1, S2, S3, S4, S5, S6 = H*Lwind, H*L_door, H*L1, H*L2, H*L3, H*L4, H*L5, H*L6  #surface of walls, windows and doors; in m^2
To, Ti = -5, 18 # outside and inside temperatures; in °C
E = 200 # solar irradiation; in W/m^2
P_int=1000+200 # power produced in the house by two people and a radiator ; in W
ho, hi = 25, 8 # outside and inside convection coefficients; in W/m2/K
lambda_wind, lambda_wall, lambda_insul, lambda_door   = 1.015, 0.22, 0.023, 0.17 #thermal conductivities of the windows, the walls, the insulation layers and the doors; in W/m/K
alphawall_o, alphawall_i = 0.55, 1 #outside and inside absorption coefficient for the walls
rho_wall, rho_isol, rho_air = 2300, 16, 1.2 # densities of the walls, the insulation material and the air ; in kg/m^3
c_wall, c_isol, c_air = 880, 1210, 1000 # specific heats of the walls, the insulation material and the air ; in J/(kg.K)
 

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
b[23]=-Ti
b[16]=b[18]=-To
 

#Matrix f
f=np.zeros(nθ)
f[11]=E*alphawall_o*S3
f[8]=E*alphawall_i*S3
f[16]=E*alphawall_o*S6
f[13]=E*alphawall_i*S6
f[7]=P_int
 

#Matrix C
C=np.zeros(nθ)
C[0]=rho_wall*c_wall*w_wall*S2
C[1]=rho_insul*c_insul*w_insul*S2
C[2]=rho_wall*c_wall*w_wall*S1
C[3]=rho_insul*c_insul*w_insul*S1
C[4]=rho_wall*c_wall*w_wall*S5
C[5]=rho_insul*c_insul*w_insul*S5
C[6]=rho_air*c_air*L1*L2*H+rho_wall*c_air*w_wall*S4
C[7]=rho_air*c_air*L6*L_wind*H
C[9]=rho_insul*c_insul*w_insul*S3
C[10]=rho_wall*c_wall*w_wall*S3
C[12]=rho_air*c_air*w_wind*S_wind
C[14]=rho_insul*c_insul*w_insul*S6
C[15]=rho_wall*c_wall*w_wall*S6
 

#Matrix y
y=np.zeros(nθ)
y[6]=y[7]=1


# Steady-state
θ = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
q = G @ (-A @ θ + b)

print(f"The temperature in the room 1 is : θ6 = {θ[6]:.2f} °C, and the temperature in the room 2 is :θ7 = {θ[7]:.2f} °C")
