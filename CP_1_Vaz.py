from unicodedata import ucd_3_2_0
from matplotlib import pyplot as plt
import numpy as np
from numpy import linalg as lg
import math


R = 100


N = R**2
U = np.zeros((N,N)) 
x = np.linspace(0,1,N)
dx = dy = 0.05
y = np.linspace(0,1,N)
F = np.zeros((R,R))
L = 3#0.005
c = N/2

# -Ad(i+1,j) + Bd(i,j) + Cd(i,j+1) + Cd(i,j-1) + Ad(i-1,j)
A = ((L**2)/dx**2) #-L²/▲x²
B = 1 + (2*(L**2)/dx**2) + (2*(L**2)/dy**2) #1 + 2L²/▲x² + 2L²/▲y²
C = ((L**2)/dy**2) #-L²/▲y²

### Laplacian
for i in range(N):
    for j in range(N):
        if i == j and i>R and i<N-R-1:
            U[i,j] = B
            if i < N and i%R!=0:
                U[i-1,i] = -A
            if i <= N and j%R!= 0:  
                U[i,i-1] = -A
            if i > 2: 
                U[i,i-R] = -C
            if i <= N-4:
                U[i,i+R] = -C

### CL borde gradient null
for i in range(N):
    if i < R: 
        U[i,:] = 0
        U[i,i] = 1
        U[i,i + R] = -1
    if i%R == 0 and i>R:
        U[i-1,:] = 0
        U[i-1,i-2] = -1
        U[i-1,i-1] = 1
    if i%R == 1 and i>=R and i<=N-R:
        U[i-1,:] = 0
        U[i-1,i-1] = 1
        U[i-1,i] = -1
    if i>= N-R:
        U[i,:] = 0
        U[i,i] = 1
        U[i,i - R] = -1

a,b = R/2, R/2

radius1 = R/5
radius2 = R/10
H = np.zeros((R,R))

### Esphere
for x in range(0,R):
    for y in range(0,R):
        if (((a-0.5)-x)**2)/radius2**2 + (((b-0.5)-y)**2)/radius1**2 - 1<=0:
            H[x,y] = 1
            # H = H.transpose()
            # F[x][y] = x * R + y + 1
            # print(x*R+y+1)
            p = (x*R+y)
            U[p,:] = 0
            U[p,p] = 1
            



Ures = H.reshape(R**2,1)
# U2 = np.linalg.inv(U)*Ures
U2 = np.linalg.solve(U,Ures)
U2 = U2.reshape(R,R)

im = plt.imshow(U2, cmap="coolwarm")
plt.colorbar(im)
plt.show()
