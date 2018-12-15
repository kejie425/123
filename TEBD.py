def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm
def Uni(tau):
    H = np.array([[h,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,-h]])
    w,v = np.linalg.eig(H)
    u = np.diag(np.exp(-w*tau))
    U = np.dot(v,np.dot(u,np.linalg.inv(v)))
    U = np.reshape(U,(2,2,2,2))
    return U
def converge1(v1,v2):
    norm = np.linalg.norm(v1-v2)
    print(norm)
    if norm <= 0.0001:
        return True
    else:
        return False

def converge2(v1,v2):
    norm = np.linalg.norm(v1-v2)
    print(norm)
    if norm <= 0.00001:
        return True
    else:
        return False

def converge3(v1,v2):
    norm = np.linalg.norm(v1-v2)
    print(norm)
    if norm <= 0.000001:
        return True
    else:
        return False

def converge4(v1,v2):
    norm = np.linalg.norm(v1-v2)
    print(norm)
    if norm <= 0.0000001:
        return True
    else:
        return False

def iTEBD_onestep1(tau,G,l):
        U = Uni(a*tau/2)
        A = 0
        B = 1
        theta = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(1,1))
        theta = np.tensordot(theta,np.diag(l[A,:]),axes=(2,0))
        theta = np.tensordot(theta,G[B,:,:,:],axes=(2,1))
        theta = np.tensordot(theta,np.diag(l[B,:]),axes=(3,0))
        theta = np.tensordot(theta,U,axes=([1,2],[0,1]))
        theta = np.reshape(np.transpose(theta,(2,0,3,1)),(d*D,d*D))
        theta = normalize(theta)
        X,Y,Z = np.linalg.svd(theta)
        Z = np.transpose(Z)
        l[A,0:D] = normalize(Y[0:D])
        X = np.reshape(X[0:d*D,0:D],(d,D,D))
        G[A,:,:,:] = np.transpose(np.tensordot(np.diag(l[B,:]**-1),X,axes=(1,1)),(1,0,2))
        Z = np.transpose(np.reshape(Z[0:d*D,0:D],(d,D,D)),(0,2,1))
        G[B,:,:,:] = np.tensordot(Z,np.diag(l[B,:]**-1),axes=(2,0))
        return G,l
def iTEBD_onestep2(tau,G,l):
        U = Uni(a*tau)
        A = 1
        B = 0
        theta = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(1,1))
        theta = np.tensordot(theta,np.diag(l[A,:]),axes=(2,0))
        theta = np.tensordot(theta,G[B,:,:,:],axes=(2,1))
        theta = np.tensordot(theta,np.diag(l[B,:]),axes=(3,0))
        theta = np.tensordot(theta,U,axes=([1,2],[0,1]))
        theta = np.reshape(np.transpose(theta,(2,0,3,1)),(d*D,d*D))
        theta = normalize(theta)
        X,Y,Z = np.linalg.svd(theta)
        Z = np.transpose(Z)
        l[A,0:D] = normalize(Y[0:D])
        X = np.reshape(X[0:d*D,0:D],(d,D,D))
        G[A,:,:,:] = np.transpose(np.tensordot(np.diag(l[B,:]**-1),X,axes=(1,1)),(1,0,2))
        Z = np.transpose(np.reshape(Z[0:d*D,0:D],(d,D,D)),(0,2,1))
        G[B,:,:,:] = np.tensordot(Z,np.diag(l[B,:]**-1),axes=(2,0))
        return G,l
def iTEBD_onestep3(tau,G,l):
        U = Uni((1-a)*tau/2)
        A = 0
        B = 1
        theta = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(1,1))
        theta = np.tensordot(theta,np.diag(l[A,:]),axes=(2,0))
        theta = np.tensordot(theta,G[B,:,:,:],axes=(2,1))
        theta = np.tensordot(theta,np.diag(l[B,:]),axes=(3,0))
        theta = np.tensordot(theta,U,axes=([1,2],[0,1]))
        theta = np.reshape(np.transpose(theta,(2,0,3,1)),(d*D,d*D))
        theta = normalize(theta)
        X,Y,Z = np.linalg.svd(theta)
        Z = np.transpose(Z)
        l[A,0:D] = normalize(Y[0:D])
        X = np.reshape(X[0:d*D,0:D],(d,D,D))
        G[A,:,:,:] = np.transpose(np.tensordot(np.diag(l[B,:]**-1),X,axes=(1,1)),(1,0,2))
        Z = np.transpose(np.reshape(Z[0:d*D,0:D],(d,D,D)),(0,2,1))
        G[B,:,:,:] = np.tensordot(Z,np.diag(l[B,:]**-1),axes=(2,0))
        return G,l
def iTEBD_onestep4(tau,G,l):
        U = Uni((1-2*a)*tau)
        A = 1
        B = 0
        theta = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(1,1))
        theta = np.tensordot(theta,np.diag(l[A,:]),axes=(2,0))
        theta = np.tensordot(theta,G[B,:,:,:],axes=(2,1))
        theta = np.tensordot(theta,np.diag(l[B,:]),axes=(3,0))
        theta = np.tensordot(theta,U,axes=([1,2],[0,1]))
        theta = np.reshape(np.transpose(theta,(2,0,3,1)),(d*D,d*D))
        theta = normalize(theta)
        X,Y,Z = np.linalg.svd(theta)
        Z = np.transpose(Z)
        l[A,0:D] = normalize(Y[0:D])
        X = np.reshape(X[0:d*D,0:D],(d,D,D))
        G[A,:,:,:] = np.transpose(np.tensordot(np.diag(l[B,:]**-1),X,axes=(1,1)),(1,0,2))
        Z = np.transpose(np.reshape(Z[0:d*D,0:D],(d,D,D)),(0,2,1))
        G[B,:,:,:] = np.tensordot(Z,np.diag(l[B,:]**-1),axes=(2,0))
        return G,l
def iTEBD_onestep5(tau,G,l):
        U = Uni((1-a)*tau/2)
        A = 0
        B = 1
        theta = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(1,1))
        theta = np.tensordot(theta,np.diag(l[A,:]),axes=(2,0))
        theta = np.tensordot(theta,G[B,:,:,:],axes=(2,1))
        theta = np.tensordot(theta,np.diag(l[B,:]),axes=(3,0))
        theta = np.tensordot(theta,U,axes=([1,2],[0,1]))
        theta = np.reshape(np.transpose(theta,(2,0,3,1)),(d*D,d*D))
        theta = normalize(theta)
        X,Y,Z = np.linalg.svd(theta)
        Z = np.transpose(Z)
        l[A,0:D] = normalize(Y[0:D])
        X = np.reshape(X[0:d*D,0:D],(d,D,D))
        G[A,:,:,:] = np.transpose(np.tensordot(np.diag(l[B,:]**-1),X,axes=(1,1)),(1,0,2))
        Z = np.transpose(np.reshape(Z[0:d*D,0:D],(d,D,D)),(0,2,1))
        G[B,:,:,:] = np.tensordot(Z,np.diag(l[B,:]**-1),axes=(2,0))
        return G,l
def iTEBD_onestep6(tau,G,l):
        U = Uni(a*tau)
        A = 1
        B = 0
        theta = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(1,1))
        theta = np.tensordot(theta,np.diag(l[A,:]),axes=(2,0))
        theta = np.tensordot(theta,G[B,:,:,:],axes=(2,1))
        theta = np.tensordot(theta,np.diag(l[B,:]),axes=(3,0))
        theta = np.tensordot(theta,U,axes=([1,2],[0,1]))
        theta = np.reshape(np.transpose(theta,(2,0,3,1)),(d*D,d*D))
        theta = normalize(theta)
        X,Y,Z = np.linalg.svd(theta)
        Z = np.transpose(Z)
        l[A,0:D] = normalize(Y[0:D])
        X = np.reshape(X[0:d*D,0:D],(d,D,D))
        G[A,:,:,:] = np.transpose(np.tensordot(np.diag(l[B,:]**-1),X,axes=(1,1)),(1,0,2))
        Z = np.transpose(np.reshape(Z[0:d*D,0:D],(d,D,D)),(0,2,1))
        G[B,:,:,:] = np.tensordot(Z,np.diag(l[B,:]**-1),axes=(2,0))
        return G,l
def iTEBD_onestep7(tau,G,l):
        U = Uni(a*tau/2)
        A = 0
        B = 1
        theta = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(1,1))
        theta = np.tensordot(theta,np.diag(l[A,:]),axes=(2,0))
        theta = np.tensordot(theta,G[B,:,:,:],axes=(2,1))
        theta = np.tensordot(theta,np.diag(l[B,:]),axes=(3,0))
        theta = np.tensordot(theta,U,axes=([1,2],[0,1]))
        theta = np.reshape(np.transpose(theta,(2,0,3,1)),(d*D,d*D))
        theta = normalize(theta)
        X,Y,Z = np.linalg.svd(theta)
        Z = np.transpose(Z)
        l[A,0:D] = normalize(Y[0:D])
        X = np.reshape(X[0:d*D,0:D],(d,D,D))
        G[A,:,:,:] = np.transpose(np.tensordot(np.diag(l[B,:]**-1),X,axes=(1,1)),(1,0,2))
        Z = np.transpose(np.reshape(Z[0:d*D,0:D],(d,D,D)),(0,2,1))
        G[B,:,:,:] = np.tensordot(Z,np.diag(l[B,:]**-1),axes=(2,0))
        return G,l

def iTEBD(tau,G,l):
    G,l = iTEBD_onestep1(tau,G,l)
    G,l = iTEBD_onestep2(tau,G,l)
    G,l = iTEBD_onestep3(tau,G,l)
    G,l = iTEBD_onestep4(tau,G,l)
    G,l = iTEBD_onestep5(tau,G,l)
    G,l = iTEBD_onestep6(tau,G,l)
    G,l = iTEBD_onestep7(tau,G,l)
    return G,l

import numpy as np
import math
h_array = np.array([3,1.1,1.01,1.001])
tau_array = np.array([0.1,0.01,0.001,0.0001])
D_array = np.array([50,100,150,200,250])
d = 2
D = D_array[2]
h = h_array[3]
G0 = np.random.rand(2,d,D,D)
l0 = np.random.rand(2,D)
a = 1/(2-pow(2,1/3))
G1 = G0.copy()
l1 = l0.copy()
G2,l2 = iTEBD(tau_array[0],G0,l0)
while not converge1(l1,l2):
    G1 = G2.copy()
    l1 = l2.copy()
    G2,l2 = iTEBD(tau_array[0],G2,l2)
G1 = G2.copy()
l1 = l2.copy()
G2,l2 = iTEBD(tau_array[1],G2,l2)
while not converge1(l1,l2):
    G1 = G2.copy()
    l1 = l2.copy()
    G2,l2 = iTEBD(tau_array[1],G2,l2)
G1 = G2.copy()
l1 = l2.copy()
G2,l2 = iTEBD(tau_array[2],G2,l2)
while not converge2(l1,l2):
    G1 = G2.copy()
    l1 = l2.copy()
    G2,l2 = iTEBD(tau_array[2],G2,l2)
G1 = G2.copy()
l1 = l2.copy()
G2,l2 = iTEBD(tau_array[3],G2,l2)
while not converge3(l1,l2):
    G1 = G2.copy()
    l1 = l2.copy()
    G2,l2 = iTEBD(tau_array[3],G2,l2)


def total_energy1(G,l):
    A = 0
    B = 1
    H = np.array([[h,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,-h]])
    H = np.reshape(H,(2,2,2,2))
    theta1 = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(0,1))
    theta1 = np.tensordot(theta1,np.diag(l[A,:]),axes=(2,0))
    theta1 = np.tensordot(theta1,G[B,:,:,:],axes=(2,1))
    theta1 = np.tensordot(theta1,np.diag(l[B,:]),axes=(3,0))
    theta1 = np.tensordot(theta1,H,axes=([1,2],[0,1]))
    theta3 = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(0,1))
    theta3 = np.tensordot(theta3,np.diag(l[A,:]),axes=(2,0))
    theta3 = np.tensordot(theta3,G[B,:,:,:],axes=(2,1))
    theta3 = np.tensordot(theta3,np.diag(l[B,:]),axes=(3,0))
    theta1 = np.tensordot(theta1,theta3,axes=([0,2,3,1],[0,1,2,3]))
    theta2 = np.tensordot(theta3,theta3,axes=([0,1,2,3],[0,1,2,3]))
    return theta1/theta2

def total_energy2(G,l):
    A = 1
    B = 0
    H = np.array([[h,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,-h]])
    H = np.reshape(H,(2,2,2,2))
    theta1 = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(0,1))
    theta1 = np.tensordot(theta1,np.diag(l[A,:]),axes=(2,0))
    theta1 = np.tensordot(theta1,G[B,:,:,:],axes=(2,1))
    theta1 = np.tensordot(theta1,np.diag(l[B,:]),axes=(3,0))
    theta1 = np.tensordot(theta1,H,axes=([1,2],[0,1]))
    theta3 = np.tensordot(np.diag(l[B,:]),G[A,:,:,:],axes=(0,1))
    theta3 = np.tensordot(theta3,np.diag(l[A,:]),axes=(2,0))
    theta3 = np.tensordot(theta3,G[B,:,:,:],axes=(2,1))
    theta3 = np.tensordot(theta3,np.diag(l[B,:]),axes=(3,0))
    theta1 = np.tensordot(theta1,theta3,axes=([0,2,3,1],[0,1,2,3]))
    theta2 = np.tensordot(theta3,theta3,axes=([0,1,2,3],[0,1,2,3]))
    return theta1/theta2
    

    

#np.savetxt('data_h1.01.txt',l2[0]**2,delimiter=',')

#left

def tensor1(Gx,l):
    n0 = np.tensordot(np.diag(l[1,:]),Gx[0,:,:,:],axes=(0,1))
    n0 = np.tensordot(n0,np.diag(l[0,:]),axes=(2,0))
    return n0

#right1
def tensor21(Gx,l):
    n0 = np.tensordot(Gx[0,:,:,:],np.diag(l[0,:]),axes = (2,0))
    return n0

#right1
def tensor22(Gx,l):
    n0 = np.tensordot(Gx[1,:,:,:],np.diag(l[1,:]),axes = (2,0))
    return n0

def tensor31(Gx,l):
    n0 = np.tensordot(Gx[0,:,:,:],np.diag(l[0,:]),axes=(2,0))
    return n0

def tensor32(Gx,l):
    n0 = np.tensordot(Gx[1,:,:,:],np.diag(l[1,:]),axes=(2,0))
    return n0

def tensor4(Gx,Gy,l):
    n0 = np.tensordot(np.diag(l[1,:]),Gx[0,:,:,:],axes=(0,1))
    n0 = np.tensordot(n0,np.diag(l[0,:]),axes=(2,0))
    n0 = np.tensordot(n0,Gy[1,:,:,:],axes=(2,1))
    n0 = np.tensordot(n0,np.diag(l[1,:]),axes=(3,0))
    return n0

def tensor5(Gx,l):
    n0 = np.tensordot(np.diag(l[0,:]),Gx[0,:,:,:],axes = (0,1))
    n0 = np.tensordot(n0,np.diag(l[1,:]),axes = (2,0))
    return(n0)
    
def correlation(N,G,l):
    Correlation = []
    G0 = G.copy()
    G[0,0,:,:] = -G0[0,0,:,:]
    G[1,0,:,:] = -G0[1,0,:,:]
    n0 = tensor5(G0,l)
    nor = np.tensordot(n0,n0,axes = ([0,1,2],[0,1,2]))
    m0 = tensor5(G,l)
    m0 = np.tensordot(n0,m0,axes = ([0,1,2],[0,1,2]))
    basic = (m0/nor)**2
    m1 = tensor4(G,G,l)
    n1 = tensor4(G0,G0,l)
    m1 = np.tensordot(n1,m1,axes = ([0,1,2,3],[0,1,2,3]))
    m1 = m1-basic
    Correlation.append(m1)
    n2 = tensor1(G0,l)
    s2 = np.tensordot(n2,n2,axes = ([0,1],[0,1]))
    m2 = tensor1(G,l)
    p2 = np.tensordot(n2,m2,axes = ([0,1],[0,1]))
    for i in range(N-1):
        if i%2 ==0:
            n3 = tensor32(G0,l)
            s3 = np.tensordot(n3,n3,axes = (0,0))
            s2 = np.tensordot(s2,s3,axes = ([0,1],[0,2]))
            m3 = tensor32(G0,l)
            p3 = np.tensordot(n3,m3,axes = (0,0))
            p2 = np.tensordot(p2,p3,axes = ([0,1],[0,2]))
        if i%2 ==1:
            n3 = tensor31(G0,l)
            s3 = np.tensordot(n3,n3,axes = (0,0))
            s2 = np.tensordot(s2,s3,axes = ([0,1],[0,2]))
            m3 = tensor31(G0,l)
            p3 = np.tensordot(n3,m3,axes = (0,0))
            p2 = np.tensordot(p2,p3,axes = ([0,1],[0,2]))
        if i%2 ==0:
            n4 = tensor21(G0,l)
            s4 = np.tensordot(n4,n4,axes = ([0,2],[0,2]))
            nor = np.tensordot(s2,s4,axes = ([0,1],[0,1]))
            m4 = tensor21(G,l)
            p4 = np.tensordot(n4,m4,axes = ([0,2],[0,2]))
            cor = np.tensordot(p2,p4,axes = ([0,1],[0,1]))           
        if i%2 ==1:
            n4 = tensor22(G0,l)
            s4 = np.tensordot(n4,n4,axes = ([0,2],[0,2]))
            nor = np.tensordot(s2,s4,axes = ([0,1],[0,1]))
            m4 = tensor22(G,l)
            p4 = np.tensordot(n4,m4,axes = ([0,2],[0,2]))
            cor = np.tensordot(p2,p4,axes = ([0,1],[0,1]))
        Correlation.append(cor/nor-basic)
        print(i,cor/nor-basic)
    return Correlation
#Corr = correlation(500,G2,l2)
#np.savetxt('data_h2_100_cor.txt',Corr,delimiter=',')
energy = (total_energy1(G2,l2)+total_energy2(G2,l2))/2
print(energy)

def magnetic(G,l):
    Correlation = []
    G0 = G.copy()
    G[0,0,:,:] = G0[0,1,:,:]
    G[0,1,:,:] = G0[0,0,:,:]
    G[1,0,:,:] = G0[1,1,:,:]
    G[1,1,:,:] = G0[1,0,:,:] 
    n0 = tensor5(G0,l)
    n1 = n0.copy()
    n1 = np.conjugate(n0)
    nor = np.tensordot(n1,n0,axes = ([0,1,2],[0,1,2]))
    m0 = tensor5(G,l)
    m = np.tensordot(n1,m0,axes = ([0,1,2],[0,1,2]))
    return m/nor

m0 = magnetic(G2,l2)
#mz.append(magnetic(G2,l2))
print(m0)

