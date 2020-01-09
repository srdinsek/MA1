import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import simplify
from random import random
from scipy.optimize import minimize
from scipy.optimize import leastsq
from scipy import linalg
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from numpy import linalg as LA
from scipy.linalg import solve
from timeit import default_timer as timer
from random import random
from secrets import randbelow
from random import SystemRandom
from scipy.stats import chisquare, norm, kstest, chi2, moment, ks_2samp
from scipy.special import erf
from mpl_toolkits.mplot3d import axes3d
from scipy.signal import gaussian
from winsound import Beep





import seaborn as sns
sns.set_palette(sns.color_palette("jet", 6))


#%% 
###########################################################################################################
#
#       1. NALOGA - VERIŽICA
#
###########################################################################################################

def verižica():
    N=input('Vnesi dolžino verižice:\t ')
    N=int(N)
    
    #predkorak
    verižica=[int(N*random()) if i>0 and i<(N-1) else 0 for i in range(N)]
    energija=[]
    print('verižice je:\t',verižica)
    
    alpha=input('Povejte svojo vrednost konstante $\alpha$:')
    kb=input('Povejte svojo vrednost konstante $k_{b}$:')
    T=input('Povejte svojo vrednost temperature T:')
    alpha=float(alpha)
    kb=float(kb)
    T=float(T)
    
    K=input('Sedaj pa mi zaupajte še število korakov')
    K=int(K)
    
    for l in range(K):
        #1.korak
        E1=sum(-alpha*verižica[i] + 0.5*(verižica[i+1]-verižica[i])*(verižica[i+1]-verižica[i]) for i in range(N-1))
        točka=int((N-2)*(random()))
        naki=random()
        gordol=np.sign(1-2*naki)
        verižica2=[verižica[i] for i in range(len(verižica))]
        verižica2[točka+1]=verižica2[točka+1]+gordol
        E2=sum(-alpha*verižica2[i] + 0.5*(verižica2[i+1]-verižica2[i])*(verižica2[i+1]-verižica2[i]) for i in range(N-1))
        
        energija.append(E1)
        
        #2.korak
        if random() <= np.exp(-(E2-E1)/(kb*T)):
            verižica=[verižica2[i] for i in range(len(verižica2))]
    
    
    
    
    plt.figure(1)
    plt.title('Burn-in za energijo alfa={0} kT={1}'.format(alpha,kb*T))
    plt.plot([i for i in range(K)],energija,alpha=0.83)
    plt.ylabel('Energija')
    plt.xlabel('Korak')
    
    plt.figure(2)
    plt.title('Verižica na koncu alfa={0} kT={1}'.format(alpha,kb*T))
    plt.plot([i for i in range(N)],[-verižica[i] for i in range(N)],'r.')
    plt.plot([i for i in range(N)],[-verižica[i] for i in range(N)],'r:')
    plt.ylabel('Višina')
    plt.xlabel('Št. kroglice')


def veriža(N,K,alpha,kb,T):
    
    #predkorak
    verižica=[int(N*random()) if i>0 and i<(N-1) else 0 for i in range(N)]
    energija=[]
    


    
    for l in range(K):
        #1.korak
        E1=sum(-alpha*verižica[i] + 0.5*(verižica[i+1]-verižica[i])*(verižica[i+1]-verižica[i]) for i in range(N-1))
        točka=int((N-2)*(random()))
        naki=random()
        gordol=np.sign(1-2*naki)
        verižica2=[verižica[i] for i in range(len(verižica))]
        verižica2[točka+1]=verižica2[točka+1]+gordol
        E2=sum(-alpha*verižica2[i] + 0.5*(verižica2[i+1]-verižica2[i])*(verižica2[i+1]-verižica2[i]) for i in range(N-1))
        
        energija.append(E1)
        
        #2.korak
        if random() <= np.exp(-(E2-E1)/(kb*T)):
            verižica=[verižica2[i] for i in range(len(verižica2))]
    
    return verižica, energija

def risar(N,K,alpha,kb,T):
    
    plt.title('Verižica na koncu alfa={0}'.format(alpha))
    for i in range(T):
        
        toto=veriža(N,K,alpha,kb,i+0.0001)[0]
        plt.plot([j for j in range(N)],[-toto[i] for i in range(N)],'.')
        plt.plot([j for j in range(N)],[-toto[i] for i in range(N)],':',label='T={}'.format(i))
    
    
    plt.ylabel('Višina')
    plt.xlabel('Št. kroglice')
    plt.legend()



#%% 
#----------------------------------------------------------------------------------------------------------
#
#       1. NALOGA - VERIŽICA - omejitev dolžine
#
#----------------------------------------------------------------------------------------------------------

def verižica():
    N=input('Vnesi dolžino verižice:\t ')
    N=int(N)
    
    #predkorak
    verižica=[int(N*random()) if i>0 and i<(N-1) else 0 for i in range(N)]
    energija=[]
    print('verižice je:\t',verižica)
    
    alpha=input('Povejte svojo vrednost konstante $\alpha$:')
    kb=input('Povejte svojo vrednost konstante $k_{b}$:')
    T=input('Povejte svojo vrednost temperature T:')
    alpha=float(alpha)
    kb=float(kb)
    T=float(T)
    
    K=input('Sedaj pa mi zaupajte še število korakov')
    K=int(K)
    
    for l in range(K):
        #1.korak
        E1=sum(-alpha*verižica[i] + 0.5*(verižica[i+1]-verižica[i])*(verižica[i+1]-verižica[i]) for i in range(N-1))
        točka=int((N-2)*(random()))
        naki=random()
        gordol=np.sign(1-2*naki)
        verižica2=[verižica[i] for i in range(len(verižica))]
        verižica2[točka+1]=verižica2[točka+1]+gordol
        if verižica2[točka+1] <= 19:
            E2=sum(-alpha*verižica2[i] + 0.5*(verižica2[i+1]-verižica2[i])*(verižica2[i+1]-verižica2[i]) for i in range(N-1))
            
            energija.append(E1)
            
            #2.korak
            if random() <= np.exp(-(E2-E1)/(kb*T)):
                verižica=[verižica2[i] for i in range(len(verižica2))]
    
    
    
    
    plt.figure(1)
    plt.title('Burn-in za energijo alfa={0} kT={1}'.format(alpha,kb*T))
    plt.plot([i for i in range(len(energija))],energija,alpha=0.83)
    plt.ylabel('Energija')
    plt.xlabel('Korak')
    
    plt.figure(2)
    plt.title('Verižica na koncu alfa={0} kT={1}'.format(alpha,kb*T))
    plt.plot([i for i in range(N)],[-verižica[i] for i in range(N)],'r.')
    plt.plot([i for i in range(N)],[-verižica[i] for i in range(N)],'r:')
    plt.ylabel('Višina')
    plt.xlabel('Št. kroglice')


def veriža(N,K,alpha,kb,T):
    
    #predkorak
    verižica=[int(N*random()) if i>0 and i<(N-1) else 0 for i in range(N)]
    energija=[]
    


    
    for l in range(K):
        #1.korak
        E1=sum(-alpha*verižica[i] + 0.5*(verižica[i+1]-verižica[i])*(verižica[i+1]-verižica[i]) for i in range(N-1))
        točka=int((N-2)*(random()))
        naki=random()
        gordol=np.sign(1-2*naki)
        verižica2=[verižica[i] for i in range(len(verižica))]
        verižica2[točka+1]=verižica2[točka+1]+gordol
        if verižica2[točka+1] <= 19:
            E2=sum(-alpha*verižica2[i] + 0.5*(verižica2[i+1]-verižica2[i])*(verižica2[i+1]-verižica2[i]) for i in range(N-1))
            
            energija.append(E1)
            
            #2.korak
            if random() <= np.exp(-(E2-E1)/(kb*T)):
                verižica=[verižica2[i] for i in range(len(verižica2))]
        
    return verižica, energija

def risar(N,K,alpha,kb,T):
    
    plt.title('Verižica na koncu alfa={0}'.format(alpha))
    for i in range(T):
        
        toto=veriža(N,K,alpha,kb,i/7+0.0001)[0]
        plt.plot([j for j in range(N)],[-toto[j] for j in range(N)],'.')
        plt.plot([j for j in range(N)],[-toto[j] for j in range(N)],':',label='T={}'.format(i/7))
    
    
    plt.ylabel('Višina')
    plt.xlabel('Št. kroglice')
    plt.legend()

def risar2(N,K,alpha,kb,T):
    
    plt.title('Verižica na koncu T={0}'.format(T))
    for i in range(alpha):
        
        toto=veriža(N,K,i/7,kb,T)[0]
        plt.plot([j for j in range(N)],[-toto[j] for j in range(N)],'.')
        plt.plot([j for j in range(N)],[-toto[j] for j in range(N)],':',label='alpha={0:.2f}'.format(i/7))
    
    
    plt.ylabel('Višina')
    plt.xlabel('Št. kroglice')
    plt.legend()


#%% 
#----------------------------------------------------------------------------------------------------------
#
#       1. NALOGA - VERIŽICA - različne dolžine
#
#----------------------------------------------------------------------------------------------------------




def verižna(N,K,alpha,kb,T,h):
    
    #predkorak
    verižica=[int(h*random()) if i>0 and i<(N-1) else 0 for i in range(N)]
    energija=[]
    


    
    for l in range(K):
        #1.korak
        E1=sum(-alpha*verižica[i] + 0.5*(verižica[i+1]-verižica[i])*(verižica[i+1]-verižica[i]) for i in range(N-1))
        točka=int((N-2)*(random()))
        naki=random()
        gordol=np.sign(1-2*naki)
        verižica2=[verižica[i] for i in range(len(verižica))]
        verižica2[točka+1]=verižica2[točka+1]+gordol
        if verižica2[točka+1] <= h:
            E2=sum(-alpha*verižica2[i] + 0.5*(verižica2[i+1]-verižica2[i])*(verižica2[i+1]-verižica2[i]) for i in range(N-1))
            
            energija.append(E1)
            
            #2.korak
            if random() <= np.exp(-(E2-E1)/(kb*T)):
                verižica=[verižica2[i] for i in range(len(verižica2))]
        
    return verižica, energija

def risarH(N,K,alpha,kb,T,h):
    
    plt.title('Verižica na koncu alfa={0} t={1}'.format(alpha,T))
    for i in range(h):
        
        toto=verižna(N,K,alpha,kb,T,i*10)[0]
        plt.plot([j for j in range(N)],[-toto[j] for j in range(N)],'.')
        plt.plot([j for j in range(N)],[-toto[j] for j in range(N)],':',label='h={}'.format(i*10))
    
    
    plt.ylabel('Višina')
    plt.xlabel('Št. kroglice')
    plt.legend()

#%% 
#----------------------------------------------------------------------------------------------------------
#
#       1. NALOGA - VERIŽICA - ohlajanje
#
#----------------------------------------------------------------------------------------------------------

def verižaN(N,K,alpha,kb,T,h):
    
    #predkorak
    verižica=[int(h*random()) if i>0 and i<(N-1) else 0 for i in range(N)]
    E1=sum(-alpha*verižica[i] + 0.5*(verižica[i+1]-verižica[i])*(verižica[i+1]-verižica[i]) for i in range(N-1))
    
    energija=[i*10 for i in range(100)]
    E2=0
    znak=[-200]
        
    for z in range(K):
        #1.korak
        E1=sum(-alpha*verižica[i] + 0.5*(verižica[i+1]-verižica[i])*(verižica[i+1]-verižica[i]) for i in range(N-1))
        točka=int((N-2)*(random()))
        naki=random()
        gordol=np.sign(1-2*naki)
        verižica2=[verižica[i] for i in range(len(verižica))]
        verižica2[točka+1]=verižica2[točka+1]+gordol
        if verižica2[točka+1] <= h:
            E2=sum(-alpha*verižica2[i] + 0.5*(verižica2[i+1]-verižica2[i])*(verižica2[i+1]-verižica2[i]) for i in range(N-1))
                
            energija.append(E1)
                
            #2.korak
            if random() <= np.exp(-(E2-E1)/(kb*T)):
                verižica=[verižica2[i] for i in range(len(verižica2))]
        if np.abs(energija[-1]-energija[-100]) > 0.1 and T>0 and z-100 > znak[-1]:
            T-=0.01
            znak.append(z)
    
    
    
    plt.figure(1)
    plt.title('Burn-in za energijo alfa={0} ohlajanje'.format(alpha))
    plt.plot([i for i in range(len(energija)-100)],[energija[i+100] for i in range(len(energija)-100)],alpha=0.83)
    plt.ylabel('Energija')
    plt.xlabel('Korak')
    
    plt.figure(2)
    plt.title('Ohlajanje alfa={}'.format(alpha))
    plt.plot([i for i in range(N)],[-verižica[i] for i in range(N)],'r.')
    plt.plot([i for i in range(N)],[-verižica[i] for i in range(N)],'r:')
    plt.ylabel('Višina')
    plt.xlabel('Št. kroglice')
    
    print(energija[-1])

#%% 
###########################################################################################################
#
#       2. NALOGA - ISINGOV MODEL
#
###########################################################################################################
#J=1
#H=1
#N=100
#kb=1
#T=1

def ising(N,K,J,H,kb,T):
    #Definicijski del
    energija=[]
    matrika=[[np.sign(1-2*random()) for i in range(N)] for j in range(N)]
    E1=sum(-J*(matrika[i%N][j%N]*matrika[(i+1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j+1)%N]+matrika[i%N][j%N]*matrika[(i-1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j-1)%N]) + H*matrika[i][j] for i in range(N) for j in range(N))
    print(E1)
    
    #algoritmični del
    for l in range(K):
        #1.korak
        točki=[int(N*random()),int(N*random())]
        
        E2= E1 + 2*J*matrika[točki[0]][točki[1]]*(matrika[točki[0]%N][(točki[1]+1)%N]+matrika[točki[0]%N][(točki[1]-1)%N]+matrika[(točki[0]+1)%N][(točki[1])%N]+matrika[(točki[0]-1)%N][(točki[1])%N]) + 2*H*matrika[točki[0]][točki[1]]
            
        energija.append(E1)
            
        #2.korak
        if random() <= np.exp(-(E2-E1)/(kb*T)):
            matrika[točki[0]][točki[1]]=-matrika[točki[0]][točki[1]]
            E1=E2
    
    
    
    #risarski del
    matrix = np.matrix(matrika)
    matrix
    fig = plt.figure(0)
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.title('Mreža z minimalno energijo J={0} H={1} kbT={2}'.format(J,H,kb*T))
    plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.YlGn,
        extent=(0.5,np.shape(matrix)[1]+0.5,0.5,np.shape(matrix)[0]+0.5))
    #plt.colorbar()
    plt.show()
    
    fig=plt.figure(1)
    plt.title('Burn-in energija za J={0} H={1} kbT={2}'.format(J,H,kb*T))
    plt.plot([i for i in range(K)],energija, alpha=0.83)
    
    
    Beep(780, 1000)

#%% 
#----------------------------------------------------------------------------------------------------------
#
#       2. NALOGA - ISINGOV MODEL - ohlajanje
#
#----------------------------------------------------------------------------------------------------------

def isingohlajanje(N,K,J,H,kb,T):
    #Definicijski del
    energija=[]
    matrika=[[np.sign(1-2*random()) for i in range(N)] for j in range(N)]
    E1=sum(-J*(matrika[i%N][j%N]*matrika[(i+1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j+1)%N]+matrika[i%N][j%N]*matrika[(i-1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j-1)%N]) + H*matrika[i][j] for i in range(N) for j in range(N))
    print(E1)
    #algoritmični del
    for l in range(K):
        #1.korak
        točki=[int(N*random()),int(N*random())]
            
        E2= E1 + 2*J*matrika[točki[0]][točki[1]]*(matrika[točki[0]%N][(točki[1]+1)%N]+matrika[točki[0]%N][(točki[1]-1)%N]+matrika[(točki[0]+1)%N][(točki[1])%N]+matrika[(točki[0]-1)%N][(točki[1])%N]) + 2*H*matrika[točki[0]][točki[1]]
            
        energija.append(E1)
            
        #2.korak
        if random() <= np.exp(-(E2-E1)/(kb*T)):
            matrika[točki[0]][točki[1]]=-matrika[točki[0]][točki[1]]
            E1=E2
            
        if T>0.1001 and l%10000==0:
            T-=0.1
            
    
    
    
    #risarski del
    matrix = np.matrix(matrika)
    matrix
    fig = plt.figure(0)
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.title('Mreža z minimalno energijo J={0} H={1} kbT={2}'.format(J,H,kb*T))
    plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.YlGn,
        extent=(0.5,np.shape(matrix)[1]+0.5,0.5,np.shape(matrix)[0]+0.5))
    #plt.colorbar()
    plt.show()
    
    fig=plt.figure(1)
    plt.title('Burn-in energija za J={0} H={1} kbT={2}'.format(J,H,kb*T))
    plt.plot([i for i in range(K)],energija, alpha=0.83)
    
    
    Beep(780, 1000)

#%% 
#----------------------------------------------------------------------------------------------------------
#
#       2. NALOGA - ISINGOV MODEL - magnetizacija in energija
#
#----------------------------------------------------------------------------------------------------------
    
def isingf(N,K,J,H,kb,T):
    #Definicijski del
    energija=[i*1000 for i in range(10000005)]
    matrika=[[np.sign(1-2*random()) for i in range(N)] for j in range(N)]
    E1=sum(-J*(matrika[i%N][j%N]*matrika[(i+1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j+1)%N]+matrika[i%N][j%N]*matrika[(i-1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j-1)%N]) + H*matrika[i][j] for i in range(N) for j in range(N))
    #algoritmični del
    while np.abs(energija[-1]-energija[-10000000])> 10:
        #1.korak
        točki=[int(N*random()),int(N*random())]
            
        E2= E1 + 2*J*matrika[točki[0]][točki[1]]*(matrika[točki[0]%N][(točki[1]+1)%N]+matrika[točki[0]%N][(točki[1]-1)%N]+matrika[(točki[0]+1)%N][(točki[1])%N]+matrika[(točki[0]-1)%N][(točki[1])%N]) + 2*H*matrika[točki[0]][točki[1]]
            
        energija.append(E1)
            
        #2.korak
        if random() <= np.exp(-(E2-E1)/(kb*T)):
            matrika[točki[0]][točki[1]]=-matrika[točki[0]][točki[1]]
            E1=E2
    
    return energija, matrika

def povprečnaE(N,K,J,H,kb,T):
    tempič=[(T-i)/20 for i in range(T)]
    sezi=[]
    for i in tempič:
        eni,mat=isingf(N,K,J,H,kb,i)
        povprečje=sum(eni[-i*100] for i in range(200))/200
        sezi.append(povprečje)
    plt.title('Povprečna energija od temperature J={0} H={1}'.format(J,H))
    plt.plot(tempič,sezi,'b:')
    plt.ylabel('Povprečna energija')
    plt.xlabel('Temperatura')
    Beep(780, 1000)
    
def specifičnatoplota(N,K,J,H,kb,T):
    tempič=[(T-i)/5 for i in range(T)]
    sezi=[]
    for i in tempič:
        eni,mat=isingf(N,K,J,H,kb,i)
        
        povprečje=sum(eni[-j*10000] for j in range(200))/200
        povprečje2=sum(eni[-j*10000]*eni[-j*10000] for j in range(200))/200
        sezi.append((-povprečje*povprečje+povprečje2)/(kb*T*T*N*N))
        
    plt.title('Specifična toplota J={0} H={1}'.format(J,H))
    plt.plot(tempič,sezi,'r:')
    plt.ylabel('Specifična toplota')
    plt.xlabel('Temperatura')
    Beep(780, 1000)
    return sezi
    

def isingM(N,K,J,H,kb,T):
    #Definicijski del
    energija=[i*100000 for i in range(10000005)]
    matrika=[[np.sign(1-2*random()) for i in range(N)] for j in range(N)]
    E1=sum(-J*(matrika[i%N][j%N]*matrika[(i+1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j+1)%N]+matrika[i%N][j%N]*matrika[(i-1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j-1)%N]) + H*matrika[i][j] for i in range(N) for j in range(N))
    #algoritmični del
    mati=[]
    while np.abs(energija[-1]-energija[-10000000])> 50:
        #1.korak
        točki=[int(N*random()),int(N*random())]
            
        E2= E1 + 2*J*matrika[točki[0]][točki[1]]*(matrika[točki[0]%N][(točki[1]+1)%N]+matrika[točki[0]%N][(točki[1]-1)%N]+matrika[(točki[0]+1)%N][(točki[1])%N]+matrika[(točki[0]-1)%N][(točki[1])%N]) + 2*H*matrika[točki[0]][točki[1]]
        
        mati.append(matrika)
        energija.append(E1)
            
        #2.korak
        if random() <= np.exp(-(E2-E1)/(kb*T)):
            matrika[točki[0]][točki[1]]=-matrika[točki[0]][točki[1]]
            E1=E2
    
    return energija, mati  
    
    
    
    
    

def magnetizacija(N,K,J,H,kb,T):
    start=timer()
    tempič=[(T-i)/25+0.001 for i in range(T)]
    sezi=[]
    for i in tempič:
        eni,mat=isingM(N,K,J,H,kb,i)
        povprečje=sum(np.abs(sum(mat[-i*100000])) for i in range(100))/100
        sezi.append(povprečje)
    plt.title('Povprečna magnetizacija od temperature J={0} H={1}'.format(J,H))
    plt.plot(tempič,sezi,'b:')
    plt.plot(tempič,sezi,'b.')
    plt.ylabel('Povprečna magnetizacija')
    plt.xlabel('Temperatura')
    end=timer()
    print(end-start)
    Beep(780, 1000)
    
def subscetibilnost(N,K,J,H,kb,T):
    start=timer()
    tempič=[(T-i)/10+0.001 for i in range(T)]
    sezi=[]
    for i in tempič:
        eni,mat=isingM(N,K,J,H,kb,i)
        povprečje=sum(np.abs(sum(mat[-j*100000])) for j in range(100))/100
        povprečje2=sum(sum(mat[-j*100000])**2 for j in range(100))/100
        sezi.append((povprečje2-povprečje*povprečje)/(N*N*kb*i))
    plt.title('Spinska sumbscetibilnost od temperature J={0} H={1}'.format(J,H))
    plt.plot(tempič,sezi,'b:')
    plt.plot(tempič,sezi,'b.')
    plt.ylabel('Spinska subscetibilnost')
    plt.xlabel('Temperatura')
    end=timer()
    print(end-start)
    Beep(780, 1000)
    
#%%
#----------------------------------------------------------------------------------------------------------
#
#       2. NALOGA - ISINGOV MODEL - magnetizacija in subscetibilnost- ciganjenje
#
#----------------------------------------------------------------------------------------------------------
    
def isingM(N,K,J,H,kb,T):
    #Definicijski del
    energija=[i*100000 for i in range(50000005)]
    matrika=[[np.sign(1-2*random()) for i in range(N)] for j in range(N)]
    E1=sum(-J*(matrika[i%N][j%N]*matrika[(i+1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j+1)%N]+matrika[i%N][j%N]*matrika[(i-1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j-1)%N]) + H*matrika[i][j] for i in range(N) for j in range(N))
    #algoritmični del
    mati=[]
    sezi=[]
    sezi2=[]
    for k in range(50000000):
        #1.korak
        točki=[int(N*random()),int(N*random())]
            
        E2= E1 + 2*J*matrika[točki[0]][točki[1]]*(matrika[točki[0]%N][(točki[1]+1)%N]+matrika[točki[0]%N][(točki[1]-1)%N]+matrika[(točki[0]+1)%N][(točki[1])%N]+matrika[(točki[0]-1)%N][(točki[1])%N]) + 2*H*matrika[točki[0]][točki[1]]
        
        if k%10000==0 and k > 20000000:
            
            
            vsota=sum(matrika)
            povprečje=np.abs(vsota)
            povprečje2=vsota*vsota
            sezi.append(povprečje)
            sezi2.append(povprečje2)
        
        
        energija.append(E1)
            
        #2.korak
        if random() <= np.exp(-(E2-E1)/(kb*T)):
            matrika[točki[0]][točki[1]]=-matrika[točki[0]][točki[1]]
            E1=E2
        vsota=sum(sezi)
        dolzina=len(sezi)
    return energija, (sum(sezi2)/len(sezi2)-vsota*vsota/(dolzina*dolzina))/(N*N*kb*T) 
    
    
    
    
    

def magnetizacija(N,K,J,H,kb,T):
    start=timer()
    tempič=[(T-i)/10+0.001 for i in range(T)]
    sezi=[]
    for i in tempič:
        eni,mat=isingM(N,K,J,H,kb,i)
        sezi.append(mat)
    plt.title('Povprečna magnetizacija od temperature J={0} H={1}'.format(J,H))
    plt.plot(tempič,sezi,'b:')
    plt.plot(tempič,sezi,'b.')
    plt.ylabel('Povprečna magnetizacija')
    plt.xlabel('Temperatura')
    end=timer()
    print(end-start)
    Beep(780, 1000)


def spinskasubscetibilnost(N,K,J,H,kb,T):
    start=timer()
    tempič=[(T-i)/10+0.001 for i in range(T)]
    sezi=[]
    for i in tempič:
        eni,mat=isingM(N,K,J,H,kb,i)
        sezi.append(mat)
    plt.title('Spinska subscetibilnost J={0} H={1}'.format(J,H))
    plt.plot(tempič,sezi,'b:')
    plt.plot(tempič,sezi,'b.')
    plt.ylabel('Spinska subscetibilnost')
    plt.xlabel('Temperatura')
    end=timer()
    print(end-start)
    Beep(780, 1000)

#%%
###########################################################################################################
#
#       2. NALOGA - ISINGOV MODEL - ANIMACIJA
#
###########################################################################################################


def isingoh(N,K,J,H,kb,T):
    #Definicijski del
    energija=[]
    matrika=[[np.sign(1-2*random()) for i in range(N)] for j in range(N)]
    E1=sum(-J*(matrika[i%N][j%N]*matrika[(i+1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j+1)%N]+matrika[i%N][j%N]*matrika[(i-1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j-1)%N]) + H*matrika[i][j] for i in range(N) for j in range(N))
    print(E1)
    #algoritmični del
    for l in range(K):
        #1.korak
        točki=[int(N*random()),int(N*random())]
            
        E2= E1 + 2*J*matrika[točki[0]][točki[1]]*(matrika[točki[0]%N][(točki[1]+1)%N]+matrika[točki[0]%N][(točki[1]-1)%N]+matrika[(točki[0]+1)%N][(točki[1])%N]+matrika[(točki[0]-1)%N][(točki[1])%N]) + 2*H*matrika[točki[0]][točki[1]]
            
        energija.append(E1)
            
        #2.korak
        if random() <= np.exp(-(E2-E1)/(kb*T)):
            matrika[točki[0]][točki[1]]=-matrika[točki[0]][točki[1]]
            E1=E2
            
        if T>0.1001 and l%10000==0:
            T-=0.1
            
    return matrika
    
    


###################################################################
    #risarski del
    matrix = np.matrix(matrika)
    matrix
    fig = plt.figure(0)
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.title('Mreža z minimalno energijo J={0} H={1} kbT={2}'.format(J,H,kb*T))
    plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.YlGn,
        extent=(0.5,np.shape(matrix)[1]+0.5,0.5,np.shape(matrix)[0]+0.5))
    #plt.colorbar()
    plt.show()
###################################################################




matrix
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.YlGn,
        extent=(0.5,np.shape(matrix)[1]+0.5,0.5,np.shape(matrix)[0]+0.5))


from matplotlib import animation






animation.FuncAnimation(fig, animate )







# First set up the figure, the axis, and the plot element we want to animate

#line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    matrix=np.matrix(isingoh(100,1000000,1,0,1,10))
    line.set_data(matrix)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()







#%%

def isingM(N,K,J,H,kb,T):
    #Definicijski del
    matrika=[[np.sign(1-2*random()) for i in range(N)] for j in range(N)]
    E1=sum(-J*(matrika[i%N][j%N]*matrika[(i+1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j+1)%N]+matrika[i%N][j%N]*matrika[(i-1)%N][j%N]+matrika[i%N][j%N]*matrika[i%N][(j-1)%N]) + H*matrika[i][j] for i in range(N) for j in range(N))
    #algoritmični del
    mati=[]
    sezi=[]
    sezi2=[]
    eni=[]
    for k in range(50000000):
        #1.korak
        točki=[int(N*random()),int(N*random())]
            
        E2= E1 + 2*J*matrika[točki[0]][točki[1]]*(matrika[točki[0]%N][(točki[1]+1)%N]+matrika[točki[0]%N][(točki[1]-1)%N]+matrika[(točki[0]+1)%N][(točki[1])%N]+matrika[(točki[0]-1)%N][(točki[1])%N]) + 2*H*matrika[točki[0]][točki[1]]
        
        if k%10000==0 and k > 30000000:
                    
            sezi.append(sum(matrika))
            eni.append(E1)

            
        #2.korak
        if random() <= np.exp(-(E2-E1)/(kb*T)):
            matrika[točki[0]][točki[1]]=-matrika[točki[0]][točki[1]]
            E1=E2
        
    d=len(sezi)
    m=sum(sezi)
    m2=sum(sezi[i]*sezi[i] for i in range(d))
        
    n=len(eni)
    s=sum(eni)
    s2=sum(eni[i]*eni[i] for i in range(n))
    return s/n, (s2/n-s*s/(n*n))/(N*N*kb*T*T), np.abs(m)/d, (m2/d-m*m/(d*d))/(N*N*kb*T) 



def risivse(N,K,J,H,kb,T):
    start=timer()
    tempič=[(T-i)/15+1 for i in range(T)]
    energija=[]
    specifična=[]
    magnetizacija=[]
    subscetbilnost=[]
    for i in tempič:
        E,c,mag,ss=isingM(100,1,1,0,1,i)
        energija.append(E)
        specifična.append(c)
        magnetizacija.append(mag)
        subscetbilnost.append(ss)
    
    plt.figure(0)
    plt.title('Povprečna energija od temperature J={0} H={1}'.format(J,H))
    plt.plot(tempič,energija,'b:')
    plt.plot(tempič,energija,'b.')
    plt.ylabel('Povprečna energija')
    plt.xlabel('Temperatura')
    plt.figure(1)
    plt.title('Specifična toplota od temperature J={0} H={1}'.format(J,H))
    plt.plot(tempič,specifična,'r:')
    plt.plot(tempič,specifična,'r.')
    plt.ylabel('Specifična toplota')
    plt.xlabel('Temperatura')
    plt.figure(2)
    plt.title('Magnetizacija od temperature J={0} H={1}'.format(J,H))
    plt.plot(tempič,magnetizacija,'y:')
    plt.plot(tempič,magnetizacija,'y.')
    plt.ylabel('Magnetizacija')
    plt.xlabel('Temperatura')
    plt.figure(3)
    plt.title('Spinska subscetibilnost od temperature J={0} H={1}'.format(J,H))
    plt.plot(tempič,subscetbilnost,'g:')
    plt.plot(tempič,subscetbilnost,'g.')
    plt.ylabel('Spinska subscetibilnost')
    plt.xlabel('Temperatura')
    end=timer()
    print(end-start)
    Beep(780, 1000)
        
    
jet= plt.get_cmap('jet')
colors = iter(jet(np.linspace(0,1,24)))
    
def risivseMag(N,K,J,H,kb,T):
    start=timer()
    tempič=[i/5+2 for i in range(T)]
    spec=[]
    for h in range(H):
        energija=[]
        specifična=[]
        magnetizacija=[]
        subscetbilnost=[]
        for j in tempič:
            E,c,mag,ss=isingM(100,1,1,h/3+0.5,1,j)
            energija.append(E)
            specifična.append(c)
            magnetizacija.append(mag)
            subscetbilnost.append(ss)
        
        plt.figure(0)
        plt.title('Povprečna energija od temperature J={0}'.format(J))
        plt.plot(tempič,energija,':',color=next(colors))
        plt.plot(tempič,energija,'.',color=next(colors),label='H={}'.format(h+1))
        plt.ylabel('Povprečna energija')
        plt.xlabel('Temperatura')
        plt.legend()
        plt.figure(1)
        plt.title('Specifična toplota od temperature J={0}'.format(J))
        plt.plot(tempič,specifična,':',color=next(colors))
        plt.plot(tempič,specifična,'.',color=next(colors),label='H={}'.format(h+1))
        plt.ylabel('Specifična toplota')
        plt.xlabel('Temperatura')
        plt.legend()
        plt.figure(2)
        plt.title('Magnetizacija od temperature J={0}'.format(J))
        plt.plot(tempič,magnetizacija,':',color=next(colors))
        plt.plot(tempič,magnetizacija,'.',color=next(colors),label='H={}'.format(h+1))
        plt.ylabel('Magnetizacija')
        plt.xlabel('Temperatura')
        plt.legend()
        plt.figure(3)
        plt.title('Spinska subscetibilnost od temperature J={0}'.format(J))
        plt.plot(tempič,subscetbilnost,':',color=next(colors))
        plt.plot(tempič,subscetbilnost,'.',color=next(colors),label='H={}'.format(h+1))
        plt.ylabel('Spinska subscetibilnost')
        plt.xlabel('Temperatura')
        spec.append(specifična)
        plt.legend()
        
    end=timer()
    print(end-start)
    Beep(780, 1000)    
    return spec


spec=risivseMag(100,1,1,3,1,25)



def magnetizacija(N,K,J,H,kb,T):
    start=timer()
    tempič=[(T-i)/10+0.001 for i in range(T)]
    sezi=[]
    for i in tempič:
        eni,mat=isingM(N,K,J,H,kb,i)
        sezi.append(mat)
    plt.title('Povprečna magnetizacija od temperature J={0} H={1}'.format(J,H))
    plt.plot(tempič,sezi,'b:')
    plt.plot(tempič,sezi,'b.')
    plt.ylabel('Povprečna magnetizacija')
    plt.xlabel('Temperatura')
    end=timer()
    print(end-start)
    Beep(780, 1000)


def spinskasubscetibilnost(N,K,J,H,kb,T):
    start=timer()
    tempič=[(T-i)/10+0.001 for i in range(T)]
    sezi=[]
    for i in tempič:
        eni,mat=isingM(N,K,J,H,kb,i)
        sezi.append(mat)
    plt.title('Spinska subscetibilnost J={0} H={1}'.format(J,H))
    plt.plot(tempič,sezi,'b:')
    plt.plot(tempič,sezi,'b.')
    plt.ylabel('Spinska subscetibilnost')
    plt.xlabel('Temperatura')
    end=timer()
    print(end-start)
    Beep(780, 1000)

#%%
#%%
###########################################################################################################
#
#       3. NALOGA - KITAJSKI TRGOVEC - brez dodatnih omejitev
#
###########################################################################################################


#4-latitude in 5-longitude sta kota
# 1 mesto
# 2 provinca
# 3 refija

mesta=[]
string='kitajska.txt'
with open(string) as f:
    podaci=[l.strip().split("\t") for l in f]
mesta=[[podaci[i][j] if j < 4 else float(podaci[i][j])*np.pi/180 for j in range(6)] for i in range(37)]
print(mesta)

def razdalja(i,j):
    return np.arccos(np.sin(mesta[i][4])*np.sin(mesta[j][4])+np.cos(mesta[i][4])*np.cos(mesta[j][4])*np.cos(mesta[i][5]-mesta[j][5]))


#plt.figure(100)
#plt.plot([mesta[i][4] for i in range(37)],[mesta[i][5] for i in range(37)],'.')






def trgovec(N,K,kb,T):
    #Definicijski del
    seznam=[i for i in range(N)]
    lsq1=sum(razdalja(seznam[i],seznam[(i+1)%N]) for i in range(N))
    L=[]
    
    #algoritmični del
    for k in range(K):
        #1.korak
        točki=[int(N*random()),int(N*random())]
        if točki[0]>točki[1]:
            točki=[točki[1],točki[0]]
       
        seznam2=[seznam[i] if i< točki[0] or i> točki[1] else seznam[točki[0]+točki[1]-i]  for i in range(N)]
        lsq2=sum(razdalja(seznam2[i],seznam2[(i+1)%N]) for i in range(N))
        
            
        #2.korak
        
        if random() <= np.exp(-(lsq2-lsq1)/(kb*(T+0.01))):
            seznam=[seznam2[i] for i in range(N)]
            lsq1=lsq2

        if k%200 and T> 0.1 :
            T-=0.1
#        if T<= 0.1 and k>=1000 and T>0.0001:
#            T-=0.0001
#        if k > 2000 and k%1000 and T>0.01:
#            T-=0.01
        L.append(lsq1)
        
    return seznam, L


start=timer()
N=20
seznam,L=trgovec(N,10000,1,10)
plt.figure(8)
plt.title('Pot popotnika med {} kitajskimi mesti'.format(N))
plt.plot([mesta[seznam[i%N]][5] for i in range(N+1)] ,[mesta[seznam[i%N]][4] for i in range(N+1)],'r')
plt.ylabel('Zemljepisna širina $[rad]$')
plt.xlabel('Zemljepisna dolžina $[rad]$')
plt.figure(9)
plt.title('Utrujenost popotnika N={}'.format(N))
plt.plot([i for i in range(len(L))] ,L,'r')
plt.ylabel('Utrujenost popotnika')
plt.xlabel('Število korakov')
end=timer()
print(end-start)
Beep(780, 1000)


#%%

def trgovecNovi(N,K,kb,T):
    #Definicijski del
    seznam=[i for i in range(N)]
    lsq1=sum(razdalja(seznam[i],seznam[(i+1)%N]) +0 if mesta[seznam[i]][3] != mesta[seznam[(i+1)%N]][3] else razdalja(seznam[i],seznam[(i+1)%N]) for i in range(N))

    L=[]
    
    #algoritmični del
    for k in range(K):
        #1.korak
        točki=[int(N*random()),int(N*random())]
        if točki[0]>točki[1]:
            točki=[točki[1],točki[0]]
       
        seznam2=[seznam[i] if i< točki[0] or i> točki[1] else seznam[točki[0]+točki[1]-i]  for i in range(N)]
        lsq2=sum(razdalja(seznam2[i],seznam2[(i+1)%N]) +0 if mesta[seznam2[i]][3] != mesta[seznam2[(i+1)%N]][3] else razdalja(seznam2[i],seznam2[(i+1)%N]) for i in range(N))

            
        #2.korak
        
        if random() <= np.exp(-(lsq2-lsq1)/(kb*(T+0.01))):
            seznam=[seznam2[i] for i in range(N)]
            lsq1=lsq2

        if k%370 and T> 0.1 :
            T-=0.1
#        if T<= 0.1 and k>=1000 and T>0.0001:
#            T-=0.0001
#        if k > 2000 and k%1000 and T>0.01:
#            T-=0.01
        L.append(lsq1)
        
    return seznam, L


start=timer()
N=37
seznam,L=trgovecNovi(N,50000,1,1)
plt.figure(8)
plt.title('Pot komunističnega popotnika med {} kitajskimi mesti'.format(N))
plt.plot([mesta[seznam[i%N]][5] for i in range(N+1)] ,[mesta[seznam[i%N]][4] for i in range(N+1)],'r')
plt.ylabel('Zemljepisna širina $[rad]$')
plt.xlabel('Zemljepisna dolžina $[rad]$')
plt.figure(9)
plt.title('Utrujenost komunističnega popotnika N={}'.format(N))
plt.plot([i for i in range(len(L))] ,L,'r')
plt.ylabel('Utrujenost popotnika')
plt.xlabel('Število korakov')
end=timer()
print(end-start)
Beep(780, 1000)




