import numpy as np
from .FFT_helpers import FQFT23,interleave2

def trotterStep(u:list,qubits:list,dt:float,eps:float,omega:float,parity:int,N:int):
    '''implements a Trotter step of system + environment, for coupling epsilon and frequency omega'''
    L=2*N
    for decal in [0,1]:
        for j in range(decal,N,2):
            u.append(['cz',qubits[2*j],qubits[2*j+1]])
            if(parity==1):
                u.append(['rxx',-dt/2,qubits[2*j],qubits[(2*j+2)%L]])
                u.append(['ryy',-dt/2,qubits[2*j],qubits[(2*j+2)%L]])
            else:
                u.append(['rxx',-dt/2*np.sign(((2*j+2)%L)-2*j),qubits[2*j],qubits[(2*j+2)%L]])
                u.append(['ryy',-dt/2*np.sign(((2*j+2)%L)-2*j),qubits[2*j],qubits[(2*j+2)%L]])
            u.append(['cz',qubits[2*j],qubits[2*j+1]])

    for j in range(N):
        u.append(['rxx',-dt*eps/4,qubits[2*j],qubits[2*j+1]])
        u.append(['ryy',-dt*eps/4,qubits[2*j],qubits[2*j+1]])

    for j in range(N):
        u.append(['rz',omega*dt/2,qubits[2*j+1]])

def initial_state_prep(u:list,qubits:list,filled:bool,mode:str,spare:bool,N:int):
    '''initial state preparation on system + environment'''
    K=[]
    for j in range(N):
        if(np.cos(2*np.pi/N*j)<0):
            K.append(j)
    for j in range(N):
        if(j in K):
            u.append(['x',qubits[2*j]])
    FQFT23(u,qubits,[2*jj for jj in range(N)],1,mode,spare,N)
    if(filled):
        for j in range(N):
            u.append(['x',qubits[2*j+1]])
            if(j%2==(N+1)%2):
                u.append(['z',qubits[2*j]])
    return len(K)%2

def create_circuit(epsilon:float,omega:float,dt:float,Ntrot:int,filled:bool,mode:str,spare:bool,N:int):
    '''circuit corresponding to one setting of omega and filled environment'''
    L=2*N
    qubits=[j for j in range(L)]

    u = []

    parity=initial_state_prep(u,qubits,filled,mode,spare,N)

    if(filled and N%2==1):
        parity=1-parity

    for t in range(Ntrot):
        trotterStep(u,qubits,dt,epsilon,omega,parity,N)

    interleave2(u,qubits,[j for j in range(L)],0,mode)

    FQFT23(u,qubits,[jj+N for jj in range(N)],-1,mode,spare,N)

    for j in range(L):
        u.append(['measure',qubits[j],j])

    return u

def create_sequence_circuit(epsilon:float,domega:float,n_omega:int,dt:float,Ntrot:int,mode:str,N:int):
    '''sequence of circuits for all values of omega and environment filling'''
    sequence=[]

    for nn in range(n_omega):
        omega=domega*(nn-(n_omega-1)/2)

        for filled in [False,True]:

            u=create_circuit(epsilon,omega,dt,Ntrot,filled,mode,True,N)

            sequence.append(u)

    return sequence
