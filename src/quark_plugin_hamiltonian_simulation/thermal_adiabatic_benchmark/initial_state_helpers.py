import numpy as np
from scipy.optimize import minimize

def convert(H):
    H_converted=[]
    list_non_trivial=[]
    for term in H:
        if(H[term]!=0):
            term_converted=[]
            for j in range(len(term)):
                if(term[j]!='I'):
                    term_converted.append([term[j],j])
                    if(j not in list_non_trivial):
                        list_non_trivial.append(j)
            H_converted.append([term_converted,H[term]])
    return H_converted,list_non_trivial

def energyTerm(term,rotations):
    res=1
    for c in term:
        theta=rotations[2*c[1]]
        alpha=rotations[2*c[1]+1]
        if(c[0]=='X'):
            res*=np.sin(2*theta)*np.cos(2*alpha)
        if(c[0]=='Y'):
            res*=np.sin(2*theta)*np.sin(2*alpha)
        if(c[0]=='Z'):
            res*=np.cos(2*theta)
    return res

def energyHamiltonian(H_converted,rotations):
    res=0
    for h in H_converted:
        res+=h[1]*energyTerm(h[0],rotations)
    return res

def computeGroundState(H):
    N=len(list(H.keys())[0])
    
    H_converted,list_non_trivial=convert(H)

    def opt(x):
        return energyHamiltonian(H_converted,np.tanh(x)*np.pi/1.9)
    
    x=minimize(opt,x0=np.ones(2*N)).x

    x=np.tanh(x)*np.pi/1.9

    U = []
    for j in list_non_trivial:
        if(x[2*j]!=0):
            U.append(['ry',x[2*j],j])
        if(x[2*j+1]!=0):
            U.append(['rz',x[2*j+1],j])

    return U
        



    
