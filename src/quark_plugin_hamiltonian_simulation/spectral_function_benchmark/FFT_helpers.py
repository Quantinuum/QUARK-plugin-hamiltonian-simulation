import numpy as np

def F2(u:list,jj:int,kk:int,sign:int):
    '''FFT on 2 sites'''
    u.append(['s',jj])
    u.append(['h',jj])
    u.append(['h',kk])
    u.append(['rzz',-np.pi/8,jj,kk])
    u.append(['h',kk])
    u.append(['h',jj])
    u.append(['sdg',jj])

    u.append(['h',jj])
    u.append(['s',kk])
    u.append(['h',kk])
    u.append(['rzz',np.pi/8,jj,kk])
    u.append(['h',kk])
    u.append(['sdg',kk])
    u.append(['h',jj])

    u.append(['z',kk])

def F3(u:list,jj:int,kk:int,ll:int,sign:int):
    '''FFT on 3 sites'''
    u.append(['rz',np.pi/4,kk])
    u.append(['tk2',-np.pi/8,-np.pi/8,0,kk,ll])
    u.append(['rz',-np.pi/4,kk])

    beta=np.arcsin(np.sqrt(2/3))
    u.append(['rz',np.pi/4,jj])
    u.append(['tk2',-beta/2,-beta/2,0,jj,kk])
    u.append(['rz',-np.pi/4,jj])

    u.append(['tk2',-np.pi/8*sign,-np.pi/8*sign,0,kk,ll])

    u.append(['rz',np.pi/2,kk])
    u.append(['rz',np.pi/4*sign,ll])

def swap(qubits:list,jj:int,kk:int):
    '''implicit swap'''
    a=qubits[jj]
    qubits[jj]=qubits[kk]
    qubits[kk]=a

def CXfanout(u:list,qubits:list,toact:list):
    '''unitary CX ladder with log depth'''
    n=len(toact)
    if(n==2):
        u.append(['cx',qubits[toact[0]],qubits[toact[1]]])
    if(n>2):
        sub_toact=[toact[1]]
        u.append(['cx',qubits[toact[n-2]],qubits[toact[n-1]]])
        for i in range(1,int(np.ceil(n/2))-1):
            u.append(['cx',qubits[toact[2*i-1]],qubits[toact[2*i]]])
            sub_toact.append(toact[2*i+1])
        if(n%2==0):
            sub_toact.append(toact[n-2])
        CXfanout(u,qubits,sub_toact)
        for i in range(1,int(np.ceil(n/2))-1):
            u.append(['cx',qubits[toact[2*i]],qubits[toact[2*i+1]]])
        u.append(['cx',qubits[toact[0]],qubits[toact[1]]])

def CXfanout_reverse(u:list,qubits:list,toact:list):
    '''inverse of CXfanout'''
    for j in toact:
        u.append(['h',qubits[j]])
    CXfanout(u,qubits,toact)
    for j in toact:
        u.append(['h',qubits[j]])

def CX_reverse(u:list,qubits:list,toact:list,mode:str):
    '''inverse of CX_direct'''
    m=len(toact)
    if(mode=='CX_ladder'):
        for j in range(m-1):
            u.append(['cx',qubits[toact[m-1-j]],qubits[toact[m-2-j]]])
    if(mode=='log2N'):
        CXfanout_reverse(u,qubits,toact)

def CX_direct(u:list,qubits:list,toact:list,mode:str):
    '''CX ladder'''
    for j in toact:
        u.append(['h',qubits[j]])
    CX_reverse(u,qubits,toact[::-1],mode)
    for j in toact:
        u.append(['h',qubits[j]])

def interleave2(u:list,qubits:list,toact:list,direction:int,mode:str):
    '''interleave operation with n=2, i.e. swaps even/odd with first half/second half'''
    list_operations=[]
    balls=[0,1]*(len(toact)//2)
    finished=False
    while(not finished):
        finished=True
        j=1
        while(j<len(toact)):
            if(np.max(balls[:j])>balls[j]):
                list_operations.append([j-1,j])
                balls[j],balls[j-1]=balls[j-1],balls[j]
                j+=2
                finished=False
            else:
                j+=1

    if(direction==0):
        for v in list_operations:
            swap(qubits,toact[v[0]],toact[v[1]])
            if(mode=='CZ_pyramid'):
                u.append(['cz',qubits[toact[v[0]]],qubits[toact[v[1]]]])

    if(mode=='CX_ladder' or mode=='log2N'):
        m=len(toact)
        CX_reverse(u,qubits,[toact[j] for j in range(1,m//2)],mode)
        for j in range(1,m//2):
            u.append(['cz',qubits[toact[j]],qubits[toact[j+m//2-1]]])
        CX_direct(u,qubits,[toact[j] for j in range(1,m//2)],mode)

    if(direction==1):
        for v in list_operations[::-1]:
            swap(qubits,toact[v[0]],toact[v[1]])
            if(mode=='CZ_pyramid'):
                u.append(['cz',qubits[toact[v[0]]],qubits[toact[v[1]]]])

def interleave3(u:list,qubits:list,toact:list,direction:int,mode:str):
    '''interleave operation with n=3, i.e. swaps lines and rows for a 3 x (L//3) lattice'''
    list_operations=[]
    balls=[0,1,2]*(len(toact)//3)
    finished=False
    while(not finished):
        finished=True
        j=1
        while(j<len(toact)):
            if(np.max(balls[:j])>balls[j]):
                list_operations.append([j-1,j])
                balls[j],balls[j-1]=balls[j-1],balls[j]
                j+=2
                finished=False
            else:
                j+=1

    if(direction==0):
        for v in list_operations:
            swap(qubits,toact[v[0]],toact[v[1]])
            if(mode=='CZ_pyramid'):
                u.append(['cz',qubits[toact[v[0]]],qubits[toact[v[1]]]])

    if(mode=='CX_ladder' or mode=='log2N'):
        m=len(toact)
        CX_reverse(u,qubits,[toact[j] for j in range(1,m//3)],mode)
        for j in range(1,m//3):
            u.append(['cz',qubits[toact[j]],qubits[toact[j+m//3-1]]])
        for j in range(1,m//3):
            u.append(['cz',qubits[toact[j]],qubits[toact[j+2*m//3-1]]])
        CX_direct(u,qubits,[toact[j] for j in range(1,m//3)],mode)

        CX_reverse(u,qubits,[toact[m//3+j] for j in range(1,m//3)],mode)
        for j in range(1,m//3):
            u.append(['cz',qubits[toact[j+m//3]],qubits[toact[j+2*m//3-1]]])
        CX_direct(u,qubits,[toact[m//3+j] for j in range(1,m//3)],mode)

    if(direction==1):
        for v in list_operations[::-1]:
            swap(qubits,toact[v[0]],toact[v[1]])
            if(mode=='CZ_pyramid'):
                u.append(['cz',qubits[toact[v[0]]],qubits[toact[v[1]]]])

def FQFT23(u:list,qubits:list,toact:list,sign:int,mode:str,spare:bool,L:int):
    '''FFT, requires that the number of qubits is only a multiple of 2 and/or 3'''
    nn=len(toact)
    if(nn==2):
        F2(u,qubits[toact[0]],qubits[toact[1]],sign)
    elif(nn==3):
        F3(u,qubits[toact[0]],qubits[toact[1]],qubits[toact[2]],sign)
    elif(nn%3==0):
        if(sign==1 and spare):
            interleave3(u,qubits,toact,0,'classical')
        else:
            interleave3(u,qubits,toact,0,mode)
        FQFT23(u,qubits,toact[:len(toact)//3],sign,mode,spare,L)
        FQFT23(u,qubits,toact[len(toact)//3:2*len(toact)//3],sign,mode,spare,L)
        FQFT23(u,qubits,toact[2*len(toact)//3:],sign,mode,spare,L)
        interleave3(u,qubits,toact,1,mode)
        for j in range(nn//3):
            u.append(['rz',np.pi*j/(nn)*sign,qubits[toact[3*j+1]]])
            u.append(['rz',np.pi*2*j/(nn)*sign,qubits[toact[3*j+2]]])
        for j in range(nn//3):
            F3(u,qubits[toact[3*j]],qubits[toact[3*j+1]],qubits[toact[3*j+2]],sign)
        if(sign==-1 and spare and len(toact)==L):
            interleave3(u,qubits,toact,0,'classical')
        else:
            interleave3(u,qubits,toact,0,mode)
    elif(nn%2==0):
        if(sign==1 and spare):
            interleave2(u,qubits,toact,0,'classical')
        else:
            interleave2(u,qubits,toact,0,mode)
        FQFT23(u,qubits,toact[:len(toact)//2],sign,mode,spare,L)
        FQFT23(u,qubits,toact[len(toact)//2:],sign,mode,spare,L)
        interleave2(u,qubits,toact,1,mode)
        for j in range(nn//2):
            u.append(['rz',np.pi*j/(nn)*sign,qubits[toact[2*j+1]]])
        for j in range(nn//2):
            F2(u,qubits[toact[2*j]],qubits[toact[2*j+1]],sign)
        if(sign==-1 and spare and len(toact)==L):
            interleave2(u,qubits,toact,0,'classical')
        else:
            interleave2(u,qubits,toact,0,mode)
    else:
        raise ValueError('system size is not a multiple of 2 and/or 3')
