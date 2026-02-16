from qiskit import QuantumCircuit

def create_circuit_qiskit(U,N):
    V=QuantumCircuit(N,N)

    for u in U:
        if(u[0]=='x'):
            V.x(u[1])
        elif(u[0]=='y'):
            V.y(u[1])
        elif(u[0]=='z'):
            V.z(u[1])
        elif(u[0]=='h'):
            V.h(u[1])
        elif(u[0]=='s'):
            V.s(u[1])
        elif(u[0]=='sdg'):
            V.sdg(u[1])
        elif(u[0]=='rx'):
            if(u[1]!=0):
                V.rx(2*u[1],u[2])
        elif(u[0]=='ry'):
            if(u[1]!=0):
                V.ry(2*u[1],u[2])
        elif(u[0]=='rz'):
            if(u[1]!=0):
                V.rz(2*u[1],u[2])
        elif(u[0]=='cx'):
            V.cx(u[1],u[2])
        elif(u[0]=='cy'):
            V.cy(u[1],u[2])
        elif(u[0]=='cz'):
            V.cz(u[1],u[2])
        elif(u[0]=='rxx'):
            if(u[1]!=0):
                V.rxx(2*u[1],u[2],u[3])
        elif(u[0]=='ryy'):
            if(u[1]!=0):
                V.ryy(2*u[1],u[2],u[3])
        elif(u[0]=='rzz'):
            if(u[1]!=0):
                V.rzz(2*u[1],u[2],u[3])
        elif(u[0]=='measure'):
            V.measure(u[1],u[2])
        elif(u[0]=='barrier'):
            V.barrier()

    return V

    
