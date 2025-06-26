import numpy as np
from itertools import combinations
from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister, ClassicalRegister


###################################################################
## Helper Functions
###################################################################


def u(r,c,k): 
    '''
    Labels unit cells in the toric code lattice at row `r` and column `c` (0-indexed)
    '''
    return (r%k)*k + (c%k);


###################################################################


def norm1(p1, p2, k):
    '''
    Computes the Manhattan norm (1-norm) between points given by the ordered pairs `p1`, `p2`
    on a k × k Toric Code, taking into account the periodic boundary 
    '''
    return sum( min( abs(p1[j]-p2[j]), k - abs(p2[j]-p1[j]) ) for j in range(len(p1)) )


###################################################################


def pair_points(points):
    '''
    Generates all possible pairings of points; 
    used to iterate over and find shortest paths to pair up defects
    '''
    def helper(remaining):
        if not remaining:
            return [[]]
        pairings = [];
        first = remaining[0];
        for i in range(1, len(remaining)):
            pair = (first, remaining[i]);
            rest = remaining[1:i] + remaining[i+1:];
            for subpairing in helper(rest):
                pairings.append([pair] + subpairing);
        return pairings

    return helper(points)
    

###################################################################
## Quantum Circuits
###################################################################


def Av_proj(vis=False):
    '''
    Creates circuit for unitary projection onto eigenstates of Av
    '''
    qr = QuantumRegister(size = 4, name = "q");
    qc = QuantumCircuit(qr, name = "U_A");

    qc.h(qr[0]);
    for j in range(1,4):
        qc.cx(qr[0], qr[j]);

    vis and print(qc.draw());
    return qc.to_gate();


###################################################################


def X_measure(n, vis=False):
    '''
    Creates the circuit to measure the product of Xs on `n` qubits
    '''
    qr = QuantumRegister(size = n, name = "q");
    ar = AncillaRegister(size = 1, name = "a");
    qc = QuantumCircuit(qr, ar, name = "A_v");

    qc.h(qr);
    for j in range(n):
        qc.cx(qr[j], ar[0]);
    qc.h(qr);

    vis and print(qc.draw());
    return qc.to_gate();


###################################################################


def Z_measure(n, vis=False):
    '''
    Creates the circuit to measure the product of Zs on `n` qubits
    '''
    qr = QuantumRegister(size = n, name = "q");
    ar = AncillaRegister(size = 1, name = "a");
    qc = QuantumCircuit(qr, ar, name = "B_p");

    for j in range(n):
        qc.cx(qr[j], ar[0]);

    vis and print(qc.draw());
    return qc.to_gate();


###################################################################


def create_prep_TC(k, vis=False):
    '''
    Creates circuit to generate the logical |00⟩ state of the k × k Toric Code
    '''
    n = 2*k**2;
    print(n, " qubits in Toric Code \n ")

    qr = QuantumRegister(size = n, name = "q");
    qc = QuantumCircuit(qr, name = "Prep TC");

    ## List out qubits for all stabilizers
    A_list = [ (2*u(r,c,k), 2*u(r,c,k)+1, 2*u(r-1,c,k), 2*u(r,c-1,k)+1) for r in range(k) for c in range(k)];
    ## reorder last row for proper application of gates
    for j in range(1,k+1):
        A_list[-j] = A_list[-j][1:] + A_list[-j][:1];

    vis and print("A_v Projector (GHZ_4) : ")
    A_proj = Av_proj(vis);

    ## Project onto +1 eigenstate of all Av
    for As in A_list[:-1]:
        qc.compose(A_proj, qubits=[qr[j] for j in As], inplace=True);

    vis and print("TC Circuit : \n", qc.draw());
    return qc.to_gate();


###################################################################
    

def create_stabilizer_measures(k, vis=False):
    '''
    Creates circuit for measurement of all stabilizers in the k × k Toric Code
    '''
    n = k**2;
    print(2*n, " stabilizer qubits ")

    data_qr = QuantumRegister(size = 2*n, name = "q");
    A_qr = AncillaRegister(size = n, name = "A");
    B_qr = AncillaRegister(size = n, name = "B");
    qc = QuantumCircuit(data_qr, A_qr, B_qr, name = "Stabilizers");
    print(qc.num_qubits, " qubits total ")

    ## List out qubits for all stabilizers
    A_list = [ (2*u(r,c,k), 2*u(r,c,k)+1, 2*u(r-1,c,k), 2*u(r,c-1,k)+1) for r in range(k) for c in range(k)];
    B_list = [ (2*u(r,c,k), 2*u(r,c,k)+1, 2*u(r+1,c,k)+1, 2*u(r,c+1,k)) for r in range(k) for c in range(k)];

    ## Measure all Av stabilizers
    vis and print("\n A_v Stabilizer Measure : ")
    A_gate = X_measure(4, vis);

    for n,As in enumerate(A_list):
        qc.compose(A_gate, qubits = (*[data_qr[j] for j in As], *[A_qr[n]]), inplace=True);

    ## Measure all Bp stabilizers
    vis and print("\n B_p Stabilizer Measure : ")
    B_gate = Z_measure(4, vis);

    for n,Bs in enumerate(B_list):
        qc.compose(B_gate, qubits = (*[data_qr[j] for j in Bs], *[B_qr[n]]), inplace=True);

    vis and print("Stabilizer Circuit : \n", qc.draw());
    return qc.to_gate();


###################################################################
    

def create_logical_prep(k, log_state, vis=False):
    '''
    Creates circuit to initialize the k × k Toric Code in logical |00⟩ to one of logical |01⟩, |10⟩, or |11⟩
    '''
    assert log_state in {"00", "01", "10", "11"}, "State must be one of the following: 00, 01, 10, 11"

    n = k**2;
    qr = QuantumRegister(size = 2*n, name = "q");
    qc = QuantumCircuit(qr, name = "Logical Prep");

    if log_state == "01":
        X1_bits = [ 2*u(r,0,k)+1 for r in range(k) ];
        
        qc.x([qr[j] for j in X1_bits]);
    elif log_state == "10":
        X2_bits = [ 2*u(0,c,k) for c in range(k) ];
        
        qc.x([qr[j] for j in X2_bits]);
    elif log_state == "11": 
        X1_bits = [ 2*u(r,0,k)+1 for r in range(k) ];
        X2_bits = [ 2*u(0,c,k) for c in range(k) ];
        
        qc.x([qr[j] for j in X1_bits]);
        qc.x([qr[j] for j in X2_bits]);

    vis and print("Logical State Prep : \n", qc.draw());
    return qc.to_gate();


###################################################################


def create_logical_measure(k, vis=False):
    '''
    Creates circuit to measure the logical Z operators on the k × k Toric Code
    '''
    n = k**2;
    data_qr = QuantumRegister(size = 2*n, name = "q");
    Z_qr = AncillaRegister(size = 2, name = "Z");
    qc = QuantumCircuit(data_qr, Z_qr, name = "Logical Measure");
    
    log_Z = Z_measure(k);

    Z1_bits = [ 2*u(0,c,k)+1 for c in range(k) ];
    Z2_bits = [ 2*u(r,0,k) for r in range(k) ];

    qc.compose(log_Z, qubits = (*[data_qr[j] for j in Z1_bits], *[Z_qr[0]]), inplace=True);
    qc.compose(log_Z, qubits = (*[data_qr[j] for j in Z2_bits], *[Z_qr[1]]), inplace=True);

    vis and print("Logical Measurement : \n", qc.draw());
    return qc.to_gate();


###################################################################


def create_single_qubit_error(k, spoiler=False, vis=False):
    '''
    Creates a simple error channel acting on a random single qubit in the k × k Toric Code 
    '''
    n = k**2;
    qr = QuantumRegister(size = 2*n, name = "q");
    qc = QuantumCircuit(qr, name = "Single Qubit Error");

    loc = np.random.randint(2*n);
    spoiler and print("Error at qubit : ", loc);

    if np.random.rand()>0.5:
        qc.z(qr[loc]);
    else:
        qc.x(qr[loc]);

    (spoiler and vis) and print("Single qubit error : \n", qc.draw());
    return qc.to_gate();


###################################################################


def create_random_Pauli_error(k, px, pz, spoiler=False, vis=False):
    '''
    Creates a random Pauli error channel acting on the k × k Toric Code, where for each qubit the probability of an X error is `px` and probability of Z error is `pz`
    '''
    n = k**2
    qr = QuantumRegister(size = 2*n, name = "q");
    qc = QuantumCircuit(qr, name = "Pauli Errors");

    Xerrs = []; Zerrs = [];
    for q in qr:
        if np.random.rand()<px: 
            qc.x(q);
            Xerrs.append(q);
        if np.random.rand()<pz: 
            qc.z(q);
            Zerrs.append(q);

    spoiler and print("\n X errors at : ", Xerrs);
    spoiler and print("\n Z errors at : ", Zerrs);
    (spoiler and vis) and print("\n Random Pauli Errors : \n", qc.draw());
    return qc.to_gate()
    

###################################################################
         

def decoder(qc, k):
    n = k**2;
    ## for each integer from 0 to 2^n-1:
    for z in range(1,2**n):
        ## get the binary digits and positions of the 1s in those digits
        ones_at = [n-(i+1) for i, bit in enumerate(format(z, f'0{n}b')) if bit == '1']
        ## Should NEVER have an odd number of errors; violates topological invariance
        if len(ones_at)%2==0:
            ## map those positions onto the rows/columns of a 2d grid
            coords = [(a//k, a%k) for a in ones_at];
            ## look at each possible pairing of digits:
            all_pairings = pair_points(coords);
            ## here we wants to select the pairings which minimize the Manhattan distance to connect the sites
            ## (assuming that for small error probability, errors won't have traveled far)
            dists = [sum(norm1(p[0], p[1], k) for p in pairs) for pairs in all_pairings];
            ## For now: pick the first pairing which works - may not be the best one, but how else to decide?
            chosen_pairing = all_pairings[dists==min(dists)]
    
            x_sites = [];
            z_sites = [];
            ## for each pair of defects:
            for p in chosen_pairing:
                ## make Pauli string for a connecting path
                ((r1,c1),(r2,c2)) = p; 

                ## First align the columns:
                if c1 == c2:
                    c0 = c2;
                elif c1<c2:
                    ## c1 on the left, c2 on the right
                    if (c2-c1)%k <= (c1-c2)%k:
                        ## move c2 left to c1
                        c0 = c1;
                        [x_sites.append(2*u(r2,c,k)) for c in range(c1+1,c2+1)];
                        [z_sites.append(2*u(r2,c,k)+1) for c in range(c1,c2)];
                    else: 
                        ## meet at the boundary
                        c0 = 0; 
                        [x_sites.append(2*u(r2,c,k)) for c in range(c2+1,k+1)];
                        [x_sites.append(2*u(r1,c,k)) for c in range(1,c1+1)];
                        [z_sites.append(2*u(r2,c,k)+1) for c in range(c2,k)];
                        [z_sites.append(2*u(r1,c,k)+1) for c in range(c1)];
                else: 
                    ## c2 on the left, c1 on the right
                    if (c1-c2)%k <= (c2-c1)%k:
                        ## move c1 left to c2
                        c0 = c2;
                        [x_sites.append(2*u(r1,c,k)) for c in range(c2+1,c1+1)];
                        [z_sites.append(2*u(r1,c,k)+1) for c in range(c2,c1)];
                    else: 
                        ## meet at the boundary
                        c0 = 0; 
                        [x_sites.append(2*u(r1,c,k)) for c in range(c1+1,k+1)];
                        [x_sites.append(2*u(r2,c,k)) for c in range(1,c2+1)];
                        [z_sites.append(2*u(r1,c,k)+1) for c in range(c1,k)];
                        [z_sites.append(2*u(r2,c,k)+1) for c in range(c2)];

                ## Next align the rows:
                if r1<r2:
                    ## r1 on top, r2 on bottom
                    if (r2-r1)%k <= (r1-r2)%k:
                        ## move r2 up to r1
                        [x_sites.append(2*u(r,c0,k)+1) for r in range(r1+1, r2+1)];
                        [z_sites.append(2*u(r,c0,k)) for r in range(r1, r2)];
                    else:
                        ## meet at the boundary
                        [x_sites.append(2*u(r,c0,k)+1) for r in range(1,r1)];
                        [x_sites.append(2*u(r,c0,k)+1) for r in range(r2+1,k+1)];
                        [z_sites.append(2*u(r,c0,k)) for r in range(r1)];
                        [z_sites.append(2*u(r,c0,k)) for r in range(r2,k)];
                else:
                    ## r2 on top, r1 on bottom
                    if (r1-r2)%k <= (r2-r1)%k:
                        ## move r1 up to r2
                        [x_sites.append(2*u(r,c0,k)+1) for r in range(r2+1, r1+1)];
                        [z_sites.append(2*u(r,c0,k)) for r in range(r2, r1)];
                    else:
                        ## meet at the boundary
                        [x_sites.append(2*u(r,c0,k)+1) for r in range(1,r2)];
                        [x_sites.append(2*u(r,c0,k)+1) for r in range(r1+1,k+1)];
                        [z_sites.append(2*u(r,c0,k)) for r in range(r2)];
                        [z_sites.append(2*u(r,c0,k)) for r in range(r1,k)];
    
            ## Apply Pauli strings to x_sites, z_sites
            with qc.if_test((qc.cregs[0], z)):
                qc.z([qc.qregs[0][j] for j in z_sites]);
            with qc.if_test((qc.cregs[1], z)):
                qc.x([qc.qregs[0][j] for j in x_sites]);

    return qc


    
            
    