# receives a tuple with L-1 elements, 
# each being an intermediate computational basis state
# along the path
# note that the first state (|x>) and the last state (|t>)
# are not included in the tuple
#
# also receives N (the number of computational basis states= 2**(num_qubits) )
#
# returns an unique integer between 0 .. N^{L-1} 
# identifying the path 
def my_hash (path, N):
    v = 0
    for l, state in enumerate(path):
        v += (state * (N**l))
    return v

# receives an integer which is the desired path ID  
# receives N (the number of computational basis states= 2**(num_qubits) )
# receives the integer hops, i.e. the number of intermediate states
#    excluding the initial and final state, for L layers (0..L-1) hops=L-1
#
# returns a list of hops intermediate states, each in the range [0 .. N-1] 
def my_hash_inverse (v, N, hops):
    path_list = []
    for i in range(hops):
        l = v % N
        path_list.append(l)
        v //= N
    path = tuple(path_list)
    return path

def verify_my_hash ():
    for i in range (16):  # make 16 comparisons
        n = random.randrange (2,6)
        N = 2**n
        hops = random.randrange (2,6)
        
        # number of possible paths
        P = N**hops
        # generate a path identifier
        v = random.randrange(P)
        
        path = my_hash_inverse (v, N, hops)
        v_rec = my_hash (path, N)
        
        print ('N={0}; hops={1}; P={2}'.format(N,hops,P))
        print ('v={0}; path={1}'.format(v,path))
        if (v==v_rec):
            print ('SUCCESS!')
        else:
            print ('ERROR: {0} != {1}'.format(v,v_rec))
        print ()

        
#verify_my_hash()


    