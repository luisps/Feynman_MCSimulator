from hash_paths import my_hash_inverse, my_hash 
from gates import gate_evaluate_w_sample, gate_evaluate_w, gate_evaluate_w_prob
import random, math
# Importing complex mean and variance functions
from complex_mean_variance import c_mean, c_var_true_mean

import csv 

def _ForwardPath (num_qubits, num_layers, layers, initial_bin_str, DEBUG_MODE = False):
    path = []
    path_w = []
    path_Pw = []
    path_prob = []
    path_Pprob = []
    
    # 1st node is the input state
    path.append(int(initial_bin_str,2))  # nodes along the path represented as integers
    path_w.append(complex(1.0, 0.0))     # at the inital state the amplitude is 1.0
    path_Pw.append(complex(1.0, 0.0))     # at the inital state the product of amplitudes is 1.0    
    path_prob.append(1.0)                # at the inital state the probability is 1.0
    path_Pprob.append(1.0)                # at the inital state the product of probabilities is 1.0
    
    for l in range(num_layers-1):    # do not connect to the endpoint (output state)
        if l==0:
            input_bin_str = initial_bin_str
        else:
            input_bin_str = output_bin_str

        output_bin_list = ['0']*num_qubits  # list with nbr_qubits '0's
                
        layer_weight = complex (1., 0.)
        layer_prob = 1.
        layer = layers[l]

        processed_qubits = []   # some gates have 2 qubits: 
                                # make sure we don't process twice

        # iterate over all qubits 
        for qubit in range(num_qubits):

            if qubit in processed_qubits:
                continue
                
            gate = layer[qubit]
            gate_name = gate[0]
            qubits = gate[1]  # list of addressed qubits
            params = gate[2]  # list of parameters to the gate
            unitary = gate[3] # unitary for parameterized gates
            pdfs = gate[4]    # pdf for parameterized branching gates
            
            for qb in qubits:
                processed_qubits.append(qb)
                
            # 2 qubit gates
            if gate_name in ['cx','cz','cp']:
                input_qbs = (int(input_bin_str[num_qubits-1-qubits[0]]),
                             int(input_bin_str[num_qubits-1-qubits[1]]))
            else: 
                input_qbs = (int(input_bin_str[num_qubits-1-qubits[0]]), )

            if gate_name in ['h','rx','ry','rz']:
                rnd = random.random()
            else:
                rnd = None
                    
            output_qbs, w, prob = gate_evaluate_w_sample (gate_name, input_qbs, rnd=rnd, unitary=unitary, pdfs=pdfs, smpl_dir="DIRECT")
            
            # update ouput_bin_str
            for qb_ndx,qb in enumerate(qubits):
                output_bin_list[num_qubits-1-qb] = str(output_qbs[qb_ndx])
                    
            layer_weight *= w
            layer_prob *= prob
            # end iterating over the layer qubits

        output_bin_str = "".join(output_bin_list)
        
        # store the state that has been reached in path,
        # together with the amplitude and probability of this transition
        
        path.append(int(output_bin_str,2))
        path_w.append(layer_weight)      
        path_Pw.append (layer_weight*path_Pw[-1])
        path_prob.append(layer_prob)               
        path_Pprob.append (layer_prob*path_Pprob[-1])
        # end iterating over all layers except the last
        
    if DEBUG_MODE:
        print ('Forward Path:')
        for i, node in enumerate(zip(path,path_w, path_prob)):
            (p, w, pr) = node
            print ('Node {0}: st={1}, w={2}, prob={3}'.format(i,p,w,pr))
        print ()
    
    return path, path_Pw, path_prob, path_Pprob
    # end ForwardPath
    
def _BackwardPath (num_qubits, num_layers, layers, final_bin_str, DEBUG_MODE = False):
    path = []
    path_w = []
    path_Pw = []
    path_prob = []
    path_Pprob = []
    
    # 1st node is the final state - the lists will be reversed at the end
    path.append(int(final_bin_str,2))  # nodes along the path represented as integers
    path_w.append(complex(1.0, 0.0))   # at the inital backward state the amplitude is 1.0
    path_Pw.append(complex(1.0, 0.0))  # at the inital backward state the product of amplitude is 1.0
    path_prob.append(1.0)              # at the inital backward state the probability is 1.0
    path_Pprob.append(1.0)             # at the inital backward state the product of probability is 1.0
    
    for l in range(num_layers-1, 0, -1):    # do not connect to the endpoint (input state)
        if l==num_layers-1:
            output_bin_str = final_bin_str
        else:
            output_bin_str = input_bin_str

        input_bin_list = ['0']*num_qubits  # list with nbr_qubits '0's
                
        layer_weight = complex (1., 0.)
        layer_prob = 1.
        layer = layers[l]

        processed_qubits = []   # some gates have 2 qubits: 
                                # make sure we don't process twice

        # iterate over all qubits 
        for qubit in range(num_qubits):

            if qubit in processed_qubits:
                continue
                
            gate = layer[qubit]
            gate_name = gate[0]
            qubits = gate[1]  # list of addressed qubits
            params = gate[2]  # list of parameters to the gate
            unitary = gate[3] # unitary for parameterized gates
            pdfs = gate[4]    # pdf for parameterized branching gates
            
            for qb in qubits:
                processed_qubits.append(qb)
                
            # 2 qubit gates
            if gate_name in ['cx','cz','cp']:
                output_qbs = (int(output_bin_str[num_qubits-1-qubits[0]]),
                             int(output_bin_str[num_qubits-1-qubits[1]]))
            else: 
                output_qbs = (int(output_bin_str[num_qubits-1-qubits[0]]), )

            if gate_name in ['h','rx','ry','rz']:
                rnd = random.random()
            else:
                rnd = None
                    
            input_qbs, w, prob = gate_evaluate_w_sample (gate_name, output_qbs, rnd=rnd, unitary=unitary, pdfs=pdfs, smpl_dir="REVERSE")
            
            # update ouput_bin_str
            for qb_ndx,qb in enumerate(qubits):
                input_bin_list[num_qubits-1-qb] = str(input_qbs[qb_ndx])
                    
            layer_weight *= w
            layer_prob *= prob
            # end iterating over the layer qubits

        input_bin_str = "".join(input_bin_list)
        
        # store the state that has been reached in path,
        # together with the prodiuct of amplitudes and probabilities along the path
        
        path.append(int(input_bin_str,2))
        path_w.append(layer_weight)      
        path_Pw.append(layer_weight*path_Pw[-1])      
        path_prob.append(layer_prob)               
        path_Pprob.append(layer_prob*path_Pprob[-1])               
        # end iterating over all layers except the last
        
    # the lists describing the path will now be reversed, 
    # such that each index actually corresponds to the vertex indexing along the path
    # An extra node has to be added to the reverse order path, 
    # which corresponds to traversing the circuit first layer (index [0] after reversal)
    # since the first layer HAS NOT BEEN SAMPLED this node will have 0 value
    path.append (0)
    path_w.append (complex(0.,0.))
    path_Pw.append (complex(0.,0.))
    path_prob.append(0.)
    path_Pprob.append(0.)
    # Now reverse in place
    path.reverse()
    path_w.reverse()
    path_Pw.reverse()
    path_prob.reverse()
    path_Pprob.reverse()
    
    if DEBUG_MODE:
        print ('Backward Path:')
        for i, node in enumerate(zip(path[1:],path_w[1:], path_prob[1:])):
            (p, w, pr) = node
            print ('Node {0}: st={1}, w={2}, prob={3}'.format(i+1,p,w,pr))
        print ()        

    return path, path_Pw, path_prob, path_Pprob
    # end BackwardPath

def _ConnectPath (num_qubits, layer, f_path, f_path_Pw, f_path_Pprob,
                  b_path, b_path_Pw, b_path_Pprob):
    
    format_str = '{{0:0{0}b}}'.format(num_qubits)

    input_bin_str =  format_str.format(f_path[-1])
    output_bin_str = format_str.format(b_path[0])

    layer_weight = complex (1., 0.)

    processed_qubits = []   # some gates have 2 qubits: 
                            # make sure we don't process twice

    # iterate over all qubits 
    for qubit in range(num_qubits):

        if qubit in processed_qubits:
            continue
                
        gate = layer[qubit]
        gate_name = gate[0]
        qubits = gate[1]  # list of addressed qubits
        unitary = gate[3] # unitary for parameterized gates
            
        for qb in qubits:
            processed_qubits.append(qb)
                
        # 2 qubit gates
        if gate_name in ['cx','cz','cp']:
            input_qbs = (int(input_bin_str[num_qubits-1-qubits[0]]),
                         int(input_bin_str[num_qubits-1-qubits[1]]))
            output_qbs = (int(output_bin_str[num_qubits-1-qubits[0]]),
                         int(output_bin_str[num_qubits-1-qubits[1]]))
        else: 
            input_qbs = (int(input_bin_str[num_qubits-1-qubits[0]]), )
            output_qbs = (int(output_bin_str[num_qubits-1-qubits[0]]), )

        w = gate_evaluate_w (gate_name, input_qbs, output_qbs, unitary=unitary)

        layer_weight *= w
        # end iterating over the layer qubits

    #weight = math.prod(f_path_w) * layer_weight * math.prod(b_path_w)
    weight = f_path_Pw * layer_weight * b_path_Pw
    prob = f_path_Pprob * b_path_Pprob 

    return weight, prob, layer_weight

def _ConnectPathMIS (num_qubits, layer, f_path, f_path_prob, f_path_Pw,
                  b_path, b_path_prob, b_path_Pw, DEBUG_MODE=False):
    
    format_str = '{{0:0{0}b}}'.format(num_qubits)
    
    connecting_layer = len(f_path) - 1
    
    if DEBUG_MODE:
        print ("Connecting layer: ", connecting_layer)

    input_bin_str =  format_str.format(f_path[-1])
    output_bin_str = format_str.format(b_path[0])
    
    # Build the new path description
    path = f_path + b_path            # num_layers+1 elements
    num_layers = len(path)-1
    new_segment_w = complex(1., 0.)   # to be evaluated
    path_w = f_path_Pw * b_path_Pw    # this will be * new_segment_w
    new_segment_prob = 1.             # to be evaluated
    path_prob = f_path_prob + [new_segment_prob] + b_path_prob
    
    if DEBUG_MODE:    
        print ("The forward path is ", f_path);
        print ("This path probs = {0} and product weight = {1}".format(f_path_prob, f_path_Pw));
        print ()
        print ("The backward path is ", b_path);
        print ("This path probs = {0} and product weight = {1}".format(b_path_prob, b_path_Pw));
        print ()
        print ("The path being evaluated is {0}, with {1} layers".format(path,num_layers))
        print ("Partial weights product is ", path_w);
        print ("Probs (with connecting layer set to 1, are ", path_prob)

    processed_qubits = []   # some gates have 2 qubits: 
                            # make sure we don't process twice

    # iterate over all qubits 
    for qubit in range(num_qubits):

        if qubit in processed_qubits:
            continue
                
        gate = layer[qubit]
        gate_name = gate[0]
        qubits = gate[1]  # list of addressed qubits
        unitary = gate[3] # unitary for parameterized gates
        pdfs = gate[4]    # pdf for parameterized branching gates
            
        for qb in qubits:
            processed_qubits.append(qb)
                
        # 2 qubit gates
        if gate_name in ['cx','cz','cp']:
            input_qbs = (int(input_bin_str[num_qubits-1-qubits[0]]),
                         int(input_bin_str[num_qubits-1-qubits[1]]))
            output_qbs = (int(output_bin_str[num_qubits-1-qubits[0]]),
                         int(output_bin_str[num_qubits-1-qubits[1]]))
        else: 
            input_qbs = (int(input_bin_str[num_qubits-1-qubits[0]]), )
            output_qbs = (int(output_bin_str[num_qubits-1-qubits[0]]), )

        w, prob = gate_evaluate_w_prob (gate_name, input_qbs, output_qbs, unitary=unitary, pdfs=pdfs)

        new_segment_w *= w
        new_segment_prob *= prob
        # end iterating over the layer qubits

    path_w *= new_segment_w
    path_prob[connecting_layer+1] = new_segment_prob
    
    if DEBUG_MODE:
        print ("Path probs= ", path_prob)
        print ("Path weight = ", path_w)

    # if the amplitude or the probability is 0 terminate
    if abs(new_segment_w)==0. or new_segment_prob==0:
        if DEBUG_MODE:
            print ("Terminating this path")
            print ()
        return complex(0.,0.), 0., complex(0.,0.)
    
    # compute MIS_weight as the sum of 
    # the probabilities of the num_layers alternative deterministic connections
    MIS_weight_reciprocal = 0.
    if DEBUG_MODE:
        print ("MIS_weight reciprocal = (", end ='')
    for det_connect in range(num_layers):
        # compute probability with deterministic transition across layer det_connect
        #prob = 1.  
        #for l in range(1,det_connect+1):
        #    prob *= path_prob[l]
        prob = math.prod (path_prob[1:det_connect+1])
        #for l in range(det_connect+2, num_layers+1):
        #    prob *= path_prob[l]
        prob *= math.prod (path_prob[det_connect+2:num_layers+1])
        MIS_weight_reciprocal += prob
        if DEBUG_MODE:
            print ("{0} + ".format(prob), end ='')

    # divide the above sum by the number of summands 
    # note: we divide because we are computing the reciprocal of the MIS_weight
    MIS_weight_reciprocal /= num_layers
        
    if DEBUG_MODE:
        print (") / {0} = {1}".format(num_layers, MIS_weight_reciprocal))
        print ()
    
    return path_w, MIS_weight_reciprocal, new_segment_w

def BD (num_qubits, num_layers, layers, initial_bin_str, final_bin_str, P, N, M, MIS=True,
               true_amplitude = None):
    
    DEBUG_MODE = False
    
    non_zero_paths = 0
    non_zero_paths_ndx = []
    total_sampled_paths = 0
    sum_amplitude = complex (0., 0.)
    # for complex variance see
    # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
    sum_amplitude_squared = 0.
    sum_true_variance = 0.
    
    format_str = '{{0:0{0}b}}'.format(num_qubits)

    # compute the reporting sample index
    report_stage = 1
    report_ndx = math.ceil(M/100)
    profile_list = []
    
    # csv output paths throughput and amplitude and variance estimates
    #csv_handle = open(r"RCS_throughputs.csv","w")
    #csv_writer = csv.writer(csv_handle, delimiter=' ',
    #                        quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    for s in range(M):
        if DEBUG_MODE:
            print ()
            print ('Sample {}'.format(s))
        
        # f_path has num_layers elements from the initial state 'f_path[0]' 
        #        to the next-to-last state 'f_path[num_layers-1]'
        #        each element is the integer to the respective state
        # f_path_Pw has the product of amplitudes associated with the transition to that state
        #        f_path_Pw[0] = 1.0 + 0.0 i
        # f_path_prob has the probability associated with the transition to that state
        #        f_path_prob[0] = 1.0
        f_path, f_path_Pw, f_path_prob, f_path_Pprob = _ForwardPath (num_qubits, num_layers, layers, initial_bin_str, DEBUG_MODE)

        # b_path has num_layers+1 elements
        #        indexes grow from the input state (b_path[0]) to the final state (b_path[num_layers])
        #.       b_path[0] has no meaning since the initial state is not sampled BACKWARDS 
        #        b_path[1]..b_path[num_layers] each contain the integer to the respective state
        # b_path_Pw has the product of amplitudes associated with the transition to that state
        #        b_path_Pw[0] has no meaning
        #        b_path_Pw[num_layers] = 1.0 + 0.0 i
        # b_path_prob has the probability associated with the transition to that state
        #        b_path_prob[0] has no meaning
        #        b_path_prob[num_layers] = 1.0
        b_path, b_path_Pw, b_path_prob, b_path_Pprob = _BackwardPath (num_qubits, num_layers, layers, final_bin_str, DEBUG_MODE)
            
        # connect 'num_layers' different paths
        for l in range(num_layers):
            if not MIS:
                w, prob, lw = _ConnectPath(num_qubits, layers[l], 
                                           f_path[:l+1], f_path_Pw[l], f_path_Pprob[l],
                                           b_path[l+1:], b_path_Pw[l+1], b_path_Pprob[l+1])
            else:
                w, prob, lw = _ConnectPathMIS(num_qubits, layers[l], 
                                              f_path[:l+1], f_path_prob[:l+1], f_path_Pw[l],
                                              b_path[l+1:], b_path_prob[l+1:], b_path_Pw[l+1], DEBUG_MODE)
                    
            if abs(w)!=0. and prob!=0.:
                path_throughput = w / prob
                sum_amplitude += path_throughput
                # for complex variance see
                # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
                sum_amplitude_squared += abs(path_throughput)**2 
                non_zero_paths += 1
                # store the path ndx
                non_zero_paths_ndx.append (my_hash(f_path[1:l+1] + b_path[l+1:-1], N))
            else:
                path_throughput = complex (0.,0.)

            total_sampled_paths += 1
                                
            if not true_amplitude is None:
                sum_true_variance += abs(path_throughput - true_amplitude)**2 

            ### write onto csv
            #csv_writer.writerow([path_throughput, sum_amplitude / total_sampled_paths, 
            #                     sum_true_variance / total_sampled_paths])

            if DEBUG_MODE:
                print ('Connected path: ',end='')
                #for i in range (l+1):
                #    print ('{} '.format(f_path[i]),end='')
                print ('{} '.format(f_path[:l+1]),end='')
                print ('| ',end='')
                #for i in range (l+1, num_layers+1):
                #    print ('{} '.format(b_path[i]),end='')
                print ('{} '.format(b_path[l+1:]),end='')
                print ()
                print ('f_w={0}; l_w={1}; b_w={2}'.format(f_path_Pw[l], lw, b_path_Pw[l+1]),end='')
                print()
                print ('f_prob={0}; l_prob=1.0; b_prob={1}'.format(f_path_Pprob[l], b_path_Pprob[l+1]),end='')
                print ()
                print ('path w={0}; prob={1}'.format(w,prob))
                print ('throughput={}'.format(path_throughput))
                print ()
                
        # end iterating over permutations of forward and backward paths
        
        # report variance and amplitude at increments of 1% of M
        #print ('Stats: report_stage ={0}, report_ndx={1}, total_sampled_paths={2}, M={3}'.format(report_stage,report_ndx, total_sampled_paths, M))
        if (report_stage==100 and (s+1)>=M) or (s+1)>=report_ndx or (s+1)>=M:
            amplitude = sum_amplitude / total_sampled_paths
            # for complex variance see
            # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
            variance = sum_amplitude_squared / total_sampled_paths - (abs(amplitude)**2)
            if not true_amplitude is None:
                true_variance = sum_true_variance / total_sampled_paths
                profile_list.append((s+1,amplitude,true_variance))
                print ('{0}%: {1} samples: estimate = {2:.4f}, true variance = {3:.4f}; est. variance = {4:.4f}'.format(report_stage, s+1, amplitude, true_variance, variance))
            else:
                profile_list.append((s+1,amplitude,variance))
                print ('{0}%: {1} samples: estimate = {2:.4f}, est. variance = {3:.4f}'.format(report_stage, s+1 ,amplitude, variance))
            report_stage += 5
            report_ndx = math.ceil(M*report_stage/100)
                                 
    # end iterating over samples
    
    amplitude = sum_amplitude / total_sampled_paths
    # for complex variance see
    # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
    variance = sum_amplitude_squared / total_sampled_paths - (abs(amplitude)**2)
    if not true_amplitude is None:
        true_variance = sum_true_variance / total_sampled_paths
        
    #csv_handle.close()

    return amplitude, non_zero_paths, total_sampled_paths, true_variance if not true_amplitude is None else variance, profile_list, non_zero_paths_ndx


    