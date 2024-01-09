from hash_paths import my_hash_inverse, my_hash 
from gates import gate_evaluate_w
import random, math
# Importing complex mean and variance functions
from complex_mean_variance import c_mean, c_var_true_mean

def Uniform (num_qubits, num_layers, layers, initial_bin_str, final_bin_str, P, N, M,
               true_amplitude = None):
    non_zero_paths = 0
    sum_amplitude = complex (0., 0.)
    # for complex variance see
    # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
    sum_amplitude_squared = 0.
    sum_true_variance = 0.
    
    format_str = '{{0:0{0}b}}'.format(num_qubits)

    # create a list of lists to hold the values of the qubits
    # at each intermediate state
    qb_v = []
    for qb in range(num_qubits):
        qb_v.append([])
        for l in range(num_layers+1):
            qb_v[qb].append(None)

    for qb in range(num_qubits):
        #inital state
        qb_v[num_qubits-qb-1][0] = int(initial_bin_str[qb])
        #final state
        qb_v[num_qubits-qb-1][num_layers] = int(final_bin_str[qb])    
    
    # compute the reporting sample index
    report_stage = 1
    report_ndx = math.ceil(M/100)
    profile_list = []
    # iterate all samples
    for s in range(M):
        path_ndx = math.floor(random.random()*P)
        path = my_hash_inverse (path_ndx, N, num_layers-1)
        #print ('{0} : {1}'.format(v, path))

        for l in range(1,num_layers):
            bstring = format_str.format(path[l-1])
            for qb in range(num_qubits):
                qb_v[qb][l] = int(bstring[num_qubits-qb-1])
            
        path_throughput = complex(1.,0)  # path initially has weight 1
        for l in range(num_layers):

            # compute the product across gates for the current path and layer
            layer_weight = complex(1., 0.)  # layer initially has weight 1
        
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
            
                for qb in qubits:
                    processed_qubits.append(qb)
                        
                # 2 qubit gates
                if gate_name in ['cx','cz','cp']:
                    input_qbs = (qb_v[qubits[0]][l], qb_v[qubits[1]][l])
                    output_qbs = (qb_v[qubits[0]][l+1],qb_v[qubits[1]][l+1])
                else: 
                    input_qbs = (qb_v[qubit][l], )
                    output_qbs = (qb_v[qubit][l+1],)
                w = gate_evaluate_w (gate_name, input_qbs, output_qbs, 
                                     params=params, unitary=unitary)

                if abs(w)==0.:
                    layer_weight = complex(0.,0.)
                    path_throughput = complex(0.,0.)
                    break      # early terination of loop over qubits
                else:
                    layer_weight *= w
            # end iterating over the layer qubits
            if abs(layer_weight) == 0.:
                path_throughput = complex(0.,0.)
                break          # early terination of loop over layers
            else:
                path_throughput *= layer_weight
        # end iterating over layers
        if not true_amplitude is None:
            sum_true_variance += abs(path_throughput - true_amplitude)**2 
        if abs(path_throughput) != 0.:
            sum_amplitude += path_throughput
            # for complex variance see
            # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
            sum_amplitude_squared += abs(path_throughput)**2
            non_zero_paths += 1

        # compute current expectyed value and estimated variance
        # NOTE: there is an overhead on ddoing this for each sample,
        #       but this approach is numerically more stable than accumulating the sums (including squared)
        # see https://cas.ee.ic.ac.uk/people/dt10/research/thomas-08-sample-mean-and-variance.pdf
        # The commented code below is wrong; check it out if required
        #delta_amplitude = (path_throughput - amplitude) / (s+1)
        #updated_amplitude = amplitude + delta_amplitude
        #delta_var = delta_amplitude * (path_throughput- updated_amplitude) / (s+1)
        #variance += delta_var
        #amplitude = updated_amplitude
            
        amplitude = sum_amplitude * P / (s+1)
        variance = sum_amplitude_squared * P / (s+1) - (abs(amplitude)**2)
        if not true_amplitude is None:
            true_variance = sum_true_variance * P / (s+1)
            profile_list.append((path_ndx, path_throughput, 1/P, amplitude,true_variance))
        else:
            profile_list.append((path_ndx, path_throughput, 1/P, amplitude, variance))
        # report variance and amplitude at increments of 1% of M
        if (report_stage==100 and (s+1)==M) or (s+1)==report_ndx or (s+1)==M:
            # for complex variance see
            # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
            if not true_amplitude is None:
                print ('{0}%: {1} samples: estimate = {2:.4f}, true variance = {3:.4f}; est. variance = {4:.4f}'.format(report_stage, s+1, amplitude, true_variance, variance))
            else:
                profile_list.append((s+1,amplitude,variance))
                print ('{0}%: {1} samples: estimate = {2:.4f}, est. variance = {3:.4f}'.format(report_stage, 
                                                                                      s+1,amplitude,
                                                                                      variance))
            report_stage += 1
            report_ndx = math.ceil(M*report_stage/100)
        
    # end iterating over samples

    amplitude = sum_amplitude * P / M
    # for complex variance see
    # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
    variance = sum_amplitude_squared * P / M - (abs(amplitude)**2)
    if not true_amplitude is None:
        true_variance = sum_true_variance * P / M

    return amplitude, non_zero_paths, true_variance if not true_amplitude is None else variance, profile_list


    