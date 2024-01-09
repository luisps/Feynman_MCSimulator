from hash_paths import my_hash_inverse, my_hash 
from gates import gate_evaluate_w_sample, gate_evaluate_w
import random, math
# Importing complex mean and variance functions
from complex_mean_variance import c_mean, c_var_true_mean

def Importance (num_qubits, num_layers, layers, initial_bin_str, final_bin_str, P, N, M,
               true_amplitude = None):
    
    DEBUG_MODE = False
    
    non_zero_paths = 0
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

    # iterate all samples
    for s in range(M):
        if DEBUG_MODE:
            print ()
            print ('Sample {}'.format(s))
        path = []
        
        path_throughput = complex(1.,0)  # path initially has weight 1

        for l in range(num_layers):    
            # the last layer is deterministically traversed
            # by connecting to the final state
            lastLayer = True if l==num_layers-1 else False
            
            if l==0:
                input_bin_str = initial_bin_str
            else:
                input_bin_str = output_bin_str

            path.append(int(input_bin_str,2))
            
            if (lastLayer):
                output_bin_str = final_bin_str
                path.append(int(final_bin_str,2))
        
            if DEBUG_MODE:
                print ('--Layer {0}: input state:{1}'.format(l,input_bin_str))

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
                
                if not lastLayer:
                    output_qbs, w, prob = gate_evaluate_w_sample (gate_name, input_qbs, rnd=rnd,
                                                              unitary=unitary, pdfs=pdfs)
                    # update ouput_bin_str
                    for qb_ndx,qb in enumerate(qubits):
                        output_bin_list[num_qubits-1-qb] = str(output_qbs[qb_ndx])
                else:
                    # 2 qubit gates
                    if gate_name in ['cx','cz','cp']:
                        output_qbs = (int(output_bin_str[num_qubits-1-qubits[0]]),
                                     int(output_bin_str[num_qubits-1-qubits[1]]))
                    else: 
                        output_qbs = (int(output_bin_str[num_qubits-1-qubits[0]]), )
                        prob = 1.

                    w = gate_evaluate_w (gate_name, input_qbs, output_qbs, 
                                         params=params, unitary=unitary)

                if DEBUG_MODE:
                    print ('gate {0}, qubits {1}, input {2}, ouput {3}, weight {4}, probability {5}'.format(gate[0], gate[1], input_qbs, output_qbs, w, prob))

                    
                if abs(w)==0. or prob==0.:
                    layer_weight = complex (0.,0.)
                    path_throughput = complex (0.,0.)
                    break     # early termination of loop over qubits
                else:
                    layer_weight *= w
                    layer_prob *= prob
            # end iterating over the layer qubits
            if abs(layer_weight) == 0.:
                #path_throughput = 0.
                path_throughput = complex (0.,0.)
                break        # early termination of loop over layers
            else:
                path_throughput *= (layer_weight/layer_prob)

            if not lastLayer:
                output_bin_str = "".join(output_bin_list)
            if DEBUG_MODE:
                print ('output bin str {2}'.format(l,input_bin_str, output_bin_str))
                print ('weight = {0}, prob = {1}'.format(layer_weight, layer_prob))
                print ('Path throughput up to now= {0}'.format(path_throughput))
                print ()
        # end iterating over layers
            
        if not true_amplitude is None:
            sum_true_variance += abs(path_throughput - true_amplitude)**2 
        if abs(path_throughput) != 0. :
            sum_amplitude += path_throughput
            # for complex variance see
            # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
            sum_amplitude_squared += (abs(path_throughput)**2) 
            non_zero_paths += 1
    
        if DEBUG_MODE:
            if abs(path_throughput) != 0.:
                print ('Non zero path: ',path)
                print ('Path throughput = {0}'.format(path_throughput))
            print ()
        
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
            
        # report variance and amplitude at increments of 1% of M
        if (report_stage==100 and (s+1)==M) or (s+1)==report_ndx or (s+1)==M:
            amplitude = sum_amplitude / (s+1)
            # for complex variance see
            # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
            variance = sum_amplitude_squared / (s+1) - (abs(amplitude)**2)
            if not true_amplitude is None:
                true_variance = sum_true_variance / (s+1)
                profile_list.append((s+1,amplitude,true_variance))
                print ('{0}%: {1} samples: estimate = {2:.4f}, true variance = {3:.4f}; est. variance = {4:.4f}'.format(report_stage, s+1, amplitude, true_variance, variance))
            else:
                profile_list.append((s+1,amplitude,variance))
                print ('{0}%: {1} samples: estimate = {2:.4f}, est. variance = {3:.4f}'.format(report_stage, 
                                                                                      s+1,amplitude,
                                                                                      variance))
            report_stage += 1
            report_ndx = math.ceil(M*report_stage/100)
        
    # end iterating over samples
    
    amplitude = sum_amplitude / M
    # for complex variance see
    # https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
    variance = sum_amplitude_squared / M - (abs(amplitude)**2)
    if not true_amplitude is None:
        true_variance = sum_true_variance / M

    return amplitude, non_zero_paths, true_variance if not true_amplitude is None else variance, profile_list


    