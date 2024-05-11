from fractions import Fraction as F

def parse_gadget(edge_weights, K):
    if K == 2:
        # 2 benicifial edge orbits
        benificial_edge_orbits = [2, 3]
        gadget_vars = [0] * 9
    elif K == 3:
        # 4 benificial edge orbits
        benificial_edge_orbits = [2, 13, 3, 5]
        gadget_vars = [0] * 33
    elif K == 4:
        # 33 benificial edge orbits
        benificial_edge_orbits = [2, 35, 122, 323, 325, 772, 3, 125, 128, 328, 332, 708, 784, 4, 131, 138, 432, 5, 6, 7, 10]
        gadget_vars = [0] * 1078
    else:
        assert False
    
    assert len(edge_weights) == len(benificial_edge_orbits)
    for w,eind in zip(edge_weights, benificial_edge_orbits):
        gadget_vars[eind] = w
    return gadget_vars

def print_gadget(gadget_vars, K):
    from verification_scripts.verify_gadget import generate_witness

    c,s1,s2,_,_ = generate_witness(K, gadget_vars)
    
    print('Information about the gadget with K=%d:' %K)

    print('Completeness:', c, float(c))
    print('Relaxed soundness:', s1, float(s1))
    print('Infinity relaxed soundness', s2, float(s2))
    print('Min-deletion ratio:', (1-s1)/(1-c), float((1-s1)/(1-c)))
    print('Min-deletion ratio (infinity):', (1 - s2)/(1 - c), float((1 - s2)/(1 - c)))


    if K == 2:
        # 2 benicifial edge orbits
        benificial_edge_orbits = [2, 3]
    elif K == 3:
        # 4 benificial edge orbits
        benificial_edge_orbits = [2, 13, 3, 5]
    elif K == 4:
        # 33 benificial edge orbits
        benificial_edge_orbits = [2, 35, 122, 323, 325, 772, 3, 125, 128, 328, 332, 708, 784, 4, 131, 138, 432, 5, 6, 7, 10]

    print('Gadget edge weights:', [gadget_vars[eind] for eind in benificial_edge_orbits])

# Here is an example of how to load a gadget into 
# the verification tool. 
# These are the gadgets shown in table 2 and table 3

def example1():
    K = 4
    edge_weights = [0] * 21
    edge_weights[2] = F(5461, 969636864)
    edge_weights[3] = F(17007,1616061440) 
    edge_weights[4] = F(437,  404015360)
    edge_weights[5] = F(19,   92346368)
    edge_weights[6] = F(13,   215360)
    gadget_vars = parse_gadget(edge_weights, K)
    print_gadget(gadget_vars, K)


def example2():
    K = 4
    edge_weights = [0] * 21
    edge_weights[2] = F(4899,   799089790)
    edge_weights[3] = F(11843,  799089790)
    edge_weights[4] = F(1427,   1917815496)
    edge_weights[5] = F(1427,   19178154960)
    edge_weights[6] = F(6094929,102283493120)
    gadget_vars = parse_gadget(edge_weights, K)
    print_gadget(gadget_vars, K)

def verify_all():
    import list_of_gadgets

    K2    = list_of_gadgets.gadgets_K2_optimized_for_relaxed_soundness 
    K2inf = list_of_gadgets.gadgets_K2_optimized_for_infinity_relaxed_soundness 
    
    K3    = list_of_gadgets.gadgets_K3_optimized_for_relaxed_soundness 
    K3inf = list_of_gadgets.gadgets_K3_optimized_for_infinity_relaxed_soundness 
    
    K4    = list_of_gadgets.gadgets_K4_optimized_for_relaxed_soundness 
    K4inf = list_of_gadgets.gadgets_K4_optimized_for_infinity_relaxed_soundness

    for K,gadget_list in [(2, K2), (2, K2inf), (3, K3), (3, K3inf), (4, K4), (4, K4inf)]:
        for edge_weights in gadget_list:
            gadget_vars = parse_gadget(edge_weights, K)
            print_gadget(gadget_vars, K)
            print()
