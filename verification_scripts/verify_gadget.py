import fractions

import verification_scripts.affine_maps as affine_maps
bfunction, btransform, generate_btransforms = affine_maps.load()


"""
Verify soundness and completeness using integer precision
"""
def verify(K, gadget_vars, infinite):
    
    # Switch to big integers
    tmp = next(x for x in gadget_vars[::-1] if x)
    if isinstance(tmp, int):
        pass
    elif isinstance(tmp, float):
        # Make
        HUGE = 1
        while True:
            if all((x * HUGE).is_integer() if x else True for x in gadget_vars):
                gadget_vars = [int(x * HUGE) for x in gadget_vars]
                break
            else:
                HUGE *= 2
    elif isinstance(tmp, fractions.Fraction):
        def gcd(a,b):
            while b:
                a,b = b,a%b
            return a
        def lcm(a,b):
            if not a:
                return b
            if not b:
                return a
            return a * (b // gcd(a,b))

        g = 0
        for x in gadget_vars[1:]:
            g = lcm(g, x.as_integer_ratio()[1])
        
        def make_rat_into_int(x):
            a,b = x.as_integer_ratio()
            return a * (g // b)
        gadget_vars = [make_rat_into_int(x) if x is not None else None for x in gadget_vars]
    else:
        # Unable to transform gadget_vars to integer
        assert False

    import verification_scripts.hypercube_precalculations as hypercube_precalculations
    g_data, edge_data, get_edge_orbit_ind, get_entire_edge_orbit = hypercube_precalculations.load(K, infinite)

    [
        g_representatives,
        g_compressed_classes,
        g_num_compressed_classes,
    ] = g_data

    [
        num_edge_orbits,
        edge_orbits_representatives,
        edge_orbits_size,
        compressed_edge_fragments,

    ] = edge_data

    
    # total_cap == 1
    total_capacity = 0
    for i in range(1, num_edge_orbits + 1):
        f1,f2 = edge_orbits_representatives[i]
        # Skip edges on the form (f,f) or (f,-f)
        if f1 == f2 or f1 == -f2:
            continue
        total_capacity += gadget_vars[i] * edge_orbits_size[i]

    
    # Completeness is the mean normalized hamming distance
    comp = 0
    for i in range(1, num_edge_orbits + 1):
        f1,f2 = edge_orbits_representatives[i]
        # Skip edges on the form (f,f) or (f,-f)
        if f1 == f2 or f1 == -f2:
            continue
        comp += gadget_vars[i] * (f1 + f2).popcount() * edge_orbits_size[i]
    comp_normalization_factor = (1 << K) * total_capacity


    import verification_scripts.dinics as dinics
    flow_graph = dinics.init()
    


    ### FLOWS
    
    # Go over all g_equivalence_classes
    for f1_compressed_class in range(1, g_num_compressed_classes + 1):
        for f2_compressed_class in range(1, g_num_compressed_classes + 1):
            
            if f1_compressed_class <= f2_compressed_class:
                continue
            
            # The capacity of edge (f1_equivalence_class, f2_equivalence_class)
            # is a sum of capacities over edges the original graph
            s = 0

            di = compressed_edge_fragments[f1_compressed_class][f2_compressed_class]
            for edge_orbit_index, c in di.items():
                sum_of_cap = c * gadget_vars[edge_orbit_index]
                s += sum_of_cap

            # Add s to the compressed flow graph
            if s:
                flow_graph.add_edge(f1_compressed_class, f2_compressed_class, s, s)

    # Soundness is max flow of the compressed graph
    sound = flow_graph.solve(1, 2)

        
    ### Calculate % usage of edge capacities 
    gadget_flows = [0] * len(gadget_vars)
    flows = flow_graph.calc_flow()[::-1]

    compressed_flow = [(u,v,fractions.Fraction(f)/total_capacity) for u,v,f in flows]

    # Go over all g_equivalence_classes
    for f1_compressed_class in range(1, g_num_compressed_classes + 1):
        for f2_compressed_class in range(1, g_num_compressed_classes + 1):
            
            if f1_compressed_class <= f2_compressed_class:
                continue
            
            # The capacity of edge (f1_equivalence_class, f2_equivalence_class)
            # is a sum of capacities over edges the original graph
            s = 0

            di = compressed_edge_fragments[f1_compressed_class][f2_compressed_class]
            for edge_orbit_index, c in di.items():
                sum_of_cap = c * gadget_vars[edge_orbit_index]
                s += sum_of_cap

            # Add s was added to the compressed flow graph
            if s:
                u,v,f = flows.pop()
                for edge_orbit_index, c in di.items():
                    total_cap = c * gadget_vars[edge_orbit_index]
                    gadget_flows[edge_orbit_index] += fractions.Fraction(abs(f) * total_cap, s)

    # Normalize values
    normalization_factor = (1 << (1 << K))
    gadget_total_cap = [cap * count if cap is not None else None for cap, count in zip(gadget_vars, edge_orbits_size)]
    gadget_flows = [fractions.Fraction(f, (cap * normalization_factor)) if cap else float('NaN') for f,cap in zip(gadget_flows, gadget_total_cap)]
   
    
    ### Calculate % of total capacity edge edge orbit is using
    gadget_relative_vars = [fractions.Fraction(cap * count, total_capacity) if cap else float('NaN') for cap, count in zip(gadget_vars, edge_orbits_size)]


    sound_normalization_factor = (1 << (1 << K)) * total_capacity 
    return 1 - fractions.Fraction(sound,sound_normalization_factor), 1 - fractions.Fraction(comp,comp_normalization_factor), compressed_flow
    #return fractions.Fraction(sound,sound_normalization_factor), fractions.Fraction(comp,comp_normalization_factor), gadget_flows, gadget_relative_vars, compressed_flow


def generate_witness(K, gadget_vars):
    s1, c1, witness1 = verify(K, gadget_vars, False)
    s2, c2, witness2 = verify(K, gadget_vars, True)

    assert c1 == c2
    assert s1 >= s2

    return c1, s1, s2, witness1, witness2
