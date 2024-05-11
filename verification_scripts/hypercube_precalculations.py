import verification_scripts.affine_maps as affine_maps
bfunction, btransform, generate_btransforms = affine_maps.load()

def load(K, infinite = False):
    
    ### START OF edge_orbit HELPER FUNCTIONS ###

    def get_edge_orbit_index(f1, f2):
        Mrepr = f_btransform_to_representative[f1]
        f1 = Mrepr(f1)
        f2 = Mrepr(f2)
        
        f1_representative_ind = f_representatives.index(f1)
        return edge_orbits[f1_representative_ind][f2]

    def get_entire_edge_orbit(f1, f2):
        found = set()
        ret = []
        for M in generate_btransforms(K, K):
            Mf1 = M(f1)
            Mf2 = M(f2)
            
            a = Mf1.__index__()
            b = Mf2.__index__()

            if a > b:
                a,b = b,a
                Mf1,Mf2 = Mf2,Mf1

            if (a,b) not in found:
                found.add((a,b))
                ret.append((Mf1, Mf2))
        return ret

    ### END OF edge_orbit HELPER FUNCTIONS ###


    if infinite:
        PICKLE_FILE_NAME = 'verification_scripts/PRECALC_FOR_K=%d_WITH_INFINITE_SOUNDNESS.pickle' % K
    else:
        PICKLE_FILE_NAME = 'verification_scripts/PRECALC_FOR_K=%d_WITH_FINITE_SOUNDNESS.pickle' % K

    from os.path import exists
    import pickle

    if not exists(PICKLE_FILE_NAME):
        ### START OF MATRIX PRECALCULATIONS ###
        
        ### END OF MATRIX PRECALCULATIONS ###

        ### START OF g ANALYSIS PRECALCULATIONS ###
        g_representatives = []
        for dim in range(K + 1):
            tmp = []
            g_representatives.append(tmp)

            found = [0] * (1 << (1 << dim))
            transforms = generate_btransforms(dim, dim)
            for g in range(1 << (1 << dim)):
                g = bfunction(g, dim)
                if found[g]:
                    continue
                tmp.append(g)
                for M in transforms:
                    g2 = M(g)
                    found[g2] += 1
        
        # The compressed classes are the nodes in the compressed graph
        g_compressed_classes = []
        for dim in range(K + 1):
            tmp = []
            g_compressed_classes.append(tmp)
            for g in g_representatives[dim]:
                tmp.append([0] * (1 << (1 << dim)))

        g_num_compressed_classes = 0
        for dim in range(K + 1):
            for g_ind, g in enumerate(g_representatives[dim]):    
                print('Working with dim and gind', dim, g_ind, flush=True)
                found = g_compressed_classes[dim][g_ind]
                
                precalculated_fixpoint_transforms = [None] * dim
                for dim2 in range(dim, K + 1):
                    tmp = []
                    precalculated_fixpoint_transforms.append(tmp)

                    transforms = generate_btransforms(dim, dim2)
                    for gprim_ind, gprim in enumerate(g_representatives[dim2]):
                        fixpoint_transforms = [M for M in transforms if g == M.transpose()(gprim)]
                        tmp.append(fixpoint_transforms)
                        
                for f in range(1 << (1 << dim)):
                    f = bfunction(f, dim)
                    if found[f]:
                        continue
                    g_num_compressed_classes += 1
                    
                    for dim2 in range(dim, dim + 1 if (not infinite) and dim != 0 else K + 1):
                        for gprim_ind, gprim in enumerate(g_representatives[dim2]):
                            fixpoint_transforms = precalculated_fixpoint_transforms[dim2][gprim_ind]
                            found2 = g_compressed_classes[dim2][gprim_ind]
                                            
                            for M in fixpoint_transforms:
                                f2 = M(f)
                                assert found2[f2] == 0 or found2[f2] == g_num_compressed_classes
#                                if found2[f2] == 0:
#                                    #print(f2, flush=True)
#                                    print('wgat§')
#                                    print('check_new', g_num_compressed_classes, [gprim(i) ^ f2(0) for i in f2.supporting_affine_subspace()])
#                                    print('check_orig', g_num_compressed_classes, [g(i) ^ f(0) for i in f.supporting_affine_subspace()])
#                                    print(gprim,g, M.Aind, M.b, M.beta, M.c)
                                found2[f2] = g_num_compressed_classes
#                    # DEBUG
#                    if g_num_compressed_classes == 3:
#                        exit()
#
#        exit()
        
        # The g_equivalence_classes of source/sink placement g
        # Are the orbits of actions M(f), over M s.t. M^T(g) = g.
        g_equivalence_classes = []
        for g in g_representatives[K]:
            g_equivalence_classes.append([0] * (1 << (1 << K)))

        g_num_equivalence_classes = 0
        g_size_equivalence_classes = []

        transforms = generate_btransforms(K, K)
        for g_ind, g in enumerate(g_representatives[K]):
            fixpoint_transforms = [M for M in transforms if g == M.transpose()(g)]
            g_size_equivalence_classes.append(len(transforms) // len(fixpoint_transforms))
            found = g_equivalence_classes[g_ind]
            
            for f in range(1 << (1 << K)):
                f = bfunction(f, K)
                if found[f]:
                    continue
                g_num_equivalence_classes += 1
                
                for M in fixpoint_transforms:
                    f2 = M(f)
                    assert found[f2] == 0 or found[f2] == g_num_equivalence_classes
                    found[f2] = g_num_equivalence_classes

        ### START OF edge orbit ANALYSIS ###

        f_representatives = []
        f_btransform_to_representative = [None] * (1 << (1 << K)) # ...[f] = tranform M such that f -> f_repr

        transforms = generate_btransforms(K, K)
        for f in range(1 << (1 << K)):
            f = bfunction(f, K)
            if f_btransform_to_representative[f] is not None:
                continue
            f_representatives.append(f)
            for M in transforms:
                f2 = M(f)
                if f_btransform_to_representative[f2] is None:
                    f_btransform_to_representative[f2] = M.inverse()


        """
        To find the edge orbits, we will only consider
        edges with at least one endpoint being a f_representative
        """
        
        num_edge_orbits = 0
        fixpoint_transforms = [[M for M in transforms if f == M(f)] for f in f_representatives]
        edge_orbits = [[0] * (1 << (1 << K)) for f in f_representatives]
        edge_orbits_representatives = [None]
        for f_representative_ind, f_representative in enumerate(f_representatives):
            found = edge_orbits[f_representative_ind]
            for f in range(1 << (1 << K)):
                f = bfunction(f, K)
                if found[f]:
                    continue
                
                """
                We want to find all edges in (f_representative, f) edge orbit
                such that at least one end point is a f_representative.

                There are 2 cases
                
                1. Edges on the form (f_representative, M(f)), over
                all M such that M(f_representative) = f_representative.
                
                2. Edges on the form (M(f_representative), M(f)) over
                all M such that M(f) = f's representative.

                """

                num_edge_orbits += 1
                edge_orbits_representatives.append((f_representative, f))

                # Case 1
                for M in fixpoint_transforms[f_representative_ind]:
                    f2 = M(f)
                    found[f2] = num_edge_orbits
                
                # Case 2
                Mrepr = f_btransform_to_representative[f]
                h_representative_ind = f_representatives.index(Mrepr(f))
                assert h_representative_ind >= 0
               
                h = Mrepr(f_representative)
                found2 = edge_orbits[h_representative_ind]
                
                for M in fixpoint_transforms[h_representative_ind]:
                    h2 = M(h)
                    found2[h2] = num_edge_orbits

        """
        All the edge orbits have been identified, using the f_representatives
        So now we will count how many edges are in each orbit,
        """

        edge_orbits_size = [0] * (num_edge_orbits + 1)
        for f_representative_ind, f_representative in enumerate(f_representatives):
            multiplicity = len(transforms) // len(fixpoint_transforms[f_representative_ind])
            found = edge_orbits[f_representative_ind]
            for f in range(1 << (1 << K)):
                f = bfunction(f, K)
                if f == f_representative:
                    edge_orbits_size[found[f]] += 2 * multiplicity 
                else:
                    edge_orbits_size[found[f]] += multiplicity 
        # Counted each edge twice, so divide by 2
        edge_orbits_size = [x//2 for x in edge_orbits_size]


        size = (1 << (1 << K))
        assert sum(edge_orbits_size) == size * (size + 1) // 2

        ### END OF edge_orbit ANALYSIS PRECALCULATIONS ###

        # Create a map from g_equivalence_class to compressed_classes (the nodes in the compressed graph)
        # In the finite case this map is identity
        #
        # The reasons why is done is that it is easier to
        # to describe edges between g_equivalence_classes
        # than edges between compressed_classes.
        g_equivalence_class_to_compressed_graph = mapper = {} 
        for g_ind,g in enumerate(g_representatives[K]):
            equivalence_classes = g_equivalence_classes[g_ind]
            compressed_classes = g_compressed_classes[K][g_ind]

            for f in range(1 << (1 << K)):
                f = bfunction(f, K)
                equivalence_class = equivalence_classes[f]
                compressed_class = compressed_classes[f]
                
                assert equivalence_class not in mapper or mapper[equivalence_class] == compressed_class
                mapper[equivalence_class] = compressed_class

        # In the compressed graph, the nodes are given "compressed" classes
        compressed_edge_fragments = [[{} for _ in range(g_num_compressed_classes + 1)] for _ in range(g_num_compressed_classes + 1)]
        for g_ind,g in enumerate(g_representatives[K]):
            print('Finding all edge fragments between equivalence classes in the compressed graph', g_ind, flush=True)
            g_count = g_size_equivalence_classes[g_ind]
            equivalence_classes = g_equivalence_classes[g_ind]

            count = [0] * (g_num_compressed_classes + 1)
            for f in range(1 << (1 << K)):
                f = bfunction(f, K)
                f_equivalence_class = equivalence_classes[f]
                f_compressed_class = g_equivalence_class_to_compressed_graph[f_equivalence_class]
                count[f_compressed_class] += 1
            
            seen = [0] * (g_num_compressed_classes + 1)
            for f1 in range(1 << (1 << K)):
                f1 = bfunction(f1, K)
                f1_equivalence_class = equivalence_classes[f1]
                f1_compressed_class = g_equivalence_class_to_compressed_graph[f1_equivalence_class]
                
                if seen[f1_compressed_class]:
                    continue
                seen[f1_compressed_class] = 1
                f1_count = count[f1_compressed_class]
                multiplicity = f1_count * g_count

                for f2 in range(1 << (1 << K)):
                    f2 = bfunction(f2, K)
                    f2_equivalence_class = equivalence_classes[f2]
                    f2_compressed_class = g_equivalence_class_to_compressed_graph[f2_equivalence_class]
                    
                    if f1 == -f2 or f1_compressed_class <= f2_compressed_class:
                        continue
                    
                    edge_orbit_index = get_edge_orbit_index(f1, f2)
                    ff1, ff2 = edge_orbits_representatives[edge_orbit_index]
                    assert ff1 != ff2

                    di = compressed_edge_fragments[f1_compressed_class][f2_compressed_class]
                    key = edge_orbit_index
                    if key in di:
                        di[key] += multiplicity
                    else:
                        di[key] = multiplicity

        



        ### START DUMP PICKLE PRECALCULATIONS ###
        
        g_data = [
            g_representatives,
            g_compressed_classes,
            g_num_compressed_classes,
        ]

        edge_data = [
            f_representatives,
            f_btransform_to_representative,
            num_edge_orbits,
            edge_orbits,
            edge_orbits_representatives,
            edge_orbits_size,
            compressed_edge_fragments,
        ]
        
        with open(PICKLE_FILE_NAME, 'wb') as my_file:
            pickle.dump([g_data, edge_data],my_file)

        ### END DUMP PICKLE PRECALCULATIONS
    else:
        ### START LOAD PICKLE PRECALCULATIONS ###
        with open(PICKLE_FILE_NAME, 'rb') as my_file:
            g_data, edge_data = pickle.load(my_file)
            
            [
                g_representatives,
                g_compressed_classes,
                g_num_compressed_classes,
            ] = g_data
            
            
            [
                f_representatives,
                f_btransform_to_representative,
                num_edge_orbits,
                edge_orbits,
                edge_orbits_representatives,
                edge_orbits_size,
                compressed_edge_fragments,
            ] = edge_data
        ### END LOAD PICKLE PRECALCULATIONS ###

    



    g_data = [
        g_representatives,
        g_compressed_classes,
        g_num_compressed_classes,
    ]

    edge_data = [
        num_edge_orbits,
        edge_orbits_representatives,
        edge_orbits_size,
        compressed_edge_fragments,
    ]

    return g_data, edge_data, get_edge_orbit_index, get_entire_edge_orbit
