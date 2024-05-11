from fractions import Fraction as F

def load_LP(K, INFINITE_SOUNDNESS=False, RESTRICTED=False, ONLY_SHORTEST_EDGES=False): 
    """
    Load LP from file using Gurobi

    Note this does not require Gurobi license.
    It only requires the Gurobi module for Python.
    """
    if K == 4 and RESTRICTED:
        MODEL_FILE = "K" + str(K) + "_compressed_restricted_" + ['finite', 'infinite'][INFINITE_SOUNDNESS] + ".mps.zip"
    elif K == 4 and ONLY_SHORTEST_EDGES:
        MODEL_FILE = "K" + str(K) + "_compressed_" + ['finite', 'infinite'][INFINITE_SOUNDNESS] + "_length_1" + ".mps.zip"
    else:
        MODEL_FILE = "K" + str(K) + "_compressed_" + ['finite', 'infinite'][INFINITE_SOUNDNESS] + ".mps.zip"
    
    import gurobipy as gp
    from gurobipy import GRB
    
    print('Loading File:', MODEL_FILE)
    m = gp.read(MODEL_FILE)
    
    soundness_var = m.getVarByName('soundness_var')
    completeness_var = m.getVarByName('completeness_var')

    gadget_vars = [None]
    for var in m.getVars():
        if var.VarName[0] == 'g':
            gadget_vars.append(var)
    
    return m, soundness_var, completeness_var, gadget_vars

def setup_LP(completeness, K, INFINITE_SOUNDNESS=False, RESTRICTED=False): 
    """
    Setup the LP used to construct the gadget by
    loading the base LP file, and then adding on
    the objective.
    """
    c = 1 - F(completeness)
    assert 2**-K <= c <= 0.5

    m, soundness_var, completeness_var, gadget_vars = load_LP(K, INFINITE_SOUNDNESS, RESTRICTED, c == 2**-K)

    import gurobipy as gp
    from gurobipy import GRB
    
    m.addConstr(completeness_var * c.denominator == c.numerator)
    m.setObjective(-soundness_var, GRB.MINIMIZE)

    return m, soundness_var, completeness_var, gadget_vars


def LP_solve_gurobi(completeness, K, INFINITE_SOUNDNESS=False, RESTRICTED=False):
    """
    Construct a gadget by solving LP using Gurobi.
    Note: Gurobi runs fast but can have significant 
          numerical errors.

    Return value: A list with the capacities of the edges of the gadget 
                  as floating point numbers
    """
    m, soundness_var, completeness_var, gadget_vars = setup_LP(completeness, K, INFINITE_SOUNDNESS, RESTRICTED)

    m.Params.Aggregate = 0
    m.Params.method = 2
    m.Params.NumericFocus = 3
    m.Params.FeasibilityTol = 1e-9
    m.Params.OptimalityTol = 1e-9
    m.Params.IntFeasTol = 1e-9
    m.Params.Threads = 1

    # Crossover turned off for better performance
    m.Params.Crossover = 0
    m.optimize()
    
    return [g.X if g is not None else None for g in gadget_vars]

def LP_solve_qsopt_ex(completeness, K, INFINITE_SOUNDNESS=False, COMPRESSED=False):
    """
    Construct a gadget by solving LP using Gurobi.
    Note: Gurobi runs fast but can have significant 
          numerical errors.

    Return value: A list with the capacities of the edges of the gadget
                  as fractions.Fraction
    """
    m, soundness_var, completeness_var, gadget_vars = setup_LP(completeness, K, INFINITE_SOUNDNESS, COMPRESSED)
    
    random_id = str(hash('a'))
    MODEL_FILE =  'tmp' + random_id + '.mps'
    SOLUTION_FILE = 'tmp' + random_id + '.sol'

    print('Creating temporary file', MODEL_FILE)
    m.write(MODEL_FILE)
    
    print('Calling QSOPT_EX and storing the output in', SOLUTION_FILE)
    import subprocess
    subprocess.run(["esolver", "-O", SOLUTION_FILE, MODEL_FILE])

    print('Parsing output file from QSOPT_Ex')
    name_to_index = {gadget_vars[i].VarName:i for i in range(len(gadget_vars)) if gadget_vars[i] is not None}
    for i in range(1, len(gadget_vars)):
        gadget_vars[i] = 0

    lines = list(open(SOLUTION_FILE))
    assert lines[0].rstrip() == "status = OPTIMAL"
    for line in lines:
        line = line.split()
        if line and line[-1][-1] == ':' and line[-1] != 'VARS:':
            break
        if line[0] in name_to_index:
            gadget_vars[name_to_index[line[0]]] = F(line[2])
    
    print('Removing temporary files')
    import os
    try:
        os.remove(MODEL_FILE)
    except OSError:
        pass
    
    try:
        os.remove(SOLUTION_FILE)
    except OSError:
        pass

    return gadget_vars


# Construct a new gadget using either the LP solver qsopt_ex (rational solver)
# or Gurobi (floating point solver)

# Note that Gurobi is used to load the LP, so Gurobipy needs to be installed

K = 4
INFINITE = True

gadget = LP_solve_qsopt_ex(1 - 2**-K, K, INFINITE)
#gadget = LP_solve_gurobi(1 - 2**-K, K, INFINITE)

import analyse_gadget
analyse_gadget.print_gadget(gadget, K)
