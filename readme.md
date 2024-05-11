This repo contains the (Python3) scripts used to construct and verify the gadgets used in the paper:

**On the NP-hardness approximation curve for 2 Max-2Lin(2)**

### list_of_gadgets.py.
Contains a list of all gadgets used in the Paper.

### analyse_gadget.py 
The verification script (that computes soundness + completeness). For perfomance reasons, it is recommended to use PyPy3 when running this script.

### generate_new_gadgets.py
Script for generating a gadget with a fixed completeness using an LP solver, Gurobi (fast and uses floating point numbers) or qsopt_ex (slower is exact).
Running this script requires installing;
* Gurobipy
* (Optionally) qsopt_ex (https://www.math.uwaterloo.ca/~bico/qsopt/ex/)
* (Optionally) License for Gurobi
