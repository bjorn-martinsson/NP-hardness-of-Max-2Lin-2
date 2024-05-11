This repo contains the (Python3) scripts used to construct and verify the gadgets used in the paper

On the NP-hardness approximation curve for 2 Max-2Lin(2)

The list of all gadgets can be found in list_of_gadgets.py.

The verification script (that computes soundness + completeness) can be found in analyse_gadget.py . For perfomance reasons, it is recommended to use PyPy3 when running this script.

Generating new gadgets can be done by running generate_new_gadgets.py
This requires;
* Gurobipy + qsopt_ex (https://www.math.uwaterloo.ca/~bico/qsopt/ex/)
* Note: License for Gurobi is not needed for running qsopt_ex
