
"""
Implements 2 classes
bfunction and btransforms
and a funtion generate_btransform(n, m)
that can be used to generatae all
btransforms fron F_n to F_m

To use this script, do

import affine_maps
bfunction, btransform, generate_btransforms = affine_maps.load()

"""

# Hardcoded max value for K
K = 4 # anything between 1 to 5 will work

def load():
    return bfunction, btransform, generate_btransforms



### START OF MATRIX HELPER FUNCTIONS ###

def mat_rank(mat, n, m):
    base = []
    for j in range(m):
        col = (mat >> n * j) % (1 << n)
        for b in base:
            col = min(b ^ col, col)
        if col:
            base.append(col)
    return len(base)

def mat_transpose(mat, n, m):
    mat2 = 0
    for i in range(n):
        for j in range(m):
            x = (mat >> (j * n + i)) & 1
            mat2 |= x << (i * m + j)
    return mat2, m, n

def mat_to_list(index):
    mat, n, m = mat_mat[index], mat_N[index], mat_M[index]
    A = []
    for i in range(n):
        row = []
        for j in range(m):
            row.append((mat >> (j * n + i)) & 1)
        A.append(row)
    return A

def mat_vec_mult(mat, n, m, vec):
    """ Multiply matrix and vector mod 2 where the vector and rows of the matrix are bit strings """
    vec2 = 0
    for bit in range(m):
        col = (mat >> n * bit) % (1 << n)
        vec2 ^= col * ((vec >> bit) & 1)
    return vec2

def mat_get_permutation(mat, n, m):
    out = [-1] * (1 << m)
    for i in range(1 << m):
        i2 = mat_vec_mult(mat, n, m, i)
        out[i] = i2
    return out


### END OF MATRIX HELPER FUNCTIONS ###


### START OF bfunction AND btransform ###
class bfunction:
    def __init__(self, A, n):
        self.A = A
        self.n = n
    
    def __add__(self, f2):
        A1 = self.A
        A2 = f2.A
        return bfunction(self.A ^ f2.A, self.n)

    def permute(self, perm):
        A = self.A
        A2 = 0
        for i in range(len(perm)):
            x = (A >> perm[i]) & 1
            A2 |= x << i
        return bfunction(A2, len(perm).bit_length() - 1)

    def __index__(self):
        return self.A
    
    def __hash__(self):
        return self.A
    
    def __str__(self):
        out = [(self.A >> i) & 1 for i in range(1 << self.n)]
        return str(out)
    
    def __eq__(self, f2):
        return self.A == f2.A and self.n == f2.n

    def __neg__(self):
        return bfunction(self.A ^ ((1 << (1 << self.n)) - 1), self.n)

    def __call__(self, i):
        return (self.A >> i) & 1

    def popcount(self):
        x = self.A
        ans = 0
        while x:
            x, xrest = divmod(x, 1 << K)
            ans += popcount[xrest]
        return ans

    def supporting_affine_subspace(self):
        if self.n == 0:
            return [0]
        non_zero = []
        half = (1 << self.n) // 2
        for i, hadi in enumerate(hadamards[self.n]):
            if (self + hadi).popcount() != half:
                non_zero.append(i)
        return non_zero

        constant = non_zero.pop()
        linear_basis = []
        for j in non_zero:
            j ^= constant
            for b in linear_basis:
                j = min(b ^ j, j)
            if j:
                linear_basis.append(j)

        affine_subspace = []
        for j in range(1 << (1 << self.n)):
            j2 = j
            j ^= constant
            for b in linear_basis:
                j = min(b ^ j, j)
            if not j:
                affine_subspace.append(j2)
        return affine_subspace


class btransform:
    def __init__(self, Aind, b, beta, c):
        self.Aind = Aind
        self.b = b
        self.beta = beta
        self.c = c

    def transpose(self):
        return btransform(mat_transpose_index[self.Aind], self.beta, self.b, self.c)

    def __call__(self, f):
        Aind = self.Aind
        b = self.b
        beta = self.beta
        c = self.c
        
        n = mat_N[Aind]
        m = mat_M[Aind]

        f = f.permute(offset_permutation[n][b])
        f = f.permute(mat_permutation[Aind])
        f += signed_hadamards[m][~beta if c else beta]

        return f

    def inverse(self):
        Aind = self.Aind
        b = self.b
        beta = self.beta
        c = self.c

        n = mat_N[Aind]
        m = mat_M[Aind]
        assert n == m

        # Weird edge case
        if n == 0:
            return self

        # Invertion formula
        # A -> A^-1
        # b -> A^-1 b
        # beta -> (A^-1)^T beta
        # c -> c xor had_((A^-1)^T beta)(b)

        Aind2 = mat_inverse_index[Aind]
        Aind2T = mat_transpose_index[Aind2]
        b2 = mat_permutation[Aind2][b]
        beta2 = mat_permutation[Aind2T][beta]
        c2 = c ^ hadamards[n][beta2](b)

        return btransform(Aind2, b2, beta2, c2)

def generate_btransforms(n, m):
    # Returns all transforms M:F_n -> F_m

    out = []
    for Aind in mat_dim_index[n][m]:
        for b in range(1 << n):
            for beta in range(1 << m):
                for c in range(2):
                    M = btransform(Aind, b, beta, c)
                    out.append(M)
    return out

### END OF bfunction AND btransform ###



### START OF PRECALCULATIONS USED BY bfunction AND btransform

popcount = [0] * (1 << K)
for i in range(1 << K):
    popcount[i] = popcount[i >> 1] + (i & 1)


# All linear functions 
hadamards = []
for dim in range(K + 1):
    tmp = []
    hadamards.append(tmp)
    for alpha in range(1 << dim):
        f = 0
        for i in range(1 << dim):
            x = popcount[alpha & i] & 1
            f |= x << i
        tmp.append(bfunction(f, dim))

# All affine functions
signed_hadamards = []
for dim in range(K + 1):
    had = hadamards[dim]
    ones = bfunction((1 << (1 << dim)) - 1, dim)

    signed_hadamards.append(had + [h + ones for h in had][::-1])


# Precalc f(x) -> f(x + b) transform to speed up future calculations
offset_permutation = []
for dim in range(K + 1):
    tmp = []
    offset_permutation.append(tmp)
    for b in range(1 << dim):
        P = [-1] * (1 << dim)
        for i in range(1 << dim):
            P[i] = b ^ i
        tmp.append(P)



# Enumerate all full rank n x m matrices (for 1 <= n,m <= K)
# Full rank meaning rank = min(n, m)
mat_N = []
mat_M = []
mat_mat = []
num_mat = 0

for n in range(K + 1):
    for m in range(K + 1):
        full_rank = min(n, m)
        for mat in range(1 << n * m):
            r = mat_rank(mat, n, m)
            if r == full_rank:
                mat_N.append(n)
                mat_M.append(m)
                mat_mat.append(mat)
                num_mat += 1
                
# Precalc f(x) -> f(A x) transform to speed up future calculations
mat_permutation = []
for index in range(num_mat):
    mat, n, m = mat_mat[index], mat_N[index], mat_M[index]
    mat_permutation.append(mat_get_permutation(mat, n, m))

# List of indices to all matrices of a certain size
mat_dim_index = [[[] for _ in range(K + 1)] for _ in range(K + 1)]
for index in range(num_mat):
    mat, n, m = mat_mat[index], mat_N[index], mat_M[index]
    mat_dim_index[n][m].append(index)

# Match each matrix with its transpose
mat_transpose_index = [-1] * num_mat # Will contain index to transpose matrix
mat_index = {}
for index in range(num_mat):
    key = mat, n, m = mat_mat[index], mat_N[index], mat_M[index]
    keyT = matT, nT, mT = mat_transpose(mat, n, m)
    if key == keyT:
        mat_transpose_index[index] = index
    elif keyT in mat_index:
        indexT = mat_index[keyT]
        del mat_index[keyT]
        mat_transpose_index[index] = indexT
        mat_transpose_index[indexT] = index
    else:
        mat_index[key] = index


# Match each invertible matrix with its transpose
mat_inverse_index = [-1] * num_mat # Will contain index to inverse matrix
mat_index = {}
for index in range(num_mat):
    key = mat, n, m = mat_mat[index], mat_N[index], mat_M[index]
    if n != m or n == 0 or m == 0:
        continue
    perm = mat_permutation[index]
    inv_perm = [0] * len(perm)
    for i in range(len(perm)):
        inv_perm[perm[i]] = i
    
    key = tuple(perm)
    key_inv = tuple(inv_perm)

    if key == key_inv:
        mat_inverse_index[index] = index
    elif key_inv in mat_index:
        index_inv = mat_index[key_inv]
        del mat_index[key_inv]
        mat_inverse_index[index] = index_inv
        mat_inverse_index[index_inv] = index
    else:
        mat_index[key] = index

### END OF PRECALCULATIONS USED BY bfunction AND btransform



# DEBUG
#for M in generate_btransforms(0,3):
#    print('M:')
#    print(M.Aind, M.b, M.beta, M.c, mat_permutation[M.Aind])
#    M2 = M.transpose()
#    print(M2.Aind, M2.b, M2.beta, M2.c, mat_permutation[M2.Aind])
#exit()
