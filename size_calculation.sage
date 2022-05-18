import matplotlib
from numpy import poly
import sage.rings
from sage.symbolic.expression_conversions import polynomial

#=============================#
# Linear Algebra Calculations #
#=============================#

# Basis for GF(2)^6
x_vec = vector([1,0,0,0,0,0])
y_vec = vector([0,1,0,0,0,0])
z_vec = vector([0,0,1,0,0,0])
u_vec = vector([0,0,0,1,0,0])
v_vec = vector([0,0,0,0,1,0])
w_vec = vector([0,0,0,0,0,1])

# Converts a Vector in GF(2)^6 into an expression in 6 variables
def vec_to_expr(V):
    e = 0
    for i in range(6):
        if (V[i] == 1):
            if (i == 0):
                e = e + x
            elif (i == 1):
                e = e + y
            elif (i == 2):
                e = e + z
            elif (i == 3):
                e = e + u
            elif (i == 4):
                e = e + v
            elif (i == 5):
                e = e + w
    return e

def act_on(M, V):
    x2 = vec_to_expr(M*x_vec)
    y2 = vec_to_expr(M*y_vec)
    z2 = vec_to_expr(M*z_vec)
    u2 = vec_to_expr(M*u_vec)
    v2 = vec_to_expr(M*v_vec)
    w2 = vec_to_expr(M*w_vec)

    V = V.subs(x == x2)
    V = V.subs(y == y2)
    V = V.subs(z == z2)
    V = V.subs(u == u2)
    V = V.subs(v == v2)
    V = V.subs(w == w2)

    return V.expand().simplify()

def convert_to_poly(v):
    p = polynomial(v, base_ring=GF(2))
    return p
    

#================================#
# Basis of Momonials Calculation #
#================================#

#List of monomials as generic sage expressions
monomial_expr = []
monomial_poly = []

# Ring of polynomials in six variables
R.<x,y,z,u,v,w> = QQ[]

#Total degree of monomial term
def total_degree(f):
    return f.degree(x) + f.degree(y) + f.degree(z) + f.degree(u) + f.degree(v) + f.degree(w)

x,y,z,u,v,w = var('x,y,z,u,v,w')

#Extract monomials of degree 3
for i in monomials([x,y,z,u,v,w], [4,4,4,4,4,4]):
    if (total_degree(i) == 3):
        monomial_expr.append(i(x=x,y=y,z=z,u=u,v=v,w=w))
        monomial_poly.append(convert_to_poly(i(x=x,y=y,z=z,u=u,v=v,w=w)))


#================================#
#      Matrix Calculations       #
#================================#

def insert_vector(M, v, i):
    M = M.insert_row(i, v)
    M = M.delete_rows([i+1])
    return M

#General linear group
G = GL(6, GF(2))

def convert_to_56_basis(v):
    vec = vector([0]*56)
    for i in v.monomials():
        vec[monomial_poly.index(i)] = 1
    return vec
        

def convert_to_56(g):
    M = matrix.identity(56)
    col = 0
    for i in monomial_expr:
        v = convert_to_56_basis(convert_to_poly(act_on(g, i)))
        M = insert_vector(M, v, col)
        col = col + 1
    return M

matrices = []



def compute_num_iso_classes():
    sum = 0
    ord = G.order()
    print(ord)
    for m in G.conjugacy_classes_representatives():
        mat = convert_to_56(m) - matrix.identity(GF(2), 56)
        sum = sum + (2^(kernel(mat).dimension()))*G.conjugacy_class(m).cardinality()
    print(floor(sum/ord))

compute_num_iso_classes()

