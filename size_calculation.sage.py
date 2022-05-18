

# This file was *autogenerated* from the file size_calculation.sage
from sage.all_cmdline import *   # import sage library

_sage_const_0 = Integer(0); _sage_const_6 = Integer(6); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_3 = Integer(3); _sage_const_4 = Integer(4); _sage_const_5 = Integer(5); _sage_const_56 = Integer(56)
import matplotlib
from numpy import poly
import sage.rings
from sage.symbolic.expression_conversions import polynomial

#================================#
# Helper/Conversion Functions #
#================================#

# Converts a Vector in GF(2)^6 to a symbolic expression in 6 variables
#
# Return type - sage.symbolic.expression.Expression
def vec_to_expr(V):
    e = _sage_const_0 
    for i in range(_sage_const_6 ):
        if (V[i] == _sage_const_1 ):
            if (i == _sage_const_0 ):
                e = e + x
            elif (i == _sage_const_1 ):
                e = e + y
            elif (i == _sage_const_2 ):
                e = e + z
            elif (i == _sage_const_3 ):
                e = e + u
            elif (i == _sage_const_4 ):
                e = e + v
            elif (i == _sage_const_5 ):
                e = e + w
    return e

# Act on a cubic homogenous poly by a matrix in GL(6, GF(2))
# 
# Return type - sage.symbolic.expression.Expression
def act_on(M, V):
    # Act M on a basis in GF(2)^6
    x2 = vec_to_expr(M*x_vec)
    y2 = vec_to_expr(M*y_vec)
    z2 = vec_to_expr(M*z_vec)
    u2 = vec_to_expr(M*u_vec)
    v2 = vec_to_expr(M*v_vec)
    w2 = vec_to_expr(M*w_vec)

    # Substitute for symbols in given expression
    V = V.subs(x == x2)
    V = V.subs(y == y2)
    V = V.subs(z == z2)
    V = V.subs(u == u2)
    V = V.subs(v == v2)
    V = V.subs(w == w2)

    # Expand and simplify final expression to get this as a lienear
    # combination of cubic monomials
    return V.expand().simplify()

# Converts a symbolic expression to a polynomial
#
# Return type - 
def convert_to_poly(v):
    p = polynomial(v, base_ring=GF(_sage_const_2 ))
    print(type(p))
    return p

#=============================#
# Linear Algebra Calculations #
#=============================#

# Standard basis for GF(2)^6
x_vec = vector([_sage_const_1 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ])
y_vec = vector([_sage_const_0 ,_sage_const_1 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ])
z_vec = vector([_sage_const_0 ,_sage_const_0 ,_sage_const_1 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ])
u_vec = vector([_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_1 ,_sage_const_0 ,_sage_const_0 ])
v_vec = vector([_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_1 ,_sage_const_0 ])
w_vec = vector([_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_1 ])

#================================#
# Basis of Momonials Calculation #
#================================#

#List of monomials as generic sage expressions
monomial_expr = []
monomial_poly = []

# Ring of polynomials in six variables
R = QQ['x, y, z, u, v, w']; (x, y, z, u, v, w,) = R._first_ngens(6)

#Total degree of monomial term
def total_degree(f):
    return f.degree(x) + f.degree(y) + f.degree(z) + f.degree(u) + f.degree(v) + f.degree(w)

x,y,z,u,v,w = var('x,y,z,u,v,w')

#Extract monomials of degree 3
for i in monomials([x,y,z,u,v,w], [_sage_const_4 ,_sage_const_4 ,_sage_const_4 ,_sage_const_4 ,_sage_const_4 ,_sage_const_4 ]):
    if (total_degree(i) == _sage_const_3 ):
        monomial_expr.append(i(x=x,y=y,z=z,u=u,v=v,w=w))
        monomial_poly.append(convert_to_poly(i(x=x,y=y,z=z,u=u,v=v,w=w)))


#================================#
#      Matrix Calculations       #
#================================#

def insert_vector(M, v, i):
    M = M.insert_row(i, v)
    M = M.delete_rows([i+_sage_const_1 ])
    return M

#General linear group
G = GL(_sage_const_6 , GF(_sage_const_2 ))

def convert_to_56_basis(v):
    vec = vector([_sage_const_0 ]*_sage_const_56 )
    for i in v.monomials():
        vec[monomial_poly.index(i)] = _sage_const_1 
    return vec
        

def convert_to_56(g):
    M = matrix.identity(_sage_const_56 )
    col = _sage_const_0 
    for i in monomial_expr:
        v = convert_to_56_basis(convert_to_poly(act_on(g, i)))
        M = insert_vector(M, v, col)
        col = col + _sage_const_1 
    return M

matrices = []



def compute_num_iso_classes():
    sum = _sage_const_0 
    ord = G.order()
    print(ord)
    for m in G.conjugacy_classes_representatives():
        mat = convert_to_56(m) - matrix.identity(GF(_sage_const_2 ), _sage_const_56 )
        sum = sum + (_sage_const_2 **(kernel(mat).dimension()))*G.conjugacy_class(m).cardinality()
    print(floor(sum/ord))

compute_num_iso_classes()


