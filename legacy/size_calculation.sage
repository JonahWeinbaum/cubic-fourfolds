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

# Converts a symbolic expression to a polynomial
#
# Return type - sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular
def convert_to_poly(v):
    p = polynomial(v, base_ring=GF(2))
    return p

#=============================#
# Linear Algebra Calculations #
#=============================#

# Standard basis for GF(2)^6
x_vec = vector([1,0,0,0,0,0])
y_vec = vector([0,1,0,0,0,0])
z_vec = vector([0,0,1,0,0,0])
u_vec = vector([0,0,0,1,0,0])
v_vec = vector([0,0,0,0,1,0])
w_vec = vector([0,0,0,0,0,1])

# Act on a cubic homogenous poly by a matrix in GL(6, GF(2))
# 
# Return type - sage.symbolic.expression.Expression
def act_on(M, vec):
    # Act M on a basis in GF(2)^6
    x2 = vec_to_expr(M*x_vec)
    y2 = vec_to_expr(M*y_vec)
    z2 = vec_to_expr(M*z_vec)
    u2 = vec_to_expr(M*u_vec)
    v2 = vec_to_expr(M*v_vec)
    w2 = vec_to_expr(M*w_vec)

    # Substitute for symbols in given expression
    vec = vec.subs(x == x2, y == y2, z == z2, u == u2, v == v2, w == w2)

    # Expand and simplify final expression to get this as a linear
    # combination of cubic monomials
    return vec.expand().simplify()

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

#===============================#
#      Matrix Calculations      #
#===============================#

# Inserts a vector into the matrix at row i
def insert_vector(M, vec, i):
    M = M.insert_row(i, vec)
    M = M.delete_rows([i+1])
    return M

# General linear group
G = GL(6, GF(2))

# Converts a cubic homogenous equation as a polynomial
# to a vector in GF(2)^56
def convert_to_vec(p):
    vec = vector([0]*56)
    for i in p.monomials():
        vec[monomial_poly.index(i)] = 1
    return vec

# Converts a matrix in GL(6, GF(2)) to its induced matrix in GL(56, GF(2))
def convert_to_56(m):
    mat = matrix.identity(56)
    row = 0
    for i in monomial_expr:
        vec = convert_to_vec(convert_to_poly(act_on(m, i)))
        mat = insert_vector(mat, vec, row)
        row = row + 1
    return mat

# Solves burnside equation to find number of orbits
def burnside_eqn():
    sum = 0
    ord = G.order()

    for m in G.conjugacy_classes_representatives():
        mat = convert_to_56(m) - matrix.identity(GF(2), 56)
        sum = sum + (2^(kernel(mat).dimension()))*G.conjugacy_class(m).cardinality()

    return(floor(sum/ord))

print(burnside_eqn())

