
from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials


# NOTE: Sutherland-Kedlaya Theorem 1 has a small typo about the normalization.

# WeilPolynomials(21, 1)


# Need code to

def extract_unit_eigvals(f):
    R = parent(f)
    x = gen(R)
    r = 0; e = 0;
    while r == 0:
        _, r = f.quo_rem((x-1)^(e+1))
        e += 1

    e -= 1
    # print(f.quo_rem((x-1)^e), e)
    return f.quo_rem((x-1)^e)[0], e


d = 21
q = 4   # Absolute value squared of Frobenius eigenvalues.

# We really can iterate through all of them.
W = WeilPolynomials(d, q)
I = iter(W)

mine = [f for f in I]
blah = [extract_unit_eigvals(f) for f in mine]
loblaw = [(g(1), g(-1)) for g,e in blah]


print({b for a,b in loblaw})
print({a for a,b in loblaw})

