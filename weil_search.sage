
from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials

# File location
# /private/var/tmp/sage-9.6-current/local/var/lib/sage/venv-python3.10.3/lib/python3.10/site-packages/sage/rings/polynomial/weil/weil_polynomials.pyx


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


###########################
# Sutherland's actualy search code. Might be relevant.

# import time
# polRing.<x> = PolynomialRing(ZZ)
# ans = [polRing(2)]
# t = time.time()
# c = 0
# for i in range(1, 11):
#     wp = WeilPolynomials(2*i, 1, sign=1, lead=2, squarefree=True, parallel=True)
#     l = [j for j in wp if not j.has_cyclotomic_factor()]
#     ans += l
#     c += wp.node_count()
#     print(len(l), "polynomials added")
#     print(c, "nodes enumerated")
#     print("time so far: ", time.time() - t, " seconds")

# with open("k3f2-lines.txt", "w") as f:
#     for i in ans:
#         f.write(str(i.list()) + "\n")


=#
