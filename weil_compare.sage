import time

DATA_PATH = "/Users/f004584/Dropbox (Dartmouth College)/Prj-cubics/cubic-F2/database/zeta_functions/"

if False:
    # Polynomial ring.
    R.<x> = PolynomialRing(Integers())

    # Load possible transcendental K3 polynomials
    with open(DATA_PATH + "k3f2-lines.txt", 'r') as F:
        k3arch = [R(eval(ll)) for ll in F.readlines()]

    # Load computed cubic zeta functions
    t = time.time()
    cubic_indices = []
    zeta_funcs = []
    with open(DATA_PATH + "zetacoeffs.csv", 'r') as F:
        for ll in F:
            tup = sage_eval(ll.replace('<', '(').replace('>', ')'))
            a, b = tup
            twob = [2*x for x in b]
            cubic_indices.append(a)
            zeta_funcs.append(R(twob))

    # Should be about 230 seconds.
    print("Load time: {}".format(time.time() - t))


cyclo = list(filter((lambda f : f.degree() <= 24),
                    (R.cyclotomic_polynomial(i) for i in range(1, 1000))))

# Remove the algebraic part.
def transcendental_part(f):
    if f == 0: return 0
    q = f
    for g in cyclo:
        while True:
            nextq, r = q.quo_rem(g)
            if r != 0: break
            q = nextq

    return q

# Conduct comparison.
transc_iter = (transcendental_part(f) for f in zeta_funcs)

flag = True
for f in transc_iter:
    if f not in k3arch and f.degree() <= 20:
        print(f)
        flag = False
        break
