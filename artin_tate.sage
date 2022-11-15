
import os, time

DATA_PATH = os.path.join("/Users/f004584/Dropbox (Dartmouth College)/Prj-cubics/cubic-F2/",
                         "database/zeta_functions/")

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



##########
# Generic functionality for removing/saturating divisors.
def remove_factors(f, fact_list):
    if f == 0: return 0
    q = f
    for g in fact_list:
        while True:
            nextq, r = q.quo_rem(g)
            if r != 0: break
            q = nextq
    return q


def remove_factor(f, g):
    return remove_factors(f, [g])


##########
# Transcendental parts and Artin-Tate

cyclo = list(filter((lambda f : f.degree() <= 24),
                    (R.cyclotomic_polynomial(i) for i in range(1, 1000))))

# Remove the algebraic part.
def transcendental_part(f):
    return remove_factors(g, cyclo)


# Given a characteristic polynomial of Frobenius, divide out by all factors
# of (1-x) and evaluate at 2. 
def zeta_special_value(f, v):
    g = remove_factor(f, 1-x)
    return g(2)
    

########################################
# Main script start.

# Conduct comparison.
import multiprocessing as mp

def check_transc_task(f):
    if f not in k3arch and f.degree() <= 20:
        print(f)
        return False
    else:
        return True
        
transc_iter = (transcendental_part(f) for f in zeta_funcs)

if __name__ == '__main__':
    with mp.Pool(40) as P:
        weil_check_lst = P.map(check_transc_task, transc_iter)


