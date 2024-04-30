#     weil_search.sage
#
# This file contains the code to call the WeilPolynomials iterator and extract the
# list of potential Weil polynomials of degrees 21 (for K3's) or degree 22
# (for cubic fourfolds) that we are interested in.
#
# As Kedlaya-Sutherland does, we only aim to enumerate the possible transcendental parts.
# The algebraic factors are just products of cyclotomic polynomials, which can be attached
# to the transcendental part after the fact.

from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
import time

datapath = "../../database/zeta_functions/"
polRing.<x> = PolynomialRing(ZZ)

#######
# K3 surfaces

if False:
    for i in range(1, 11):
        ans = [polRing(2)]
        t = time.time()
        c = 0
        
        # NOTE: Kedlaya-Sutherland set the squarefree option to true.
        wp = WeilPolynomials(2*i, 1, sign=1, lead=2, squarefree=False, parallel=True)
        l = [j for j in wp if not j.has_cyclotomic_factor()]
        ans += l
        c += wp.node_count()
        print(len(l), "polynomials added")
        print(c, "nodes enumerated")
        print("time so far: ", time.time() - t, " seconds")

    with open(datapath + "k3f2-lines.txt", "w") as f:
        for i in ans:
            f.write(str(i.list()) + "\n")


#######
# Cubic fourfolds.

if True:
    for i in range(1, 3): #range(1, 12):
        ans = [polRing(2)]
        t = time.time()
        c = 0
        
        # NOTE: Kedlaya-Sutherland set the squarefree option to true.
        # NOTE: On my machine, I have to set parallel to false..
        wp = WeilPolynomials(2*i, 1, sign=1, lead=2, squarefree=False, parallel=False)
        l = [j for j in wp if not j.has_cyclotomic_factor()]
        ans += l
        c += wp.node_count()
        print(len(l), "polynomials added")
        print(c, "nodes enumerated")
        print("time so far: ", time.time() - t, " seconds")

    with open(datapath + "transcendental_weil_polynomials_up_to_deg_22.csv", "a") as f:
        for i in ans:
            f.write(str(i.list()) + "\n")

