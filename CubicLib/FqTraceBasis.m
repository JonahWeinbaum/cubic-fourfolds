
// Functions to go to and from our C++ representation and Magma's representation.
// Also, contains some code to generate various tables.

function CppRep(x : Bits:=false)
    // x -- Finite field element.
    s := ChangeUniverse(Eltseq(x), Integers());

    if Bits then
        return s;
    else
        return Seqint(s, 2);
    end if;
end function;


function MagmaRep(k, x)
    // k -- Finite field GF(2^m)
    // x -- positive integer less than q.
    s := Intseq(x, 2);
    if #s eq 0 then return Zero(k); end if;
    
    _<t> := PolynomialRing(GF(2));
    xrep := &+[t^(i-1) * s[i] : i in [1..#s]];
    return k ! xrep;
end function;

// Create basis for trace zero subspace. Also return the z such that z^ + z = tr
// TODO: Decide whether to keep the placeholder.
function TraceBasis(k)
    // NOTE: We reverse to compute the Echelon form in the right order.
    A := Matrix([Reverse(Eltseq(b^2+b)) : b in Basis(k)]);
    
    E, U := EchelonForm(A);
    newBasis := [k ! Eltseq(row) : row in Rows(U)];

    return [z^2 + z : z in newBasis], newBasis;
end function;


function generateCppTraceBases()

    ptvs := []; // Pretrace values
    tbs  := [];
    for i in [1..22] do
        tb, b := TraceBasis(GF(2^i));
        cppb := [CppRep(x) : x in b];
        ptvs[i] := cppb;

        cpptb := [CppRep(x) : x in tb];
        tbs[i] := cpptb;
    end for;

    // Note the placeholder. C is 0-indexed.
    return [[0]] cat tbs, [[0]] cat ptvs;
end function;
