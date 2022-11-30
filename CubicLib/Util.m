

intrinsic RandomPolynomial(R, d) -> RngMPolElt
{Return a random polynomial of degree d from R}
    k := BaseRing(R);
    require IsFinite(k) : "Base ring is not finite.";

    monsd := MonomialsOfDegree(R, d);
    return &+[Random(k) * m : m in monsd];    
end intrinsic;

/////////////////////////////////////////////////
//
// Dictionary manipulations.
//
/////////////////////////////////////////////////

intrinsic CacheMap(func, A) -> Any
{Basically implement a `map` idiom, via caching.}
    Cache := AssociativeArray();
    valset := {x : x in A};
    
    for x in valset do
        Cache[x] := func(x);
    end for;

    case Type(A):        
    when SeqEnum: return [Cache[x] : x in A];
    when List: return [* Cache[x] : x in A *];
    when Assoc:
        B := AssociativeArray();
        for k in Keys(A) do
            B[k] := Cache[A[k]];
        end for;
        return B;
    else:
    error "Not implemented.";
    end case;        
end intrinsic;


intrinsic CacheFilter(func, A) -> Any
{Basically implement a `filter` idiom, via caching.}
    Cache := AssociativeArray();
    valset := {x : x in A};
    
    for x in valset do
        Cache[x] := func(x);
    end for;

    case Type(A):        
    when SeqEnum: return [x : x in A | Cache[x]];
    when List: return [* x : x in A | Cache[x] *];
    when Assoc:
        B := AssociativeArray();
        for k in Keys(A) do
            v := A[k];
            if Cache[v] then
                B[k] := v;
            end if;
        end for;
        return B;
    else:
    error "Not implemented.";
    end case;
end intrinsic;


intrinsic InverseAssociativeArray(A::Assoc) -> Assoc
{Given an associative array, return the associative array with the keys and values reversed.
(Keys are placed in a set.)}

    B := AssociativeArray();
    for k in Keys(A) do
        v := A[k];
        if v in Keys(B) then
            S := B[v];
            Include(~S, v);
            B[v] := S;
        else
            B[v] := {k};
        end if;
    end for;
    
    return B;
end intrinsic;

