

intrinsic RandomPolynomial(R, d) -> RngMPolElt
{Return a random polynomial of degree d from R}
    k := BaseRing(R);
    require IsFinite(k) : "Base ring is not finite.";

    monsd := MonomialsOfDegree(R, d);
    return &+[Random(k) * m : m in monsd];    
end intrinsic;
