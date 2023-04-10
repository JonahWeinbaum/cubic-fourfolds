AttachSpec("CubicLib.spec");

// Fetch values Jack computed.
signs := eval Read("../../database/zeta_functions/signsequence");

// Get the first 100 examples with each sign.
firstmin1 := [];
i := 1;
while #firstmin1 lt 100 do
    if signs[i] eq -1 then
        Append(~firstmin1, i);
    end if;
    i := i + 1;
end while;

first1 := [];
i := 1;
while #first1 lt 100 do
    if signs[i] eq 1 then
        Append(~first1, i);
    end if;
    i := i + 1;
end while;

// Load cubics and compare signs.
cubics := LoadCubicOrbitData(: Flat, OnlySmooth);

for i in first1[1..5] do
    sign := FunctionalEquationSign(cubics[i]);
    assert sign eq 1;
end for;
