Attach("../cubic-fourfolds/CubicLib.m");
load "data-processing/read-data-lines-index.m";
k := FiniteField(2);
F4 := FiniteField(4);
basF4 := Basis(F4);
k6 := RSpace(k, 6);
k4 := RSpace(k, 4);

G_4 := GL(6, F4);
P5_4<[x]> := ProjectiveSpace(F4, 5);
R_4 := CoordinateRing(P5_4);

V_4, Bit_4 := GModule(G_4, R_4, 3);

ptsonlines := AssociativeArray();
monoevaluated := AssociativeArray();

//Loop over all lines
for form in lines do
    nullsp := NullspaceOfTranspose(form);
    bas := Basis(nullsp);
    nullsp := ExtendField(nullsp, F4);
    pt1 := bas[1];
    pt2 := bas[2];
    pt3 := (nullsp!pt1)*1 + (nullsp!pt2)*F4.1;

    ptsonlines[form] := <pt1, pt2, pt3>;
    
    monoevaluated[form] := <0,0,0,0>;

    for mono in Basis(V_4) do 

        eval1 := Evaluate(mono @@ Bit_4, Eltseq(pt1));
        eval2 := Evaluate(mono @@ Bit_4, Eltseq(pt2));
        eval3 := Evaluate(mono @@ Bit_4, Eltseq(pt3));

        eval3bas1 := Eltseq(eval3)[1];
        eval3bas2 := Eltseq(eval3)[2];

        monoevaluated[form][1] :=BitwiseOr (ShiftLeft(monoevaluated[form][1], 1), Integers()!eval1);
        monoevaluated[form][2] :=BitwiseOr (ShiftLeft(monoevaluated[form][2], 1), Integers()!eval2);
        monoevaluated[form][3] :=BitwiseOr (ShiftLeft(monoevaluated[form][3], 1), Integers()!eval3bas1);
         monoevaluated[form][4] :=BitwiseOr (ShiftLeft(monoevaluated[form][4], 1), Integers()!eval3bas2);
          end for;
end for;

function paritycalc(bitstr)
    parity := false;
    while bitstr ne 0 do
        parity := not parity;
        bitstr := BitwiseAnd(bitstr, (bitstr - 1));
    end while;  
    return parity;
end function;

function LinesThrough(cubic)
    f := CubicToInt(cubic);
    lines := [];
    for form in Keys(monoevaluated) do 
        evals := monoevaluated[form];

        eval1 := BitwiseAnd(monoevaluated[form][1], f);

        if  paritycalc(eval1) then 
            continue;
        end if;

        eval2 := BitwiseAnd(monoevaluated[form][2], f);

        if  paritycalc(eval2) then 
            continue;
        end if;

        eval3 := BitwiseAnd(monoevaluated[form][3], f);

        if  paritycalc(eval3) then 
            continue;
        end if;

        eval4 := BitwiseAnd(monoevaluated[form][4], f);

        if  paritycalc(eval4) then 
            continue;
        end if;

        lines cat:= [form];

    end for;

    return lines;
end function;

orbdata := LoadCubicOrbitData(: Flat:=true);


k := FiniteField(2);
G := GL(6, k);
P5<[x]> := ProjectiveSpace(k, 5);
R := CoordinateRing(P5);

V, Bit := GModule(G, R, 3);

fname := Sprintf("lines-data.data");
file := Open(fname, "w");
i := 1;
for f  in orbdata do
    seq := serialize_line_info(LinesThrough(f@@Bit));
    WriteBytes(file, seq);
    if i/1000 in Integers() then
        i;
    end if;
    i := i+1;
end for;

/*// Compute lines for all cubics in database
linesdata := [];

for i in [1..85] do
o := orbdata[i];
time linesdata cat:= [[LinesThrough(f) : f in o ] ];
 i;
end for;

 Ar := AssociativeArray();
 
 linesseq := Setseq(lines);
 
for i in [1..#linesseq] do
Ar[linesseq[i]] := i;   
end for;

for i in [1..85] do                              
fname := Sprintf("../data/lines-%o.data", i);
 PrintFile(fname, linesdata[i]);
end for;
*/


