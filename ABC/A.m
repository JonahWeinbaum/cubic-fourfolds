// Task (A) from Jack's email.

// (A) Which singular cubics have a line contained entirely in the
// smooth locus?? Obviously no such lines for the most singular cubics
// out there, e.g. the triple hyperplane x^3=0. But for cubics with
// isolated singularities, is there such a line?

AttachSpec("../CubicLib/CubicLib.spec");
cubics := LoadCubicOrbitData(: Quick, Flat);

P6 := Proj(Parent(cubics[1]));
mons1 := [P6.i : i in [1..Dimension(P6)+1]];

function LineToScheme(PP, mat)
    mons1 := [PP.i : i in [1..Dimension(PP)+1]];
    return Scheme(PP, [Polynomial(Eltseq(row), mons1) : row in Rows(mat)]);
end function;

function CheckJackA(f)    
    lines := LinesThrough(f);
    for ll in lines do
        L := LineToScheme(P6, ll);
        if IsEmpty(L meet sing) then
            return true;
        end if;
    end for;

    // All lines meet the singular locus.
    return false;
end function;

fname := DatabaseDirectory() * "linear_subspaces/smooth_lines_singular_cubics/jackA.data";
N := #cubics;

for i in [1..N] do
    f := cubics[i];
    if not IsSmooth(f) and CheckJackA(f) then
        Write(fname, Sprintf("%o", i));
    end if;
end for;

// time alist := [CheckJackA(f) : f in cubics | not IsSmooth(f)];
