
load "Resultants.m2"

for i from 59 to 62 do
(

inputfilenames = ("x00", "x01", "x02", "x03", "x04", "x05", "x06", "x07", "x08", "x09", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30", "x31", "x32", "x33", "x34", "x35", "x36", "x37", "x38", "x39", "x40", "x41", "x42", "x43", "x44", "x45", "x46", "x47", "x48", "x49", "x50", "x51", "x52", "x53", "x54", "x55", "x56", "x57", "x58", "x59", "x60", "x61", "x62" );


outputfilenames = ("out00", "out01", "out02", "out03", "out04", "out05", "out06", "out07", "out08", "out09", "out10", "out11", "out12", "out13", "out14", "out15", "out16", "out17", "out18", "out19", "out20", "out21", "out22", "out23", "out24", "out25", "out26", "out27", "out28", "out29", "out30", "out31", "out32", "out33", "out34", "out35", "out36", "out37", "out38", "out39", "out40", "out41", "out42", "out43", "out44", "out45", "out46", "out47", "out48", "out49", "out50", "out51", "out52", "out53", "out54", "out55", "out56", "out57", "out58", "out59", "out60", "out61", "out62" );



R = ZZ[u,v,w,x,y,z];

theseq = (); 

scanLines(l -> theseq = append(theseq, discriminant value l % 8), "disctests/" | inputfilenames#i);

outputfilenames#i << theseq << endl << close;
print i;
clearAll;
)


