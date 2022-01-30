(* ::Package:: *)

(* ::Text:: *)
(*Output table parameters*)


koutmax=1; (* koutmax will be the upper decade (1,10,100, etc) *)
koutmin=0.001;
OutputLogInk=1; (*Spaced even in log-k*)


(* ::Text:: *)
(*Input table*)


InpkT=PSLT[[All,1]];
kinMinInput=InpkT[[1]];
kinMaxInput=InpkT[[Length@InpkT]];
Print["Input kmin = ",kinMinInput, ",  input kmax = ",kinMaxInput]



(* ::Text:: *)
(*Output Table*)


Proposedkoutmax =koutmax;
Proposedkoutmin =koutmin;


(*Proposedkoutmax =1000;
Proposedkoutmin =0.0000001;*)

koutmax = Min[Proposedkoutmax,kinMaxInput];
koutmin = Max[Proposedkoutmin,kinMinInput];


If[Proposedkoutmin< kinMinInput, Print["Warning: The output kmin should be considerably larger than input linear pk kmin"];
koutmin= kinMinInput/0.06];

If[Proposedkoutmax > kinMaxInput, Print["Warning: The output kmax should be considerably smaller than the k max of input linear pk"];
koutmax= kinMaxInput/18.];


If[OutputLogInk==1,
decadesout = Log10[koutmax/koutmin];
pperdecadeout = 40;
(*pperdecadeout = 60;*)
sizekT = Ceiling[decadesout]*pperdecadeout+1;
delta = Ceiling[decadesout]/(sizekT - 1);
kT = Table[10^(Log10[koutmin] + (ii - 1)*delta), {ii, 1, sizekT}];
sizekT = Dimensions[kT][[1]];
strkT="Output files log-k spaced: kmin = "<>ToString[First[kT]]<>",  kmax = "<> ToString[Last[kT]]<>",  number of ks: "<> ToString[sizekT];
Print[strkT];
];



(* ::Text:: *)
(*Inner integration table*)


(*pperdecade = 26;*)  (*Precision*)
(*pperdecade = 50;*)
pperdecade = 26;

pmin = Max[InpkT[[1]], 0.01 kT[[1]]];
pmax = Min[InpkT[[Length@InpkT]], 16 kT[[sizekT]]];
decades = Log10[pmax/pmin];
sizepT = pperdecade Ceiling[decades];
deltap = decades/sizepT;
pT = Table[10^(Log10[pmin] + pii*deltap), {pii, 0, sizepT}];
sizepT = Dimensions[pT][[1]];
strpT="Inner integration pmin = "<> ToString[pT[[1]]]<> ",  pmax = "<> ToString[pT[[sizepT]]]<>",  number of ps: "<> ToString[sizepT];    
Print[strpT];


If[ttft==True,faopen = OpenAppend[outlogfile];WriteLine[faopen, "..."];WriteLine[faopen, strkT];WriteLine[faopen, strpT];Close[faopen]];

