(* ::Package:: *)

(* ::Section:: *)
(*P=PNW + PW*)


AllT = Import[FileIni];
AllTheaders = AllT[[1]];
AllT = Drop[AllT, 1];
Dimensions[AllT];

kT = AllT[[All, 1]];
PSL = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 2]]}]];
fk = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 3]]}]];
P22dd = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 4]]}]];
P22dt = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 5]]}]];
P22tt = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 6]]}]];
P13dd = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 7]]}]];
P13dt = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 8]]}]];
P13tt = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 9]]}]];
I1udd1A = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 10]]}]];
I2uud1A = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 11]]}]];
I2uud2A = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 12]]}]];
I3uuu2A = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 13]]}]];
I3uuu3A = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 14]]}]];
I2uudd1B = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 15]]}]];
I2uudd2B = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 16]]}]];
I3uuud1B = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 17]]}]];
I3uuud2B = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 18]]}]];
I3uuud3B = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 19]]}]];
I4uuuu1B = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 20]]}]];
I4uuuu2B = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 21]]}]];
I4uuuu3B = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 22]]}]];
I4uuuu4B = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 23]]}]];
I2uudd1C = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 24]]}]];
I2uudd2C = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 25]]}]];
I3uuud1C = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 26]]}]];
I3uuud2C = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 27]]}]];
I3uuud3C = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 28]]}]];
I4uuuu1C = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 29]]}]];
I4uuuu2C = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 30]]}]];
I4uuuu3C = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 31]]}]];
I4uuuu4C = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 32]]}]];

(*Bias terms*)

Pb2b1 = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 36]]}]];
Pbs2b1 = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 37]]}]];
Pb22 = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 38]]}]];
Pb2s2 = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 39]]}]];
Pbs22 = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 40]]}]];
Pb2t = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 41]]}]];
Pbs2t = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 42]]}]];
sigma32pk = Interpolation[Transpose[{AllT[[All, 1]], AllT[[All, 43]]}]];


sigma2 = AllT[[1, 33]];
sigma2v = AllT[[1, 34]];
f0 = AllT[[1, 35]];
Print["f(k=0) = ", f0];


PddL[k_] := PSL[k];
Pdd[k_] := PddL[k] + P22dd[k] + P13dd[k];
Pddloop[k_] := P22dd[k] + P13dd[k];
PdtL[k_] := (fk[k]/f0) PSL[k];
Pdt[k_] := PdtL[k] + P22dt[k] + P13dt[k];
Pdtloop[k_] := P22dt[k] + P13dt[k];
PttL[k_] := (fk[k]/f0)^2 PSL[k];
Ptt[k_] := PttL[k] + P22tt[k] + P13tt[k];
Pttloop[k_] := P22tt[k] + P13tt[k];




PddXloop[k_, mu_, b1_, b2_, bs2_, b3nl_] := 
  b1^2 Pddloop[k] + 2 b1 b2 Pb2b1[k] + 2 b1 bs2 Pbs2b1[k] + 
   b2^2 Pb22[k] + bs2^2 Pbs22[k] + 2 b2 bs2 Pb2s2[k] + 
   2 b1 b3nl sigma32pk[k];

PdtXloop[k_, mu_, b1_, b2_, bs2_, b3nl_] :=  
 b1 Pdtloop[k] + b2 Pb2t[k] + bs2 Pbs2t[k] + b3nl fk[k]/f0 sigma32pk[k];
PttXloop[k_, mu_, b1_, b2_, bs2_, b3nl_] :=  Pttloop[k] ;


Af[k_, mu_, f0_] := 
  f0 mu^2 I1udd1A[k] + f0^2 ( mu^2 I2uud1A[k] + mu^4 I2uud2A[k]) + 
   f0^3 (mu^4 I3uuu2A[k] + mu^6 I3uuu3A[k]);
Bf[k_, mu_, f0_] := 
  f0^2 ( mu^2 I2uudd1B[k] + mu^4 I2uudd2B[k]) + 
   f0^3 (mu^2 I3uuud1B[k] + mu^4 I3uuud2B[k] + mu^6 I3uuud3B[k]) + 
   f0^4 (mu^2 I4uuuu1B[k] + mu^4 I4uuuu2B[k] + mu^6 I4uuuu3B[k] + 
      mu^8 I4uuuu4B[k]);
Cf[k_, mu_, f0_] := 
  f0^2 ( mu^2 I2uudd1C[k] + mu^4 I2uudd2C[k]) + 
   f0^3 (mu^2 I3uuud1C[k] + mu^4 I3uuud2C[k] + mu^6 I3uuud3C[k]) + 
   f0^4 (mu^2 I4uuuu1C[k] + mu^4 I4uuuu2C[k] + mu^6 I4uuuu3C[k] + 
      mu^8 I4uuuu4C[k]);


ATNS[k_, mu_, b1_] := b1^3 Af[k, mu, f0/b1];
BTNS[k_, mu_, b1_] := b1^4 Bf[k, mu, f0/b1];
CTNS[k_, mu_, b1_] := b1^4 Cf[k, mu, f0/b1];

GTNS[k_, mu_, b1_] := -f0^2 mu^2 k^2 sigma2v b1^2 PddL[k] - 2 f0^3 mu^4 k^2 sigma2v b1 PdtL[k] - f0^4 mu^6 k^2 sigma2v  PttL[k];

pPMEloop[k_, mu_, b1_, b2_, bs2_, b3nl_] := 
  PddXloop[k, mu, b1, b2, bs2, b3nl] + 
   2 mu^2 f0 PdtXloop[k, mu, b1, b2, bs2, b3nl]  + 
   mu^4 f0^2 PttXloop[k, mu, b1, b2, bs2, b3nl] + ATNS[k, mu, b1] + 
   BTNS[k, mu, b1] + CTNS[k, mu, b1] + GTNS[k, mu, b1];

PKaiser[k_,mu_,b1_]:=(b1+mu^2 fk[k])^2 PddL[k];
Pctilde[k_,mu_,b1_,ctilde_]:= ctilde (mu k f0)^4 sigma2v^2 PKaiser[k,mu,b1];

Peft[k_,mu_,b1_,alpha0_,alpha2_,alpha4_,ctilde_]:=(alpha0+alpha2 mu^2+alpha4 mu^4)k^2 PddL[k]+Pctilde[k,mu,b1,ctilde];
Pshot[k_,mu_,alphashot0_,alphashot2_,PshotP_]:=PshotP(alphashot0+alphashot2 k^2 mu^2);

Ploop[k_, mu_, b1_, b2_, bs2_, b3nl_,alpha0_,alpha2_,alpha4_,ctilde_,alphashot0_,alphashot2_,PshotP_]:=pPMEloop[k, mu, b1, b2, bs2, b3nl]+Peft[k,mu,b1,alpha0,alpha2,alpha4,ctilde]+Pshot[k,mu,alphashot0,alphashot2,PshotP];


Print["test Ploop(k=0.2,mu=0.8)=",Ploop[0.2, 0.8, 1.2, 0.2, 0.01, 0.3,0.2, 0.8, 1.2, 0.2, 0.01, 0.3,0]];


(* ::Section:: *)
(*P=PNW*)


AllTNW = Import[FileIniNW];
AllTheadersNW = AllTNW[[1]];
AllTNW = Drop[AllTNW, 1];
Dimensions[AllTNW];


kTNW = AllTNW[[All, 1]];
PSLNW = Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 2]]}]];
fkNW = Interpolation[
   Transpose[{AllTNW[[All, 1]], AllTNW[[All, 3]]}]];
P22ddNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 4]]}]];
P22dtNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 5]]}]];
P22ttNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 6]]}]];
P13ddNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 7]]}]];
P13dtNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 8]]}]];
P13ttNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 9]]}]];
I1udd1ANW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 10]]}]];
I2uud1ANW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 11]]}]];
I2uud2ANW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 12]]}]];
I3uuu2ANW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 13]]}]];
I3uuu3ANW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 14]]}]];
I2uudd1BNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 15]]}]];
I2uudd2BNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 16]]}]];
I3uuud1BNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 17]]}]];
I3uuud2BNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 18]]}]];
I3uuud3BNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 19]]}]];
I4uuuu1BNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 20]]}]];
I4uuuu2BNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 21]]}]];
I4uuuu3BNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 22]]}]];
I4uuuu4BNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 23]]}]];
I2uudd1CNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 24]]}]];
I2uudd2CNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 25]]}]];
I3uuud1CNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 26]]}]];
I3uuud2CNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 27]]}]];
I3uuud3CNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 28]]}]];
I4uuuu1CNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 29]]}]];
I4uuuu2CNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 30]]}]];
I4uuuu3CNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 31]]}]];
I4uuuu4CNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 32]]}]];


(*Bias terms*)

Pb2b1NW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 36]]}]];
Pbs2b1NW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 37]]}]];
Pb22NW = Interpolation[
   Transpose[{AllTNW[[All, 1]], AllTNW[[All, 38]]}]];
Pb2s2NW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 39]]}]];
Pbs22NW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 40]]}]];
Pb2tNW = Interpolation[
   Transpose[{AllTNW[[All, 1]], AllTNW[[All, 41]]}]];
Pbs2tNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 42]]}]];
sigma32pkNW = 
  Interpolation[Transpose[{AllTNW[[All, 1]], AllTNW[[All, 43]]}]];


sigma2NW = AllTNW[[1, 33]];
sigma2vNW = AllTNW[[1, 34]];
f0NW = AllTNW[[1, 35]];
Print["f(k=0) = ", f0NW];

PddLNW[k_] := PSLNW[k];
PddNW[k_] := PddLNW[k] + P22ddNW[k] + P13ddNW[k];
PddloopNW[k_] := P22ddNW[k] + P13ddNW[k];
PdtLNW[k_] := (fkNW[k]/f0NW) PSLNW[k];
PdtNW[k_] := PdtLNW[k] + P22dtNW[k] + P13dtNW[k];
PdtloopNW[k_] := P22dtNW[k] + P13dtNW[k];
PttLNW[k_] := (fkNW[k]/f0NW)^2 PSLNW[k];
PttNW[k_] := PttLNW[k] + P22ttNW[k] + P13ttNW[k];
PttloopNW[k_] := P22ttNW[k] + P13ttNW[k];


PddXloopNW[k_, mu_, b1_, b2_, bs2_, b3nl_] := 
  b1^2 PddloopNW[k] + 2 b1 b2 Pb2b1NW[k] + 2 b1 bs2 Pbs2b1NW[k] + 
   b2^2 Pb22NW[k] + bs2^2 Pbs22NW[k] + 2 b2 bs2 Pb2s2NW[k] + 
   2 b1 b3nl sigma32pkNW[k];

PdtXloopNW[k_, mu_, b1_, b2_, bs2_, b3nl_] :=  
 b1 PdtloopNW[k] + b2 Pb2tNW[k] + bs2 Pbs2tNW[k] + 
  b3nl fkNW[k]/f0NW sigma32pkNW[k];
PttXloopNW[k_, mu_, b1_, b2_, bs2_, b3nl_] :=  PttloopNW[k] ;


AfNW[k_, mu_, f0_] := 
  f0 mu^2 I1udd1ANW[k] + 
   f0^2 ( mu^2 I2uud1ANW[k] + mu^4 I2uud2ANW[k]) + 
   f0^3 (mu^4 I3uuu2ANW[k] + mu^6 I3uuu3ANW[k]);
BfNW[k_, mu_, f0_] := 
  f0^2 ( mu^2 I2uudd1BNW[k] + mu^4 I2uudd2BNW[k]) + 
   f0^3 (mu^2 I3uuud1BNW[k] + mu^4 I3uuud2BNW[k] + 
      mu^6 I3uuud3BNW[k]) + 
   f0^4 (mu^2 I4uuuu1BNW[k] + mu^4 I4uuuu2BNW[k] + 
      mu^6 I4uuuu3BNW[k] + mu^8 I4uuuu4BNW[k]);
CfNW[k_, mu_, f0_] := 
  f0^2 ( mu^2 I2uudd1CNW[k] + mu^4 I2uudd2CNW[k]) + 
   f0^3 (mu^2 I3uuud1CNW[k] + mu^4 I3uuud2CNW[k] + 
      mu^6 I3uuud3CNW[k]) + 
   f0^4 (mu^2 I4uuuu1CNW[k] + mu^4 I4uuuu2CNW[k] + 
      mu^6 I4uuuu3CNW[k] + mu^8 I4uuuu4CNW[k]);


ATNSNW[k_, mu_, b1_] := b1^3 AfNW[k, mu, f0NW/b1];
BTNSNW[k_, mu_, b1_] := b1^4 BfNW[k, mu, f0NW/b1];
CTNSNW[k_, mu_, b1_] := b1^4 CfNW[k, mu, f0NW/b1];

GTNSNW[k_, mu_, b1_] := -f0NW^2 mu^2 k^2 sigma2vNW b1^2 PddLNW[k] - 
   2 f0NW^3 mu^4 k^2 sigma2vNW b1 PdtLNW[k] - 
   f0NW^4 mu^6 k^2 sigma2vNW  PttLNW[k];

pPMEloopNW[k_, mu_, b1_, b2_, bs2_, b3nl_] := 
  PddXloopNW[k, mu, b1, b2, bs2, b3nl] + 
   2 mu^2 f0NW PdtXloopNW[k, mu, b1, b2, bs2, b3nl]  + 
   mu^4 f0NW^2 PttXloopNW[k, mu, b1, b2, bs2, b3nl] + 
   ATNSNW[k, mu, b1] + BTNSNW[k, mu, b1] + CTNSNW[k, mu, b1] + 
   GTNSNW[k, mu, b1];


PKaiserNW[k_,mu_,b1_]:=(b1+mu^2 fkNW[k])^2 PddLNW[k];
PctildeNW[k_,mu_,b1_,ctilde_]:= ctilde (mu k f0NW)^4 sigma2vNW^2 PKaiserNW[k,mu,b1];

PeftNW[k_,mu_,b1_,alpha0_,alpha2_,alpha4_,ctilde_]:=(alpha0+alpha2 mu^2+alpha4 mu^4)k^2 PddLNW[k]+PctildeNW[k,mu,b1,ctilde];
PshotNW[k_,mu_,alphashot0_,alphashot2_,PshotP_]:=PshotP(alphashot0+alphashot2 k^2 mu^2);

PloopNW[k_, mu_, b1_, b2_, bs2_, b3nl_,alpha0_,alpha2_,alpha4_,ctilde_,alphashot0_,alphashot2_,PshotP_]:=pPMEloopNW[k, mu, b1, b2, bs2, b3nl]+PeftNW[k,mu,b1,alpha0,alpha2,alpha4,ctilde]+PshotNW[k,mu,alphashot0,alphashot2,PshotP];


Print["test PloopNW(k=0.2,mu=0.8)=",PloopNW[0.2, 0.8, 1.2, 0.2, 0.01, 0.3,0.2, 0.8, 1.2, 0.2, 0.01, 0.3,0]];


(* ::Section:: *)
(*Damping*)


kosc = 1/104.
ks = 0.4;
Sigma2 = 1/(6 Pi^2)NIntegrate[PSLNW[k] (1 - SphericalBesselJ[0, k/kosc] + 
      2 SphericalBesselJ[2, k/kosc]), {k, 0.0001, ks}];
deltaSigma2 = 
 1/(2 Pi^2)NIntegrate[PSLNW[k] SphericalBesselJ[2, k/kosc], {k, 0.0001, ks}];

Sigma2Total[mu_] := (1 + f0 mu^2 (2 + f0NW)) Sigma2 + f0^2 mu^2 (mu^2 - 1) deltaSigma2;


(* ::Section:: *)
(*NLO EFT*)


DNLOell0[k_, b1_] := b1^2/5 + (2 b1 fk[k])/7 + fk[k]^2/9;
DNLOell2[k_, b1_] := (4 b1^2)/7 + (20 b1 fk[k])/21 + (40 fk[k]^2)/99;
DNLOell4[k_, b1_] := (8 b1^2)/35 + (48 b1 fk[k])/77 + (48 fk[k]^2)/
   143;


(* ::Section:: *)
(*IR-Resummed pk*)


pk[k_, mu_, b1_, b2_, bs2_, b3nl_, alpha0_, alpha2_, alpha4_, ctilde_,
    alphashot0_, alphashot2_, 
   PshotP_] := (b1 + fk[k] mu^2)^2 (PSLNW[k] + 
      Exp[-k^2 Sigma2Total[mu]] (PSL[k] - PSLNW[k]) (1 + 
         k^2 Sigma2Total[mu])) + 
   Exp[-k^2 Sigma2Total[mu]] Ploop[k, mu, b1, b2, bs2, b3nl, alpha0, 
     alpha2, alpha4, ctilde, alphashot0, alphashot2, 
     PshotP] + (1 - Exp[-k^2 Sigma2Total[mu]]) PloopNW[k, mu, b1, b2, 
     bs2, b3nl, alpha0, alpha2, alpha4, ctilde, alphashot0, 
     alphashot2, PshotP];
