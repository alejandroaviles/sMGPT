(* ::Package:: *)

(*PT Models*)
(* 1. PTkernels=1  (FullKernels) *)
(* 2. PTkernels=2  (MG no Screenings) *)
(* 3. PTkernels=3  (LCDM kernels) *)
(* 4. PTkernels=4  (LCDM kernels but keeping f(k) *)
(* 5. PTkernels=5  (EdS kernels) *)
(* 6. PTkernels=6  (EdS kernels but keeping f(k) *)

PTkernels=1;





fR0=10^(-5.);
presuffix="Fz05_NW" (*suffix for output files*)

zev=0.5; (*output redshift*)
h=0.697;om=0.281;(*BLi simulations*)
(*h=0.6774;om=0.3089;*)(*lightcone simulations*)


inputpk="../input/pkl_nw_WMAP9y.dat" 
InputpkIsLCDM = True;
zinputpk=0;


If[PTkernels==1,suffix=presuffix<>"_FullMG"];
If[PTkernels==2,suffix=presuffix<>"_NoScreen"];
If[PTkernels==3,suffix=presuffix<>"_LCDM"];
If[PTkernels==4,suffix=presuffix<>"_LCDMfk"];
If[PTkernels==5,suffix=presuffix<>"_EdS"];
If[PTkernels==6,suffix=presuffix<>"_EdSfk"];

If[step==0 && PTkernels==1,Print["PTkernels=1: Using full MG kernels"]];
If[step==0 && PTkernels==2,Print["PTkernels=2: Using MG kernels without screenings"]];
If[step==0 && PTkernels==3,Print["PTkernels=3: Using LCDM kernels"]];
If[step==0 && PTkernels==4,Print["PTkernels=4: Using LCDM kernels but keeping f(k)"]];
If[step==0 && PTkernels==5,Print["PTkernels=5: Using EdS kernels"]];
If[step==0 && PTkernels==6,Print["PTkernels=6: Using EdS kernels but keeping f(k)"]];




(*If[bLCDM==1,fR0=10^-50.;Print["Using LCDM kernels"]];*)

etaev=Log[1/(1.+zev)];
If[step==0,Print["suffix = ", suffix]];
If[step==0,Print["fR0 = -",fR0]];
If[step==0,Print["OmegaM0 = ",om,", h = ", h]];
If[step==0,Print["redshift z = ",zev," , eta = ", etaev]];

(*with screening Sc= 1, no screening Sc= 0*)
If[PTkernels==1,Sc=1;If[step==0,Print["Sc=1: Screenings on"]]];
If[PTkernels==2||PTkernels==3||PTkernels==4||PTkernels==5||PTkernels==6,Sc=0;If[step==0,Print["Sc=0: Screenings off"]]];

