(* ::Package:: *)

(*PT Models*)
(* 1. PTkernels=1  (FullKernels) *)
(* 2. PTkernels=2  (MG no Screenings) *)
(* 3. PTkernels=3  (LCDM kernels) *)
(* 4. PTkernels=4  (LCDM kernels but keeping f(k) *)
(* 5. PTkernels=5  (EdS kernels) *)

outputbasedir="outputs/"<>presuffix;

If[! DirectoryQ[outputbasedir],
 CreateDirectory[outputbasedir];
 ];



If[PTkernels==1,outputdir=outputbasedir<>"/FullMG"];
If[PTkernels==2,outputdir=outputbasedir<>"/NoScreen"];
If[PTkernels==3,outputdir=outputbasedir<>"/LCDM"];
If[PTkernels==4,outputdir=outputbasedir<>"/LCDMfk"];
If[PTkernels==5,outputdir=outputbasedir<>"/EdS"];
If[PTkernels==6,outputdir=outputbasedir<>"/EdSfk"];

If[! DirectoryQ[outputdir],CreateDirectory[outputdir]];
chat="Files will be stored in the directory ./"<>outputdir;
Print[chat];


 


outlogfile=outputdir<>"/outlog"<>suffix<>".txt";
(*If[! FileExistsQ[outlogfile],CopyFile["void.txt",outlogfile]];*)
If[! FileExistsQ[outlogfile],Put[outlogfile]];


fopen = OpenWrite[outlogfile];

WriteLine[fopen,DateString[]];
WriteLine[fopen,"..."];

If[PTkernels==1,WriteLine[fopen,"PTkernels=1: Using full MG kernels"]];
If[PTkernels==2,WriteLine[fopen,"PTkernels=2: Using MG kernels without screenings"]];
If[PTkernels==3,WriteLine[fopen,"PTkernels=3: Using LCDM kernels"]];
If[PTkernels==4,WriteLine[fopen,"PTkernels=4: Using LCDM kernels but keeping f(k)"]];
If[PTkernels==5,WriteLine[fopen,"PTkernels=5: Using EdS kernels"]];
If[PTkernels==6,WriteLine[fopen,"PTkernels=6: Using EdS kernels but keeping f(k)"]];


WriteLine[fopen,"..."];
WriteLine[fopen,"name: "<> suffix];
WriteLine[fopen,"fR0 = -"<> ToString[fR0]];
WriteLine[fopen,"OmegaM0 = "<>ToString[om]];
WriteLine[fopen,"h = "<>ToString[h]];
WriteLine[fopen,"output redshift = "<>ToString[zev]];
WriteLine[fopen,"..."];
WriteLine[fopen,"Input pk: "<>inputpk];
WriteLine[fopen,"Input pk is LCDM: "<> ToString[InputpkIsLCDM]];
WriteLine[fopen,"Input pk redshift:  zin="<>ToString[zinputpk]];

Close[fopen];
