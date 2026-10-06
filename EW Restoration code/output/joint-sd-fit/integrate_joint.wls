(* Joint linear rate response. Columns: SM, coefficient of x=n7/1e-4, coefficient of y=C_Hq3/(0.01 TeV^-2).
   pT bins 100-200 and 200-1500 GeV. Common upper mass cut in both channels.
   Interference order only: no ad hoc addition of separate nonlinear scans. *)
base=Directory[];dir=base<>"/output/n7-masscut";
qOrder=ToExpression[Environment["SCAN_ORDER"]];If[!IntegerQ[qOrder],qOrder=48];
gauss[n_]:=Module[{v,w,b},b=Table[k/Sqrt[4.k^2-1],{k,n-1}];{v,w}=Eigensystem[DiagonalMatrix[b,1]+DiagonalMatrix[b,-1]];SortBy[Transpose[{v,2 w[[All,1]]^2}],First]];
gA=gauss[qOrder];gM=gauss[If[qOrder>=128,24,16]];
raw=Association[];
Do[
 an=Get[dir<>"/"<>channel<>"_amplitudes.wl"][[1,2]];
 ac=Get[base<>"/output/chq3-masscut/"<>channel<>"_amplitudes.wl"][[1,2]];
 smn=an/.nn->0;smc=ac/.cc->0;
 smcheck=Max[Flatten[Table[Abs[(smn/smc-1)/.{E0->e,\[Theta]->t}],{e,{250.,500.,1500.}},{t,{0.4,1.2,2.3}}]]];
 Print["SM amplitude relative agreement ",channel," ",smcheck];
 If[smcheck>10^-7,Print["SM amplitude mismatch"];Quit[1]];
 coeff={smn^2,2 smn (D[an,nn]/.nn->0) 10^-4,2 smc (D[ac,cc]/.cc->0) .01};
 Do[
 expr=coeff/.{E0->mass/2,\[Theta]->ArcCos[z]};
 fun=With[{ex=N[expr]},Compile[{{mass,_Real},{z,_Real}},ex,RuntimeOptions->"Speed"]];
 Print["compiled ",channel," ",pol," test ",fun[1000.,0.3]];
 m1=80.379;m2=If[channel=="WZ",91.1876,125.];shad=14000.^2;
 pbar[m_]:=Sqrt[(m^2-(m1+m2)^2)(m^2-(m1-m2)^2)]/(2m);
 threshold[p_]:=Sqrt[m1^2+p^2]+Sqrt[m2^2+p^2];
 Do[
 ld=Import[dir<>"/luminosity_"<>ToString[iset]<>".dat","Table"];
 If[!MatrixQ[ld,NumericQ],Print["bad luminosity"];Quit[1]];
 lum=Interpolation[Log[ld],InterpolationOrder->3];
 Do[
 yields=Table[
 low=bin[[1]];high=bin[[2]];
 If[mcut<=threshold[low],ConstantArray[0.,3],
splits=Sort[DeleteDuplicates[Join[{Log[threshold[low]],Log[mcut]},If[threshold[high]<mcut,{Log[threshold[high]]},{}]]]];
 mesh=Flatten[Table[Most[Subdivide[splits[[j]],splits[[j+1]],Ceiling[(splits[[j+1]]-splits[[j]])/0.08]]],{j,Length[splits]-1}]];AppendTo[mesh,Last[splits]];
 Total[Table[Module[{lmid=(mesh[[j]]+mesh[[j+1]])/2,lhalf=(mesh[[j+1]]-mesh[[j]])/2},
 Total[Table[Module[{lm=lmid+lhalf gM[[k,1]],mm,pp,cl,ch,cm,cr,pref,angles},
 mm=Exp[lm];pp=pbar[mm];cl=Sqrt[Max[0.,1-(high/pp)^2]];ch=Sqrt[Max[0.,1-(low/pp)^2]];cm=(cl+ch)/2;cr=(ch-cl)/2;
 pref=2 mm^2/shad Exp[lum[lm]] pp/(256 Pi^2 (mm/2)^3*12);
 angles=Total[Table[gA[[a,2]](fun[mm,cm+cr gA[[a,1]]]+fun[mm,-cm-cr gA[[a,1]]]),{a,Length[gA]}]]cr;
 lhalf gM[[k,2]]pref angles],{k,Length[gM]}]]],{j,Length[mesh]-1}]]],
 {bin,{{100.,200.},{200.,1500.}}}];
 key=StringRiffle[{channel,pol,ToString[iset],ToString[Round[mcut]]},"_"];
 AssociateTo[raw,key->Re[yields]];Print[key," ",yields];
 Put[raw,base<>"/output/joint-sd-fit/response_integrals_order"<>ToString[qOrder]<>".wl"];
 ,{mcut,{500.,1000.,2000.,3000.}}];
 ,{iset,{1,2}}];
 ,{pol,{"long"}}];
 ,{channel,{"Wh","WZ"}}];
Export[base<>"/output/joint-sd-fit/response_integrals_order"<>ToString[qOrder]<>".json",raw,"RawJSON"];Quit[];
