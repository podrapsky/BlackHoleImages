(* ::Package:: *)

(* ::Title:: *)
(*KerrNullGeodesics Package*)


(* ::Section::Closed:: *)
(*Define Usage for Public Functions*)


BeginPackage["KerrNullGeodesics`"]

KerrNullGeoDistant::usage = "KerrNullGeoDistant[a,\[Theta]o,\[Alpha],\[Beta]] returns a KerrNullGeoFunction which stores information about the trajectory. a, \[Alpha], \[Beta] given in units of the BH mass";
KerrNullGeoDistantFunction::usage = "KerrNullGeoDistantFunction[a,\[Theta]o,\[Alpha],\[Beta],assoc] an object for storing the trajectory and its parameters in the assoc Association.";

Begin["`Private`"];


KerrNullGeoDistant::OutOfBounds = "Out of bounds error: `1`"


(* ::Section::Closed:: *)
(*Constants of Motion*)


DistantNullConstantsOfMotion[a_, \[Theta]o_, \[Alpha]_, \[Beta]_] := <|
	"\[ScriptL]" ->  -\[Alpha] Sin[\[Theta]o],
	"\[Eta]" -> \[Beta]^2+(\[Alpha]^2-a^2) (Cos[\[Theta]o])^2
|>


(* ::Section::Closed:: *)
(*Polar Motion*)


(* ::Text:: *)
(*Taken from Gralla & Lupsasca, arXiv:1910.12881v3*)


Options[OrdinaryPolarMotion] = {"ReturnValues" -> "All"}
Options[VorticalPolarMotion] = {"ReturnValues" -> "All"}


OrdinaryPolarMotion[a_, \[Eta]_, \[ScriptL]_, \[Theta]o_, \[Nu]\[Theta]_, \[Lambda]x_, OptionsPattern[]] := Module[{\[CapitalDelta]\[Theta], u1, u2, EPrime, G\[Theta]o, Gto, G\[Phi]o, \[CapitalPsi], \[Theta], Gt, G\[Phi], return}, 
\[CapitalDelta]\[Theta]= 1/2 (1-(\[Eta]+\[ScriptL]^2)/a^2);
u1=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; u2=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
G\[Theta]o=-1/Sqrt[-u1 a^2] EllipticF[Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1];
\[CapitalPsi][\[Lambda]_] := JacobiAmplitude[Sqrt[-u1 a^2] (\[Lambda] + \[Nu]\[Theta] G\[Theta]o), u2/u1];

G\[Phi]o=-1/Sqrt[-u1 a^2] EllipticPi[u2, Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1];
G\[Phi]= Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[1/Sqrt[-u1 a^2] EllipticPi[u2, \[CapitalPsi][Global`\[Lambda]], u2/u1]-\[Nu]\[Theta] G\[Phi]o]], Listable]; 

\[Theta] = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[ArcCos[-\[Nu]\[Theta] Sqrt[u2] Sin[\[CapitalPsi][Global`\[Lambda]]]]]], Listable];

If[OptionValue["ReturnValues"]!="OmitT", 

  EPrime[\[Phi]_, k_] := (EllipticE[\[Phi], k]-EllipticF[\[Phi], k])/(2 k);
  Gto=2 u2/Sqrt[-u1 a^2] EPrime[Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1];
  Gt=Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[-2 u2/Sqrt[-u1 a^2] EPrime[\[CapitalPsi][Global`\[Lambda]], u2/u1] - \[Nu]\[Theta] Gto]], Listable];
  return = <|"\[Theta]" -> \[Theta], "Gt" -> Gt, "G\[Phi]" -> G\[Phi]|>,

  return = <|"\[Theta]" -> \[Theta], "G\[Phi]" -> G\[Phi]|>
];
  
return
]


VorticalPolarMotion[a_, \[Eta]_, \[ScriptL]_, \[Theta]o_, \[Nu]\[Theta]_, \[Lambda]x_, OptionsPattern[]] := Module[{h, \[CapitalDelta]\[Theta], u1, u2, \[CapitalUpsilon], \[CapitalUpsilon]\[Lambda], G\[Theta]o, Gto, G\[Phi]o,\[Theta], Gt, G\[Phi], return},
h = Sign[Cos[\[Theta]o]];
\[CapitalDelta]\[Theta]= 1/2 (1-(\[Eta]+\[ScriptL]^2)/a^2);
u1=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; u2=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2];
\[CapitalUpsilon][\[CapitalTheta]_] := ArcSin[Sqrt[((Cos[\[CapitalTheta]])^2-u1)/(u2-u1)]];
G\[Theta]o = -(h/Sqrt[u1 a^2]) EllipticF[\[CapitalUpsilon][\[Theta]o], 1-u2/u1];
\[CapitalUpsilon]\[Lambda][\[Lambda]_] := JacobiAmplitude[Sqrt[u1 a^2] (\[Lambda] + \[Nu]\[Theta] G\[Theta]o), 1-u2/u1];

G\[Phi]o = -(h/((1-u1) Sqrt[u1 a^2])) EllipticPi[(u2-u1)/(1-u1), \[CapitalUpsilon][\[Theta]o], 1-u2/u1];
G\[Phi] = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[1/((1-u1) Sqrt[u1 a^2]) EllipticPi[(u2-u1)/(1-u1), \[CapitalUpsilon]\[Lambda][Global`\[Lambda]], 1-u2/u1]-\[Nu]\[Theta] G\[Phi]o]], Listable];

\[Theta] = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[ArcCos[h Sqrt[u1+(u2-u1) (Sin[\[CapitalUpsilon]\[Lambda][Global`\[Lambda]]])^2]]]], Listable];

If[OptionValue["ReturnValues"]!="OmitT",
  
  Gto = -h Sqrt[u1/a^2] EllipticE[\[CapitalUpsilon][\[Theta]o], 1-u2/u1];
  Gt = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Sqrt[u1/a^2] EllipticE[\[CapitalUpsilon]\[Lambda][Global`\[Lambda]], 1-u2/u1]-\[Nu]\[Theta] Gto]], Listable];
  return = <|"\[Theta]" -> \[Theta], "Gt" -> Gt, "G\[Phi]" -> G\[Phi]|>,

  return = <|"\[Theta]" -> \[Theta], "G\[Phi]" -> G\[Phi]|>
];

return
]


(* ::Section::Closed:: *)
(*Radial Roots*)


(* ::Text:: *)
(*Taken from Gralla & Lupsasca, arXiv:1910.12881v3*)


RadialRoots[a_, \[Eta]_, \[ScriptL]_] := Module[{A, B, C, P, Q, \[CapitalOmega]1, \[CapitalOmega]2, \[Omega]1, \[Omega]2, z},
A=a^2-\[Eta]-\[ScriptL]^2; B=2(\[Eta]+(\[ScriptL]-a)^2); C=-a^2 \[Eta];
P= -A^2/12 - C; Q=-A/3 ((A/6)^2-C) - B^2/8;
\[CapitalOmega]1 = -Q/2 - Sqrt[(P/3)^3+(Q/2)^2];
\[CapitalOmega]2 = -Q/2 + Sqrt[(P/3)^3+(Q/2)^2];
\[Omega]1=If[Element[\[CapitalOmega]1, Reals], Surd[\[CapitalOmega]1, 3], Power[\[CapitalOmega]1, 1/3]];
\[Omega]2=If[Element[\[CapitalOmega]2, Reals], Surd[\[CapitalOmega]2, 3] , Power[\[CapitalOmega]2, 1/3]];
z=Sqrt[(\[Omega]2+\[Omega]1-A/3)/2];

<|
  "r1" -> -z-Sqrt[-A/2-z^2+B/(4 z)],
  "r2" -> -z+Sqrt[-A/2-z^2+B/(4 z)],
  "r3" -> +z-Sqrt[-A/2-z^2-B/(4 z)],
  "r4" -> +z+Sqrt[-A/2-z^2-B/(4 z)]
|>
]


(* ::Section::Closed:: *)
(*Radial Motion*)


(* ::Text:: *)
(*Taken from Gralla & Lupsasca, arXiv:1910.12881v3*)


(* ::Text:: *)
(*For a distant observer, Case 1 is irrelevant because r is bounded there.*)
(*We introduce the reduced integral RedIt defined as the limit of It[ro] - ro as ro -> \[Infinity].*)


DistantRadialMotionCase2[roots_, a_, \[Eta]_, \[ScriptL]_] := Module[{r1, r2, r3, r4, rp, rm, xo, k, I0o, X, r, Pip, Pim, RedPi1, E\[Lambda], Ip, Iminus, RedI1, deriv, RedI2, RedIt, I\[Phi], \[Lambda]x},
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;
xo=Sqrt[(r3-r1)/(r4-r1)];  k=((r3-r2) (r4-r1))/((r3-r1) (r4-r2)); rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];
I0o = 2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[xo], k];
If[r4<rm, 
  \[Lambda]x = I0o - 2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[Sqrt[((rp-r4) (r3-r1))/((rp-r3) (r4-r1))]], k],
  \[Lambda]x = 2 I0o;
];
X[\[Lambda]_] := Sqrt[(r3-r1) (r4-r2)]/2 (\[Lambda]-I0o);

r = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(r4 (r3-r1) - r3 (r4-r1) (JacobiSN[X[Global`\[Lambda]], k])^2)/((r3-r1) - (r4-r1) (JacobiSN[X[Global`\[Lambda]], k])^2)]]], Listable];

RedPi1[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)](EllipticPi[(r4-r1)/(r3-r1), JacobiAmplitude[X[\[Lambda]], k], k] + EllipticF[ArcSin[Sqrt[(r3-r1)/(r4-r1)]], k] - EllipticPi[(r3-r2)/(r4-r2), ArcSin[Sqrt[(r3-r1)/(r4-r1)]], k] + 1/(2 Sqrt[(1-(r3-r2)/(r4-r2))((r4-r1)/(r3-r1)-1)]) Log[4/(r3-r1+r4-r2)]);
E\[Lambda][\[Lambda]_] := Sqrt[(r3-r1)(r4-r2)](EllipticE[JacobiAmplitude[X[\[Lambda]], k], k] + EllipticE[ArcSin[xo], k]);
Pip[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)] (r4-r3)/((rp-r3) (rp-r4)) (EllipticPi[((rp-r3) (r4-r1))/((rp-r4) (r3-r1)), JacobiAmplitude[X[\[Lambda]], k], k] + EllipticPi[((rp-r3) (r4-r1))/((rp-r4) (r3-r1)), ArcSin[xo], k]);
Pim[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)] (r4-r3)/((rm-r3) (rm-r4)) (EllipticPi[((rm-r3) (r4-r1))/((rm-r4) (r3-r1)), JacobiAmplitude[X[\[Lambda]], k], k] + EllipticPi[((rm-r3) (r4-r1))/((rm-r4) (r3-r1)), ArcSin[xo], k]);
Ip[\[Lambda]_] := -\[Lambda]/(rp-r3)-Pip[\[Lambda]]; Iminus[\[Lambda]_] := -\[Lambda]/(rm-r3)-Pim[\[Lambda]];
RedI1[\[Lambda]_] := r3 \[Lambda] + (r4-r3) RedPi1[\[Lambda]];
deriv[expr_, var_] := expr'[var] /. Re'[e_]:> 1;
RedI2[\[Lambda]_] := deriv[r, \[Lambda]]/(r[\[Lambda]]-r3) + r3 - (r1 r4 + r2 r3)/2 \[Lambda] - E\[Lambda][\[Lambda]];

I\[Phi] = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[ScriptL]/2) Ip[Global`\[Lambda]] - (rm- a \[ScriptL]/2) Iminus[Global`\[Lambda]])]]], Listable];
RedIt = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 RedI1[Global`\[Lambda]] + RedI2[Global`\[Lambda]] + 2 Log[2]]]], Listable];

<|"r" -> r, "I\[Phi]" -> I\[Phi], "RedIt" -> RedIt, "\[Lambda]x" -> Re[\[Lambda]x]|>
]


DistantRadialMotionCase3[roots_, a_, \[Eta]_, \[ScriptL]_] := Module[{r1, r2, r3, r4, rp, rm, A, B, xo, k, p1, f1, RedR1, RedR2, R1, R2, I0o, \[Lambda]x, X, r, \[Alpha]0, \[Alpha]p, \[Alpha]m, Pip, Pim, RedPi1, RedPi2, Ip, Iminus, RedI1, RedI2, RedIt, I\[Phi]},
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;
A = Abs[Sqrt[(r3-r2) (r4-r2)]]; B = Abs[Sqrt[(r3-r1) (r4-r1)]];
xo = (A-B)/(A+B); k = ((A+B)^2-(r2-r1)^2)/(4 A B); rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];
p1[\[Gamma]_, j_] := Sqrt[(\[Gamma]^2-1)/(j+(1-j) \[Gamma]^2)];
f1[\[Gamma]_, \[Xi]_, j_] := p1[\[Gamma], j]/2 Log[Abs[(p1[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2] + Sin[\[Xi]])/(p1[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2] - Sin[\[Xi]])]];
RedR1[\[Gamma]_, \[Xi]_, j_] := 1/(1-\[Gamma]^2) (EllipticF[\[Xi], j] - EllipticPi[(j(\[Gamma]^2-1))/\[Gamma]^2, \[Xi], j] - (\[Gamma] Sqrt[A B])/(r2-r1) Log[(4(r2-r1))/(B^2-A^2)] + (\[Gamma] Sqrt[A B])/(r2-r1) Log[(B^2-A^2)/(4(r2-r1))+(A B (r2-r1))/(B^2-A^2)]); 
RedR2[\[Gamma]_, \[Xi]_, j_] := 1/(\[Gamma]^2-1) (EllipticF[\[Xi], j] -\[Gamma]^2/(j + (1-j)\[Gamma]^2) EllipticE[\[Xi], j]);
R1[\[Gamma]_, \[Xi]_, j_] := 1/(1-\[Gamma]^2) (EllipticPi[\[Gamma]^2/(\[Gamma]^2-1), \[Xi], j]-\[Gamma] f1[\[Gamma], \[Xi], j]);
R2[\[Gamma]_, \[Xi]_, j_] := 1/(\[Gamma]^2-1) (EllipticF[\[Xi], j] -\[Gamma]^2/(j + (1-j)\[Gamma]^2) (EllipticE[\[Xi], j] - (\[Gamma] Sin[\[Xi]]Sqrt[1-j (Sin[\[Xi]])^2])/(1+\[Gamma] Cos[\[Xi]]))) + 1/(j + (1-j) \[Gamma]^2) (2j - \[Gamma]^2/(\[Gamma]^2-1))R1[\[Gamma], \[Xi], j];
I0o = 1/Sqrt[A B] EllipticF[ArcCos[xo], k];
\[Lambda]x = I0o - 1/Sqrt[A B] EllipticF[ArcCos[(A (rp-r1) - B (rp-r2))/(A (rp-r1) + B (rp-r2))], k];
X[\[Lambda]_] := Sqrt[A B] (\[Lambda] - I0o);

r = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[((B r2 - A r1) + (B r2 + A r1) JacobiCN[X[Global`\[Lambda]], k])/((B-A) + (B+A) JacobiCN[X[Global`\[Lambda]], k])]]], Listable];

\[Alpha]0 = (B+A)/(B-A);
\[Alpha]p = (B (rp - r2) + A (rp - r1))/(B (rp - r2) - A (rp - r1)); \[Alpha]m = (B (rm - r2) + A (rm - r1))/(B (rm - r2) - A (rm - r1));
Pip[\[Lambda]_] := ((2 (r2-r1) Sqrt[A B])/(B (rp-r2) - A (rp-r1)))(R1[\[Alpha]p, JacobiAmplitude[X[\[Lambda]], k], k]+R1[\[Alpha]p, ArcCos[xo], k]); 
Pim[\[Lambda]_] := ((2 (r2-r1) Sqrt[A B])/(B (rm-r2) - A (rm-r1)))(R1[\[Alpha]m, JacobiAmplitude[X[\[Lambda]], k], k]+R1[\[Alpha]m, ArcCos[xo], k]);
RedPi1[\[Lambda]_] := (2(r2-r1)Sqrt[A B])/(B^2-A^2) (R1[\[Alpha]0, JacobiAmplitude[X[\[Lambda]], k], k] + RedR1[\[Alpha]0, ArcCos[xo], k]);
RedPi2[\[Lambda]_] := ((2(r2-r1)Sqrt[A B])/(B^2-A^2))^2 (R2[\[Alpha]0, JacobiAmplitude[X[\[Lambda]], k], k] + RedR2[\[Alpha]0, ArcCos[xo], k]);
Ip[\[Lambda]_] := -((B+A) \[Lambda] + Pip[\[Lambda]])/(B (rp-r2) + A (rp-r1)); Iminus[\[Lambda]_] := -((B+A) \[Lambda] + Pim[\[Lambda]])/(B (rm-r2) + A (rm-r1));
RedI1[\[Lambda]_] := ((B r2 + A r1)/(B + A))\[Lambda] + RedPi1[\[Lambda]];
RedI2[\[Lambda]_] := ((B r2 + A r1)/(B + A))^2 \[Lambda] + 2((B r2 + A r1)/(B + A)) (2(r2-r1)Sqrt[A B])/(B^2-A^2) R1[\[Alpha]0, JacobiAmplitude[X[\[Lambda]], k], k] + Sqrt[A B] RedPi2[\[Lambda]] + (A^2-B^2)/(2 (r2-r1))-(r1+r2)+(B r2 + A r1)/(B + A);

I\[Phi] = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[ScriptL]/2) Ip[Global`\[Lambda]] - (rm- a \[ScriptL]/2) Iminus[Global`\[Lambda]])]]], Listable];
RedIt = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 RedI1[Global`\[Lambda]] + RedI2[Global`\[Lambda]] + 2 Log[2]]]], Listable];

<|"r" -> r, "I\[Phi]" -> I\[Phi], "RedIt" -> RedIt, "\[Lambda]x" -> Re[\[Lambda]x]|>
]


DistantRadialMotionCase4[roots_, a_, \[Eta]_, \[ScriptL]_] := Module[{r1, r2, r3, r4, rp, rm, C, D, k, a2, b1, g0, I0o, \[Lambda]x, x, X, p2, f2, RedS1, S1, RedS2, S2, r, Pip, Pim, RedPi1, RedPi2, Ip, Iminus, RedI1, RedI2, I\[Phi], RedIt},
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;
C = Sqrt[(r3-r1) (r4-r2)]; D = Sqrt[(r3-r2) (r4-r1)];
k=(4 C D)/(C + D)^2; a2=Sqrt[-((r2-r1)^2/4)]; b1= (r3+r4)/2; g0=Sqrt[(4 a2^2-(C-D)^2)/((C+D)^2-4 a2^2)]; 
rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];
I0o = 2/(C+D) EllipticF[\[Pi]/2+ArcTan[g0], k];
\[Lambda]x = I0o - 2/(C+D) EllipticF[ArcTan[(rp+b1)/a2]+ArcTan[g0], k];
x[\[Rho]_] = (\[Rho]+b1)/a2; X[\[Lambda]_] = (C + D)/2 (-\[Lambda]+I0o);
p2[\[Gamma]_, j_] := Sqrt[(\[Gamma]^2+1)/(1-j + \[Gamma]^2)];
f2[\[Gamma]_, \[Xi]_, j_] := p2[\[Gamma], j]/2 Log[Abs[(1-p2[\[Gamma], j])/(1+p2[\[Gamma], j]) (1+p2[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2] )/(1-p2[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2])]];
RedS1[\[Gamma]_, \[Xi]_, j_] := EllipticF[\[Xi], j] - 1/(1+\[Gamma]^2) (\[Gamma]^2 EllipticPi[j/(1+\[Gamma]^2), \[Xi], j] + (\[Gamma] (C+D))/(4 a2) Log[(64 a2^2)/((2 a2 + C + D)^2 (\[Gamma]^2 (C+D)^2 + 4 a2^2))]);
S1[\[Gamma]_, \[Xi]_, j_] := 1/(1+\[Gamma]^2) (EllipticF[\[Xi],j] +\[Gamma]^2 EllipticPi[1+\[Gamma]^2, \[Xi], j]-\[Gamma] f2[\[Gamma], \[Xi], j]);
RedS2[\[Gamma]_, \[Xi]_, j_] := -1/((1+\[Gamma]^2)(1-j+\[Gamma]^2)) ((1-j)EllipticF[\[Xi], j] + \[Gamma]^2 EllipticE[\[Xi], j] - \[Gamma]^3);
S2[\[Gamma]_, \[Xi]_, j_] := -1/((1+\[Gamma]^2)(1-j+\[Gamma]^2)) ((1-j)EllipticF[\[Xi], j] + \[Gamma]^2 EllipticE[\[Xi], j] + (\[Gamma]^2 Sqrt[1-j (Sin[\[Xi]])^2](\[Gamma]-Tan[\[Xi]]))/(1+\[Gamma] Tan[\[Xi]]) - \[Gamma]^3) + (1/(1+\[Gamma]^2)+(1-j)/(1-j+\[Gamma]^2))S1[\[Gamma], \[Xi], j];

r = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[-a2 ((g0-JacobiSC[X[Global`\[Lambda]], k])/(1+g0 JacobiSC[X[Global`\[Lambda]], k])) - b1]]], Listable];

Pip[\[Lambda]_] := -(2/(C+D)) ((1+g0^2)/(g0 (g0+x[rp]))) (S1[(g0 x[rp] - 1)/(g0 + x[rp]), JacobiAmplitude[X[\[Lambda]], k], k] - S1[(g0 x[rp] - 1)/(g0 + x[rp]), \[Pi]/2+ArcTan[g0], k]);
Pim[\[Lambda]_] := -(2/(C+D)) ((1+g0^2)/(g0 (g0+x[rm]))) (S1[(g0 x[rm] - 1)/(g0 + x[rm]), JacobiAmplitude[X[\[Lambda]], k], k] - S1[(g0 x[rm] - 1)/(g0 + x[rm]), \[Pi]/2+ArcTan[g0], k]);
RedPi1[\[Lambda]_] := -2/(C+D) (a2/g0 (1+g0^2))(S1[g0, JacobiAmplitude[X[\[Lambda]], k], k] - RedS1[g0, \[Pi]/2+ArcTan[g0], k]);
RedPi2[\[Lambda]_] := -2/(C+D) (a2/g0 (1+g0^2))^2 (S2[g0, JacobiAmplitude[X[\[Lambda]], k], k] - RedS2[g0, \[Pi]/2+ArcTan[g0], k]);
Ip[\[Lambda]_] := g0/(a2 (1-g0 x[rp])) (\[Lambda]-Pip[\[Lambda]]);
Iminus[\[Lambda]_] := g0/(a2 (1-g0 x[rm])) (\[Lambda]-Pim[\[Lambda]]);
RedI1[\[Lambda]_] := (a2/g0-b1)\[Lambda] - RedPi1[\[Lambda]];
RedI2[\[Lambda]_] := (a2/g0-b1)^2 \[Lambda] + 4 (a2/g0-b1) 1/(C+D) a2/g0 (1+g0^2) S1[g0, JacobiAmplitude[X[\[Lambda]], k], k] + RedPi2[\[Lambda]]+b1-(2 g0 C D)/(Sqrt[g0^2+1] Sqrt[(C-D)^2+g0^2 (C+D)^2]);


I\[Phi] = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[ScriptL]/2) Ip[Global`\[Lambda]] - (rm- a \[ScriptL]/2) Iminus[Global`\[Lambda]])]]], Listable];
RedIt = Function[{Global`\[Lambda]}, If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 RedI1[Global`\[Lambda]] + RedI2[Global`\[Lambda]] + 2 Log[2]]]], Listable];

<|"r" -> r, "I\[Phi]" -> I\[Phi], "RedIt" -> RedIt, "\[Lambda]x" -> Re[\[Lambda]x]|>
]


(* ::Section::Closed:: *)
(*Equator Intersections*)


EquatorIntersectionMinoTimes[a_, \[Eta]_, \[ScriptL]_, \[Theta]o_, \[Nu]\[Theta]_, \[Lambda]x_] := Module[{\[CapitalDelta]\[Theta], u1, u2, G\[Theta]o, equator\[Lambda], j0, j, t},
\[CapitalDelta]\[Theta]= 1/2 (1-(\[Eta]+\[ScriptL]^2)/a^2);
u1=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; u2=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
G\[Theta]o=-1/Sqrt[-u1 a^2] EllipticF[Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1]; (*Sometimes the ArcSin argument is slightly over 1 probably due to numerical errors*)
equator\[Lambda] = {};
If[(-\[Nu]\[Theta] G\[Theta]o) < 0, j0=1, j0=0];
K = EllipticK[u2/u1];
For[j=j0, j<=10, j++,
t=2 j K/Sqrt[-u1 a^2]-\[Nu]\[Theta] G\[Theta]o;
If[t>\[Lambda]x, Break[]];
equator\[Lambda] = Append[equator\[Lambda], t]];
equator\[Lambda]
]


(* ::Section:: *)
(*Emission Parameters*)


(* ::Text:: *)
(*Only applicable for the thin disk approximation with matter orbiting on the direct circular equatorial orbits.*)
(*The parameter \[Kappa] represents the ratio between the frequency of observed and emitted photons. \[Theta]loc, \[Phi]loc are the angles of emitted photon in emitter's local spherical coordinate system.*)


EmissionParameters[a_, \[Eta]_, \[ScriptL]_, \[Theta]o_, rem_, \[Theta]em_] := Module[{Z1, Z2, rms, A, B, \[CapitalDelta], \[Omega], \[CapitalOmega], utcg, R, \[CapitalTheta], ploct, plocr, ploc\[Theta], \[Kappa], \[Theta]loc, \[Phi]loc},
Z1=1+(1-a^2)^(1/3) ((1+a)^(1/3) + (1-a)^(1/3)); Z2=(3 a^2 + Z1^2)^(1/2); rms=3+Z2-((3-Z1) (3+Z1+2 Z2))^(1/2);
If[rem<rms, <|"\[Kappa]" -> -1, "\[Theta]loc" -> -1, "\[Phi]loc" -> -1|>,
A = rem^2; \[CapitalDelta] = rem^2-2 rem + a^2; B = (rem^2 + a^2)^2 - \[CapitalDelta] a^2; \[Omega] = (2 a rem)/B;
\[CapitalOmega] =1/(Sqrt[rem^3]+ a); utcg=((\[CapitalDelta] A)/B - (\[Omega]-\[CapitalOmega])^2 B/A)^(-1/2);
R = (rem^2 + a^2 - a \[ScriptL])^2 - \[CapitalDelta] (\[Eta] + (\[ScriptL]-a)^2);
\[CapitalTheta] = \[Eta] + a^2 (Cos[\[Theta]em])^2 - \[ScriptL]^2/(Tan[\[Theta]em])^2;

ploct = -utcg (1- \[CapitalOmega] \[ScriptL]);
plocr = -Sqrt[(R/(\[CapitalDelta] A))];
ploc\[Theta] = Sign[\[Pi]/2 - \[Theta]o] Sqrt[\[CapitalTheta]/A];

\[Kappa] = -1/ploct;
\[Theta]loc = ArcCos[ploc\[Theta]/-ploct];
\[Phi]loc = ArcCos[plocr/Sqrt[ploct^2-ploc\[Theta]^2]];

<|"\[Kappa]" -> \[Kappa], "\[Theta]loc" -> \[Theta]loc, "\[Phi]loc" -> \[Phi]loc|>
 ]]


(* ::Section::Closed:: *)
(*Public Functions*)


Options[KerrNullGeoDistant] = {"Rotation" -> "Counterclockwise", "PhiRange" -> {-\[Infinity], \[Infinity]}}
SyntaxInformation[KerrNullGeoDistant] = {"ArgumentsPattern"->{_,_,_,_}};

KerrNullGeoDistant[a_, \[Theta]o_, \[Alpha]_, \[Beta]_, OptionsPattern[]] := Module[ {consts, \[Eta], \[ScriptL], roots, r1, r2, r3, r4, rp, rm, k, r, \[Theta], \[Phi], G\[Phi], I\[Phi], \[CapitalDelta]v, Gt, RedIt, \[Lambda]x, \[Phi]x, \[Theta]x, assoc, equator\[Lambda], type, \[Kappa], \[Theta]loc, \[Phi]loc},

If[a<=0 || a>=1, Message[KerrNullGeoDistantDistant::OutOfBounds, "Parameter a must be between 0 and 1."]; Return[];];
If[\[Theta]o<0 || \[Theta]o>\[Pi], Message[KerrNullGeoDistantDistant::OutOfBounds, "Parameter \[Theta]o must be between 0 and \[Pi]."]; Return[];];

If[OptionValue["Rotation"]=="Clockwise", \[Alpha]=-\[Alpha]];

consts = DistantNullConstantsOfMotion[a, \[Theta]o, \[Alpha], \[Beta]];
{\[ScriptL], \[Eta]} = {"\[ScriptL]", "\[Eta]"} /. consts;
If[Abs[\[Eta]]< 10 $MachineEpsilon || Abs[\[Eta]+(Abs[\[ScriptL]]-a)^2]< 10 $MachineEpsilon, \[Eta]=\[Eta] + 10 $MachineEpsilon]; (*To avoid undefined expresions in polar motion*)

roots = RadialRoots[a, \[Eta], \[ScriptL]];
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;

type = "PhotonCapture";
If[Im[r2] != 0,
  {r, I\[Phi], \[Lambda]x, RedIt} = {"r", "I\[Phi]", "\[Lambda]x", "RedIt"} /. DistantRadialMotionCase4[roots, a, \[Eta], \[ScriptL]],
  If[Im[r4] != 0, 
    {r, I\[Phi], \[Lambda]x, RedIt} = {"r", "I\[Phi]", "\[Lambda]x", "RedIt"} /. DistantRadialMotionCase3[roots, a, \[Eta], \[ScriptL]],
    {r, I\[Phi], \[Lambda]x, RedIt} = {"r", "I\[Phi]", "\[Lambda]x", "RedIt"} /. DistantRadialMotionCase2[roots, a, \[Eta], \[ScriptL]];
    If[r3>1-Sqrt[1-a^2], type = "PhotonEscape"]
  ]
];

ssign[x_] := If[x<0, -1, 1];

If[\[Eta]>0, 
  {\[Theta], G\[Phi], Gt} = {"\[Theta]", "G\[Phi]", "Gt"} /. OrdinaryPolarMotion[a, \[Eta], \[ScriptL], \[Theta]o, -ssign[\[Beta]], \[Lambda]x],
  {\[Theta], G\[Phi], Gt} = {"\[Theta]", "G\[Phi]", "Gt"} /. VorticalPolarMotion[a, \[Eta], \[ScriptL], \[Theta]o, -ssign[\[Beta]], \[Lambda]x]; (*Minus Sign[\[Beta]] because we use negative Mino time*)
];


If[\[Eta]>0, 
equator\[Lambda] = EquatorIntersectionMinoTimes[a, \[Eta], \[ScriptL], \[Theta]o, -Sign[\[Beta] + $MachineEpsilon], \[Lambda]x],
equator\[Lambda] = {}
];

If[OptionValue["PhiRange"][[2]]==\[Infinity], 
  \[Phi]=Function[{Global`\[Lambda]}, Evaluate[I\[Phi][Global`\[Lambda]]+\[ScriptL] G\[Phi][Global`\[Lambda]]], Listable],
  \[Phi]=Function[{Global`\[Lambda]}, Evaluate[Mod[I\[Phi][Global`\[Lambda]]+\[ScriptL] G\[Phi][Global`\[Lambda]], OptionValue["PhiRange"][[2]]-OptionValue["PhiRange"][[1]], OptionValue["PhiRange"][[1]]]], Listable]
];
\[CapitalDelta]v=Function[{Global`\[Lambda]}, Evaluate[RedIt[Global`\[Lambda]]+a^2 Gt[Global`\[Lambda]] + r[Global`\[Lambda]] + 2 Log[r[Global`\[Lambda]]/2]], Listable];

If[type == "PhotonEscape", \[Theta]x=\[Theta][\[Lambda]x]; \[Phi]x=\[Phi][\[Lambda]x], \[Theta]x=-1; \[Phi]x=-1];

If[Length[equator\[Lambda]] == 0, 
  {\[Kappa], \[Theta]loc, \[Phi]loc} = {-1, -1, -1}, 
  {\[Kappa], \[Theta]loc, \[Phi]loc} = {"\[Kappa]", "\[Theta]loc", "\[Phi]loc"} /. EmissionParameters[a, \[Eta], \[ScriptL], \[Theta]o, r[equator\[Lambda][[1]]], \[Theta][equator\[Lambda][[1]]]]
];

assoc = <|
	"Trajectory" -> {\[CapitalDelta]v, r, \[Theta], \[Phi]},
	"ConstantsOfMotion" -> consts,
	"RadialRoots" -> {r1, r2, r3, r4},
	"EquatorIntersectionMinoTimes" -> equator\[Lambda],
	"TrajectoryType" -> type,
	"MinoTimeOfCapture" -> \[Lambda]x,
	"EscapeCoordinates" -> {\[Theta]x, \[Phi]x},
	"EmissionParameters" -> {\[Kappa], \[Theta]loc, \[Phi]loc}
|>;

KerrNullGeoDistantFunction[a, \[Theta]o, \[Alpha], \[Beta], assoc]
]


KerrNullGeoDistantFunction[a_, \[Theta]o_, \[Alpha]_, \[Beta]_, assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]];
KerrNullGeoDistantFunction[a_, \[Theta]o_, \[Alpha]_, \[Beta]_, assoc_][y_?StringQ] := assoc[y];
Keys[g_KerrNullGeoDistantFunction]^:=Keys[g[[5]]];


(* ::Section::Closed:: *)
(*Close the Package*)


End[];

EndPackage[];
