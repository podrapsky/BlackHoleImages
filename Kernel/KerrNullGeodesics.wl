(* ::Package:: *)

(* ::Title:: *)
(*KerrNullGeodesics Package*)


(* ::Section:: *)
(*Define Usage for Public Functions*)


BeginPackage["KerrNullGeodesics`"]

KerrNullGeoDistant::usage = "KerrNullGeoDistant[a,\[Theta]o,\[Alpha],\[Beta]] returns a KerrNullGeoFunction which stores information about the trajectory. a, \[Alpha], \[Beta] given in units of the BH mass";
KerrNullGeoDistantFunction::usage = "KerrNullGeoDistantFunction[a,\[Theta]o,\[Alpha],\[Beta],assoc] an object for storing the trajectory and its parameters in the assoc Association.";

Begin["`Private`"];


KerrNullGeoDistant::OutOfBounds = "Out of bounds error: `1`"


(* ::Section:: *)
(*Constants of Motion*)


DistantNullConstantsOfMotion[a_, \[Theta]o_, \[Alpha]_, \[Beta]_] := <|
	"\[Lambda]" ->  -\[Alpha] Sin[\[Theta]o],
	"\[Eta]" -> \[Beta]^2+(\[Alpha]^2-a^2) (Cos[\[Theta]o])^2
|>


(* ::Section:: *)
(*Polar Motion*)


(* ::Text:: *)
(*Taken from Gralla & Lupsasca, arXiv:1910.12881v3*)


Options[OrdinaryPolarMotion] = {"ReturnValues" -> "All"}
Options[VorticalPolarMotion] = {"ReturnValues" -> "All"}


OrdinaryPolarMotion[a_, \[Eta]_, \[Lambda]_, \[Theta]o_, \[Nu]\[Theta]_, OptionsPattern[]] := Module[{\[CapitalDelta]\[Theta], u1, u2, EPrime, G\[Theta]o, Gto, G\[Phi]o, \[CapitalPsi], \[Theta], Gt, G\[Phi], return}, 
\[CapitalDelta]\[Theta]= 1/2 (1-(\[Eta]+\[Lambda]^2)/a^2);
u1=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; u2=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
G\[Theta]o=-1/Sqrt[-u1 a^2] EllipticF[Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1];
\[CapitalPsi][\[Tau]_] := JacobiAmplitude[Sqrt[-u1 a^2] (\[Tau] + \[Nu]\[Theta] G\[Theta]o), u2/u1];

G\[Phi]o=-1/Sqrt[-u1 a^2] EllipticPi[u2, Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1];
G\[Phi]= Function[{Global`\[Tau]}, Evaluate[1/Sqrt[-u1 a^2] EllipticPi[u2, \[CapitalPsi][Global`\[Tau]], u2/u1]-\[Nu]\[Theta] G\[Phi]o], Listable]; 

\[Theta] = Function[{Global`\[Tau]}, Evaluate[ArcCos[-\[Nu]\[Theta] Sqrt[u2] Sin[\[CapitalPsi][Global`\[Tau]]]]], Listable];

If[OptionValue["ReturnValues"]!="OmmitT", 

  EPrime[\[Phi]_, k_] := (EllipticE[\[Phi], k]-EllipticF[\[Phi], k])/(2 k);
  Gto=2 u2/Sqrt[-u1 a^2] EPrime[Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1];
  Gt=Function[{Global`\[Tau]}, Evaluate[-2 u2/Sqrt[-u1 a^2] EPrime[\[CapitalPsi][Global`\[Tau]], u2/u1] - \[Nu]\[Theta] Gto], Listable];
  return = <|"\[Theta]" -> \[Theta], "Gt" -> Gt, "G\[Phi]" -> G\[Phi]|>,

  return = <|"\[Theta]" -> \[Theta], "G\[Phi]" -> G\[Phi]|>
];
  
return
]


VorticalPolarMotion[a_, \[Eta]_, \[Lambda]_, \[Theta]o_, \[Nu]\[Theta]_, OptionsPattern[]] := Module[{h, \[CapitalDelta]\[Theta], u1, u2, \[CapitalUpsilon], \[CapitalUpsilon]\[Tau], G\[Theta]o, Gto, G\[Phi]o,\[Theta], Gt, G\[Phi], return},
h = Sign[Cos[\[Theta]o]];
\[CapitalDelta]\[Theta]= 1/2 (1-(\[Eta]+\[Lambda]^2)/a^2);
u1=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; u2=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2];
\[CapitalUpsilon][\[CapitalTheta]_] := ArcSin[Sqrt[((Cos[\[CapitalTheta]])^2-u1)/(u2-u1)]];
G\[Theta]o = -(h/Sqrt[u1 a^2]) EllipticF[\[CapitalUpsilon][\[Theta]o], 1-u2/u1];
\[CapitalUpsilon]\[Tau][\[Tau]_] := JacobiAmplitude[Sqrt[u1 a^2] (\[Tau] + \[Nu]\[Theta] G\[Theta]o), 1-u2/u1];

G\[Phi]o = -(h/((1-u1) Sqrt[u1 a^2])) EllipticPi[(u2-u1)/(1-u1), \[CapitalUpsilon][\[Theta]o], 1-u2/u1];
G\[Phi] = Function[{Global`\[Tau]}, Evaluate[1/((1-u1) Sqrt[u1 a^2]) EllipticPi[(u2-u1)/(1-u1), \[CapitalUpsilon]\[Tau][Global`\[Tau]], 1-u2/u1]-\[Nu]\[Theta] G\[Phi]o], Listable];

\[Theta] = Function[{Global`\[Tau]}, Evaluate[ArcCos[h Sqrt[u1+(u2-u1) (Sin[\[CapitalUpsilon]\[Tau][Global`\[Tau]]])^2]]], Listable];

If[OptionValue["ReturnValues"]!="OmmitT",
  
  Gto = -h Sqrt[u1/a^2] EllipticE[\[CapitalUpsilon][\[Theta]o], 1-u2/u1];
  Gt = Function[{Global`\[Tau]}, Evaluate[Sqrt[u1/a^2] EllipticE[\[CapitalUpsilon]\[Tau][Global`\[Tau]], 1-u2/u1]-\[Nu]\[Theta] Gto], Listable];
  return = <|"\[Theta]" -> \[Theta], "Gt" -> Gt, "G\[Phi]" -> G\[Phi]|>,

  return = <|"\[Theta]" -> \[Theta], "G\[Phi]" -> G\[Phi]|>
];

return
]


(* ::Section::Closed:: *)
(*Radial Roots*)


(* ::Text:: *)
(*Taken from Gralla & Lupsasca, arXiv:1910.12881v3*)


RadialRoots[a_, \[Eta]_, \[Lambda]_] := Module[{A, B, C, P, Q, \[CapitalOmega]1, \[CapitalOmega]2, \[Omega]1, \[Omega]2, z},
A=a^2-\[Eta]-\[Lambda]^2; B=2(\[Eta]+(\[Lambda]-a)^2); C=-a^2 \[Eta];
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


DistantRadialMotionCase2[roots_, a_, \[Eta]_, \[Lambda]_] := Module[{r1, r2, r3, r4, rp, rm, xo, k, I0o, X, r, Pip, Pim, Ip, Iminus, I\[Phi], \[Tau]x},
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;
xo=Sqrt[(r3-r1)/(r4-r1)];  k=((r3-r2) (r4-r1))/((r3-r1) (r4-r2)); rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];
I0o = 2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[xo], k];
If[r4<rm, 
  \[Tau]x = I0o - 2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[Sqrt[((rp-r4) (r3-r1))/((rp-r3) (r4-r1))]], k],
  \[Tau]x = 2 I0o;
];
X[\[Tau]_] := Sqrt[(r3-r1) (r4-r2)]/2 (\[Tau]-I0o);

r = Function[{Global`\[Tau]}, Evaluate[Re[(r4 (r3-r1) - r3 (r4-r1) (JacobiSN[X[Global`\[Tau]], k])^2)/((r3-r1) - (r4-r1) (JacobiSN[X[Global`\[Tau]], k])^2)]], Listable];

Pip[\[Tau]_] := 2/Sqrt[(r3-r1) (r4-r2)] (r4-r3)/((rp-r3) (rp-r4)) (EllipticPi[((rp-r3) (r4-r1))/((rp-r4) (r3-r1)), JacobiAmplitude[X[\[Tau]], k], k] + EllipticPi[((rp-r3) (r4-r1))/((rp-r4) (r3-r1)), ArcSin[xo], k]);
Pim[\[Tau]_] := 2/Sqrt[(r3-r1) (r4-r2)] (r4-r3)/((rm-r3) (rm-r4)) (EllipticPi[((rm-r3) (r4-r1))/((rm-r4) (r3-r1)), JacobiAmplitude[X[\[Tau]], k], k] + EllipticPi[((rm-r3) (r4-r1))/((rm-r4) (r3-r1)), ArcSin[xo], k]);
Ip[\[Tau]_] := -\[Tau]/(rp-r3)-Pip[\[Tau]]; Iminus[\[Tau]_] := -\[Tau]/(rm-r3)-Pim[\[Tau]];

I\[Phi] = Function[{Global`\[Tau]}, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[Lambda]/2) Ip[Global`\[Tau]] - (rm- a \[Lambda]/2) Iminus[Global`\[Tau]])]], Listable];

<|"r" -> r, "I\[Phi]" -> I\[Phi], "\[Tau]x" -> Re[\[Tau]x]|>
]


DistantRadialMotionCase3[roots_, a_, \[Eta]_, \[Lambda]_] := Module[{r1, r2, r3, r4, rp, rm, A, B, xo, k, p1, f1, R1, I0o, \[Tau]x, X, r, \[Alpha]p, \[Alpha]m, Pip, Pim, Ip, Iminus, I\[Phi]},
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;
A = Re[Sqrt[(r3-r2) (r4-r2)]]; B = Re[Sqrt[(r3-r1) (r4-r1)]];
xo = (A-B)/(A+B); k = ((A+B)^2-(r2-r1)^2)/(4 A B); rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];
p1[\[Gamma]_, j_] := Sqrt[(\[Gamma]^2-1)/(j+(1-j) \[Gamma]^2)];
f1[\[Gamma]_, \[Xi]_, j_] := p1[\[Gamma], j]/2 Log[Abs[(p1[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2] + Sin[\[Xi]])/(p1[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2] - Sin[\[Xi]])]];
R1[\[Gamma]_, \[Xi]_, j_] := 1/(1-\[Gamma]^2) (EllipticPi[\[Gamma]^2/(\[Gamma]^2-1), \[Xi], j]-\[Gamma] f1[\[Gamma], \[Xi], j]);
I0o = 1/Sqrt[A B] EllipticF[ArcCos[xo], k];
\[Tau]x = I0o - 1/Sqrt[A B] EllipticF[ArcCos[(A (rp-r1) - B (rp-r2))/(A (rp-r1) + B (rp-r2))], k];
X[\[Tau]_] := Sqrt[A B] (\[Tau] - I0o);

r = Function[{Global`\[Tau]}, Evaluate[Re[((B r2 - A r1) + (B r2 + A r1) JacobiCN[X[Global`\[Tau]], k])/((B-A) + (B+A) JacobiCN[X[Global`\[Tau]], k])]], Listable];

\[Alpha]p = (B (rp - r2) + A (rp - r1))/(B (rp - r2) - A (rp - r1)); \[Alpha]m = (B (rm - r2) + A (rm - r1))/(B (rm - r2) - A (rm - r1));
Pip[\[Tau]_] := ((2 (r2-r1) Sqrt[A B])/(B (rp-r2) - A (rp-r1)))(R1[\[Alpha]p, JacobiAmplitude[X[\[Tau]], k], k]+R1[\[Alpha]p, ArcCos[xo], k]); 
Pim[\[Tau]_] := ((2 (r2-r1) Sqrt[A B])/(B (rm-r2) - A (rm-r1)))(R1[\[Alpha]m, JacobiAmplitude[X[\[Tau]], k], k]+R1[\[Alpha]m, ArcCos[xo], k]); 
Ip[\[Tau]_] := -((B+A) \[Tau] + Pip[\[Tau]])/(B (rp-r2) + A (rp-r1)); Iminus[\[Tau]_] := (-(((B+A) \[Tau] + Pim[\[Tau]])/(B (rm-r2) + A (rm-r1))));

I\[Phi] = Function[{Global`\[Tau]}, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[Lambda]/2) Ip[Global`\[Tau]] - (rm- a \[Lambda]/2) Iminus[Global`\[Tau]])]], Listable];

<|"r" -> r, "I\[Phi]" -> I\[Phi], "\[Tau]x" -> Re[\[Tau]x]|>
]


DistantRadialMotionCase4[roots_, a_, \[Eta]_, \[Lambda]_] := Module[{r1, r2, r3, r4, rp, rm, C, D, k, a2, b1, g0, I0o, \[Tau]x, x, X, p2, f2, S1, r, Pip, Pim, Ip, Iminus, I\[Phi]},
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;
C = Sqrt[(r3-r1) (r4-r2)]; D = Sqrt[(r3-r2) (r4-r1)];
k=(4 C D)/(C + D)^2; a2=Sqrt[-((r2-r1)^2/4)]; b1= (r3+r4)/2; g0=Sqrt[(4 a2^2-(C-D)^2)/((C+D)^2-4 a2^2)]; 
rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];
I0o = 2/(C+D) EllipticF[\[Pi]/2+ArcTan[g0], k];
\[Tau]x = I0o - 2/(C+D) EllipticF[ArcTan[(rp+b1)/a2]+ArcTan[g0], k];
x[\[Rho]_] = (\[Rho]+b1)/a2; X[\[Tau]_] = (C + D)/2 (-\[Tau]+I0o);
p2[\[Gamma]_, j_] := Sqrt[(\[Gamma]^2+1)/(1-j + \[Gamma]^2)];
f2[\[Gamma]_, \[Xi]_, j_] := p2[\[Gamma], j]/2 Log[Abs[(1-p2[\[Gamma], j])/(1+p2[\[Gamma], j]) (1+p2[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2] )/(1-p2[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2])]];
S1[\[Gamma]_, \[Xi]_, j_] := 1/(1+\[Gamma]^2) (EllipticF[\[Xi],j] +\[Gamma]^2 EllipticPi[1+\[Gamma]^2, \[Xi], j]-\[Gamma] f2[\[Gamma], \[Xi], j]);

r = Function[{Global`\[Tau]}, Evaluate[Re[-a2 ((g0-JacobiSC[X[Global`\[Tau]], k])/(1+g0 JacobiSC[X[Global`\[Tau]], k])) - b1]], Listable];

Pip[\[Tau]_] := -(2/(C+D)) ((1+g0^2)/(g0 (g0+x[rp]))) (S1[(g0 x[rp] - 1)/(g0 + x[rp]), JacobiAmplitude[X[\[Tau]], k], k] - S1[(g0 x[rp] - 1)/(g0 + x[rp]), \[Pi]/2+ArcTan[g0], k]);
Pim[\[Tau]_] := -(2/(C+D)) ((1+g0^2)/(g0 (g0+x[rm]))) (S1[(g0 x[rm] - 1)/(g0 + x[rm]), JacobiAmplitude[X[\[Tau]], k], k] - S1[(g0 x[rm] - 1)/(g0 + x[rm]), \[Pi]/2+ArcTan[g0], k]);
Ip[\[Tau]_] := g0/(a2 (1-g0 x[rp])) (\[Tau]-Pip[\[Tau]]);
Iminus[\[Tau]_] := g0/(a2 (1-g0 x[rm])) (\[Tau]-Pim[\[Tau]]);

I\[Phi] = Function[{Global`\[Tau]}, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[Lambda]/2) Ip[Global`\[Tau]] - (rm- a \[Lambda]/2) Iminus[Global`\[Tau]])]], Listable];

<|"r" -> r, "I\[Phi]" -> I\[Phi], "\[Tau]x" -> Re[\[Tau]x]|>
]


(* ::Section::Closed:: *)
(*Equator Intersections*)


EquatorIntersectionMinoTimes[a_, \[Eta]_, \[Lambda]_, \[Theta]o_, \[Nu]\[Theta]_, \[Tau]x_] := Module[{\[CapitalDelta]\[Theta], u1, u2, G\[Theta]o, equator\[Tau], j0, j, t},
\[CapitalDelta]\[Theta]= 1/2 (1-(\[Eta]+\[Lambda]^2)/a^2);
u1=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; u2=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
G\[Theta]o=-1/Sqrt[-u1 a^2] EllipticF[Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1]; (*Sometimes the ArcSin argument is slightly over 1 probably due to numerical errors*)
equator\[Tau] = {};
If[(-\[Nu]\[Theta] G\[Theta]o) < 0, j0=1, j0=0];
For[j=j0, True, j++,
t=1/Sqrt[-u1 a^2] EllipticF[j \[Pi], u2/u1]-\[Nu]\[Theta] G\[Theta]o;
If[t>\[Tau]x, Break[]];
equator\[Tau] = Append[equator\[Tau], t]];
equator\[Tau]
]


(* ::Section:: *)
(*Public Functions*)


Options[KerrNullGeoDistant] = {"Rotation" -> "Counterclockwise"}
SyntaxInformation[KerrNullGeoDistant] = {"ArgumentsPattern"->{_,_,_,_}};

KerrNullGeoDistant[a_, \[Theta]o_, \[Alpha]_, \[Beta]_, OptionsPattern[]] := Module[ {consts, \[Eta], \[Lambda], roots, r1, r2, r3, r4, rp, rm, k, r, \[Theta], \[Phi], G\[Phi], I\[Phi], \[Tau]x, assoc, equator\[Tau], type},

If[a<=0 || a>=1, Message[KerrNullGeoDistantDistant::OutOfBounds, "Parameter a must be between 0 and 1."]; Return[];];
If[\[Theta]o<0 || \[Theta]o>\[Pi], Message[KerrNullGeoDistantDistant::OutOfBounds, "Parameter \[Theta]o must be between 0 and \[Pi]."]; Return[];];

If[OptionValue["Rotation"]=="Clockwise", \[Alpha]=-\[Alpha]];

consts = DistantNullConstantsOfMotion[a, \[Theta]o, \[Alpha], \[Beta]];
{\[Lambda], \[Eta]} = {"\[Lambda]", "\[Eta]"} /. consts;
If[Abs[\[Eta]]< 10 $MachineEpsilon || Abs[\[Eta]+(Abs[\[Lambda]]-a)^2]< 10 $MachineEpsilon, \[Eta]=\[Eta] + 10 $MachineEpsilon]; (*To avoid undefined expresions in polar motion*)
If[\[Eta]>0, 
  {\[Theta], G\[Phi]} = {"\[Theta]", "G\[Phi]"} /. OrdinaryPolarMotion[a, \[Eta], \[Lambda], \[Theta]o, -Sign[\[Beta]], ReturnValues -> "OmmitT"],
  {\[Theta], G\[Phi]} = {"\[Theta]", "G\[Phi]"} /. VorticalPolarMotion[a, \[Eta], \[Lambda], \[Theta]o, -Sign[\[Beta]], ReturnValues -> "OmmitT"]; (*Minus Sign[\[Beta]] because we use negative Mino time*)
];

roots = RadialRoots[a, \[Eta], \[Lambda]];
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;

type = "PhotonCapture";
If[Im[r2] != 0,
  {r, I\[Phi], \[Tau]x} = {"r", "I\[Phi]", "\[Tau]x"} /. DistantRadialMotionCase4[roots, a, \[Eta], \[Lambda]],
  If[Im[r4] != 0, 
    {r, I\[Phi], \[Tau]x} = {"r", "I\[Phi]", "\[Tau]x"} /. DistantRadialMotionCase3[roots, a, \[Eta], \[Lambda]],
    {r, I\[Phi], \[Tau]x} = {"r", "I\[Phi]", "\[Tau]x"} /. DistantRadialMotionCase2[roots, a, \[Eta], \[Lambda]];
    If[r3>1-Sqrt[1-a^2], type = "PhotonEscape"]
  ]
];

If[\[Eta]>0, 
equator\[Tau] = EquatorIntersectionMinoTimes[a, \[Eta], \[Lambda], \[Theta]o, -Sign[\[Beta] + $MachineEpsilon], \[Tau]x],
equator\[Tau] = {}
];

\[Phi]=Function[{Global`\[Tau]}, Evaluate[I\[Phi][Global`\[Tau]]+\[Lambda] G\[Phi][Global`\[Tau]]], Listable];

assoc = <|
	"Trajectory" -> {r, \[Theta], \[Phi]},
	"ConstantsOfMotion" -> consts,
	"RadialRoots" -> {r1, r2, r3, r4},
	"EquatorIntersectionMinoTimes" -> equator\[Tau],
	"TrajectoryType" -> type,
	"MinoTimeOfCapture" -> \[Tau]x
|>;

KerrNullGeoDistantFunction[a, \[Theta]o, \[Alpha], \[Beta], assoc]
]


KerrNullGeoDistantFunction[a_, \[Theta]o_, \[Alpha]_, \[Beta]_, assoc_][\[Tau]_/;StringQ[\[Tau]] == False] := Through[assoc["Trajectory"][\[Tau]]];
KerrNullGeoDistantFunction[a_, \[Theta]o_, \[Alpha]_, \[Beta]_, assoc_][y_?StringQ] := assoc[y];
Keys[g_KerrNullGeoDistantFunction]^:=Keys[g[[5]]];


(* ::Section::Closed:: *)
(*Close the Package*)


End[];

EndPackage[];
