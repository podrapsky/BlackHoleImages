(* ::Package:: *)

(* ::Title:: *)
(*KerrNullGeodesics Package*)


(* ::Text:: *)
(*The  KerrNullGeodesics  package  generates  null  geodesics in the exterior of a Kerr black hole (BH) based on the expression of the analytical solution as given by Gralla & Lupsasca (2020). Null geodesics of the Kerr exterior. Physical Review D, 101(4), 044032,  arXiv:1910.12881. *)
(**)
(*At intermediate steps the package uses the \[Alpha],\[Beta] parametrization of the null geodesic parameters as introduced by Bardeen, J. M. (1973). Timelike and null geodesics in the Kerr metric. Black holes, 215.*)
(**)
(*The primary use case for the package is the imaging of accretion disks near black holes. As such, the public functions return objects that include lists of intersections through the equatorial plane for convenience. However, the package is written as flexible enough to find many other use cases.*)


(* ::Section:: *)
(*Define Usage for Public Functions*)


BeginPackage["BlackHoleImages`KerrNullGeodesics`"]

KerrNullGeo::usage = "KerrNullGeo[a,xs,ps] returns a KerrNullGeoFunction which stores information about the trajectory of a light-ray starting from specified initial conditions. The black hole spin a, position xs, and wavevector ps  are assumed to be given in units of the BH mass, unless mass M is specified (optional argument).";
KerrNullGeoFunction::usage = "KerrNullGeoFunction[a,xs,ps,M,assoc] is an object for storing the trajectory and its parameters in the assoc Association.";
KerrNullGeoDistant::usage = "KerrNullGeoDistant[a,\[Theta]o,\[Alpha],\[Beta]] returns a KerrNullGeoDistantFunction which stores information about the trajectory of a light-ray scattering off the black hole from infinity. The spin a, and Bardeen's impact parameters \[Alpha], \[Beta] are assumed to be given in units of the BH mass";
KerrNullGeoDistantFunction::usage = "KerrNullGeoDistantFunction[a,\[Theta]o,\[Alpha],\[Beta],assoc] is an object for storing the trajectory and its parameters in the assoc Association.";

Begin["`Private`"];


KerrNullGeo::OutOfBounds = "Out of bounds error: `1`"
KerrNullGeo::ListSize = "Parameters `1` or `2` is not a list of length `3`."


(* ::Section:: *)
(*Constants of Motion*)


NullConstantsOfMotion[a_, \[Theta]s_, pts_, p\[Theta]s_, p\[Phi]s_] := <|
	(*Specific angular momentum wrt spin axis:*)
	"\[ScriptL]" -> -p\[Phi]s/pts,
	(*"Vertical" part of specific angular momentum squared, also Carter constant per energy squared:*)
	"\[Eta]" -> (p\[Theta]s^2 - (Cos[\[Theta]s])^2 (a^2 pts^2 - p\[Phi]s^2/(Sin[\[Theta]s])^2))/pts^2 
|>


DistantNullConstantsOfMotion[a_, \[Theta]o_, \[Alpha]_, \[Beta]_] := <|
	(*Bardeen's \[Alpha],\[Beta], \[Theta]o is the angle at which the observer is viewing the BH or the \[Theta] angle to which the lightray asymptotes*)
	"\[ScriptL]" ->  -\[Alpha] Sin[\[Theta]o],
	"\[Eta]" -> \[Beta]^2+(\[Alpha]^2-a^2) (Cos[\[Theta]o])^2
|>


(* ::Section::Closed:: *)
(*Polar Motion*)


(* ::Text:: *)
(*Polar motion takes two forms, ordinary motion, and vortical motion as described in Section III of  Gralla & Lupsasca, arXiv:1910.12881v3 (we follow their notation closely with equation numbers in comments).*)


Options[OrdinaryPolarMotion] = {"ReturnValues" -> "All"}
Options[VorticalPolarMotion] = {"ReturnValues" -> "All"}


OrdinaryPolarMotion[a_, \[Eta]_, \[ScriptL]_, \[Theta]o_, \[Nu]\[Theta]_, \[Lambda]x_, OptionsPattern[]] := Module[{\[CapitalDelta]\[Theta], u1, u2, EPrime, G\[Theta]o, Gto, G\[Phi]o, \[CapitalPsi], \[Theta], Gt, G\[Phi], return}, 
(*Eq. (19):*)
\[CapitalDelta]\[Theta]= 1/2 (1-(\[Eta]+\[ScriptL]^2)/a^2);
u1=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; u2=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
(*Eq. (29):*)
G\[Theta]o=-1/Sqrt[-u1 a^2] EllipticF[Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1];
(*Eq. (46):*)
\[CapitalPsi][\[Lambda]_] := JacobiAmplitude[Sqrt[-u1 a^2] (\[Lambda] + \[Nu]\[Theta] G\[Theta]o), u2/u1];
(*Eq.(30):*)
G\[Phi]o=-1/Sqrt[-u1 a^2] EllipticPi[u2, Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1];
G\[Phi]= Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[1/Sqrt[-u1 a^2] EllipticPi[u2, \[CapitalPsi][Global`\[Lambda]], u2/u1]-\[Nu]\[Theta] G\[Phi]o]]], Listable]; 

\[Theta] = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[ArcCos[-\[Nu]\[Theta] Sqrt[u2] Sin[\[CapitalPsi][Global`\[Lambda]]]]]]], Listable];

If[OptionValue["ReturnValues"]!="OmitT", 

  EPrime[\[Phi]_, k_] := (EllipticE[\[Phi], k]-EllipticF[\[Phi], k])/(2 k);
  (*Eq. (31):*)
  Gto=2 u2/Sqrt[-u1 a^2] EPrime[Re[ArcSin[Cos[\[Theta]o]/Sqrt[u2]]], u2/u1];
  Gt=Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[-2 u2/Sqrt[-u1 a^2] EPrime[\[CapitalPsi][Global`\[Lambda]], u2/u1] - \[Nu]\[Theta] Gto]]], Listable];
  return = <|"\[Theta]" -> \[Theta], "Gt" -> Gt, "G\[Phi]" -> G\[Phi]|>,

  return = <|"\[Theta]" -> \[Theta], "G\[Phi]" -> G\[Phi]|>
];
  
return
]


VorticalPolarMotion[a_, \[Eta]_, \[ScriptL]_, \[Theta]o_, \[Nu]\[Theta]_, \[Lambda]x_, OptionsPattern[]] := Module[{h, \[CapitalDelta]\[Theta], u1, u2, \[CapitalUpsilon], \[CapitalUpsilon]\[Lambda], G\[Theta]o, Gto, G\[Phi]o,\[Theta], Gt, G\[Phi], return},
(*Eq. (54)*)
h = Sign[Cos[\[Theta]o]];
(*Eq. (19):*)
\[CapitalDelta]\[Theta]= 1/2 (1-(\[Eta]+\[ScriptL]^2)/a^2);
u1=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; u2=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2];
(*Eq. (59):*)
\[CapitalUpsilon][\[CapitalTheta]_] := ArcSin[Sqrt[((Cos[\[CapitalTheta]])^2-u1)/(u2-u1)]];
(*Eq. (56):*)
G\[Theta]o = -(h/Sqrt[u1 a^2]) EllipticF[\[CapitalUpsilon][\[Theta]o], 1-u2/u1];
(*Eq. (66):*)
\[CapitalUpsilon]\[Lambda][\[Lambda]_] := JacobiAmplitude[Sqrt[u1 a^2] (\[Lambda] + \[Nu]\[Theta] G\[Theta]o), 1-u2/u1];

(*Eq. (57):*)
G\[Phi]o = -(h/((1-u1) Sqrt[u1 a^2])) EllipticPi[(u2-u1)/(1-u1), \[CapitalUpsilon][\[Theta]o], 1-u2/u1];
G\[Phi] = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[1/((1-u1) Sqrt[u1 a^2]) EllipticPi[(u2-u1)/(1-u1), \[CapitalUpsilon]\[Lambda][Global`\[Lambda]], 1-u2/u1]-\[Nu]\[Theta] G\[Phi]o]]], Listable];

\[Theta] = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[ArcCos[h Sqrt[u1+(u2-u1) (Sin[\[CapitalUpsilon]\[Lambda][Global`\[Lambda]]])^2]]]]], Listable];

If[OptionValue["ReturnValues"]!="OmitT",
  (*Eq. (58):*)
  Gto = -h Sqrt[u1/a^2] EllipticE[\[CapitalUpsilon][\[Theta]o], 1-u2/u1];
  Gt = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Sqrt[u1/a^2] EllipticE[\[CapitalUpsilon]\[Lambda][Global`\[Lambda]], 1-u2/u1]-\[Nu]\[Theta] Gto]]], Listable];
  return = <|"\[Theta]" -> \[Theta], "Gt" -> Gt, "G\[Phi]" -> G\[Phi]|>,

  return = <|"\[Theta]" -> \[Theta], "G\[Phi]" -> G\[Phi]|>
];

return
]


(* ::Section::Closed:: *)
(*Radial Roots*)


(* ::Text:: *)
(*Taken from Section IV-A of Gralla & Lupsasca, arXiv:1910.12881v3. The computation uses Ferrari's method to find the roots of the radial potential. This yields a more transparent structure of real and imaginary parts of the roots than e.g. the brute-force Mathematica implementation.*)


RadialRoots[a_, \[Eta]_, \[ScriptL]_] := Module[{A, B, C, P, Q, \[CapitalOmega]1, \[CapitalOmega]2, \[Omega]1, \[Omega]2, z},
(*Eq. (79)-(81):*)
A=a^2-\[Eta]-\[ScriptL]^2; B=2(\[Eta]+(\[ScriptL]-a)^2); C=-a^2 \[Eta];
(*Eq. (85)-(86):*)
P= -A^2/12 - C; Q=-A/3 ((A/6)^2-C) - B^2/8;
(*Eq. (90):*)
\[CapitalOmega]1 = -Q/2 - Sqrt[(P/3)^3+(Q/2)^2];
\[CapitalOmega]2 = -Q/2 + Sqrt[(P/3)^3+(Q/2)^2];
\[Omega]1=If[Element[\[CapitalOmega]1, Reals], Surd[\[CapitalOmega]1, 3], Power[\[CapitalOmega]1, 1/3]];
\[Omega]2=If[Element[\[CapitalOmega]2, Reals], Surd[\[CapitalOmega]2, 3] , Power[\[CapitalOmega]2, 1/3]];
(*Eq. (94):*)
z=Sqrt[(\[Omega]2+\[Omega]1-A/3)/2];

<|(*Eq. (95):*)
  "r1" -> -z-Sqrt[-A/2-z^2+B/(4 z)],
  "r2" -> -z+Sqrt[-A/2-z^2+B/(4 z)],
  "r3" -> +z-Sqrt[-A/2-z^2-B/(4 z)],
  "r4" -> +z+Sqrt[-A/2-z^2-B/(4 z)]
|>
]


(* ::Section:: *)
(*Radial Motion*)


(* ::Text:: *)
(*There are various cases for the values of the four roots of the radial potential for null geodesics as discussed in Section IV B of Gralla & Lupsasca, arXiv:1910.12881v3 (Equation numbers below). The various cases correspond to light-rays "oscillating" between a white-hole and black-hole horizon (Case I), approaching from past null infinity and scattering back to future null infinity (Case II),  ?? (Case III), and a light-ray approaching from past null infinity, plunging into the BH, and emerging from the white-hole horizon and approaching future null infinity.  TODO*)


Options[RadialMotionCase2] = {"Observer" -> "Regular"}
Options[RadialMotionCase3] = {"Observer" -> "Regular"}
Options[RadialMotionCase4] = {"Observer" -> "Regular"}


RadialMotionCase1[roots_, a_, \[Eta]_, \[ScriptL]_, rs_, \[Nu]r_] := Module[{r1, r2, r3, r4, rp, rm, xs, k,  I0s, \[Lambda]x, X, r, Pi1, E\[Lambda], Pip, Pim, Ip, Iminus, I1, deriv, I2, I\[Phi], It},
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;
xs=Sqrt[(rs-r2)/(rs-r1) (r3-r1)/(r3-r2)]; k=((r3-r2) (r4-r1))/((r3-r1) (r4-r2)); rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];

I0s = 2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[xs], k];
If[\[Nu]r<0,
  \[Lambda]x = -(2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[Sqrt[((rp-r2) (r3-r1))/((rp-r1) (r3-r2))]], k] - I0s),
  \[Lambda]x = 4/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[Sqrt[((r3-r4) (r3-r1))/((r3-r1) (r3-r2))]], k] - 2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[Sqrt[((rp-r2) (r3-r1))/((rp-r1) (r3-r2))]], k] - I0s;
];
X[\[Lambda]_] := Sqrt[(r3-r1) (r4-r2)]/2 (\[Lambda] + \[Nu]r I0s);

r = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(r2 (r3-r1) - r1 (r3-r2) (JacobiSN[X[Global`\[Lambda]], k])^2)/((r3-r1) - (r3-r2) (JacobiSN[X[Global`\[Lambda]], k])^2)]]]], Listable];

Pi1[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)](EllipticPi[(r3-r2)/(r3-r1), JacobiAmplitude[X[\[Lambda]], k], k] -\[Nu]r EllipticPi[(r3-r2)/(r3-r1), ArcSin[xs], k]);
E\[Lambda][\[Lambda]_] := Sqrt[(r3-r1)(r4-r2)](EllipticE[JacobiAmplitude[X[\[Lambda]], k], k] - \[Nu]r EllipticE[ArcSin[xs], k]);
Pip[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)] (r2-r1)/((rp-r1) (rp-r2)) (EllipticPi[((rp-r1) (r3-r2))/((rp-r2) (r3-r1)), JacobiAmplitude[X[\[Lambda]], k], k] - \[Nu]r EllipticPi[((rp-r1) (r3-r2))/((rp-r2) (r3-r1)), ArcSin[xs], k]);
Pim[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)] (r2-r1)/((rm-r1) (rm-r2)) (EllipticPi[((rm-r1) (r3-r2))/((rm-r2) (r3-r1)), JacobiAmplitude[X[\[Lambda]], k], k] - \[Nu]r EllipticPi[((rm-r1) (r3-r2))/((rm-r2) (r3-r1)), ArcSin[xs], k]);
Ip[\[Lambda]_] := -\[Lambda]/(rp-r1)-Pip[\[Lambda]]; Iminus[\[Lambda]_] := -\[Lambda]/(rm-r1)-Pim[\[Lambda]];
I1[\[Lambda]_] := r1 \[Lambda] + (r2-r1) Pi1[\[Lambda]];
deriv[expr_, var_] := expr'[var] /. Re'[e_]:> 1;
I2[\[Lambda]_] := deriv[r, \[Lambda]]/(r[\[Lambda]]-r1) - \[Nu]r Sqrt[(rs^2+a^2-a \[ScriptL])^2-(rs^2-2 rs + a^2) (\[Eta]+(\[ScriptL]-a)^2)]/(rs-r1)- (r1 r4 + r2 r3)/2 \[Lambda] - E\[Lambda][\[Lambda]];

I\[Phi] = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[ScriptL]/2) Ip[Global`\[Lambda]] - (rm- a \[ScriptL]/2) Iminus[Global`\[Lambda]])]]]], Listable];
It = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 I1[Global`\[Lambda]] + I2[Global`\[Lambda]]]]]], Listable];

<|"r" -> r, "I\[Phi]" -> I\[Phi], "It" -> It, "\[Lambda]x" -> Re[\[Lambda]x]|>
]


RadialMotionCase2[roots_, a_, \[Eta]_, \[ScriptL]_, Optional[rs_/; NumberQ[rs], 1], Optional[\[Nu]r_/; NumberQ[\[Nu]r], -1], OptionsPattern[]] := Module[{Distant, r1, r2, r3, r4, rp, rm, xo, k, I0o, \[Lambda]x, X, r, Pi1, E\[Lambda], Pip, Pim, Ip, Iminus, I1, deriv, I2, I\[Phi], It},
Distant=If[OptionValue["Observer"]=="Distant", True, False];
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"}/. roots;
k=((r3-r2) (r4-r1))/((r3-r1) (r4-r2)); rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];

xo= If[Distant, Sqrt[(r3-r1)/(r4-r1)],Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)]]; 

I0o = 2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[xo], k];
\[Lambda]x = Re[
If[\[Nu]r<0, 
  If[r4<rm,
    -(2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[Sqrt[((rp-r4) (r3-r1))/((rp-r3) (r4-r1))]], k] - I0o),
    If[Distant,
    2 I0o,
    2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[Sqrt[(r3-r1)/(r4-r1)]], k] + I0o (*EllipticF[0, k] == 0 at turning point r4*)
    ]
  ],
  2/Sqrt[(r3-r1) (r4-r2)] EllipticF[ArcSin[Sqrt[(r3-r1)/(r4-r1)]], k] - I0o;
]
];
X[\[Lambda]_] := Sqrt[(r3-r1) (r4-r2)]/2 (\[Lambda] + \[Nu]r I0o);

r = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(r4 (r3-r1) - r3 (r4-r1) (JacobiSN[X[Global`\[Lambda]], k])^2)/((r3-r1) - (r4-r1) (JacobiSN[X[Global`\[Lambda]], k])^2)]]]], Listable];

E\[Lambda][\[Lambda]_] := Sqrt[(r3-r1)(r4-r2)](EllipticE[JacobiAmplitude[X[\[Lambda]], k], k] - \[Nu]r EllipticE[ArcSin[xo], k]);
Pip[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)] (r4-r3)/((rp-r3) (rp-r4)) (EllipticPi[((rp-r3) (r4-r1))/((rp-r4) (r3-r1)), JacobiAmplitude[X[\[Lambda]], k], k] - \[Nu]r EllipticPi[((rp-r3) (r4-r1))/((rp-r4) (r3-r1)), ArcSin[xo], k]);
Pim[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)] (r4-r3)/((rm-r3) (rm-r4)) (EllipticPi[((rm-r3) (r4-r1))/((rm-r4) (r3-r1)), JacobiAmplitude[X[\[Lambda]], k], k] - \[Nu]r EllipticPi[((rm-r3) (r4-r1))/((rm-r4) (r3-r1)), ArcSin[xo], k]);
Ip[\[Lambda]_] := -\[Lambda]/(rp-r3)-Pip[\[Lambda]]; Iminus[\[Lambda]_] := -\[Lambda]/(rm-r3)-Pim[\[Lambda]];
deriv[expr_, var_] := expr'[var] /. Re'[e_]:> 1;

I\[Phi] = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[ScriptL]/2) Ip[Global`\[Lambda]] - (rm- a \[ScriptL]/2) Iminus[Global`\[Lambda]])]]]], Listable];

If[Distant,
  Pi1[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)](EllipticPi[(r4-r1)/(r3-r1), JacobiAmplitude[X[\[Lambda]], k], k] + EllipticF[ArcSin[Sqrt[(r3-r1)/(r4-r1)]], k] - EllipticPi[(r3-r2)/(r4-r2), ArcSin[Sqrt[(r3-r1)/(r4-r1)]], k] + 1/(2 Sqrt[(1-(r3-r2)/(r4-r2))((r4-r1)/(r3-r1)-1)]) Log[4/(r3-r1+r4-r2)]);
  I1[\[Lambda]_] := r3 \[Lambda] + (r4-r3) Pi1[\[Lambda]];
  I2[\[Lambda]_] := deriv[r, \[Lambda]]/(r[\[Lambda]]-r3) + r3 - (r1 r4 + r2 r3)/2 \[Lambda] - E\[Lambda][\[Lambda]];
  It = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 I1[Global`\[Lambda]] + I2[Global`\[Lambda]] + 2 Log[2]]]]], Listable],
  (*else*)
  Pi1[\[Lambda]_] := 2/Sqrt[(r3-r1) (r4-r2)](EllipticPi[(r4-r1)/(r3-r1), JacobiAmplitude[X[\[Lambda]], k], k] -\[Nu]r EllipticPi[(r4-r1)/(r3-r1), ArcSin[xo], k]);
  I1[\[Lambda]_] := r3 \[Lambda] + (r4-r3) Pi1[\[Lambda]];
  I2[\[Lambda]_] := deriv[r, \[Lambda]]/(r[\[Lambda]]-r3) - \[Nu]r Sqrt[(rs^2+a^2-a \[ScriptL])^2-(rs^2-2 rs + a^2) (\[Eta]+(\[ScriptL]-a)^2)]/(rs-r3)- (r1 r4 + r2 r3)/2 \[Lambda] - E\[Lambda][\[Lambda]];
  It = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 I1[Global`\[Lambda]] + I2[Global`\[Lambda]]]]]], Listable];
];



<|"r" -> r, "I\[Phi]" -> I\[Phi], "It" -> It, "\[Lambda]x" -> Re[\[Lambda]x]|>
]


RadialMotionCase3[roots_, a_, \[Eta]_, \[ScriptL]_, Optional[rs_/; NumberQ[rs], 1], Optional[\[Nu]r_/; NumberQ[\[Nu]r], -1], OptionsPattern[]] := Module[{Distant, r1, r2, r3, r4, rp, rm, A, B, xo, k, p1, f1, RedR1, RedR2, R1, R2, I0o, \[Lambda]x, X, r, \[Alpha]0, \[Alpha]p, \[Alpha]m, Pip, Pim, Ip, Iminus, I\[Phi], Pi1, Pi2, I1, I2, It},
If[OptionValue["Observer"]=="Distant", Distant=True, Distant=False];
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;
A = Abs[Sqrt[(r3-r2) (r4-r2)]]; B = Abs[Sqrt[(r3-r1) (r4-r1)]];
k = ((A+B)^2-(r2-r1)^2)/(4 A B); rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];

\[Alpha]0 = (B+A)/(B-A);
\[Alpha]p = (B (rp - r2) + A (rp - r1))/(B (rp - r2) - A (rp - r1)); \[Alpha]m = (B (rm - r2) + A (rm - r1))/(B (rm - r2) - A (rm - r1));

p1[\[Gamma]_, j_] := Sqrt[(\[Gamma]^2-1)/(j+(1-j) \[Gamma]^2)];
f1[\[Gamma]_, \[Xi]_, j_] := p1[\[Gamma], j]/2 Log[Abs[(p1[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2] + Sin[\[Xi]])/(p1[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2] - Sin[\[Xi]])]];
RedR1[\[Gamma]_, \[Xi]_, j_] := 1/(1-\[Gamma]^2) (EllipticF[\[Xi], j] - EllipticPi[(j(\[Gamma]^2-1))/\[Gamma]^2, \[Xi], j] - (\[Gamma] Sqrt[A B])/(r2-r1) Log[(4(r2-r1))/(B^2-A^2)] + (\[Gamma] Sqrt[A B])/(r2-r1) Log[(B^2-A^2)/(4(r2-r1))+(A B (r2-r1))/(B^2-A^2)]); 
RedR2[\[Gamma]_, \[Xi]_, j_] := 1/(\[Gamma]^2-1) (EllipticF[\[Xi], j] -\[Gamma]^2/(j + (1-j)\[Gamma]^2) EllipticE[\[Xi], j]);
R1[\[Gamma]_, \[Xi]_, j_] := 1/(1-\[Gamma]^2) (EllipticPi[\[Gamma]^2/(\[Gamma]^2-1), \[Xi], j]-\[Gamma] f1[\[Gamma], \[Xi], j]);
R2[\[Gamma]_, \[Xi]_, j_] := 1/(\[Gamma]^2-1) (EllipticF[\[Xi], j] -\[Gamma]^2/(j + (1-j)\[Gamma]^2) (EllipticE[\[Xi], j] - (\[Gamma] Sin[\[Xi]]Sqrt[1-j (Sin[\[Xi]])^2])/(1+\[Gamma] Cos[\[Xi]]))) + 1/(j + (1-j) \[Gamma]^2) (2j - \[Gamma]^2/(\[Gamma]^2-1))R1[\[Gamma], \[Xi], j];

If[Distant,xo = (A-B)/(A+B), xo = (1 - (B (rs-r2))/(A (r-r1)))/(1 + (B (rs-r2))/(A (r-r1)))];
I0o = 1/Sqrt[A B] EllipticF[ArcCos[xo], k];
\[Lambda]x = \[Nu]r (1/Sqrt[A B] EllipticF[ArcCos[(A (rp-r1) - B (rp-r2))/(A (rp-r1) + B (rp-r2))], k] - I0o);
X[\[Lambda]_] := Sqrt[A B] (\[Lambda] +\[Nu]r I0o);

r = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[((B r2 - A r1) + (B r2 + A r1) JacobiCN[X[Global`\[Lambda]], k])/((B-A) + (B+A) JacobiCN[X[Global`\[Lambda]], k])]]]], Listable];
  
Pip[\[Lambda]_] := ((2 (r2-r1) Sqrt[A B])/(B (rp-r2) - A (rp-r1)))(R1[\[Alpha]p, JacobiAmplitude[X[\[Lambda]], k], k]-\[Nu]r R1[\[Alpha]p, ArcCos[xo], k]); 
Pim[\[Lambda]_] := ((2 (r2-r1) Sqrt[A B])/(B (rm-r2) - A (rm-r1)))(R1[\[Alpha]m, JacobiAmplitude[X[\[Lambda]], k], k]-\[Nu]r R1[\[Alpha]m, ArcCos[xo], k]);
Ip[\[Lambda]_] := -((B+A) \[Lambda] + Pip[\[Lambda]])/(B (rp-r2) + A (rp-r1)); Iminus[\[Lambda]_] := -((B+A) \[Lambda] + Pim[\[Lambda]])/(B (rm-r2) + A (rm-r1));

I\[Phi] = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[ScriptL]/2) Ip[Global`\[Lambda]] - (rm- a \[ScriptL]/2) Iminus[Global`\[Lambda]])]]]], Listable];

If[Distant,
  Pi1[\[Lambda]_] := (2(r2-r1)Sqrt[A B])/(B^2-A^2) (R1[\[Alpha]0, JacobiAmplitude[X[\[Lambda]], k], k] + RedR1[\[Alpha]0, ArcCos[xo], k]);
  Pi2[\[Lambda]_] := ((2(r2-r1)Sqrt[A B])/(B^2-A^2))^2 (R2[\[Alpha]0, JacobiAmplitude[X[\[Lambda]], k], k] + RedR2[\[Alpha]0, ArcCos[xo], k]);
  I1[\[Lambda]_] := ((B r2 + A r1)/(B + A))\[Lambda] + Pi1[\[Lambda]];
  I2[\[Lambda]_] := ((B r2 + A r1)/(B + A))^2 \[Lambda] + 2((B r2 + A r1)/(B + A)) (2(r2-r1)Sqrt[A B])/(B^2-A^2) R1[\[Alpha]0, JacobiAmplitude[X[\[Lambda]], k], k] + Sqrt[A B] Pi2[\[Lambda]] + (A^2-B^2)/(2 (r2-r1))-(r1+r2)+(B r2 + A r1)/(B + A);
  It = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 I1[Global`\[Lambda]] + I2[Global`\[Lambda]] + 2 Log[2]]]]], Listable],
  (*else*)
  Pi1[\[Lambda]_] := (2(r2-r1)Sqrt[A B])/(B^2-A^2) (R1[\[Alpha]0, JacobiAmplitude[X[\[Lambda]], k], k] - \[Nu]r R1[\[Alpha]0, ArcCos[xo], k]);
  Pi2[\[Lambda]_] := ((2(r2-r1)Sqrt[A B])/(B^2-A^2))^2 (R2[\[Alpha]0, JacobiAmplitude[X[\[Lambda]], k], k] - \[Nu]r R2[\[Alpha]0, ArcCos[xo], k]);
  I1[\[Lambda]_] := ((B r2 + A r1)/(B + A))\[Lambda] + Pi1[\[Lambda]];
  I2[\[Lambda]_] := ((B r2 + A r1)/(B + A))^2 \[Lambda] + 2((B r2 + A r1)/(B + A)) Pi1[\[Lambda]] + Sqrt[A B] Pi2[\[Lambda]];
  It = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 I1[Global`\[Lambda]] + I2[Global`\[Lambda]]]]]], Listable];
];

<|"r" -> r, "I\[Phi]" -> I\[Phi], "It" -> It, "\[Lambda]x" -> Re[\[Lambda]x]|>
]


RadialMotionCase4[roots_, a_, \[Eta]_, \[ScriptL]_, Optional[rs_/; NumberQ[rs], 1], Optional[\[Nu]r_/; NumberQ[\[Nu]r], -1], OptionsPattern[]] := Module[{Distant, r1, r2, r3, r4, rp, rm, C, D, k, a2, b1, g0, I0o, \[Lambda]x, x, X, p2, f2, RedS1, S1, RedS2, S2, r, Pip, Pim, Pi1, Pi2, Ip, Iminus, I1, I2, It, I\[Phi]},
If[OptionValue["Observer"]=="Distant", Distant=True, Distant=False];
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;
C = Sqrt[(r3-r1) (r4-r2)]; D = Sqrt[(r3-r2) (r4-r1)];
k=(4 C D)/(C + D)^2; a2=Sqrt[-((r2-r1)^2/4)]; b1= (r3+r4)/2; g0=Sqrt[(4 a2^2-(C-D)^2)/((C+D)^2-4 a2^2)]; 
rp = 1+Sqrt[1-a^2]; rm = 1-Sqrt[1-a^2];

x[\[Rho]_] = (\[Rho]+b1)/a2;
p2[\[Gamma]_, j_] := Sqrt[(\[Gamma]^2+1)/(1-j + \[Gamma]^2)];
f2[\[Gamma]_, \[Xi]_, j_] := p2[\[Gamma], j]/2 Log[Abs[(1-p2[\[Gamma], j])/(1+p2[\[Gamma], j]) (1+p2[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2] )/(1-p2[\[Gamma], j] Sqrt[1- j (Sin[\[Xi]])^2])]];
RedS1[\[Gamma]_, \[Xi]_, j_] := EllipticF[\[Xi], j] - 1/(1+\[Gamma]^2) (\[Gamma]^2 EllipticPi[j/(1+\[Gamma]^2), \[Xi], j] + (\[Gamma] (C+D))/(4 a2) Log[(64 a2^2)/((2 a2 + C + D)^2 (\[Gamma]^2 (C+D)^2 + 4 a2^2))]);
S1[\[Gamma]_, \[Xi]_, j_] := 1/(1+\[Gamma]^2) (EllipticF[\[Xi],j] +\[Gamma]^2 EllipticPi[1+\[Gamma]^2, \[Xi], j]-\[Gamma] f2[\[Gamma], \[Xi], j]);
RedS2[\[Gamma]_, \[Xi]_, j_] := -1/((1+\[Gamma]^2)(1-j+\[Gamma]^2)) ((1-j)EllipticF[\[Xi], j] + \[Gamma]^2 EllipticE[\[Xi], j] - \[Gamma]^3);
S2[\[Gamma]_, \[Xi]_, j_] := -1/((1+\[Gamma]^2)(1-j+\[Gamma]^2)) ((1-j)EllipticF[\[Xi], j] + \[Gamma]^2 EllipticE[\[Xi], j] + (\[Gamma]^2 Sqrt[1-j (Sin[\[Xi]])^2](\[Gamma]-Tan[\[Xi]]))/(1+\[Gamma] Tan[\[Xi]]) - \[Gamma]^3) + (1/(1+\[Gamma]^2)+(1-j)/(1-j+\[Gamma]^2))S1[\[Gamma], \[Xi], j];

If[Distant, I0o = 2/(C+D) EllipticF[\[Pi]/2+ArcTan[g0], k], I0o = 2/(C+D) EllipticF[ArcTan[x[rs]]+ArcTan[g0], k]];
\[Lambda]x = \[Nu]r (2/(C+D) EllipticF[ArcTan[(rp+b1)/a2]+ArcTan[g0], k] - I0o);
X[\[Lambda]_] = (C + D)/2 (\[Nu]r \[Lambda]+I0o);

r = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[-a2 ((g0-JacobiSC[X[Global`\[Lambda]], k])/(1+g0 JacobiSC[X[Global`\[Lambda]], k])) - b1]]]], Listable];

If[Distant,
  Pip[\[Lambda]_] := -(2/(C+D)) ((1+g0^2)/(g0 (g0+x[rp]))) (S1[(g0 x[rp] - 1)/(g0 + x[rp]), JacobiAmplitude[X[\[Lambda]], k], k] - S1[(g0 x[rp] - 1)/(g0 + x[rp]), \[Pi]/2+ArcTan[g0], k]);
  Pim[\[Lambda]_] := -(2/(C+D)) ((1+g0^2)/(g0 (g0+x[rm]))) (S1[(g0 x[rm] - 1)/(g0 + x[rm]), JacobiAmplitude[X[\[Lambda]], k], k] - S1[(g0 x[rm] - 1)/(g0 + x[rm]), \[Pi]/2+ArcTan[g0], k]);
  Pi1[\[Lambda]_] := -2/(C+D) (a2/g0 (1+g0^2))(S1[g0, JacobiAmplitude[X[\[Lambda]], k], k] - RedS1[g0, \[Pi]/2+ArcTan[g0], k]);
  Pi2[\[Lambda]_] := -2/(C+D) (a2/g0 (1+g0^2))^2 (S2[g0, JacobiAmplitude[X[\[Lambda]], k], k] - RedS2[g0, \[Pi]/2+ArcTan[g0], k]);
  Ip[\[Lambda]_] := g0/(a2 (1-g0 x[rp])) (\[Lambda]-Pip[\[Lambda]]);
  Iminus[\[Lambda]_] := g0/(a2 (1-g0 x[rm])) (\[Lambda]-Pim[\[Lambda]]);
  I1[\[Lambda]_] := (a2/g0-b1)\[Lambda] - Pi1[\[Lambda]];
  I2[\[Lambda]_] := (a2/g0-b1)^2 \[Lambda] + 4 (a2/g0-b1) 1/(C+D) a2/g0 (1+g0^2) S1[g0, JacobiAmplitude[X[\[Lambda]], k], k] + Pi2[\[Lambda]]+b1-(2 g0 C D)/(Sqrt[g0^2+1] Sqrt[(C-D)^2+g0^2 (C+D)^2]);
  It = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 I1[Global`\[Lambda]] + I2[Global`\[Lambda]] + 2 Log[2]]]]], Listable],
  (*else*)
  Pip[\[Lambda]_] := \[Nu]r (2/(C+D)) ((1+g0^2)/(g0 (g0+x[rp]))) (S1[(g0 x[rp] - 1)/(g0 + x[rp]), JacobiAmplitude[X[\[Lambda]], k], k] - S1[(g0 x[rp] - 1)/(g0 + x[rp]), ArcTan[x[rs]]+ArcTan[g0], k]);
  Pim[\[Lambda]_] := \[Nu]r (2/(C+D)) ((1+g0^2)/(g0 (g0+x[rm]))) (S1[(g0 x[rm] - 1)/(g0 + x[rm]), JacobiAmplitude[X[\[Lambda]], k], k] - S1[(g0 x[rm] - 1)/(g0 + x[rm]), ArcTan[x[rs]]+ArcTan[g0], k]);
  Pi1[\[Lambda]_] := \[Nu]r 2/(C+D) (a2/g0 (1+g0^2))(S1[g0, JacobiAmplitude[X[\[Lambda]], k], k] - S1[g0, ArcTan[x[rs]]+ArcTan[g0], k]);
  Pi2[\[Lambda]_] := \[Nu]r 2/(C+D) (a2/g0 (1+g0^2))^2 (S2[g0, JacobiAmplitude[X[\[Lambda]], k], k] - S2[g0, ArcTan[x[rs]]+ArcTan[g0], k]);
  Ip[\[Lambda]_] := g0/(a2 (1-g0 x[rp])) (\[Lambda]-Pip[\[Lambda]]);
  Iminus[\[Lambda]_] := g0/(a2 (1-g0 x[rm])) (\[Lambda]-Pim[\[Lambda]]);
  I1[\[Lambda]_] := (a2/g0-b1)\[Lambda] - Pi1[\[Lambda]];
  I2[\[Lambda]_] := (a2/g0-b1)^2 \[Lambda] - 2(a2/g0-b1) Pi1[\[Lambda]] + Pi2[\[Lambda]];
  It = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[4/(rp-rm) (rp(rp - (a \[ScriptL])/2) Ip[Global`\[Lambda]] - rm(rm - (a \[ScriptL])/2)Iminus[Global`\[Lambda]]) +4 Global`\[Lambda] + 2 I1[Global`\[Lambda]] + I2[Global`\[Lambda]]]]]], Listable]
];

I\[Phi] = Function[{Global`\[Lambda]}, Evaluate[If[Global`\[Lambda]>\[Lambda]x || Global`\[Lambda]<0, Undefined, Evaluate[Re[(2 a)/(rp-rm) ((rp - a \[ScriptL]/2) Ip[Global`\[Lambda]] - (rm- a \[ScriptL]/2) Iminus[Global`\[Lambda]])]]]], Listable];
 
<|"r" -> r, "I\[Phi]" -> I\[Phi], "It" -> It, "\[Lambda]x" -> Re[\[Lambda]x]|>
]


(* ::Section:: *)
(*Equator Intersections*)


(* ::Text:: *)
(*This is a convenience function for rendering images of thin disks in the equatorial plane of the Kerr BH. Thanks to the symmetries of the Elliptic integrals, one is able to find all the values of Mino time at which a null geodesic intersects the equatorial plane \[Theta]=\[Pi]/2. *)


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


Options[EmissionParameters] = {"PhiRange" -> {-\[Pi], \[Pi]}}

EmissionParameters[a_, \[Eta]_, \[ScriptL]_, \[Theta]o_, rem_, \[Theta]em_, OptionsPattern[]] := Module[{Z1, Z2, rms, A, B, \[CapitalDelta], \[Omega], \[CapitalOmega], utcg, R, \[CapitalTheta], ploct, plocr, ploc\[Theta], \[Kappa], \[Theta]loc, \[Phi]loc},
Z1=1+Surd[1-a^2,3] (Surd[1+a,3] + Surd[1-a,3]); Z2=Sqrt[3 a^2 + Z1^2]; rms=3+Z2-Sqrt[(3-Z1) (3+Z1+2 Z2)];
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
\[Phi]loc = Mod[Sign[-1+(\[CapitalDelta] A^2/(B^2 (\[Omega]-\[CapitalOmega]))+\[Omega]) \[ScriptL]] ArcCos[plocr/Sqrt[ploct^2-ploc\[Theta]^2]], OptionValue["PhiRange"][[2]]-OptionValue["PhiRange"][[1]], OptionValue["PhiRange"][[1]]];

<|"\[Kappa]" -> \[Kappa], "\[Theta]loc" -> \[Theta]loc, "\[Phi]loc" -> \[Phi]loc|>
 ]]


MomentumFromParameters[rem_, a_ \[Omega]loc_, \[Theta]loc_, \[Phi]loc_, M_:1] := Module[{ploct, plocr, ploc\[Theta], ploc\[Phi], A, B, \[CapitalDelta], \[Omega], \[CapitalOmega], utcg, et\[Phi], pt, pr, p\[Theta], p\[Phi], vecpt, vecp\[Phi]},
ploct = \[Omega]_loc Quantity["ReducedPlanckConstant"] Quantity["GravitationalConstant"]  / (Quantity["SpeedOfLight"])^3/ M;
plocr = ploct Sin[\[Theta]loc] Cos[\[Phi]loc];
ploc\[Phi] = ploct Sin[\[Theta]loc] Sin[\[Phi]loc];
ploc\[Theta] = ploct Cos[\[Theta]loc];

A = rem^2; \[CapitalDelta] = rem^2-2 rem + a^2; B = (rem^2 + a^2)^2 - \[CapitalDelta] a^2; \[Omega] = (2 a rem)/B;
\[CapitalOmega] =1/(Sqrt[rem^3]+ a); utcg=((\[CapitalDelta] A)/B - (\[Omega]-\[CapitalOmega])^2 B/A)^(-1/2);
et\[Phi] = ((\[CapitalDelta]^2 A^3)/(B^3 (\[Omega]-\[CapitalOmega])) - (\[CapitalDelta] A)/B)^(-1/2);

vecpt = ploct utcg + ploc\[Phi] et\[Phi];
pr = plocr Sqrt[\[CapitalDelta]/A];
p\[Theta] = ploc\[Theta] Sqrt[1/A];
vecp\[Phi] = ploct utcg \[CapitalOmega] + ploc\[Phi] et\[Phi] ((\[CapitalDelta] A^2)/(B^2 (\[Omega]-\[CapitalOmega]))+\[Omega]);

pt = ((-\[CapitalDelta] A)/B + (\[Omega]^2 B)/A) vecpt - (\[Omega] B)/A vecp\[Phi];
p\[Phi] = - ((\[Omega] B)/A) vecpt + B/A vecp\[Phi];

{pt, pr, p\[Theta], p\[Phi]}  
]


(* ::Section:: *)
(*Public Functions*)


Options[KerrNullGeo] = {"Momentum" -> "Momentum", "PhiRange" -> {-\[Infinity], \[Infinity]}}
SyntaxInformation[KerrNullGeoDistant] = {"ArgumentsPattern"->{_,_,_}};

KerrNullGeo[a_, xs_, ps_, M_:1, OptionsPattern[]] := Module[{ts, rs, \[Theta]s, \[Phi]s, pts, prs, p\[Theta]s, p\[Phi]s, consts, \[ScriptL], \[Eta], roots, r1, r2, r3, r4, type, r, I\[Phi], \[Lambda]x, It, \[Theta], G\[Phi], Gt, equator\[Lambda], \[Phi], t, \[Kappa], \[Theta]loc, \[Phi]loc, \[Theta]x, \[Phi]x, assoc},
If[a<=0 || a>=1, Message[KerrNullGeo::OutOfBounds, "Parameter a must be between 0 and 1."]; Return[];];
If[OptionValue["Momentum"]=="WaveVector", ps = ps Quantity["ReducedPlanckConstant"] Quantity["GravitationalConstant"]  / (Quantity["SpeedOfLight"])^3/ M];

If[ Length[xs]==Length[ps] == 4,
  {{ts, rs, \[Theta]s, \[Phi]s}, {pts, prs, p\[Theta]s, p\[Phi]s}} = {xs, ps};
  If[\[Theta]s<0 || \[Theta]s>\[Pi], Message[KerrNullGeo::OutOfBounds, "Parameter \[Theta]o must be between 0 and \[Pi]."]; Return[]],
  If[ Length[xs]==Length[ps] == 3,
    {ts, rs, \[Phi]s} = xs; \[Theta]s = \[Pi]/2;
    {pts, prs, p\[Theta]s, p\[Phi]s} = MomentumFromParameters[rs, a, ps[0], ps[1], ps[2], M],
    Message[KerrNullGeo::ListSize, xs, ps, "3 or 4"]; Return[]
  ]
];

consts = NullConstantsOfMotion[a, \[Theta]s, pts, p\[Theta]s, p\[Phi]s];
{\[ScriptL], \[Eta]} = {"\[ScriptL]", "\[Eta]"} /. consts;
If[Abs[\[Eta]]< 10 $MachineEpsilon || Abs[\[Eta]+(Abs[\[ScriptL]]-a)^2]< 10 $MachineEpsilon, \[Eta]=\[Eta] + 10 $MachineEpsilon]; (*To avoid undefined expresions in polar motion*)

roots = RadialRoots[a, \[Eta], \[ScriptL]];
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;

type = "PhotonCapture";
If[Im[r2] != 0,
  {r, I\[Phi], \[Lambda]x, It} = {"r", "I\[Phi]", "\[Lambda]x", "It"} /. RadialMotionCase4[roots, a, \[Eta], \[ScriptL], rs, Sign[prs]],
  If[Im[r4] != 0, 
    {r, I\[Phi], \[Lambda]x, It} = {"r", "I\[Phi]", "\[Lambda]x", "It"} /. RadialMotionCase3[roots, a, \[Eta], \[ScriptL], rs, Sign[prs]],
    If[r4<1+Sqrt[1-a^2],
      {r, I\[Phi], \[Lambda]x, It} = {"r", "I\[Phi]", "\[Lambda]x", "It"} /. RadialMotionCase2[roots, a, \[Eta], \[ScriptL], rs, Sign[prs]];
      If[r3>1-Sqrt[1-a^2], type = "PhotonEscape"],
      {r, I\[Phi], \[Lambda]x, It} = {"r", "I\[Phi]", "\[Lambda]x", "It"} /. RadialMotionCase1[roots, a, \[Eta], \[ScriptL], rs, Sign[prs]];
      ];
  ]
];


If[\[Eta]>0, 
  {\[Theta], G\[Phi], Gt} = {"\[Theta]", "G\[Phi]", "Gt"} /. OrdinaryPolarMotion[a, \[Eta], \[ScriptL], \[Theta]s, Sign[p\[Theta]s], \[Lambda]x];
  equator\[Lambda] = EquatorIntersectionMinoTimes[a, \[Eta], \[ScriptL], \[Theta]s, Sign[p\[Theta]s], \[Lambda]x],
  {\[Theta], G\[Phi], Gt} = {"\[Theta]", "G\[Phi]", "Gt"} /. VorticalPolarMotion[a, \[Eta], \[ScriptL], \[Theta]s, Sign[p\[Theta]s], \[Lambda]x]; 
  equator\[Lambda] = {};
];


If[OptionValue["PhiRange"][[2]]==\[Infinity], 
  \[Phi]=Function[{Global`\[Lambda]}, Evaluate[I\[Phi][Global`\[Lambda]]+\[ScriptL] G\[Phi][Global`\[Lambda]] +\[Phi]s], Listable],
  \[Phi]=Function[{Global`\[Lambda]}, Evaluate[Mod[I\[Phi][Global`\[Lambda]]+\[ScriptL] G\[Phi][Global`\[Lambda]] + \[Phi]s, OptionValue["PhiRange"][[2]]-OptionValue["PhiRange"][[1]], OptionValue["PhiRange"][[1]]]], Listable]
];

t = Function[{Global`\[Lambda]}, Evaluate[It[Global`\[Lambda]] + a^2 Gt[Global`\[Lambda]] + ts]];
If[type == "PhotonEscape", \[Theta]x=\[Theta][\[Lambda]x]; \[Phi]x=\[Phi][\[Lambda]x], \[Theta]x=-1; \[Phi]x=-1];

If[Length[equator\[Lambda]] == 0, 
  {\[Kappa], \[Theta]loc, \[Phi]loc} = {-1, -1, -1}, 
  {\[Kappa], \[Theta]loc, \[Phi]loc} = {"\[Kappa]", "\[Theta]loc", "\[Phi]loc"} /. EmissionParameters[a, \[Eta], \[ScriptL], \[Theta]o, r[equator\[Lambda][[1]]], \[Theta][equator\[Lambda][[1]]], "PhiRange" -> If[OptionValue["PhiRange"][[2]]!=\[Infinity], OptionValue["PhiRange"]]]
];

assoc = <|
	"Trajectory" -> {t, r, \[Theta], \[Phi]},
	"ConstantsOfMotion" -> consts,
	"RadialRoots" -> {r1, r2, r3, r4},
	"EquatorIntersectionMinoTimes" -> equator\[Lambda],
	"EquatorIntersectionCoordinates" -> {t[equator\[Lambda]], r[equator\[Lambda]], \[Phi][equator\[Lambda]]},
	"TrajectoryType" -> type,
	"MinoTimeOfCapture" -> \[Lambda]x,
	"EscapeCoordinates" -> {\[Theta]x, \[Phi]x},
	"EmissionParameters" -> {\[Kappa], \[Theta]loc, \[Phi]loc}
|>;

KerrNullGeoFunction[a, xs, ps, M, assoc]
]


KerrNullGeoFunction[a_, xs_, ps_, M_, assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]];
KerrNullGeoFunction[a_, xs_, ps_, M_, assoc_][y_?StringQ] := assoc[y];
Keys[g_KerrNullGeoFunction]^:=Keys[g[[5]]];


Options[KerrNullGeoDistant] = {"Rotation" -> "Counterclockwise", "PhiRange" -> {-\[Infinity], \[Infinity]}}
SyntaxInformation[KerrNullGeoDistant] = {"ArgumentsPattern"->{_,_,_,_}};

KerrNullGeoDistant[a_, \[Theta]o_, \[Alpha]_, \[Beta]_, OptionsPattern[]] := Module[ {consts, \[Eta], \[ScriptL], roots, r1, r2, r3, r4, rp, rm, k, r, \[Theta], \[Phi], G\[Phi], I\[Phi], \[CapitalDelta]v, Gt, RedIt, \[Lambda]x, \[Phi]x, \[Theta]x, assoc, equator\[Lambda], type, \[Kappa], \[Theta]loc, \[Phi]loc},

If[a<=0 || a>=1, Message[KerrNullGeo::OutOfBounds, "Parameter a must be between 0 and 1."]; Return[];];
If[\[Theta]o<0 || \[Theta]o>\[Pi], Message[KerrNullGeo::OutOfBounds, "Parameter \[Theta]o must be between 0 and \[Pi]."]; Return[];];

If[OptionValue["Rotation"]=="Clockwise", \[Alpha]=-\[Alpha]];

consts = DistantNullConstantsOfMotion[a, \[Theta]o, \[Alpha], \[Beta]];
{\[ScriptL], \[Eta]} = {"\[ScriptL]", "\[Eta]"} /. consts;
If[Abs[\[Eta]]< 10 $MachineEpsilon || Abs[\[Eta]+(Abs[\[ScriptL]]-a)^2]< 10 $MachineEpsilon, \[Eta]=\[Eta] + 10 $MachineEpsilon]; (*To avoid undefined expresions in polar motion*)

roots = RadialRoots[a, \[Eta], \[ScriptL]];
{r1, r2, r3, r4} = {"r1", "r2", "r3", "r4"} /. roots;

type = "PhotonCapture";
If[Im[r2] != 0,
  {r, I\[Phi], \[Lambda]x, RedIt} = {"r", "I\[Phi]", "\[Lambda]x", "It"} /. RadialMotionCase4[roots, a, \[Eta], \[ScriptL], "Observer" -> "Distant"],
  If[Im[r4] != 0, 
    {r, I\[Phi], \[Lambda]x, RedIt} = {"r", "I\[Phi]", "\[Lambda]x", "It"} /. RadialMotionCase3[roots, a, \[Eta], \[ScriptL], "Observer" -> "Distant"],
    {r, I\[Phi], \[Lambda]x, RedIt} = {"r", "I\[Phi]", "\[Lambda]x", "It"} /. RadialMotionCase2[roots, a, \[Eta], \[ScriptL], "Observer" -> "Distant"];
    If[r3>1-Sqrt[1-a^2], type = "PhotonEscape"]
  ]
];

ssign[x_] := If[x<0, -1, 1];

If[\[Eta]>0, 
  {\[Theta], G\[Phi], Gt} = {"\[Theta]", "G\[Phi]", "Gt"} /. OrdinaryPolarMotion[a, \[Eta], \[ScriptL], \[Theta]o, -ssign[\[Beta]], \[Lambda]x];
  equator\[Lambda] = EquatorIntersectionMinoTimes[a, \[Eta], \[ScriptL], \[Theta]o, -ssign[\[Beta]], \[Lambda]x],
  {\[Theta], G\[Phi], Gt} = {"\[Theta]", "G\[Phi]", "Gt"} /. VorticalPolarMotion[a, \[Eta], \[ScriptL], \[Theta]o, -ssign[\[Beta]], \[Lambda]x]; (*Minus Sign[\[Beta]] because we use negative Mino time*)
  equator\[Lambda] = {}; 
];


If[OptionValue["PhiRange"][[2]]==\[Infinity], 
  \[Phi]=Function[{Global`\[Lambda]}, Evaluate[I\[Phi][Global`\[Lambda]]+\[ScriptL] G\[Phi][Global`\[Lambda]]], Listable],
  \[Phi]=Function[{Global`\[Lambda]}, Evaluate[Mod[I\[Phi][Global`\[Lambda]]+\[ScriptL] G\[Phi][Global`\[Lambda]], OptionValue["PhiRange"][[2]]-OptionValue["PhiRange"][[1]], OptionValue["PhiRange"][[1]]]], Listable]
];
\[CapitalDelta]v=Function[{Global`\[Lambda]}, Evaluate[RedIt[Global`\[Lambda]]+a^2 Gt[Global`\[Lambda]] + r[Global`\[Lambda]] + 2 Log[r[Global`\[Lambda]]/2]], Listable];

If[type == "PhotonEscape", \[Theta]x=\[Theta][\[Lambda]x]; \[Phi]x=\[Phi][\[Lambda]x], \[Theta]x=-1; \[Phi]x=-1];

If[Length[equator\[Lambda]] == 0, 
  {\[Kappa], \[Theta]loc, \[Phi]loc} = {-1, -1, -1}, 
  {\[Kappa], \[Theta]loc, \[Phi]loc} = {"\[Kappa]", "\[Theta]loc", "\[Phi]loc"} /. EmissionParameters[a, \[Eta], \[ScriptL], \[Theta]o, r[equator\[Lambda][[1]]], \[Theta][equator\[Lambda][[1]]], "PhiRange" -> If[OptionValue["PhiRange"][[2]]!=\[Infinity], OptionValue["PhiRange"], {-\[Pi], \[Pi]}]]
];

assoc = <|
	"Trajectory" -> {\[CapitalDelta]v, r, \[Theta], \[Phi]},
	"ConstantsOfMotion" -> consts,
	"RadialRoots" -> {r1, r2, r3, r4},
	"EquatorIntersectionMinoTimes" -> equator\[Lambda],
	"EquatorIntersectionCoordinates" -> {\[CapitalDelta]v[equator\[Lambda]], r[equator\[Lambda]], \[Phi][equator\[Lambda]]},
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
