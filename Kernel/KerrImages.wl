(* ::Package:: *)

(* ::Title:: *)
(*KerrImages Package*)


(* ::Section:: *)
(*Define Usage for Public Functions*)


BeginPackage["BlackHoleImages`KerrImages`",
	{"BlackHoleImages`KerrNullGeodesics`", "BlackHoleImages`AlphaDiskModel`"}
]

GenerateTemplate::usage = "GenerateTemplate[directory, name, a, \[Theta]o, imageSize, maxBardeenCoordinate, Options] generates a template containing information about null geodesics specified by (dimensionless) a, \[Theta]o (distant observer's \[Theta]) and maximal Bardeen coordinate that is shown. The number of geodesics to be generated is specified in imageSize in the form {xsize, ysize}. The template is saved in directory/name.mx file, which can be used later with DiskImageFromTemplate.";
DiskImage::usage = "DiskImage[a, \[Theta]o, \[Alpha], m, mdot, imageSize, maxBardeenCoordinate, Options] returns a list containing tables with information specified by the option 'Output'->{'information1', 'information2', ...}. Geodesics are specified by (dimensionless) a, \[Theta]o (distant observer's \[Theta]) and maximal Bardeen coordinate that is shown. The number of geodesics to be generated is specified in imageSize in the form {xsize, ysize}. Disk is specified by BH mass m (in solar mass by default), matter inflow mdot and parameter \[Alpha] (both in Shakura & Sunyaev definition by default).";
DiskImageFromTemplate::usage = "DiskImageFromTemplate[file, \[Alpha], m, mdot, Options] returns a list containing tables with information specified by the option 'Output'->{'information1', 'information2', ...}. Template is passed in file. Disk is specified by BH mass m (in solar mass by default), matter inflow mdot and parameter \[Alpha] (both in Shakura & Sunyaev definition by default).";
StellarBackgroundFromTemplate::usage = "StellarBackgroundFromTemplate[templateFile, \[Theta]o, imageFile, angle, Ratio_:0.8, bgColor_:{0.,0.,0.}] generates an image of stellar background given by imageFile distorted by geometry given by the template stored in templateFile and \[Theta]o. The image's part Ratio (default is 0.8) spans angle on the celestial sphere. The background color can be specified in bgColor as RGB list of size 3 (default is black)."

Begin["`Private`"];


(* ::Section:: *)
(*Public Functions*)


Options[GenerateTemplate] = {"Rotation" -> "Counterclockwise", "PhiRange" -> {-\[Pi], \[Pi]} (*BlackHoleImages`KerrNullGeodesics`KerrNullGeoDistant options*)}

GenerateTemplate[directory_, name_, a_, \[Theta]o_, imageSize_, maxBardeenCoordinate_, OptionsPattern[]] := Module[{x, y, i, j, imax, jmax, step, geod, row, template, file},
(*set the dimensions*)
If[Length[imageSize]==2,
x = imageSize[[1]]; y = imageSize[[2]]; imax = maxBardeenCoordinate; jmax = imax y/x,
Message[imageSize, "is not a list of size 2."]
];

(*how much two adjacent pixels differ*)
step = 2 imax/x;

(*temporary containers*)
row={};
template={};

(*show progress*)
j=-jmax;
Print[Dynamic[N[100 (j+jmax)/(2 jmax)]], " %"];

(*generate the template*)
For[j=-jmax, j<jmax,j+=step,
	For[i=-imax, i<imax ,i+=step,
		geod = BlackHoleImages`KerrNullGeodesics`KerrNullGeoDistant[a, \[Theta]o, i, j, "Rotation" -> OptionValue["Rotation"],"PhiRange" -> OptionValue["PhiRange"]];
		AppendTo[row, {geod["EmissionCoordinates"], geod["EmissionParameters"], geod["EscapeCoordinates"]}];
	];
	AppendTo[template, row];
	row = {};
];

(*export the template to given destination*)
file = FileNameJoin[{directory, name<>".mx"}];
Export[file, template];


]


Options[DiskImage] = {"InputUnits" -> "NovikovThorne", "OutputUnits" -> "SI", "rUnits" -> "BHMass",  (*BlackHoleImages`AlphaDiskModel`DiskParams options*)
					"Grid"->True, (*BlackHoleImages`AlphaDiskModel`ObservedDiskElement option*)
					"Rotation" -> "Counterclockwise", "PhiRange" -> {-\[Pi], \[Pi]}, (*BlackHoleImages`KerrNullGeodesics`KerrNullGeoDistant options*)
					"Output" -> {"MaximalFrequency"} (*Specifies what information is requested in the form {"information1", "information 2",...},
													where "informationN" must be element of the association returned by BlackHoleImages`AlphaDiskModel`ObservedDiskElement
													("PhysicalTemperature", "EffectiveTemperature", "SpectralFluxDensity", "FluxDensity" or "PeakFrequency") *)
}		

DiskImage[a_, \[Theta]o_, \[Alpha]_, m_, mdot_, imageSize_, maxBardeenCoordinate_, OptionsPattern[]] := Module[{x, y, i, j, imax, jmax, step, disk, length, row, matrix, geod, element},
(*set the dimensions*)
If[Length[imageSize]==2,
x = imageSize[[1]]; y = imageSize[[2]]; imax = maxBardeenCoordinate; jmax = imax y/x,
Message[imageSize, "is not a list of size 2."]
];

(*how much two adjacent pixels differ*)
step = 2 imax/x;

(*how many tables user wants to generate = length of the "Output" option*)
length = Length[OptionValue["Output"]];
matrix = Table[{}, length];
row = Table[{}, length];

(*generate parameters of the disk*)
disk = BlackHoleImages`AlphaDiskModel`DiskParams[a, \[Alpha], m, mdot, "InputUnits" -> OptionValue["InputUnits"], "OutputUnits" -> OptionValue["OutputUnits"], "rUnits" -> OptionValue["rUnits"]];

(*generate the requested data*)
For[j=-jmax, j<jmax,j+=step,
	For[i=-imax, i<imax ,i+=step,
		geod = BlackHoleImages`KerrNullGeodesics`KerrNullGeoDistant[a, \[Theta]o, i, j, "Rotation" -> OptionValue["Rotation"],"PhiRange" -> OptionValue["PhiRange"]];
		element = BlackHoleImages`AlphaDiskModel`ObservedDiskElement[disk, geod, "Grid"->OptionValue["Grid"]];
		row = MapThread[Append, {row, Table[element[OptionValue["Output"][[k]]], {k, 1, length}]}];
	];
	matrix = MapThread[Append, {matrix, row}];
	row = Table[{}, length];
];

(*return a list containing tables with requested output*)
matrix
]


Options[DiskImageFromTemplate] = {"InputUnits" -> "NovikovThorne", "OutputUnits" -> "SI", "rUnits" -> "BHMass",  (*BlackHoleImages`AlphaDiskModel`DiskParams options*)
								  "Grid"->True, (*BlackHoleImages`AlphaDiskModel`ObservedDiskElement option*)
								  "Output" -> {"MaximalFrequency"} (*Specifies what information is requested in the form {"information1", "information 2",...},
																	where "informationN" must be element of the association returned by BlackHoleImages`AlphaDiskModel`ObservedDiskElement
																	("PhysicalTemperature", "EffectiveTemperature", "SpectralFluxDensity", "FluxDensity" or "PeakFrequency") *)
}		

DiskImageFromTemplate[file_, a_, \[Alpha]_, m_, mdot_, OptionsPattern[]] := Module[{template, i, j, imax, jmax, disk, length, matrix, row, shard, element},
(*import the template*)
template = Import[file];

(*get the dimensions*)
jmax = Length[template];
imax = Length[template[[1]]];

(*how many tables user wants to generate = length of the "Output" option*)
length = Length[OptionValue["Output"]];
matrix = Table[{}, length];
row = Table[{}, length];

(*generate parameters of the disk*)
disk = BlackHoleImages`AlphaDiskModel`DiskParams[a, \[Alpha], m, mdot, "InputUnits" -> OptionValue["InputUnits"], "OutputUnits" -> OptionValue["OutputUnits"], "rUnits" -> OptionValue["rUnits"]];

(*generate the requested data*)
For[j=1, j<=jmax,j+=1,
	For[i=1, i<=imax ,i+=1,
		shard = template[[j, i]];
		element = BlackHoleImages`AlphaDiskModel`ObservedDiskElement[disk, <|"EmissionCoordinates"-> shard[[1]], "EmissionParameters" -> shard[[2]]|>, "Grid"->OptionValue["Grid"]];
		row = MapThread[Append, {row, Table[element[OptionValue["Output"][[k]]], {k, 1, length}]}];
	];
	matrix = MapThread[Append, {matrix, row}];
	row = Table[{}, length];
];

(*return a list containing tables with requested output*)
matrix
]


StellarBackgroundFromTemplate[templateFile_, \[Theta]o_,imageFile_, angle_, Ratio_:0.8, bgColor_:{0.,0.,0.}] := Module[{template, image, data, i, j, imax, jmax, Jmax, Imax,ratio, disk, length, matrix, row, shard, element, \[Theta], \[Phi], \[CapitalDelta]\[Phi], \[CapitalDelta]\[Theta], append},
(*import the template and image*)
template = Import[templateFile];
image = Import[imageFile];
data = ImageData[image];

(*get the dimensions*)
jmax = Length[template];
imax = Length[template[[1]]];
{Imax, Jmax} = ImageDimensions[image];

(*how many pixels of the original image correspond to a pixel of the template*)
ratio = Ratio Min[Imax, Jmax]/({imax, jmax}[[Position[{Imax, Jmax}, Min[Imax, Jmax]][[1]][[1]]]]);


matrix = {};
row = {};

(*generate the requested data*)
For[j=1, j<=jmax,j+=1,
	For[i=1, i<=imax ,i+=1,
		shard = template[[j, i]];
		\[Theta] = shard[[3, 1]];
		\[Phi] = shard[[3, 2]];
		\[CapitalDelta]\[Theta] = \[Theta]-(\[Pi]-\[Theta]o);
		\[CapitalDelta]\[Phi] = Mod[\[Phi] , 2\[Pi]]-\[Pi];

		(*if the geodesic "points" outside of the original image, the background color is used instead*)
		If[Abs[Round[ratio (2j -jmax)/2 + Jmax \[CapitalDelta]\[Phi]/angle]] > Jmax/2 - 1 || Abs[Round[ratio (2i -imax)/2 + Imax \[CapitalDelta]\[Theta]/angle]]>Imax/2 - 1,
			append = bgColor,
			(*else*)
			append = data[[Round[Jmax/2+ratio (2j -jmax)/2 + Jmax \[CapitalDelta]\[Phi]/angle],
			Round[Imax/2+ratio (2i -imax)/2 + Imax \[CapitalDelta]\[Theta]/angle]]]
		];
		
		If[\[Theta]==-1, append={0,0,0}];
		
		row = Append[row, append];
	];
	
	matrix = Append[matrix, row];
	row = {};
];

(*return a list containing tables with requested output*)
matrix
]


(* ::Section:: *)
(*Close the Package*)


End[];

EndPackage[];
