(* Wolfram Language Package *)
(* Created by Claire Guerin June 30, 2020*)

BeginPackage["policing-demography`"]
(*Needs[ "chiefstitution`plots`"]*)
(*Needs[ "chiefstitution`onegeneration`"]*)
(* Exported symbols added here with SymbolName::usage *) 

gradient::usage="gradient[]"

Begin["`Private`"]
(* Implementation of the package *)
resources[x_,xn_,n_] := Block[{baseres,cost,benef,police,tot},
	baseres = rb / n;
	cost = x * c1 + x * x * c2;
	benef = b * ((x + (n - 1) * xn) * (1 - p))^gamma;
	police = d * (1 - x) * (p * (x + (n - 1) * xn) / n)^eta;
	tot = baseres * (1 - cost + benef / n - police);
	Return[tot]]

fitness[x_,xn_,n_] := Block[{r,w},
	r = resources[x,xn,n];
	w = r / (1 + th * r);
	Return[w]]

dynamics[xm_,n_,nr_] := (1 - m) * fitness[xm,xm,n] * n + m * nr

monopopulation[] := Block[{names, values, monopop},
	names = {x,xn,xm,nr,r,rbar};
	values = {y,y,y,n,((1 - m)^2) / (1 + (1 - (1 - m)^2) * (n - 1)),(1 + r * (n - 1)) / n};
	monopop = Map[Rule[names[[#]], values[[#]]] &, Range[Length[names]]];
	Return[monopop]]

gradient[] := Block[{monopop,dn,se,sw,s},
	monopop = monopopulation[];
	dn = D[dynamics[xm,n,nr], xm] * (1 - (1 - m) * D[dynamics[xm,n,nr], n])^(-1) * (1 - m) * rbar //. monopop // Simplify;
	sw = D[fitness[x,xn,x], x] + D[w[x, xn, n], xn] * r //. monopop // Simplify;
	se = dn * D[w[x,xn,n], n] //. monopop // Simplify;
	s = sw + se;
	Return[s]]

End[]

EndPackage[]
