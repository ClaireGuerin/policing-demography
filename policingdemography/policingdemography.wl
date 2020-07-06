(* ::Package:: *)
(* Created by Claire Guerin June 30, 2020*)

BeginPackage["policingdemography`"]

(* Exported symbols added here with SymbolName::usage *) 
search::usage="search[possible parameter values to be explored] needs to be fed {th, rb, c1, c2, b, d, p, \[Eta], \[Gamma], m} in that order";

Begin["`Private`"]
(* Implementation of the package *)

(* FUNCTIONS AND EXPRESSIONS *)
resources[x_,xn_,n_] := Block[{baseres,cost,benef,police,tot},
	(*VALIDATED*)
	(*pars = {th -> 0.1, rb -> 20, c1 -> 0, c2 -> 0.1, b -> 2, d -> 0.3, p -> 0.5, eta -> 1, gamma -> 1, m->0.5};*)
	baseres = rb / n;
	cost = x * c1 + x * x * c2;
	benef = b * ((x + (n - 1) * xn) * (1 - p))^gamma;
	police = d * (1 - x) * (p * (x + (n - 1) * xn) / n)^eta;
	tot = baseres * (1 - cost + benef / n - police);
	Return[tot]]

fitness[x_,xn_,n_] := Block[{r,w},
	(*VALIDATED*)
	(*pars = {th -> 0.1, rb -> 20, c1 -> 0, c2 -> 0.1, b -> 2, d -> 0.3, p -> 0.5, eta -> 1, gamma -> 1, m->0.5};*)
	r = resources[x,xn,n];
	w = r / (1 + th * r);
	Return[w]]

dynamics[xm_,n_,nr_] := (1 - m) * fitness[xm,xm,n] * n + m * nr

(* Monomorphic population & Selection Gradient*)
monopopulation = {x -> y, xn -> y, xm -> y, nr -> n, r -> ((1 - m)^2)/(1 + (1 - (1 - m)^2) * (n - 1)), rbar -> 1/n + (1 - 1/n) * r}
gradient = ((1 - m)*(m - 1)*n^4*rb^2*(-(b*gamma*((-n)*(p - 1)*y)^gamma) + c1*n*y + 2*c2*n*y^2 + d*eta*n*(p*y)^eta - d*n*y*(p*y)^eta - d*eta*n*y*(p*y)^eta)*(b*(gamma - 2)*((-n)*(p - 1)*y)^gamma + n*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1)))/(y*(m^2*(-(n - 1)) + 2*m*(n - 1) + 1)*(b*rb*th*((-n)*(p - 1)*y)^gamma - n*rb*th*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1) + n^2)^4*(1 - ((m - 1)^2*rb*(b^2*rb*th*((-n)*(p - 1)*y)^(2*gamma) + n^2*(b*(gamma - 1)*((-n)*(p - 1)*y)^gamma + rb*th*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1)^2) - 2*b*n*rb*th*((-n)*(p - 1)*y)^gamma*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1)))/(b*rb*th*((-n)*(p - 1)*y)^gamma - n*rb*th*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1) + n^2)^2)) - (n^2*rb*(b*gamma*((-n)*(p - 1)*y)^gamma + c1*n*y*(m^2*(n - 1) - 2*m*(n - 1) - 1) + 2*c2*n*y^2*(m^2*(n - 1) - 2*m*(n - 1) - 1) - d*m^2*n^2*y*(p*y)^eta + d*m^2*n*y*(p*y)^eta + 2*d*m*n^2*y*(p*y)^eta - 2*d*m*n*y*(p*y)^eta - d*eta*n*(p*y)^eta + d*n*y*(p*y)^eta + d*eta*n*y*(p*y)^eta))/(y*(m^2*(n - 1) - 2*m*(n - 1) - 1)*(b*rb*th*((-n)*(p - 1)*y)^gamma - n*rb*th*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1) + n^2)^2);
(*VALIDATED*)

(* PARAMETER SEARCH *)

numerics[parSet_] := Block[{ecoequi,nZero,nOne,nTrend},

	(**********************************************************************************)
	
	(* CALCULATE THE ECOLOGICAL EQUILIBRIUM WITH GIVEN PARAMETER SET *)
	ecoequi = DeleteCases[Flatten[NSolve[(n == dynamics[xm,n,nr]) /. monopopulation /. parSet, n]], n -> 0.`];
	(** ->[nZero] Get the size of the ecological equilibrium when y=0 **)
	nZero = n /. ecoequi /. y->0;
	(** ->[nOne] Get the size of the ecological equilibrium when y=1 **)
	nOne = n /. ecoequi /. y->1;
	(** ->[nTrend] Check whether n equilibrium increases with y (get the sign of the derivative for y in ]0,1[) **)
	nTrend = Refine[Sign[D[n /. ecoequi, y]], 0 < y < 1];
	(**********************************************************************************)
	
	(* CALCULATE SELECTION GRADIENT AT ECOLOGICAL EQUILIBRIUM FOR GIVEN PARAMETER SET *)
	s = gradient /. parSet;
	(** ->[sEqui] Get the selection gradient sign for a range of y between 0 and 1 (included)  **)
	gradientZeroOut = s /. ecoequi /. y -> # & /@ Range[0.1, 1, 0.1];
	gradientZeroIn = s /. ecoequi /. y -> 10^(-10);
	(* errorTest = Check[s /. {n->nOne}, StringForm["Error for n*(y=1)=``", nOne]];  *)
	sEqui = Prepend[gradientZeroOut,gradientZeroIn];
	(** ->[sZero] Get the sign / value of the selection gradient for y=0 and n=neq(1) **)
	sZero = s /. {n -> nOne, y -> 10^(-10)};
	(* sZero = Sign[s /. {n -> nOne, y -> 10^(-10)}] *)
	(**********************************************************************************)

	Return[<|"nZero"->nZero,"nOne"->nOne,"nTrend"->nTrend,"sEqui"->sEqui,"sZero"->sZero|>]]
	(* Return[errorTest]] *)

search[values_List] := Block[{names, combinationsList, combinationsRule,results},
	(*Assign the pars values to the private set of parameter rules*)
	names = {th, rb, c1, c2, b, d, p, eta, gamma, m};
	If[Length[values]!=Length[names],Message[search::len, Length[values]],
	combinationsList = Tuples[values];
	combinationsRule = Map[MapThread[Rule, {names, combinationsList[[#]]}] &, Range[Length[combinationsList]]];

	(*Calculate all the things we need*)
	results = numerics[#] &/@ combinationsRule;
	Return[<|"comb"->combinationsList, "res"->results|>]]]

search::len = "You gave `1` parameters, this model needs exactly 10";

End[]

EndPackage[]

(* 
selectionGradient[] := Block[{dn,se,sw,s},
	(*!!!! WRONG !!!!*)
	(*pars = {th -> 0.1, rb -> 20, c1 -> 0, c2 -> 0.1, b -> 2, d -> 0.3, p -> 0.5, eta -> 1, gamma -> 1, m->0.5};*)
	dn = D[dynamics[xm,n,nr], xm] * (1 - (1 - m) * D[dynamics[xm,n,nr], n])^(-1) * (1 - m) * rbar //. monopopulation // Simplify;
	sw = D[fitness[x,xn,x], x] + D[fitness[x, xn, n], xn] * r //. monopopulation // Simplify;
	se = dn * D[fitness[x,xn,n], n] //. monopopulation // Simplify;
	s = sw + se;
	Return[s]] *)