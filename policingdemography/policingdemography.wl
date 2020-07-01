(* ::Package:: *)
(* Created by Claire Guerin June 30, 2020*)

BeginPackage["policingdemography`"]

(* Exported symbols added here with SymbolName::usage *) 
resources::usage="";
fitness::usage="";
dynamics::usage="";
search::usage="";

Begin["`Private`"]
(* Implementation of the package *)

(*doTheCoolStuff[s_,f_,pars_] := Block[{sNum,sSign},
	ecologicalEq = DeleteCases[Flatten[NSolve[(n == f) /. pars, n]], n -> 0.`];
	sNum = s /. pars /. m -> 0.5 /. ecologialEq;
	sZero = sNum /. y->0;
	sSign = Refine[Reduce[sNum > 0], 0 < y < 1];
	Return[ecologicalEq]]*)
(*FUNCTIONS AND EXPRESSIONS*)
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

monopopulation = {x -> y, xn -> y, xm -> y, nr -> n, r -> ((1 - m)^2)/(1 + (1 - (1 - m)^2) * (n - 1)), rbar -> 1/n + (1 - 1/n) * r}
gradient = ((1 - m)*(m - 1)*n^4*rb^2*(-(b*gamma*((-n)*(p - 1)*y)^gamma) + c1*n*y + 2*c2*n*y^2 + d*eta*n*(p*y)^eta - d*n*y*(p*y)^eta - d*eta*n*y*(p*y)^eta)*(b*(gamma - 2)*((-n)*(p - 1)*y)^gamma + n*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1)))/(y*(m^2*(-(n - 1)) + 2*m*(n - 1) + 1)*(b*rb*th*((-n)*(p - 1)*y)^gamma - n*rb*th*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1) + n^2)^4*(1 - ((m - 1)^2*rb*(b^2*rb*th*((-n)*(p - 1)*y)^(2*gamma) + n^2*(b*(gamma - 1)*((-n)*(p - 1)*y)^gamma + rb*th*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1)^2) - 2*b*n*rb*th*((-n)*(p - 1)*y)^gamma*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1)))/(b*rb*th*((-n)*(p - 1)*y)^gamma - n*rb*th*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1) + n^2)^2)) - (n^2*rb*(b*gamma*((-n)*(p - 1)*y)^gamma + c1*n*y*(m^2*(n - 1) - 2*m*(n - 1) - 1) + 2*c2*n*y^2*(m^2*(n - 1) - 2*m*(n - 1) - 1) - d*m^2*n^2*y*(p*y)^eta + d*m^2*n*y*(p*y)^eta + 2*d*m*n^2*y*(p*y)^eta - 2*d*m*n*y*(p*y)^eta - d*eta*n*(p*y)^eta + d*n*y*(p*y)^eta + d*eta*n*y*(p*y)^eta))/(y*(m^2*(n - 1) - 2*m*(n - 1) - 1)*(b*rb*th*((-n)*(p - 1)*y)^gamma - n*rb*th*(c1*y + c2*y^2 + d*(p*y)^eta - d*y*(p*y)^eta - 1) + n^2)^2);
(*VALIDATED*)

(*PARAMETER SEARCH*)

numerics[parSet_,mig_] := Block[{},
	(* CALCULATE THE ECOLOGICAL EQUILIBRIUM WITH GIVEN PARAMETER SET *)
	ecoequi = DeleteCases[Flatten[NSolve[(n == dynamics[xm,n,nr]) /. parSet, n]], n -> 0.`];

	(** ->[nZero] Get the size of the ecological equilibrium when y=0 **)

	(** ->[nOne] Get the size of the ecological equilibrium when y=1 **)

	(** ->[nTrend] Check whether n equilibrium increases with y (get the sign of the derivative for y in ]0,1[) **)

	(* CALCULATE SELECTION GRADIENT AT ECOLOGICAL EQUILIBRIUM FOR GIVEN PARAMETER SET *)
	s = gradient /. parSet /. mig -> 0.5 /. ecoequi;

	(** ->[sEqui] Get the selection gradient sign for a range of y between 0 and 1 (included)  **)

	(** ->[sZero] Get the sign of the selection gradient for y=0 and n=neq(1) **)

	Return[nZero,nOne,nTrend,sEqui,sZero]]

search[values_List,mig_] := Block[{names, combinationsList, combinationsRule,results},
	(*Assign the pars values to the private set of parameter rules*)
	names = {th, rb, c1, c2, b, d, p, eta, gamma};
	combinationsList = Tuples[values];
	combinationsRule = Map[MapThread[Rule, {names, combinationsList[[#]]}] &, Range[Length[combinationsList]]];

	(*Calculate all the things we need*)
	results = numerics[#,mig] &/@ combinationsRule;
	Return[results]]

End[]

EndPackage[]


selectionGradient[] := Block[{dn,se,sw,s},
	(*!!!! WRONG !!!!*)
	(*pars = {th -> 0.1, rb -> 20, c1 -> 0, c2 -> 0.1, b -> 2, d -> 0.3, p -> 0.5, eta -> 1, gamma -> 1, m->0.5};*)
	dn = D[dynamics[xm,n,nr], xm] * (1 - (1 - m) * D[dynamics[xm,n,nr], n])^(-1) * (1 - m) * rbar //. monopopulation // Simplify;
	sw = D[fitness[x,xn,x], x] + D[fitness[x, xn, n], xn] * r //. monopopulation // Simplify;
	se = dn * D[fitness[x,xn,n], n] //. monopopulation // Simplify;
	s = sw + se;
	Return[s]]