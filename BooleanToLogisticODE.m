(* ::Package:: *)

(* ::Title:: *)
(*Boolean-to-Logistic ODE Translation (The recursive De Morgan formula)*)


(* ::Author:: *)
(*Ismail Belgacem 2026*)


(* ::Affiliation:: *)
(*Paris*)


(* ============================================================================ *)
(*     * BooleanToOdeSystem  -- maps a Boolean network to a product-of-logistics   *)
(*                           ODE system via the recursive De Morgan formula;    *)
(*     * logisticp        -- increasing (activation) logistic kernel;           *)
(*     * logisticm        -- decreasing (repression) logistic kernel.           *)
(*                                                                              *)
(* ============================================================================ *)


BeginPackage["BooleanToLogisticODE`"];

(* Allow reloading in the same kernel: lift any protection from a previous load. *)
Unprotect["BooleanToLogisticODE`*"];


(* ----------------------------------------------------------------------------- *)
(* ----------------------------------------------------------------------------- *)

BooleanToOdeSystem::usage="BooleanToOdeSystem[F,v, \[Kappa],\[Gamma],\[Theta], n, var0, fstep] converts a Boolean network into ODE system. The structures are associations the keys of which are variables of the network.
PARAMETERS:
\[FilledSmallSquare] F : Boolean network, given as a list of rules variable \[Rule] Boolean formula.
\[FilledSmallSquare] v : variable of the system, default t.
\[FilledSmallSquare] \[Kappa] : Expression rate (Association var \[Rule] Real).
\[FilledSmallSquare] \[Gamma] : Decay rate (Association var \[Rule] Real).
\[FilledSmallSquare] \[Theta] : Threshold of expression (Association var \[Rule] Real).
\[FilledSmallSquare] n :  Steepenss of Hill function (The Hill coeficeint) 
\[FilledSmallSquare] var0: initial condition (Association). 
\[FilledSmallSquare] fstep: pair of step functions {function inhibition, function activation}
";

logisticp::usage="logisticp[x, thr, n] is the increasing (activation) logistic kernel 1/(1+Exp[-(n/thr)(x-thr)]). The steepness is matched to a Hill function of cooperativity n at the threshold via \[Lambda]=n/thr.";

logisticm::usage="logisticm[x, thr, n] is the decreasing (repression) logistic kernel 1/(1+Exp[(n/thr)(x-thr)]). The steepness is matched to a Hill function of cooperativity n at the threshold via \[Lambda]=n/thr.";

(* Message issued by BooleanToOdeSystem on an unrecognised Boolean formula. *)
b2o::nnarg="Syntax error formula not recognized aborted: `1`  ";


Begin["`Private`"];

(* ----------------------------------------------------------------------------- *)
(*  Logistic kernels (lambda_i = n / thr).                                       *)
(* ----------------------------------------------------------------------------- *)

logisticp[x_?NumericQ, thr_?NumericQ, n_?NumericQ] := 1/(1 + Exp[-(n/thr) (x - thr)]);

logisticm[x_?NumericQ, thr_?NumericQ, n_?NumericQ] := 1/(1 + Exp[ (n/thr) (x - thr)]);


(* ----------------------------------------------------------------------------- *)
(*  Boolean network  ->  product-of-logistics ODE system  (PURE De Morgan map).   *)
(*                                                                              *)
(*  Each Boolean rule is mapped, through the recursive De Morgan correspondence,  *)
(*  to a product of increasing/decreasing logistic factors; every gene obeys      *)
(*       x_i'[v] == kappa_i * Phi_i(x) - gamma_i * x_i[v],   x_i[0] == var0_i.    *)
(*                                                                              *)
(*  The rules F are assumed to be in minimal disjunctive normal form. Since the   *)
(*  soft map is neither idempotent nor distributive, that canonical form is what  *)
(*  makes Phi a function of the Boolean FUNCTION rather than of how it is written. *)

Clear[BooleanToOdeSystem]
BooleanToOdeSystem[F_, v_Symbol, \[Kappa]_Association, \[Gamma]_Association,
                \[Theta]_Association, n_, var0_Association,
                fstep:{_Symbol, _Symbol}] :=
  Module[{B2O, fneg = fstep[[1]], fpos = fstep[[2]]},

    B2O[network_List]   := Flatten[B2O /@ network];
    B2O[x_Symbol -> f_] := {x'[v] == \[Kappa][x]*B2O[f] - \[Gamma][x]*x[v],
                             x[0]  == var0[x]};

    B2O[var_Symbol]       := fpos[var[v], \[Theta][var], n];

    (* Negation of a single variable *)
    B2O[\[Not] var_Symbol]     := fneg[var[v], \[Theta][var], n];

    (* FIX 2: negation of compound expression \[LongDash] stays in [0,1] *)
    B2O[\[Not] f_]             := 1 - B2O[f];

    (* AND: cooperative \[LongDash] product stays in [0,1] *)
    B2O[f1_ \[And] f2_]        := B2O[f1] * B2O[f2];

    (* OR: De Morgan \[LongDash] stays in [0,1], biologically correct *)
    B2O[f1_ \[Or] f2_]        := 1 - (B2O[\[Not] f1] * B2O[\[Not] f2]);

    B2O[True]  := 1;
    B2O[False] := 0;
    B2O[x_]    := (Message[b2o::nnarg, x]; Abort[]);

    B2O[F]  (* translate the supplied rules directly; they are assumed minimal-DNF
              (the driver applies BooleanMinimize before calling this routine) *)
  ]

End[];


(* ----------------------------------------------------------------------------- *)
SetAttributes[#, {Protected}] & /@ {BooleanToOdeSystem, logisticp, logisticm};

EndPackage[]
