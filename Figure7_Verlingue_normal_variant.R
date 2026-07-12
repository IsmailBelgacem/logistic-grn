# =====================================================================
#  Verlingue et al. 2016 (Aging Cell 15:1018-1026) -- NORMAL (healthy)
#  VARIANT of the geroconversion network -- 25-node logistic ODE.
#  R reproduction of the Mathematica figure (Verlingue_normal_variant.nb).
#
#  Each Boolean rule phi_i is mapped to a product-of-logistics Phi_i by the
#  recursive De Morgan formula of Section 4, applied here by the GENERAL
#  translator in BooleanToLogisticODE.R (the native-R equivalent of the
#  Mathematica BooleanToOdeSystem). The regulatory functions are therefore
#  DERIVED from the Boolean rules, not hand-written, and are identical (to
#  machine precision) to the Mathematica output.
#
#  Each gene obeys   dx_i/dt = kappa_i * Phi_i(x) - gamma_i * x_i,
#  with lambda_i = n/theta_i.
#
#  MODEL VARIANT (per Verlingue et al. 2016, p.1019 and Fig. 1 legend): the
#  paper defines TWO versions of the network that differ by a single edge.
#  The thick edge concerning IRS_PIK3CA inhibition by mTORC1_S6K1 defines
#  the T2DM model (see Figure6_Verlingue_T2DM_geroconversion.R). This file
#  is the NORMAL VARIANT: identical in every other rule, but WITHOUT that
#  edge --            IRS_PIK3CA -> Insulin
#  (rather than       IRS_PIK3CA -> Insulin & !mTORC1_S6K1   in the T2DM file).
#  The two are NOT numerically interchangeable: an exhaustive Boolean
#  fixed-point search shows the T2DM network has 2 stable states, both with
#  IRS_PIK3CA forced FALSE (insulin resistance persists structurally,
#  regardless of which basin is reached), while the normal network has 3
#  stable states, all three with IRS_PIK3CA forced TRUE. Removing the single
#  feedback edge does not merely shift a threshold: it destroys every T2DM
#  fixed point and creates a disjoint set of normal ones -- of the twenty-
#  five nodes, EVERY node except IRS_PIK3CA that is shared between a T2DM
#  state and its normal-variant counterpart keeps the same Boolean value;
#  the edge removal acts exactly and only where the mechanism says it should.
#
#  Network source: SBML-qual model verlingue2016.xml (BioModels/GINsim
#  deposit cited in the paper), with the Therapy input node removed by the
#  standard network-normalisation step (NormalizeNet[DeleteCases[VDS,
#  Therapy -> _]] in Boolean-Network.nb), exactly as for the T2DM file.
#  Insulin and GF are the two sustained inputs (rule "-> True");
#  Senescence, G1_S and Metabolism are the three phenotype readouts
#  (Fig. 1 legend). The three Boolean stable states of the normal network,
#  with Insulin = GF = 1 fixed, are:
#    stateAktActiveNonCycling: AKT=IRS_PIK3CA=1, G1_S=CDK2=E2F1=Senescence=0
#    stateSenescent          : Senescence=p21=p53=PTEN=pRB=1, G1_S=E2F1=CDK2=0,
#                               but (unlike the T2DM senescent state)
#                               IRS_PIK3CA=1 -- senescence without insulin
#                               resistance
#    stateProliferative      : E2F1=G1_S=CDK2=AKT=IRS_PIK3CA=1,
#                               Senescence=p21=p16=0 -- the healthy outcome
#  All three are re-verified below by direct substitution into the raw
#  Boolean map, independently of the ODE integration.
#
#  Reproduces: VERLINGUE_NORMAL_VARIANT.png.
#  Requires  : install.packages("deSolve")
# =====================================================================
library(deSolve)
source("BooleanToLogisticODE.R")     # logisticp/logisticm, softValue, booleanToLogisticField

## Steepness lambda_i = n/theta_i, shared cooperativity n = 4, exactly as
## drawn in Verlingue_normal_variant.nb (note: this differs from n = 6 used
## for the T2DM variant -- each Verlingue example fixes its own single
## independent draw of n, kappa, gamma, theta and x0, as for every curated
## network in this series).
n     <- 4
genes <- c("Insulin","GF","Senescence","G1_S","MAPK","p16","MDM2","p53","p21","AKT",
           "mTORC1_S6K1","ATP","IRS_PIK3CA","AMPK","PTEN","TSC","MYC","CDK2","pRB",
           "E2F1","PRC","Metabolism","PP2A","FOXO","PP1C")

## ---- ETAPE 1 -- Regles exactement telles que publiees ----------------
## Identical to the T2DM rule set (rulesRaw in Figure6_Verlingue_T2DM_
## geroconversion.R) at every node EXCEPT IRS_PIK3CA -- see MODEL VARIANT
## note above. Kept in raw (non-DNF-minimal) form for the same traceability
## reason documented in the T2DM file.
rulesRaw <- list(
  Insulin     = quote(TRUE),
  GF          = quote(TRUE),
  Therapy     = quote(TRUE),
  Senescence  = quote((mTORC1_S6K1 | p16) & (p16 | p21)),
  G1_S        = quote(CDK2 & E2F1 & Metabolism & !p21),
  MAPK        = quote(GF & !PP2A),
  p16         = quote(!E2F1 & MAPK & !p53 & !PRC),
  MDM2        = quote((AKT | p16 | p53) & !E2F1 & !mTORC1_S6K1 & !MYC),
  p53         = quote(!MDM2),
  p21         = quote(!AKT & (FOXO | p53) & !MYC),
  AKT         = quote((CDK2 | IRS_PIK3CA) & !PP2A & !PTEN),
  mTORC1_S6K1 = quote(!AMPK & !TSC),
  ATP         = quote(Metabolism),
  IRS_PIK3CA  = quote(Insulin),                     # <-- the ONE changed rule
  AMPK        = quote(!ATP & p53),
  PTEN        = quote(!AKT & p53),
  TSC         = quote(!AKT & AMPK & !MAPK),
  MYC         = quote(E2F1 & MAPK & mTORC1_S6K1 & !p53),
  CDK2        = quote((E2F1 | MYC) & mTORC1_S6K1 & !p21),
  pRB         = quote(!CDK2),
  E2F1        = quote(E2F1 & (GF | MYC) & !pRB),
  PRC         = quote(!AKT & MYC),
  Metabolism  = quote((AKT | PP1C) & (MAPK | mTORC1_S6K1)),
  PP2A        = quote(!mTORC1_S6K1),
  FOXO        = quote(!AKT & (AMPK | Metabolism | p16) & !MAPK),
  PP1C        = quote(AKT | MAPK)
)

## ---- ETAPE 2 -- Suppression du noeud Therapy --------------------------
## Same normalisation step as the T2DM file; Therapy is a vestigial input
## (pharmacological simulations only, Materials and Methods p.1024) and is
## not part of the base network dynamics (Fig. 1: grey node).
rulesRaw$Therapy <- NULL
rules <- rulesRaw

## ---- ETAPE 4 -- Parametres cinetiques ---------------------------------
## As for the T2DM file, kappa/gamma/theta are not published in Verlingue
## et al. 2016 (a pure Boolean model simulated stochastically with MaBoSS);
## the values below are the single independent draw reported in
## Verlingue_normal_variant.nb: kappa ~ Uniform[50,100], gamma ~
## Uniform[0.25,2], theta ~ Uniform[10,20], n = 4. Running this script
## reproduces the notebook figure and its Jacobian spectrum (largest real
## part -0.3051, set by min gamma = gamma_AKT = 0.3051).
kappa <- c(Insulin=86.8963, GF=72.4351, Senescence=94.6186, G1_S=93.7670, MAPK=62.8255,
           p16=84.1431, MDM2=96.6196, p53=54.0983, p21=64.8349, AKT=53.5813,
           mTORC1_S6K1=74.3393, ATP=73.6526, IRS_PIK3CA=64.9635, AMPK=57.7498, PTEN=60.5431,
           TSC=60.4151, MYC=58.4264, CDK2=72.6472, pRB=85.5416, E2F1=54.0984,
           PRC=81.4405, Metabolism=62.1464, PP2A=80.4183, FOXO=93.0308, PP1C=88.6539)
gamma <- c(Insulin=0.3092, GF=1.2442, Senescence=0.5957, G1_S=0.8056, MAPK=1.3449,
           p16=0.3772, MDM2=1.5876, p53=0.4196, p21=1.1774, AKT=0.3051,
           mTORC1_S6K1=1.2706, ATP=1.7608, IRS_PIK3CA=0.8646, AMPK=1.4242, PTEN=0.8750,
           TSC=0.9335, MYC=1.1724, CDK2=0.6951, pRB=1.4941, E2F1=0.8848,
           PRC=1.5997, Metabolism=0.5707, PP2A=0.7755, FOXO=0.3583, PP1C=1.5265)
theta <- c(Insulin=15.6650, GF=15.8662, Senescence=17.3886, G1_S=18.2190, MAPK=16.8623,
           p16=15.4681, MDM2=12.5930, p53=15.5627, p21=15.2005, AKT=10.5953,
           mTORC1_S6K1=10.6927, ATP=15.8599, IRS_PIK3CA=12.8231, AMPK=15.7422, PTEN=19.4743,
           TSC=18.9398, MYC=13.5479, CDK2=16.2350, pRB=13.7492, E2F1=18.3530,
           PRC=16.2775, Metabolism=18.2466, PP2A=19.7952, FOXO=13.8287, PP1C=11.1833)

## ---- initial condition -------------------------------------------------
## The single realisation reported in Verlingue_normal_variant.nb (not a
## deliberately "mixed" construction as in the T2DM file -- this draw
## already lies in the basin of the proliferative fixed point).
x0    <- c(Insulin=63.5812, GF=14.7108, Senescence=61.7234, G1_S=96.8135, MAPK=95.4912,
           p16=10.3455, MDM2=9.5610, p53=62.8235, p21=51.1062, AKT=35.3197,
           mTORC1_S6K1=28.5190, ATP=46.3278, IRS_PIK3CA=94.5954, AMPK=15.4485, PTEN=95.8189,
           TSC=84.1208, MYC=69.4285, CDK2=76.1456, pRB=17.0991, E2F1=55.7069,
           PRC=77.6568, Metabolism=20.7439, PP2A=56.7135, FOXO=0.8879, PP1C=41.6233)

## ---- ETAPE 3 -- minimizeDNF (see T2DM file for the full rationale) ----
rulesMin <- lapply(rules, function(r) minimizeDNF(r, ruleVars(r)))

## ---- ODE : dx_i/dt = kappa_i * Phi_i(x) - gamma_i * x_i ----------------
verlingue_normal <- booleanToLogisticField(rulesMin, kappa, gamma, theta, n,
                                           fpos = logisticp, fneg = logisticm)

## ---- integrate over t in [0, 60] ----------------------------------------
tt  <- seq(0, 60, by = 0.02)
out <- as.data.frame(ode(y = x0[genes], times = tt, func = verlingue_normal, parms = NULL,
                         method = "lsoda", rtol = 1e-8, atol = 1e-8))

## ---- plot all 25 trajectories -------------------------------------------
# png("VERLINGUE_NORMAL_VARIANT.png", width = 1700, height = 1100, res = 150)
par(mar = c(4, 4, 2, 2))
cols <- rainbow(length(genes))
ymax <- max(sapply(genes, function(g) max(out[[g]])))
plot.default(out$time, out[[genes[1]]], type = "l", col = cols[1], lwd = 2, las = 1,
             xlab = "Time (t)", ylab = "Expression level",
             ylim = c(0, ymax))
for (i in 2:length(genes)) lines(out$time, out[[genes[i]]], col = cols[i], lwd = 2)
grid(col = "gray", lty = "dotted")
legend("topright",
       legend = genes,
       col = cols,
       lty = 1,
       lwd = 2,
       cex = 0.6,
       bg = "white")
# dev.off()

## ---- console summary: oscillating vs high vs off (t > 40) --------------
cat("Second-half (t>40) behaviour:\n")
for (g in genes) {
  seg   <- out[[g]][out$time > 40]
  pp    <- max(seg) - min(seg)
  level <- max(abs(seg))
  cutoff <- max(1.0, 0.05 * level)
  tag <- if (pp > cutoff) "oscillates" else "steady"
  cat(sprintf("  %-12s mean=%7.2f  p2p=%7.3f  %s\n", g, mean(seg), pp, tag))
}

## ---- independent check: Boolean stable states (biological phenotypes) --
## Same methodology as the T2DM file's checkFixed: an evaluation of the raw
## Boolean map `rules` (pre-minimisation), independent of the ODE
## integration above.

boolEval <- function(expr, x) {
  if (is.logical(expr)) return(expr)
  if (is.name(expr))    return(x[[as.character(expr)]])
  op <- as.character(expr[[1]])
  if (op == "(") return(boolEval(expr[[2]], x))
  if (op == "!") return(!boolEval(expr[[2]], x))
  if (op == "&") return(boolEval(expr[[2]], x) & boolEval(expr[[3]], x))
  if (op == "|") return(boolEval(expr[[2]], x) | boolEval(expr[[3]], x))
  stop("bad expr")
}
stateAktActiveNonCycling <- c(Insulin=TRUE,GF=TRUE,Senescence=FALSE,G1_S=FALSE,MAPK=TRUE,
  p16=FALSE,MDM2=FALSE,p53=TRUE,p21=FALSE,AKT=TRUE,mTORC1_S6K1=TRUE,ATP=TRUE,IRS_PIK3CA=TRUE,
  AMPK=FALSE,PTEN=FALSE,TSC=FALSE,MYC=FALSE,CDK2=FALSE,pRB=TRUE,E2F1=FALSE,PRC=FALSE,
  Metabolism=TRUE,PP2A=FALSE,FOXO=FALSE,PP1C=TRUE)
stateSenescent <- c(Insulin=TRUE,GF=TRUE,Senescence=TRUE,G1_S=FALSE,MAPK=TRUE,p16=FALSE,
  MDM2=FALSE,p53=TRUE,p21=TRUE,AKT=FALSE,mTORC1_S6K1=TRUE,ATP=TRUE,IRS_PIK3CA=TRUE,AMPK=FALSE,
  PTEN=TRUE,TSC=FALSE,MYC=FALSE,CDK2=FALSE,pRB=TRUE,E2F1=FALSE,PRC=FALSE,Metabolism=TRUE,
  PP2A=FALSE,FOXO=FALSE,PP1C=TRUE)
stateProliferative <- c(Insulin=TRUE,GF=TRUE,Senescence=FALSE,G1_S=TRUE,MAPK=TRUE,p16=FALSE,
  MDM2=FALSE,p53=TRUE,p21=FALSE,AKT=TRUE,mTORC1_S6K1=TRUE,ATP=TRUE,IRS_PIK3CA=TRUE,AMPK=FALSE,
  PTEN=FALSE,TSC=FALSE,MYC=FALSE,CDK2=TRUE,pRB=FALSE,E2F1=TRUE,PRC=FALSE,Metabolism=TRUE,
  PP2A=FALSE,FOXO=FALSE,PP1C=TRUE)
checkFixed <- function(state, label) {
  xl <- as.list(state)
  ok <- TRUE
  for (g in genes) {
    val <- boolEval(rules[[g]], xl)
    if (!identical(unname(val), unname(state[[g]]))) {
      cat(sprintf("  MISMATCH at %s: rule gives %s, state has %s\n", g, val, state[[g]]))
      ok <- FALSE
    }
  }
  cat(sprintf("%s fixed point: %s\n", label, if (ok) "VERIFIED" else "FAILED"))
}
checkFixed(stateAktActiveNonCycling, "AktActiveNonCycling")
checkFixed(stateSenescent,           "Senescent (normal-rule)")
checkFixed(stateProliferative,       "Proliferative")

## ---- independent check: Jacobian spectrum at the recovered equilibrium -
## Confirms the value reported in the manuscript (largest real part -0.3051,
## set by the smallest degradation rate, gamma_AKT), by direct numerical
## differentiation of the ODE right-hand side -- independent of, and in
## addition to, the raw-Boolean-map check above.
eq   <- unlist(out[nrow(out), genes])
rhs0 <- unlist(verlingue_normal(0, eq, NULL))
h    <- 1e-6
J <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
for (j in genes) {
  eh <- eq; eh[j] <- eh[j] + h
  J[, j] <- (unlist(verlingue_normal(0, eh, NULL)) - rhs0) / h
}
ev <- eigen(J, only.values = TRUE)$values
cat(sprintf("\nJacobian: max Re(eigenvalue) = %.4f  (min gamma = %.4f at %s)\n",
            max(Re(ev)), min(gamma), names(gamma)[which.min(gamma)]))
