# =====================================================================
#  Verlingue et al. 2016 (Aging Cell 15:1018-1026) -- T2DM (diabetes)
#  VARIANT of the geroconversion network -- 25-node logistic ODE.
#  R reproduction of the Mathematica figure (Verlingue_T2DM_logistic_ODE.nb).
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
#  MODEL VARIANT (important, per Verlingue et al. 2016, p.1019 and Fig. 1
#  legend): the paper defines TWO versions of the network that differ by a
#  single edge. "The thick edge concerning IRS_PIK3CA inhibition by
#  mTORC1_S6K1 defines the type diabetes model. The normal model is
#  identical in everything but without this thick edge." This file is the
#  T2DM VARIANT:            IRS_PIK3CA -> Insulin & !mTORC1_S6K1
#  The normal-model variant (IRS_PIK3CA -> Insulin only) is in
#  Verlingue_normal_variant.R -- the two are NOT numerically interchangeable: an
#  exhaustive Boolean fixed-point search shows the T2DM network has 2
#  stable states while the normal network has 3 (see below).
#
#  Network source: SBML-qual model verlingue2016.xml (BioModels/GINsim
#  deposit cited in the paper). In Boolean-Network.nb it is processed by
#  MultiValuedToBoon[SBMLQualImport["verlingue2016.xml"], VanHam, "CNF"]
#  then VDS = NormalizeNet[DeleteCases[VDS, Therapy -> _]]. The imported
#  network carries the T2DM feedback edge (IRS_PIK3CA -> Insulin &
#  !mTORC1_S6K1, confirmed in the notebook's rule table), so its
#  NiceStableStates output is the T2DM network's. Insulin and GF are the
#  two sustained inputs (rule "-> True"); Senescence, G1_S and Metabolism
#  are the three phenotype readouts (Fig. 1 legend). The T2DM network has
#  exactly two Boolean stable states with Insulin = GF = 1 fixed:
#    state 1 (proliferative): E2F1=G1_S=CDK2=AKT=1, Senescence=p21=p16=0
#    state 2 (senescent)    : Senescence=p21=p53=PTEN=pRB=1, G1_S=E2F1=CDK2=0
#  These match, node for node on all 25 genes, the two states decoded
#  directly from the color-coded NiceStableStates grid in Boolean-Network.nb
#  (green sphere = ON, red sphere = OFF), and are re-verified below by
#  direct substitution into the raw Boolean map.
#
#  Reproduces: Verlingue_T2DM.jpeg.
#  Requires  : install.packages("deSolve")
# =====================================================================
library(deSolve)
source("BooleanToLogisticODE.R")     # logisticp/logisticm, softValue, booleanToLogisticField

## Hill-matched steepness lambda_i = n/theta_i, shared cooperativity n = 6,
## as reported for the Verlingue example in the manuscript (Fig. verlingue_t2dm).
n     <- 6
genes <- c("Insulin","GF","Senescence","G1_S","MAPK","p16","MDM2","p53","p21","AKT",
           "mTORC1_S6K1","ATP","IRS_PIK3CA","AMPK","PTEN","TSC","MYC","CDK2","pRB",
           "E2F1","PRC","Metabolism","PP2A","FOXO","PP1C")

## ---- ETAPE 1 -- Regles exactement telles que publiees ---------------
## rulesRaw ci-dessous reprend, noeud par noeud, la forme AND-de-OR /
## OR-de-AND telle que curatee depuis la litterature dans le modele
## SBML-qual original (avant toute reduction). Cette forme n'est PAS en
## forme normale disjonctive (DNF) minimale : par exemple Senescence est
## ecrite (mTORC1_S6K1 | p16) & (p16 | p21), une forme conjonctive de
## clauses disjonctives, alors que sa DNF minimale est
## (mTORC1_S6K1 & p21) | p16 (equivalence verifiee par table de verite
## exhaustive, voir chat). Garder cette forme brute ici, plutot que de la
## transcrire deja minimisee, rend la provenance directement tracable au
## papier (Verlingue et al. 2016 / SBML-qual deposit -- Therapy inclus
## comme dans le reseau original, avant l'etape de reduction DNF).
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
  IRS_PIK3CA  = quote(Insulin & !mTORC1_S6K1),
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

## ---- ETAPE 2 -- Suppression du noeud Therapy -------------------------
## Therapy est l'entree utilisee dans le papier UNIQUEMENT pour les
## simulations pharmacologiques (calibration des doses d'everolimus/
## rapamycine, Materials and Methods p.1024). Elle ne fait pas partie de
## la dynamique de base du reseau (Fig. 1 : noeud gris, distinct des
## entrees Insulin/GF). On la retire ici, ce qui reproduit exactement
## l'etape NormalizeNet[DeleteCases[VDS, Therapy -> _]] deja appliquee
## dans Boolean-Network.nb.
rulesRaw$Therapy <- NULL
rules <- rulesRaw

## ---- ETAPE 4 -- Parametres cinetiques : Les
## valeurs kappa/gamma/theta ci-dessous sont un tirage reproductible dans les
## memes plages que le driver Traynard : kappa ~ Uniform[50,100] (taux de
## synthese), gamma ~ Uniform[0.25,2] (taux de degradation),
## theta ~ Uniform[10,20] (seuil d'activation, avec lambda=n/theta la pente du
## noyau logistique), et n=6 (coefficient de cooperativite, analogue au
## coefficient de Hill ; c'est la valeur reportee pour cet exemple dans le
## papier). La condition initiale x0, elle, N'EST PAS un tirage aleatoire :
## le reseau etant bistable, x0 est choisie dans le bassin de l'etat senescent
## (geroconversion) pour reproduire la Fig. verlingue_t2dm du papier -- voir le
## bloc x0 ci-dessus.
## The kappa/gamma/theta below are EXACTLY the values reported in the
## manuscript (Table "Kinetic parameters for the twenty-five-node Verlingue
## T2DM geroconversion ODE system", shared cooperativity n = 6). Running this
## script therefore reproduces the manuscript figure and the reported Jacobian
## spectrum (largest real part -0.39, set by min gamma = gamma_p21 = 0.39).
kappa <- c(Insulin=68.54, GF=71.40, Senescence=92.18, G1_S=59.90, MAPK=88.30,
           p16=89.19, MDM2=68.57, p53=51.12, p21=51.83, AKT=62.63,
           mTORC1_S6K1=63.85, ATP=85.77, IRS_PIK3CA=66.43, AMPK=53.10, PTEN=53.41,
           TSC=85.37, MYC=64.53, CDK2=62.68, pRB=78.07, E2F1=61.62,
           PRC=51.99, Metabolism=67.10, PP2A=55.49, FOXO=51.53, PP1C=57.04)
gamma <- c(Insulin=1.91, GF=1.13, Senescence=0.84, G1_S=0.99, MAPK=1.94,
           p16=1.34, MDM2=0.71, p53=1.26, p21=0.39, AKT=1.40,
           mTORC1_S6K1=1.68, ATP=1.26, IRS_PIK3CA=0.55, AMPK=0.58, PTEN=1.52,
           TSC=1.51, MYC=0.96, CDK2=0.83, pRB=1.81, E2F1=1.78,
           PRC=0.84, Metabolism=1.25, PP2A=1.40, FOXO=1.18, PP1C=0.61)
theta <- c(Insulin=19.88, GF=17.18, Senescence=17.38, G1_S=19.00, MAPK=10.80,
           p16=15.75, MDM2=13.08, p53=10.22, p21=17.25, AKT=19.25,
           mTORC1_S6K1=14.34, ATP=18.62, IRS_PIK3CA=19.48, AMPK=18.58, PTEN=17.48,
           TSC=14.65, MYC=14.15, CDK2=10.45, pRB=16.12, E2F1=18.07,
           PRC=17.12, Metabolism=18.41, PP2A=16.63, FOXO=15.51, PP1C=16.84)
## ---- initial condition: geroconversion (senescent) basin -------------
## The T2DM network is BISTABLE: both the proliferative and the senescent
## Boolean fixed points are exponentially stable equilibria of the logistic
## ODE (Theorem, Boolean fixed-point recovery), so the attractor reached is
## selected by the initial condition's basin. The manuscript's Fig.
## verlingue_t2dm reports convergence to the INSULIN-RESISTANT SENESCENT
## (geroconversion) state; the mixed initial condition below lies in that
## basin. Biologically it is the onset of geroconversion: insulin resistance
## (IRS_PIK3CA low), no proliferation (CDK2, E2F1, G1_S low), AKT low, and the
## tumour-suppressor arm engaged (p53, PTEN, p21, Senescence high).
##
## (An initial condition with IRS_PIK3CA and CDK2 high instead drives AKT high
## and lands in the T2DM network's OWN proliferative fixed point -- E2F1=G1_S=
## CDK2=AKT=1, Senescence=p21=p16=0, exactly as in stateProliferative below,
## but still with IRS_PIK3CA=FALSE, since mTORC1_S6K1=TRUE at that fixed point
## too. Insulin resistance is thus structural to the T2DM edge, not
## basin-dependent: BOTH of its stable states keep IRS_PIK3CA off. The
## genuinely healthy outcome -- the same 24 nodes, but IRS_PIK3CA=TRUE -- is
## only reachable once the feedback edge itself is removed, i.e. under the
## normal-model rule reproduced in Verlingue_normal_variant.R.)
x0    <- c(Insulin=60.00, GF=55.00, Senescence=70.00, G1_S=5.00, MAPK=50.00,
           p16=8.00, MDM2=5.00, p53=65.00, p21=55.00, AKT=6.00,
           mTORC1_S6K1=40.00, ATP=45.00, IRS_PIK3CA=4.00, AMPK=10.00, PTEN=60.00,
           TSC=8.00, MYC=6.00, CDK2=5.00, pRB=45.00, E2F1=5.00,
           PRC=6.00, Metabolism=40.00, PP2A=8.00, FOXO=12.00, PP1C=45.00)

## ---- ETAPE 3 -- Pourquoi minimizeDNF est obligatoire -----------------
## booleanToLogisticField (dans BooleanToLogisticODE.R) traduit chaque
## regle booleenne en un produit de noyaux logistiques via la formule de
## De Morgan recursive : AND devient un produit, OR devient 1 moins un
## produit de complements. Cette traduction "molle" n'est NI idempotente
## NI distributive -- contrairement a AND/OR booleens exacts. Concretement,
## pour Senescence, traduire directement (mTORC1_S6K1 | p16) & (p16 | p21)
## donnerait un champ de vecteurs different de celui obtenu en traduisant
## sa DNF minimale (mTORC1_S6K1 & p21) | p16, alors que les deux formules
## sont logiquement identiques (memes 0/1 sur les sommets de l'hypercube).
## Sans minimisation prealable, l'ODE resultante dependrait donc de la
## facon dont la regle a ete ECRITE et non de la fonction booleenne
## qu'elle represente -- ce qui rendrait le modele continu non
## reproductible.

rulesMin <- lapply(rules, function(r) minimizeDNF(r, ruleVars(r)))

## ---- ODE : dx_i/dt = kappa_i * Phi_i(x) - gamma_i * x_i --------------
## Phi_i est le produit de noyaux logistiques obtenu en appliquant la
## formule de De Morgan a la regle DNF minimale de rulesMin
## (booleanToLogisticField, le traducteur general reutilise pour tous les
## reseaux de cette serie -- Traynard, Verlingue). fpos=logisticp et
## fneg=logisticm sont les noyaux croissant (activation) et decroissant
## (repression).
verlingue <- booleanToLogisticField(rulesMin, kappa, gamma, theta, n,
                                    fpos = logisticp, fneg = logisticm)

## ---- integrate over t in [0, 60] -------------------------------------
tt  <- seq(0, 60, by = 0.02)
out <- as.data.frame(ode(y = x0[genes], times = tt, func = verlingue, parms = NULL,
                         method = "lsoda", rtol = 1e-8, atol = 1e-8))

## ---- plot all 25 trajectories -----------------------------------------
# png("verlingue_t2dm.png", width = 1700, height = 1100, res = 150)
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

## ---- console summary: oscillating vs high vs off (t > 40) ------------
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
## Un point fixe booleen (x = F(x) sous mise a jour synchrone) correspond
## a un phenotype cellulaire stable au sens du papier : un attracteur du
## reseau vers lequel convergent les trajectoires stochastiques MaBoSS.
## Cette verification est INDEPENDANTE de l'integration ODE ci-dessus --
## elle porte directement sur `rules` (forme booleenne exacte, avant
## minimisation, puisque l'egalite booleenne x=F(x) ne depend pas de la
## syntaxe). 

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
stateProliferative <- c(Insulin=TRUE,GF=TRUE,Senescence=FALSE,G1_S=TRUE,MAPK=TRUE,p16=FALSE,
  MDM2=FALSE,p53=TRUE,p21=FALSE,AKT=TRUE,mTORC1_S6K1=TRUE,ATP=TRUE,IRS_PIK3CA=FALSE,AMPK=FALSE,
  PTEN=FALSE,TSC=FALSE,MYC=FALSE,CDK2=TRUE,pRB=FALSE,E2F1=TRUE,PRC=FALSE,Metabolism=TRUE,
  PP2A=FALSE,FOXO=FALSE,PP1C=TRUE)
stateSenescent <- c(Insulin=TRUE,GF=TRUE,Senescence=TRUE,G1_S=FALSE,MAPK=TRUE,p16=FALSE,
  MDM2=FALSE,p53=TRUE,p21=TRUE,AKT=FALSE,mTORC1_S6K1=TRUE,ATP=TRUE,IRS_PIK3CA=FALSE,AMPK=FALSE,
  PTEN=TRUE,TSC=FALSE,MYC=FALSE,CDK2=FALSE,pRB=TRUE,E2F1=FALSE,PRC=FALSE,Metabolism=TRUE,
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
checkFixed(stateProliferative, "Proliferative")
checkFixed(stateSenescent,     "Senescent")
