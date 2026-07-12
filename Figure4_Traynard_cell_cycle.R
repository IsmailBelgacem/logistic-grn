# =====================================================================
#  Section 5  --  Traynard 11-gene mammalian cell-cycle logistic ODE.
#                 R reproduction of the Mathematica figure.
#
#  Each Boolean rule phi_i is mapped to a product-of-logistics Phi_i by the
#  recursive De Morgan formula of Section 4, applied here by the GENERAL
#  translator in BooleanToLogisticODE.R (the native-R equivalent of the
#  Mathematica BooleanToODESys). The regulatory functions are therefore
#  DERIVED from the Boolean rules, not hand-written.
#
#  Each gene obeys   dx_i/dt = kappa_i * Phi_i(x) - gamma_i * x_i,
#  with lambda_i = n/theta_i.  With CycD active (proliferative regime) the
#  Boolean network has NO fixed point but a cyclic attractor, and the ODE
#  settles onto a SUSTAINED LIMIT CYCLE (period ~ 5): Cdc20, Cdh1, CycA,
#  CycB, E2F, UbcH10 oscillate; CycD, CycE, Skp2 stay high; p27, Rb stay off.
#
#  Reproduces: Traynard.jpeg.
#  Requires  : install.packages("deSolve")
# =====================================================================
library(deSolve)
source("BooleanToLogisticODE.R")     # logisticp/logisticm, softValue, booleanToLogisticField

n     <- 4
genes <- c("Cdc20","Cdh1","CycA","CycB","CycE","E2F","p27","Rb","Skp2","UbcH10","CycD")

## ---- Boolean network (Traynard et al.); already in minimal DNF -------
rules <- list(
  Cdc20  = quote(CycB),
  Cdh1   = quote((!CycA & !CycB) | p27),
  CycA   = quote((!Cdc20 & !Cdh1 & CycA) | (!Cdc20 & !Cdh1 & E2F & !Rb) |
                 (CycA & !UbcH10) | (E2F & !Rb & !UbcH10)),
  CycB   = quote((!Cdc20 & !Cdh1) | (!Cdh1 & !UbcH10)),
  CycE   = quote(E2F & !Rb),
  E2F    = quote((!Cdc20 & !CycA & !Rb) | (!Cdc20 & p27 & !Rb) |
                 (!CycA & !CycB & !Rb) | (!CycB & p27 & !Rb)),
  p27    = quote((!CycA & !CycB & !CycD & !CycE) | (!CycA & !CycB & !CycD & p27) |
                 (!CycB & !CycD & !CycE & p27) | (!CycD & !Skp2)),
  Rb     = quote((!CycA & !CycB & !CycD & !CycE) | (!CycA & !CycD & p27) |
                 (!CycB & !CycD & p27) | (!CycD & !CycE & p27)),
  Skp2   = quote(!Cdh1 | !Rb),
  UbcH10 = quote((Cdc20 & UbcH10) | !Cdh1 | (CycA & UbcH10) | (CycB & UbcH10)),
  CycD   = quote(CycD)
)

## ---- kinetic parameters (Table 2) -----------------------------------
kappa <- c(Cdc20=74.33, Cdh1=61.00, CycA=79.17, CycB=76.70, CycE=91.06,
           E2F=56.50, p27=68.79, Rb=53.20, Skp2=66.65, UbcH10=73.51, CycD=64.75)
gamma <- c(Cdc20=0.70, Cdh1=1.45, CycA=1.74, CycB=0.94, CycE=0.58,
           E2F=0.58, p27=0.68, Rb=1.24, Skp2=0.43, UbcH10=1.95, CycD=0.76)
theta <- c(Cdc20=19.18, Cdh1=11.84, CycA=19.29, CycB=18.84, CycE=18.90,
           E2F=17.30, p27=14.69, Rb=11.73, Skp2=12.95, UbcH10=12.05, CycD=19.89)
x0    <- c(Cdc20=0.54, Cdh1=38.68, CycA=96.61, CycB=64.56, CycE=40.97,
           E2F=69.56, p27=5.55, Rb=32.23, Skp2=63.64, UbcH10=45.80, CycD=58.36)

## ---- reduce each rule to minimal DNF BEFORE translating -------------
## R equivalent of the BooleanMinimize step in the Mathematica driver. It is
## kept EXPLICIT here (not hidden inside the translator) so the canonicalisation
## is visible and the reduced rules cannot be mistaken for a translation error.
## booleanToLogisticField is then a pure De Morgan translator. For this network
## the rules are already in minimal DNF, so this step is a transparent no-op.
rulesMin <- lapply(rules, function(r) minimizeDNF(r, ruleVars(r)))

## ---- vector field built by the De Morgan translator -----------------
traynard <- booleanToLogisticField(rulesMin, kappa, gamma, theta, n,
                                   fpos = logisticp, fneg = logisticm)

## ---- integrate over t in [0, 60] ------------------------------------
tt  <- seq(0, 60, by = 0.02)
out <- as.data.frame(ode(y = x0[genes], times = tt, func = traynard, parms = NULL,
                         method = "lsoda", rtol = 1e-8, atol = 1e-8))

## ---- plot all 11 trajectories ---------------------------------------
# png("Traynard.png", width = 1500, height = 1000, res = 150)
par(mar = c(4.5, 5, 2, 8), xpd = NA)
cols <- rainbow(11)
ymax <- max(sapply(genes, function(g) max(out[[g]])))
plot.default(out$time, out[[genes[1]]], type = "l", col = cols[1], lwd = 2, las = 1,
             xlab = "Time (t)", ylab = "Expression level",
             ylim = c(0, ymax))
for (i in 2:11) lines(out$time, out[[genes[i]]], col = cols[i], lwd = 2)
grid(col = "gray", lty = "dotted")
legend(par("usr")[2], par("usr")[4], legend = genes, col = cols, lty = 1, lwd = 2,
       bg = "white", cex = 0.85)
# dev.off()

## ---- console summary: oscillating vs high vs off (t > 40) -----------
cat("Second-half (t>40) behaviour:\n")
for (g in genes) {
  seg   <- out[[g]][out$time > 40]
  pp    <- max(seg) - min(seg)
  level <- max(abs(seg))
  cutoff <- max(1.0, 0.05 * level)
  tag <- if (pp > cutoff) "oscillates" else "steady"
  cat(sprintf("  %-7s mean=%7.2f  p2p=%7.3f  cutoff=%6.2f  %s\n", g, mean(seg), pp, cutoff, tag))
}
