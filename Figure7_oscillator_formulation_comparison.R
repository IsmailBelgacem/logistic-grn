# =====================================================================
#  Two-gene oscillator (Section 3.1) under THREE modelling formulations
#    (a) Hill functions          (slope-matched: n = lambda*theta)
#    (b) Product-of-logistics    (this work, exactly Section 3.1)
#    (c) Samuilik weighted sum   (shared threshold theta_i = sum(w)/2)
#  + the biological-threshold Samuilik collapse (Remark 6.1)
#
#  Reproduces Figure 7 (the 2x2 four-formulation oscillator comparison).
#  (Figure 8, the static repression-function comparison, is in Figure8_repression_function_comparison.R.)
# =====================================================================
library(deSolve)

## ---- shared parameters (identical to Section 3.1 / Fig. 2) ----------
lambda <- 3
kappa1 <- 3 ; gamma1 <- 0.25
kappa2 <- 4 ; gamma2 <- 0.5
theta1 <- 4 ; theta2 <- 3          # biological thresholds (EC50 / Kd)
state  <- c(x1 = 1, x2 = 1)
time   <- seq(0, 25, by = 0.01)

## ---- regulatory response functions ----------------------------------
# logistic (this work)
f_plus  <- function(x, theta, lambda) 1 / (1 + exp(-lambda * (x - theta)))
f_minus <- function(x, theta, lambda) 1 - f_plus(x, theta, lambda)
# Hill
h_plus  <- function(x, theta, n) x^n / (x^n + theta^n)
h_minus <- function(x, theta, n) theta^n / (x^n + theta^n)
# Samuilik single increasing sigmoid of a signed weighted sum
sam <- function(arg, mu) 1 / (1 + exp(-mu * arg))

## ---- (a) HILL : matched at each threshold via lambda = n/theta -------
n_act <- lambda * theta1     # = 12  (gene 2 activated by gene 1)
n_rep <- lambda * theta2     # =  9  (gene 1 repressed by gene 2)
sys_hill <- function(t, s, p) with(as.list(s), {
  dx1 <- kappa1 * h_minus(x2, theta2, n_rep) - gamma1 * x1
  dx2 <- kappa2 * h_plus (x1, theta1, n_act) - gamma2 * x2
  list(c(dx1, dx2))
})

## ---- (b) LOGISTIC : exactly Section 3.1 ------------------------------
sys_log <- function(t, s, p) with(as.list(s), {
  dx1 <- kappa1 * f_minus(x2, theta2, lambda) - gamma1 * x1
  dx2 <- kappa2 * f_plus (x1, theta1, lambda) - gamma2 * x2
  list(c(dx1, dx2))
})

## ---- (c) SAMUILIK : prescription theta_i = sum_j w_ij / 2 ------------
mu  <- lambda                # matched steepness
w12 <- -1 ; w21 <- 1         # repression of x1 by x2 ; activation of x2 by x1
th1S <- w12 / 2              # = -0.5
th2S <- w21 / 2              # = +0.5
sys_sam <- function(t, s, p) with(as.list(s), {
  dx1 <- kappa1 * sam(w12 * x2 - th1S, mu) - gamma1 * x1
  dx2 <- kappa2 * sam(w21 * x1 - th2S, mu) - gamma2 * x2
  list(c(dx1, dx2))
})

## ---- Samuilik with MEASURED thresholds (Remark 6.1, collapse) --------
# repression of gene1 uses theta_12 = 3, activation of gene2 uses theta_21 = 4
# (exactly the logistic theta2, theta1). Repressor midpoint x_c = theta_12/w12 = -3 < 0.
sys_sam_bio <- function(t, s, p) with(as.list(s), {
  dx1 <- kappa1 * sam(w12 * x2 - theta2, mu) - gamma1 * x1   # pinned ~ 0
  dx2 <- kappa2 * sam(w21 * x1 - theta1, mu) - gamma2 * x2
  list(c(dx1, dx2))
})

## ---- integrate -------------------------------------------------------
H <- as.data.frame(ode(state, time, sys_hill,    NULL))
L <- as.data.frame(ode(state, time, sys_log,     NULL))
S <- as.data.frame(ode(state, time, sys_sam,     NULL))
B <- as.data.frame(ode(state, time, sys_sam_bio, NULL))

cat(sprintf("equilibria  Hill (%.3f,%.3f) | Logistic (%.3f,%.3f) | Samuilik (%.3f,%.3f) | bio-collapse (%.4f,%.4f)\n",
            tail(H$x1,1),tail(H$x2,1), tail(L$x1,1),tail(L$x2,1),
            tail(S$x1,1),tail(S$x2,1), tail(B$x1,1),tail(B$x2,1)))

## ---- helper: one panel in the Fig. 2 style ---------------------------
panel <- function(df, ttl, ymax = 6) {
  plot.default(df$time, df$x1, type = "l", col = "red", lwd = 2,
               xlab = "Time (t)",
               ylab = expression(paste(x[1], "(t), ", x[2], "(t)")),
               xlim = c(0, 24), ylim = c(0, ymax), las = 1, main = ttl)
  lines(df$time, df$x2, col = "blue", lwd = 2)
  grid(col = "gray", lty = "dotted", lwd = 1.5)
  legend("topright", legend = c(expression(x[1](t)), expression(x[2](t))),
         col = c("red", "blue"), lty = 1, lwd = 2, bg = "white")
}

## ---- FIGURE: four-formulation comparison (2 x 2, common 0-6 scale) ----
# png("oscillator_comparison.png", width = 1600, height = 1150, res = 160)
par(mfrow = c(2, 2), mar = c(4, 5, 3, 1))
panel(H,  "(a) Hill")
panel(L,  "(b) Logistic (this work)")
panel(S,  expression(paste("(c) Samuilik, prescription  ",
                           theta[i] == frac(1,2) * Sigma * w[ij])))
panel(B,  expression(paste("(d) Samuilik, measured  ",
                           theta[12] == 3, ", ", theta[21] == 4)))
# dev.off()
