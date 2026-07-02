# =====================================================================
#  Figure 2  --  Two-gene logistic oscillator  (Section 3.1)
#
#  Gene 1 is repressed by gene 2 (decreasing logistic f^-),
#  gene 2 is activated by gene 1 (increasing logistic f^+):
#
#     dx1/dt = kappa1 * f^-(x2, theta2, lambda) - gamma1 * x1
#     dx2/dt = kappa2 * f^+(x1, theta1, lambda) - gamma2 * x2
#
#  The trajectory shows damped oscillations spiralling into the unique
#  equilibrium (x1*, x2*) ~ (3.87, 3.25)  (Theorem 3.1: globally
#  asymptotically stable, no Hopf bifurcation).
#
#  Reproduces: Oscillateur_original.png
#  Requires  : install.packages("deSolve")
# =====================================================================
library(deSolve)

## ---- parameters (identical to the manuscript) -----------------------
lambda <- 3
kappa1 <- 3 ; gamma1 <- 0.25
kappa2 <- 4 ; gamma2 <- 0.5
theta1 <- 4 ; theta2 <- 3
state  <- c(x1 = 1, x2 = 1)          # initial condition x1(0)=x2(0)=1
time   <- seq(0, 40, by = 0.01)

## ---- logistic response functions ------------------------------------
f_plus  <- function(x, theta, lambda) 1 / (1 + exp(-lambda * (x - theta)))  # increasing
f_minus <- function(x, theta, lambda) 1 / (1 + exp( lambda * (x - theta)))  # decreasing

## ---- vector field ----------------------------------------------------
oscillator <- function(t, s, p) with(as.list(s), {
  dx1 <- kappa1 * f_minus(x2, theta2, lambda) - gamma1 * x1
  dx2 <- kappa2 * f_plus (x1, theta1, lambda) - gamma2 * x2
  list(c(dx1, dx2))
})

## ---- integrate -------------------------------------------------------
out <- as.data.frame(ode(y = state, times = time, func = oscillator, parms = NULL))
cat(sprintf("equilibrium ~ (%.4f, %.4f)\n", tail(out$x1, 1), tail(out$x2, 1)))

## ---- plot (Fig. 2 style: x1 red, x2 blue, dotted grey grid) ----------
# png("Oscillateur_original.png", width = 1400, height = 950, res = 160)
par(mar = c(4.5, 5, 2, 1))
plot.default(out$time, out$x1, type = "l", col = "red", lwd = 2,
             xlab = "Time (t)",
             ylab = expression(paste(x[1], "(t),  ", x[2], "(t)")),
             ylim = range(0, out$x1, out$x2), las = 1)
lines(out$time, out$x2, col = "blue", lwd = 2)
grid(col = "gray", lty = "dotted", lwd = 1.5)
legend("topright", legend = c(expression(x[1](t)), expression(x[2](t))),
       col = c("red", "blue"), lty = 1, lwd = 2, bg = "white")
# dev.off()
