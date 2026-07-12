# =====================================================================
#  Figure 3  --  Logistic genetic toggle switch  (Section 3.2)
#
#     dx1/dt = kappa * f^-(x2, theta, lambda) - gamma * x1
#     dx2/dt = kappa * f^-(x1, theta, lambda) - gamma * x2
#
#  Symmetric parameters kappa=10, gamma=1, theta=5  (so rho = kappa/gamma = 10,
#  theta = rho/2, and the pitchfork threshold is lambda* = 4/rho = 0.4).
#
#  (a) Bistable phase portrait at lambda = 2: two nullclines, two stable
#      nodes ~ (9.9995, 0.0005) and (0.0005, 9.9995) with det J = +1,
#      eigenvalues {-1.001,-0.999}; central saddle (5,5) with det J = -24,
#      eigenvalues {4,-6}. Trajectories from a grid are coloured by basin.
#  (b) Bifurcation diagram x1* vs lambda: a single branch for lambda < 0.4
#      splits at lambda* = 0.4 into a saddle (red) and two stable nodes
#      (green) -- a supercritical pitchfork.
#
#  Reproduces: toggle_switch.png
#  Requires  : install.packages("deSolve")
# =====================================================================
library(deSolve)

## ---- parameters ------------------------------------------------------
kappa <- 10 ; gamma <- 1 ; theta <- 5
rho   <- kappa / gamma                 # = 10
lam_a <- 2                             # steepness for panel (a)
lam_star <- 4 / rho                    # pitchfork threshold = 0.4

f_minus <- function(x, theta, lambda) 1 / (1 + exp(lambda * (x - theta)))

## ---- toggle vector field (steepness lam passed through parms) --------
toggle <- function(t, s, p) with(as.list(s), {
  dx1 <- kappa * f_minus(x2, theta, p$lam) - gamma * x1
  dx2 <- kappa * f_minus(x1, theta, p$lam) - gamma * x2
  list(c(dx1, dx2))
})

## =====================================================================
##  PANEL (a): bistable phase portrait at lambda = 2
## =====================================================================
## return map T and its derivative (for equilibria / classification)
g_fac <- function(x, lam) lam * f_minus(x, theta, lam) * (1 - f_minus(x, theta, lam))
Tmap  <- function(x1, lam) rho * f_minus(rho * f_minus(x1, theta, lam), theta, lam)
Tder  <- function(x1, lam) {                       # T'(x1) = rho^2 g1 g2
  u  <- rho * f_minus(x1, theta, lam)
  rho * rho * g_fac(x1, lam) * g_fac(u, lam)
}
equilibria <- function(lam, nscan = 4000) {
  xs <- seq(1e-6, rho - 1e-6, length.out = nscan)
  gg <- Tmap(xs, lam) - xs
  rts <- numeric(0)
  for (i in seq_len(nscan - 1)) {
    if (gg[i] == 0 || gg[i] * gg[i + 1] < 0)
      rts <- c(rts, uniroot(function(z) Tmap(z, lam) - z,
                            lower = xs[i], upper = xs[i + 1])$root)
  }
  sort(unique(round(rts, 6)))
}

eqx1 <- equilibria(lam_a)              # {~0.0005, 5, ~9.9995}
eqs  <- data.frame(x1 = eqx1, x2 = rho * f_minus(eqx1, theta, lam_a))

# classify by determinant det J = gamma^2 (1 - T'(x1*))
detJ <- gamma^2 * (1 - Tder(eqx1, lam_a))
node_idx   <- which(detJ > 0)          # stable nodes
saddle_idx <- which(detJ < 0)          # saddle

# nullclines: x1 = rho f^-(x2)  and  x2 = rho f^-(x1)
grid_x <- seq(0, rho, length.out = 400)
nc1_x2 <- grid_x                       # for dx1=0: x1 = rho f^-(x2)
nc1_x1 <- rho * f_minus(grid_x, theta, lam_a)
nc2_x1 <- grid_x                       # for dx2=0: x2 = rho f^-(x1)
nc2_x2 <- rho * f_minus(grid_x, theta, lam_a)

# basins: integrate from a grid, colour by nearest stable node reached
tt <- seq(0, 60, by = 0.1)
ic <- expand.grid(x1 = seq(0.5, 9.5, by = 1.0), x2 = seq(0.5, 9.5, by = 1.0))
node_hi <- eqs[node_idx[which.max(eqs$x1[node_idx])], ]   # ~ (9.9995, 0.0005)
basin_col <- character(nrow(ic))
traj <- vector("list", nrow(ic))
for (k in seq_len(nrow(ic))) {
  z <- as.data.frame(ode(c(x1 = ic$x1[k], x2 = ic$x2[k]), tt, toggle,
                         parms = list(lam = lam_a)))
  traj[[k]] <- z
  endpt <- c(tail(z$x1, 1), tail(z$x2, 1))
  basin_col[k] <- if (endpt[1] > endpt[2]) "#d73027" else "#4575b4"  # red / blue basin
}

## =====================================================================
##  PANEL (b): bifurcation diagram x1* vs lambda
## =====================================================================
lams <- seq(0.05, 2.0, by = 0.01)
bp_lam <- numeric(0); bp_x1 <- numeric(0); bp_stab <- logical(0)
for (lam in lams) {
  for (xr in equilibria(lam)) {
    bp_lam  <- c(bp_lam, lam)
    bp_x1   <- c(bp_x1, xr)
    bp_stab <- c(bp_stab, (gamma^2 * (1 - Tder(xr, lam))) > 0)  # TRUE = stable node
  }
}

## ---- draw both panels -----------------------------------------------
# png("toggle_switch.png", width = 1700, height = 800, res = 160)
par(mfrow = c(1, 2), mar = c(4.5, 5, 3, 1))

## (a) phase portrait
plot.default(NA, NA, xlim = c(0, rho), ylim = c(0, rho), las = 1,
             xlab = expression(x[1]), ylab = expression(x[2]),
             main = expression(paste("(a) Bistable phase portrait, ", lambda == 2)))
for (k in seq_len(nrow(ic)))
  lines(traj[[k]]$x1, traj[[k]]$x2, col = adjustcolor(basin_col[k], 0.5), lwd = 1)
lines(nc1_x1, nc1_x2, col = "blue", lwd = 2.5)   # dx1 = 0
lines(nc2_x1, nc2_x2, col = "red",  lwd = 2.5)   # dx2 = 0
points(eqs$x1[node_idx],   eqs$x2[node_idx],   pch = 19, cex = 1.6, col = "black")
points(eqs$x1[saddle_idx], eqs$x2[saddle_idx], pch =  4, cex = 1.8, lwd = 3, col = "black")
legend("top", legend = c(expression(dot(x)[1] == 0), expression(dot(x)[2] == 0),
                         "stable node", "saddle"),
       col = c("blue", "red", "black", "black"),
       lty = c(1, 1, NA, NA), pch = c(NA, NA, 19, 4), lwd = 2, bg = "white", cex = 0.9)

## (b) bifurcation diagram
plot.default(bp_lam[bp_stab], bp_x1[bp_stab], pch = 19, cex = 0.4, col = "forestgreen",
             xlim = range(lams), ylim = c(0, rho), las = 1,
             xlab = expression(lambda), ylab = expression(x[1]^"*"),
             main = "(b) Pitchfork bifurcation")
points(bp_lam[!bp_stab], bp_x1[!bp_stab], pch = 19, cex = 0.4, col = "red")
abline(v = lam_star, lty = "dotted", lwd = 1.5, col = "gray40")
text(lam_star, rho * 0.92, expression(lambda^"*" == 0.4), pos = 4, col = "gray40")
legend("right", legend = c("stable node", "saddle"),
       col = c("forestgreen", "red"), pch = 19, bg = "white", cex = 0.9)
# dev.off()

## ---- console summary -------------------------------------------------
cat("lambda = 2 equilibria (x1*, x2*):\n"); print(round(eqs, 4))
cat(sprintf("det J: %s\n", paste(sprintf("%.2f", detJ), collapse = ", ")))
cat(sprintf("pitchfork threshold lambda* = 4/rho = %.2f\n", lam_star))
