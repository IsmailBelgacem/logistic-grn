# =====================================================================
#  Figure 8  --  Repression: logistic vs Samuilik vs Hill  (Section 6.3)
#
#  Three decreasing regulatory functions of a single repressor x >= 0:
#
#   * our logistic  f^-(w x, theta, lambda) = 1/(1 + e^{-lambda(theta - w x)}),
#       with w = +1  ->  midpoint x_c = theta/w = +3  (biologically meaningful)
#   * Samuilik      f2^-(x, theta, mu)      = 1/(1 + e^{-mu(w x - theta)}),
#       with w = -1  ->  midpoint x_c = theta/w = -3  (outside x >= 0, stays ~0)
#   * Hill          h^-(x, theta, n)        = theta^n / (x^n + theta^n)
#
#  Parameters: n = 4, theta = 3, lambda = mu = n/theta ~ 1.3333.
#  Our logistic and Hill both pass through 1/2 at x = theta = 3; the
#  Samuilik form is pinned near 0 throughout the admissible domain.
#
#  Reproduces: Coparaison_3.png
#  Requires  : base R only.
# =====================================================================

## ---- parameters ------------------------------------------------------
n      <- 4
theta  <- 3
lambda <- n / theta          # ~ 1.3333
mu     <- lambda
w_our  <- 1                  # positive weight, decreasing logistic
w_sam  <- -1                 # negative weight inside an increasing sigmoid

## ---- the three functions --------------------------------------------
f_our <- function(x) 1 / (1 + exp(-lambda * (theta - w_our * x)))  # our f^-
f_sam <- function(x) 1 / (1 + exp(-mu     * (w_sam * x - theta)))  # Samuilik f2^-
h_dec <- function(x) theta^n / (x^n + theta^n)                     # Hill h^-

x <- seq(0, 10, length.out = 600)

## ---- plot ------------------------------------------------------------
# png("Coparaison_3.png", width = 1400, height = 950, res = 160)
par(mar = c(4.5, 5, 2, 1))
plot.default(x, f_our(x), type = "l", col = "blue", lwd = 2.5, las = 1,
             xlab = "Repressor concentration  x",
             ylab = "Repression response",
             ylim = c(0, 1))
lines(x, h_dec(x), col = "forestgreen", lwd = 2.5, lty = 2)
lines(x, f_sam(x), col = "red", lwd = 2.5)
abline(h = 0.5, col = "gray", lty = "dotted", lwd = 1.2)
abline(v = theta, col = "gray", lty = "dotted", lwd = 1.2)
points(theta, 0.5, pch = 19, cex = 1.1)             # shared midpoint of f_our and Hill
legend("topright",
       legend = c(expression(paste("our  ", f^"-", "(x, ", theta, ", ", lambda, "),  w = +1")),
                  expression(paste("Hill  ", h^"-", "(x, ", theta, ", n)")),
                  expression(paste("Samuilik  ", f[2]^"-", "(x, ", theta, ", ", mu, "),  w = -1"))),
       col = c("blue", "forestgreen", "red"), lty = c(1, 2, 1), lwd = 2.5, bg = "white")
# dev.off()

## ---- console summary -------------------------------------------------
cat(sprintf("x_c (our)  = theta/w = %+.0f   f^-(3)   = %.3f\n", theta / w_our, f_our(3)))
cat(sprintf("x_c (Sam.) = theta/w = %+.0f   f2^-(0)  = %.4f (decreasing, pinned near 0)\n",
            theta / w_sam, f_sam(0)))
cat(sprintf("Hill h^-(0) = %.3f,  h^-(3) = %.3f\n", h_dec(0), h_dec(3)))
