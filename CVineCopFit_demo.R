## Example: 10‑dimensional C‑vine with variable 10 as the central node
## Load necessary packages
library(rvinecopulib)
library(MASS)        # For multivariate normal data
library(VineCopula)  # Optional: for tau matrix & converting to R‑vine matrix

## Simulate data and convert to pseudo observations
set.seed(123)
n  <- 1000        # number of observations
d  <- 10          # dimension

# simple equicorrelated covariance matrix for the example
Sigma <- matrix(0.5, nrow = d, ncol = d)
diag(Sigma) <- 1

# simulate multivariate normal and convert to uniform margins
x <- MASS::mvrnorm(n, rep(0, d), Sigma)
u <- apply(x, 2, rank) / (n + 1)   # pseudo‑observations

## Specify C‑vine structure with the 10rd variable as central node manually 
# In cvine_structure() the first element is the root of the first tree
order <- 1:10     
c_vine_structure <- cvine_structure(order)
c_vine_structure

## Optional: Visualise the C‑vine structure
plot(c_vine_structure)   # star‑shaped tree with node 10 in the centre ?

## Fit the C‑vine copula model
# family_set = "par" chooses among many parametric copulas; you can restrict it.
c_vine_fit <- vinecop(u, structure = c_vine_structure)

## Inspect the fitted model
print(c_vine_fit)
print(summary(c_vine_fit))     # shows pair‑copulas and estimated parameters

## (Optional) Cross‑check with VineCopula:
# Compute the empirical Kendall’s τ matrix and check that variable 3 has high dependence
tau_mat <- TauMatrix(u) 
colSums(abs(tau_mat))    # sum of |τ| per variable; see central node selection:contentReference[oaicite:4]{index=4}

# Convert the cvine structure to VineCopula’s R‑vine matrix notation
rvine_matrix <- as_rvine_matrix(c_vine_structure)
rvine_matrix            # lower‑triangular matrix representation

# Observing different trees based on fitted Cvine model good to explore them in detail for the verification
plot(c_vine_fit, tree = 1, var_names = "legend", edge_labels = "family_tau")
plot(c_vine_fit, tree = 2, var_names = "legend", edge_labels = "family_tau")
plot(c_vine_fit, tree = 3, var_names = "legend", edge_labels = "family_tau")

# SUGGESTED ONE IS BELOW FOR OUR PROJECT
## Example: This part calculates the order based on the empirical kendall's tau matrix instead 

# Determine the C‑vine order yourself: The root node in a C‑vine is the variable with the strongest overall 
# dependencies to the other variables. In VineCopula this is defined as the variable with the largest column 
# sum in the empirical Kendall’s‑tau matrix. You can compute this order and pass it to cvine_structure()

# same simulated u: n×10 matrix of pseudo-observations
# Calculation of tau matrix based on data 
tau_mat <- VineCopula::TauMatrix(u)
# Looking at the colsums for making an order
sum_abs_tau <- colSums(abs(tau_mat))
order <- order(sum_abs_tau, decreasing = TRUE)   # highest sum first
# Creating the order of cvine structure based on above tau based ordering
c_vine_structure <- cvine_structure(order)
plot(c_vine_structure) # I presume, we can still get similar variables in our SIMD data set as star node!

# fit the C‑vine automatically
c_vine_fit <- vinecop(u, structure = c_vine_structure,
                      family_set = "all", selcrit = "aic")

# Looking at the fitted Cvine summary
summary(c_vine_fit)



