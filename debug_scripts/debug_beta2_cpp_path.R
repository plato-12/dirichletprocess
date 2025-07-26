#!/usr/bin/env Rscript

library(dirichletprocess)

# Test beta2 C++ path issue
beta2Obj <- BetaMixture2Create()
betaObj <- BetaMixtureCreate()

testTheta <- list()
testTheta[[1]] <- array(0.5, dim=c(1,1,1))
testTheta[[2]] <- array(0.5, dim=c(1,1,1))
names(testTheta) <- c("mu", "nu")

cat("=== DEBUGGING BETA2 C++ PATH ===\n")

# Test 1: Regular beta object 
cat("1. Regular beta object:\n")
cat("   Class:", paste(class(betaObj), collapse=", "), "\n")

# Extract distribution type like the main function does
class_vals <- class(betaObj)
filter_vals <- class_vals != "list" & class_vals != "MixingDistribution"
dist_type <- class_vals[filter_vals][1]
cat("   Distribution type:", dist_type, "\n")

oldLik <- Likelihood(betaObj, c(0.1, 0.2), testTheta)
cat("   Result class:", class(oldLik), "\n")
cat("   Result dim:", dim(oldLik), "\n")
cat("   Result structure: "); str(oldLik)

# Test 2: beta2 object directly
cat("\n2. beta2 object directly:\n")
cat("   Class:", paste(class(beta2Obj), collapse=", "), "\n")

class_vals2 <- class(beta2Obj)
filter_vals2 <- class_vals2 != "list" & class_vals2 != "MixingDistribution"
dist_type2 <- class_vals2[filter_vals2][1]
cat("   Distribution type:", dist_type2, "\n")

newLik <- Likelihood(beta2Obj, c(0.1, 0.2), testTheta)
cat("   Result class:", class(newLik), "\n")
cat("   Result dim:", dim(newLik), "\n")
cat("   Result structure: "); str(newLik)

# Test 3: Manually converted object
cat("\n3. Manually converted beta2->beta object:\n")
temp_mdObj <- beta2Obj
class(temp_mdObj) <- c("list", "beta", "nonconjugate")
cat("   Class:", paste(class(temp_mdObj), collapse=", "), "\n")

class_vals3 <- class(temp_mdObj)
filter_vals3 <- class_vals3 != "list" & class_vals3 != "MixingDistribution"
dist_type3 <- class_vals3[filter_vals3][1]
cat("   Distribution type:", dist_type3, "\n")

testLik <- Likelihood(temp_mdObj, c(0.1, 0.2), testTheta)
cat("   Result class:", class(testLik), "\n")
cat("   Result dim:", dim(testLik), "\n")
cat("   Result structure: "); str(testLik)

cat("\n=== COMPARISON ===\n")
cat("beta == beta2:", isTRUE(all.equal(oldLik, newLik)), "\n")
cat("beta == converted:", isTRUE(all.equal(oldLik, testLik)), "\n")