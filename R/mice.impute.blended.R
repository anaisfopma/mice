#####################
# MODIFIED FUNCTION #
#####################

# removed use.matcher argument, added weight (p) argument and xname argument
mice.impute.blended <- function(y, ry, x, wy = NULL, donors = 5L,
                            matchtype = 1L, ridge = 1e-05, blend = TRUE, p, xname, ...) {
  
  if (is.null(wy)) {
    wy <- !ry
  }
  x <- cbind(1, as.matrix(x))
  ynum <- y
  if (is.factor(y)) {
    ynum <- as.integer(y)
  }
  parm <- .norm.draw(ynum, ry, x, ridge = ridge, ...)
  
  if (matchtype == 0L) {
    
    # select the predictors for the donors
    x_donors <- x[ry, , drop = FALSE] 
    
    # obtain predicted values for the donors
    yhatobs <- x_donors %*% parm$coef
    
    # make x_donors a data frame so we can later add the distance measures
    x_donors <- as.data.frame(x_donors)
    
    # add the predicted values for the donors to this data frame
    x_donors$yhatobs <- yhatobs
    
    # select the predictors for the target
    x_target <- x[wy, , drop = FALSE]
    
    # obtain predicted value for the target
    yhatmis <- as.numeric(x_target %*% parm$coef)
    
  }
  
  # do the same for the other matchtypes
  if (matchtype == 1L) {
    x_donors <- x[ry, , drop = FALSE] 
    yhatobs <- x_donors %*% parm$coef
    x_donors <- as.data.frame(x_donors)
    x_donors$yhatobs <- yhatobs
    x_target <- x[wy, , drop = FALSE]
    yhatmis <- as.numeric(x_target %*% parm$beta)
  }
  
  if (matchtype == 2L) {
    x_donors <- x[ry, , drop = FALSE] 
    yhatobs <- x_donors %*% parm$beta
    x_donors <- as.data.frame(x_donors)
    x_donors$yhatobs <- yhatobs
    x_target <- x[wy, , drop = FALSE]
    yhatmis <- as.numeric(x_target %*% parm$beta)
  }
  
  # here we use code adapted from the .pmm.match function instead of matcher and matchindex
  if (blend == TRUE) {
    
    # calculate predictive distance and add it to the x_donors data frame
    x_donors$pd <- x_donors$yhatobs - yhatmis
    f <- x_donors$pd > 0
    a1 <- ifelse(any(f), min(x_donors$pd[f]), 1)
    x_donors$pd <- x_donors$pd + runif(length(x_donors$pd), 0, a1 / 10^10)
    
    # assign a rank to each donor for the predictive distance
    x_donors$pd_rank <- rank(x_donors$pd, ties.method = "random")
    
    # calculate Mahalanobis distance and add it to the x_donors data frame
    x_donors$md <- sqrt(mahalanobis(x_donors[ , xname], as.numeric(x_target[ , xname]), cov(x_donors[ , xname])))
    
    # assign a rank to each donor for the Mahalanobis distance
    x_donors$md_rank <- rank(x_donors$md, ties.method = "random")
    
    # calculate blended distance with weight = p
    x_donors$bd <- p * x_donors$pd_rank + (1-p) * x_donors$md_rank
    
    # select k closest matches 
    # if k = 1
    if (donors == 1) {
      return(y[which.min(x_donors$bd)])
    }
    # if k > 1
    donors <- min(x_donors$bd, donors)
    donors <- max(x_donors$bd, 1) 
    # I don't understand exactly what this code does: first we take the 5 donors with the smallest blended distance,
    # but then the largest value?
    ds <- sort.int(x_donors$bd, partial = donors)
    m <- sample(y[x_donors$bd <= ds[donors]], 1)
    # Should we sample one or take the mean instead? 
    
    # return imputed/predicted value for the target
    return(m)
    
  }
  
}

############################################
# TEST MODIFIED FUNCTION ON SIMULATED DATA #
############################################

# load packages
library(faux)
library(mice)

# simulate 1 dataset with 1000 observations of standardised height measurements at 
# birth, 1, 2, 3, 6, and 14 months, mean = 0, sd = 1, correlation = 0.7
set.seed(847)
dat <- rnorm_multi(1000, 6, 0, 1, 0.8, varnames = c("hgt_z_birth", "hgt_z_1mo", "hgt_z_2mo", "hgt_z_3mo", "hgt_z_6mo", "hgt_z_14mo"))

# set the height of the target at 14 months to NA, this is the missing value we want to impute/predict
dat[1,6] <- NA

# specify the predictors
xname <- c("hgt_z_birth", "hgt_z_1mo", "hgt_z_2mo", "hgt_z_3mo", "hgt_z_6mo")

# logical vector indicating whether a case is complete or not
r <- stats::complete.cases(dat[, xname])

# matrix with predictors for y
x <- dat[r, xname]

# vector with the missing value we want to impute/predict
y <- dat[r, "hgt_z_14mo"]

# logical vector indicating for each case if hgt_z_14mo is missing
ry <- !is.na(y)



# The mice.impute.pmm gives (relatively) similar predicted values every time:
set.seed(864)
mice.impute.pmm(y, ry, x)

set.seed(198)
mice.impute.pmm(y, ry, x)

set.seed(449)
mice.impute.pmm(y, ry, x)

# But the mice.impute.blended function gives values that are very different from each other:
set.seed(864)
mice.impute.blended(y, ry, x, p = 0.5, xname = xname)

set.seed(198)
mice.impute.blended(y, ry, x, p = 0.5, xname = xname)

mice.impute.blended(y, ry, x, p = 0.5, xname = xname)
set.seed(449)

# Even when p = 1, which means that it should be the same as pmm:
set.seed(864)
mice.impute.blended(y, ry, x, p = 1, xname = xname)

set.seed(198)
mice.impute.blended(y, ry, x, p = 1, xname = xname)

mice.impute.blended(y, ry, x, p = 1, xname = xname)
set.seed(449)

# For the last two seeds, the result is also exactly the same for p = 0.5 and p = 1


