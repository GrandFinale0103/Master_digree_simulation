##---
## 0
## load packages

pacman::p_load(tidyverse, psych, mirt)
source('./necessary user functions.R')

## 0
##---

##---
## 1
## setting theta

theta <- rnorm(500) # !IV option (select distribution)

## 1
##---

##---
## 2
## setting item parameter

# a(discrimination) parameters
a_para <- rep(1, 20) 

# scoring parameters (reparameterized NRM model)
scoring_first_vec <- c(1,2,3,3.5) # !IV option (select parameter of category 4)
scoring_vec <- c(scoring_first_vec, rep(1:4, 19))
scoring <- matrix(scoring_vec, nrow = 20, byrow = TRUE) 

# b(boundary) parameters
boundary_first_vec <- c(-3, -1, 1, 0) # !IV option (select parameter of category 4)
boundary_vec <- c(boundary_first_vec, rep(1:4, 19))
boundary <- matrix(boundary_vec, nrow = 20, byrow = TRUE) 

# transform boundaries matrix to intercepts matrix
c_mat <- boundary_to_c(
  a = a_para,
  scoring = scoring,
  boundary = boundary,
  return_df = FALSE
)
c_mat

# 응답행렬 생성
resp <- simulate_nrm_response(
  theta = theta,
  a = a_para,
  scoring = scoring,
  c = c_mat,
  seed = sample(1:1000, 1)
)

resp