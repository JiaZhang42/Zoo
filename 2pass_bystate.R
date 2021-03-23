rm(list = ls())
library(tidyverse)
library(lubridate)
library(locpol)
library(abind)
library(glmnet)


# load state variables, factors, and stock returns
state1 <- read_csv('~/T10Y2Y.csv')
pub_factors <- read_csv('~/pubfactors19592017.csv', col_types = cols(DATE = col_date('%Y%m%d'))) %>%
  mutate(DATE = floor_date(DATE, 'month'))
num_factors <- read_csv('~/numfactors.csv', col_types = cols(DATE = col_date('%Y%m%d'))) %>%
  mutate(DATE = floor_date(DATE, 'month'))
dum_factors <- read_csv('~/dumfactors.csv', col_types = cols(DATE = col_date('%Y%m%d'))) %>%
  mutate(DATE = floor_date(DATE, 'month'))
stock_returns <- read_csv('~/stockreturnsm19572019.csv', col_types = cols(DATE = col_date('%Y%m%d'))) 


# state variables
state2 <- state1 %>% 
  transmute(DATE = floor_date(DATE, 'month'), STATE = T10Y2Y - mean(T10Y2Y)) %>% 
  mutate(STATE = lag(STATE)) %>% na.omit()
# state - monthly average, stock/factor returns - whole month, should lag state by 1

# factors
factors <- state2 %>% 
  inner_join(pub_factors, by = 'DATE') %>% 
  inner_join(num_factors, by = 'DATE') %>% 
  inner_join(dum_factors, by = 'DATE')


# NA values in factors
# num_of_na <- function(x){
#   sum(is.na(x))
# }
# apply(factors, 2, num_of_na)
colSums(is.na(factors))
# HMLill, HMLdivo include NAs

#delete rows in state without factor values
state <- state2 %>% 
  semi_join(factors, by = 'DATE')

#sorted state value without date
sorted_state <- state %>% 
  arrange(STATE) %>% 
  select(STATE) %>% 
  distinct()

# stock returns
returns <- stock_returns %>%
  filter(!(RET %in% c('C', 'B', NA))) %>% 
  mutate(RET = as.numeric(RET), DATE = floor_date(DATE, 'month')) %>% 
  group_by(PERMNO) %>% 
  mutate(RET = RET - mean(RET, na.rm = T)) %>% 
  pivot_wider(names_from = PERMNO, values_from = RET)

# remove raw datasets
rm(state1, state2, stock_returns, dum_factors, num_factors, pub_factors)



# local linear estimation using `locpol` package (plug-in bindwidth):
ll_estimation <- function(x, y, xval = x, kernel = gaussK){
  plug_bw <- pluginBw(x, y, deg = 1, kernel)
  locPolSmootherC(x, y, xval, plug_bw, deg = 1, kernel)$beta0
}


# estimate rhat_stocks: T by n
permno <- colnames(returns)[-1]
rhat_stocks <- sorted_state
for (i in 1:length(permno)) {
  stock_i <- returns[,c(1, i+1)] %>% 
    na.omit()
  dat <- inner_join(state, stock_i, by = 'DATE') %>% 
    rename(RET = 3) %>% 
    group_by(STATE) %>% 
    summarise(RET = mean(RET))
  tryCatch({
    rhat_stock_i <- tibble(STATE = dat$STATE, !!quo_name(permno[i]) := ll_estimation(dat$STATE, dat$RET))
    rhat_stocks <- rhat_stocks %>% 
      left_join(rhat_stock_i, by = 'STATE')
  }, error = function(e) cat("ERROR :",conditionMessage(e), "\n")
  )
}



# estimate vhat_factors: T by (p+d)
factor_name <- colnames(factors)[-c(1,2)]
vhat_factors <- sorted_state
for (j in 1:length(factor_name)) {
  factor_j <- factors[,c(2, j+2)] %>%
    na.omit() 
  dat <- factor_j %>% 
    transmute(STATE, RET = .[[2]] - mean(.[[2]])) %>% 
    group_by(STATE) %>% 
    summarise(RET = mean(RET))
  vhat_factor_j <- tibble(STATE = dat$STATE, !!quo_name(factor_name[j]) := ll_estimation(dat$STATE, dat$RET))
  vhat_factors <- vhat_factors %>% left_join(vhat_factor_j, by = 'STATE')
}

saveRDS(rhat_stocks, "~/new_bystate/rhat_stocks.rds")
saveRDS(vhat_factors, "~/new_bystate/vhat_factors.rds")

###################################
# estimate rv(s_t): n by (p+d) by T
###################################
# rhat_stocks <- read_csv('D:/zoo/Data/rhat_stocks.csv')
# vhat_factors <- read_csv('D:/zoo/Data/vhat_factors.csv')
rhat_stocks <- readRDS("D:/zoo/Data/results0221/rhat_stocks.rds")
vhat_factors <- readRDS("D:/zoo/Data/results0221/vhat_factors.rds")
permno_new <- colnames(rhat_stocks)[-1]
returns_new1 <- returns[,c('DATE', permno_new)]

returns_new <- state %>% 
  left_join(returns_new1, by = 'DATE') %>% 
  group_by(STATE) %>% 
  summarise_at(vars(1:length(permno_new)+1), mean, na.rm = T)
factors_new <- state %>% 
  left_join(factors[-2], by = 'DATE') %>% 
  group_by(STATE) %>% 
  summarise_at(vars(1:length(factor_name)+1), mean, na.rm = T)
rm(returns_new1)

rv <- NULL
for (s in 1:dim(sorted_state)[1]) {
  mat_s <- t(as.matrix(returns_new[s,-1])) %*% as.matrix(factors_new[s,-1])
  rv <- abind(rv, mat_s, along = 3)
}
#dimnames(rv)[[3]] <- as.character(state$DATE)


##############################################
# estimate rvhat(s_t): n by (p+d) by T
##############################################
rvhat <- array(NA, dim = dim(rv))
for (i in 1:length(permno_new)) {
  for (j in 1:length(factor_name)) {
    rv_ij <- tibble(STATE = sorted_state$STATE, RET = rv[i, j, ]) %>% na.omit()
    tryCatch({
      rvhat_ij <- rv_ij %>% 
        mutate(RET = ll_estimation(rv_ij$STATE, rv_ij$RET))
      rvhat_ij_full <- sorted_state %>% left_join(rvhat_ij, by = 'STATE')
      rvhat[i, j, ] <- rvhat_ij_full$RET
    }, error = function(e) cat("ERROR :",conditionMessage(e), "\n")
    )
  }
}

##############################################
# estimate rhat(s_t)vhat(s_t): n by (p+d) by T
##############################################
rhatvhat <- NULL
for (s in 1:dim(sorted_state)[1]) {
  mat_s <- t(as.matrix(rhat_stocks[s,-1])) %*% as.matrix(vhat_factors[s,-1])
  rhatvhat <- abind(rhatvhat, mat_s, along = 3)
}
# dimnames(rhatvhat)[[3]] <- as.character(state$DATE)

# compute covhat(s_t): n by (p+d) by T
covhat <- rvhat - rhatvhat

# saveRDS(rvhat, "~/new_bystate/rvhat.rds")
# saveRDS(rhatvhat, "~/new_bystate/rhatvhat.rds")
# saveRDS(covhat, "~/new_bystate/covhat.rds")

###check
object.size(rv)
object.size(rvhat)
object.size(rhatvhat)
object.size(covhat)

rv[1,1,100:200]
rhatvhat[1,1,100:200]
as_tibble(rhatvhat[1,1,]) %>% na.omit()
as_tibble(rvhat[1,1,]) %>% na.omit()
as_tibble(covhat[1,1,]) %>% na.omit()

############### the 2nd pass #################
## DS LASSO single-point
rhat_stocks <- read_csv('D:/zoo/Data/rhat_stocks.csv')
covhat <- readRDS("D:/zoo/Data/covhat.rds")
# factor_name <- dimnames(covhat)[[2]]

grid <- 10^seq(0, -2, length=100)

selection <- function(x, y, lambda = NULL, seednum = 1){
  set.seed(seednum)
  cv.out <- cv.glmnet(x, y, alpha = 1, lambda = lambda)
  #  plot(cv.out)
  out <- glmnet(x, y, alpha = 1, lambda = lambda)
  lasso.coef <- predict(out, s = cv.out$lambda.min, type = 'coefficients')
  dimnames(lasso.coef)[[1]][lasso.coef@i+1][-1]
}


# factor_d of interest
d <- 2 
d <- 3
d <- 4
d <- 7
d <- 8
d <- 9
d <- 10
d <- 11
d <- 12

ds_selection <- matrix(NA, dim(covhat)[2], dim(covhat)[3])
ss_selection <- ds_selection
ps_ols <- vector("list", dim(covhat)[3])
# names(ps_ols) <- dimnames(covhat)[[3]]
factor_d_est <- matrix(NA, 6, dim(covhat)[3]) #estimate, se, t, p, CI_left, CI_right


for (s in 1:dim(covhat)[3]) {
  rhat_s <- t(rhat_stocks[s, -1]) #vector of r
  covhhat_s <- covhat[ , -d, s] #matrix of h
  covghat_s <- covhat[ , d, s] #vector of g
  dat_rs <- cbind(rhat_s, covhhat_s) %>%
    .[,colSums(is.na(.)) != nrow(.)] %>% 
    na.omit() # for 1st selection
  dat_gs <- covhat[, , s] %>% 
    .[,colSums(is.na(.)) != nrow(.)] %>% 
    na.omit() # for 2nd selection
  # whole col/row NAs are generated by NA in factor/stock returns at state s
  # other NAs are generated by local polynomial estimation (maybe)
  if (dim(dat_rs)[1] > 0) {
    selection1 <- selection(dat_rs[,-1], dat_rs[,1], lambda = grid)
    selection2 <- selection(dat_gs[,-d], dat_gs[,d], lambda = grid)
    selection_full <- as_tibble(c(selection1, selection2)) %>% distinct()
    ss_selection[,s] <- factor_name %in% selection1
    ds_selection[,s] <- factor_name %in% pull(selection_full) #output
    dat_s <- as_tibble(cbind(dat_rs[,1], dat_gs)) %>% 
      rename(rhat_s = V1) %>% 
      select(rhat_s, d+1, pull(selection_full))
    ps_ols[[s]] <- lm(rhat_s~., data = dat_s) #output
    estimate_d <- as_tibble(summary(ps_ols[[s]])$coefficients[2,])
    ci_d <- as_tibble(t(confint(ps_ols[[s]], factor_name[d]))) %>% 
      rename(value = 1) #should equal estimate[1]-/+estimate[2]*qnorm(0.975)
    factor_d_est[,s] <- pull(bind_rows(estimate_d, ci_d)) #output
  }
}

saveRDS(ds_selection, "~/new_bystate/ds_selection.rds")
saveRDS(ss_selection, "~/new_bystate/ss_selection.rds")
saveRDS(ps_ols, "~/new_bystate/ps_ols.rds")
saveRDS(factor_d_est, "~/new_bystate/factor_d_est.rds")


# num of selected factors over state 改date
sselect_num <- tibble(STATE = sorted_state$STATE, selnum = apply(ss_selection, 2, sum))
sselect_num %>% 
  ggplot(aes(x = STATE, y = selnum)) + 
  geom_line() + 
  labs(title = 'Number of factors selected by single-selection lasso')

dselect_num <- tibble(STATE = sorted_state$STATE, selnum = apply(ds_selection, 2, sum))
dselect_num %>% 
  ggplot(aes(x = STATE, y = selnum)) + 
  geom_line() + 
  labs(title = 'Number of factors selected by double-selection lasso')

# times to be selected each factor
persistency <- tibble(factor = factor_name, sel_rate = apply(ds_selection, 1, sum)/dim(ds_selection)[2])
persistency %>% arrange(desc(sel_rate))

ds_selection_neg <- ds_selection[,sorted_state<0]
persistency_neg <- tibble(factor = factor_name, sel_rate = apply(ds_selection_neg, 1, sum)/dim(ds_selection_neg)[2])
persistency_neg %>% arrange(desc(sel_rate))

ds_selection_pos <- ds_selection[,sorted_state>0]
persistency_pos <- tibble(factor = factor_name, sel_rate = apply(ds_selection_pos, 1, sum)/dim(ds_selection_pos)[2])
persistency_pos %>% arrange(desc(sel_rate))
# to save:
# - selection results at each t (ds_selection)
# - post-selection estimation results (ps_ols)
# - CI of factor d (factor_d_est)


#to solve:
# - num of selection (change cross-sectional CV to time-series): 10^(-5)
# - plot CI over time: confint(lm.object)
# - t = 25 error, covhhat_t %>% na.omit() has 0 entry
# HMLdivo (No.92 factor) NA value in covhat[, , t=25] because NA value in factors
# whole cols of NA in covhat are removed.


factor_d_state <- as_tibble(t(factor_d_est)) %>% 
  rename(estimate = V1, se = V2, t_value = V3, p_value = V4, CI_L = V5, CI_R = V6) %>% 
  mutate(STATE = sorted_state$STATE)

factor_d_state %>% 
  ggplot(aes(x = STATE, y = estimate)) + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_L, ymax = CI_R))

factor_d_state %>% 
  ggplot(aes(x = STATE, y = estimate)) + 
  geom_line() + 
  labs(title = paste('Estimation of SDF loadings of', factor_name[d]))

factor_d_state %>% 
  ggplot(aes(x = STATE, y = t_value)) + 
  geom_line() + 
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dashed', colour = 'blue') + 
  labs(title = paste('t-value of SDF loadings of', factor_name[d]))


