#------
library(devtools)
load_all()
X<- generate_X(n=50, p=10, p_b=0, cov_type="cov_equi", rho=0.5)
#Xk=create.seq(X)
Xk=create.pc(X,n_ko = 2,ncomp=2)
Xk=create.sparse_seq(X)
Xk
colnames(Xk[[1]])


#-----
set.seed(2024)
X = generate_X(n=80,p=100,rho=0.3)
y <- generate_y(X, p_nn=10, a=3)
Xk = create.shrink_Gaussian(X = X, n_ko = 10)
res1 = knockoff.filter(X, y, Xk, statistic = stat.glmnet_coefdiff,
                       offset = 1, fdr = 0.1)
res1
perf_eval(res1$shat,1:10,11:100)

# Logistic Regression
lp <- generate_lp(X, p_nn=10, a=3)
pis <- plogis(lp)
Y <- factor(rbinom(n, 1, pis))
res2 = knockoff.filter(X, Y, Xk, statistic = stat.glmnet_coefdiff,
                       family = 'binomial', offset = 0, fdr = 0.2)
res2
perf_eval(res2$shat,1:10,11:100)
