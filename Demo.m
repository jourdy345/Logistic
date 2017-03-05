n = 200;
p = 15;
beta_true = randn(p,1);
X = randn(n,p)*2.5;
logisticNoise = log(exp(exprnd(1,[n,1]))-1);
Y = (X*beta_true+logisticNoise) > 0;
[beta_save,omega_save] = LogisticPolyaGamma(Y,X,5,6000,10000);
[beta2_save,ystar_save,nu_save] = LogisticKolmogorovSmirnov(Y,X,5,6000,10000);
