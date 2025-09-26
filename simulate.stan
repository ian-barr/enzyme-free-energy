
// The input data is a series of vectors of length 2.
data {
  vector[2] kf;
  vector[2] kr;
  vector[2] Kmf;
  vector[2] Kmr;
  vector[2] Dkf;
  vector[2] DKmr;
  vector[2] DKmf;
  vector[2] Dkr;
  vector[2] slope;
  vector[2] theta;
  vector[2] kez;
  vector[2] Keq;
}

transformed data{
  real<lower=0> kfp;
  real<lower=0> krp;
  real<lower=0> Spf;
  real<lower=0> Spr;
  kfp = kf[1];
  krp = kr[1];
  Spf = kf[1]/Kmf[1];
  Spr = kr[1]/Kmr[1];
}
// The parameters accepted by the model. 
// 
parameters {
  real<lower=0, upper=1e3> K1;
  real<lower=kfp, upper=1e12> k2;
  real<lower=kfp, upper=1e12> k3;
  real<lower=kfp, upper=1e12> k4;
  real<lower=krp, upper=1e12> k_1; // reverse constants
  real<lower=krp, upper=1e12> k_2;
  real<lower=krp, upper=1e12> k_3;
  real<lower=0, upper=1e3> K4;
  real<lower=0, upper=500> Dk2;// intrinsic KIE For
  real<lower=0, upper=500> Dk_3;// intrinsic KIE Rev
  real<lower=0> bet;
  // real<lower=0> Mu2;
  // real<lower=0> SD;
}
transformed parameters{
  real<lower=0> mukf;
  real<lower=0> mukr;
  real<lower=0> muKmf;
  real<lower=0> muKmr;
  real<lower=0> muDkf;
  real<lower=0> muDkr;
  real<lower=0> muDKmf;
  real<lower=0> muDKmr;
  real<lower=0> mutheta;
  real<lower=0> muslope;
  real<lower=0> mukez;
  real<lower=0> muKeq;
  mutheta = (k_2*(1+k_3/k4))/(k3*(1+k2/k_1));
  mukf = k2/(1 + k2/k4 + ((k2+k_2)/k3)*(1 + k_3/k4)) ;
  mukr = k_3/(1 + k_3/k_1 + ((k3+k_3)/k_2)*(1 + k2/k_1)) ;
  mukez = 1.0/(k_2/k2 + k3/k_3 +1);
  muslope = ( k_1*k_2*k_3 + k2*k3*k4)/(k_1*k_2*k_3 + k2*k3*k4 + k_1*k3*k4 + k_1*k_2*k4);
  muKmf = (1/K1)*(1 + k2/k_1 +(k_2/k3)*(1+k_3/k4))/(1 + k2/k4 + ((k2 + k_2)/k3)*(1+k_3/k4));//
  muKmr = (K4)*(1 + k_3/k4 +(k3/k_2)*(1+k2/k_1))/(1 + k_3/k_1 + ((k3 + k_3)/k_2)*(1+k2/k_1));//
  muDkf = ( Dk2 + (k2/k3)*(1+k_3/k4) + k2/k4 + Dk2*(k_2/k3)*(1+k_3/k4)) / ( 1 + (k2/k3)*(1+k_3/k4) + k2/k4 + (k_2/k3)*(1+k_3/k4));//
  muDkr = ( Dk_3 + (k_3/k_2)*(1+k2/k_1) + k_3/k_1 + Dk_3*(k3/k_2)*(1+k2/k_1))/( 1 + (k_3/k_2)*(1+k2/k_1) + k_3/k_1 + (k3/k_2)*(1+k2/k_1));//
  muDKmf = ( Dk2 + (k2/k_1) + Dk2*(k_2/k3)*(1+k_3/k4)) / ( 1 + (k2/k_1) + (k_2/k3)*(1+k_3/k4));//
  muDKmr = ( Dk_3 + (k_3/k4) + Dk_3*(k3/k_2)*(1+k2/k_1))/(  1 + (k_3/k4) + (k3/k_2)*(1+k2/k_1));//
  muKeq = (K1*k2*k3*K4)/(k_2*k_3);
}
// The model to be estimated. 
// 
// 
model {
  kf[1] ~ normal(mukf, kf[2]);
  kr[1] ~ normal(mukr, kr[2]);
  Kmf[1] ~ normal(muKmf, Kmf[2]);
  Kmr[1] ~ normal(muKmr, Kmr[2]);
  slope[1] ~ normal(muslope, slope[2]);
  theta[1] ~ normal(mutheta, theta[2]);
  kez[1]  ~ normal(mukez, kez[2]);
  Dkf[1] ~  normal(muDkf, Dkf[2]);
  Dkr[1] ~ normal(muDkr, Dkr[2]);
  DKmf[1] ~  normal(muDKmf, DKmf[2]);
  DKmr[1] ~ normal(muDKmr, DKmr[2]);
  Keq[1] ~ normal(muKeq, Keq[2]);
  // priors for ks

Dk2 ~ lognormal(1,0.5);
Dk_3 ~ lognormal(1,0.5);
K4 ~ exponential(.025);
k2 ~ exponential(bet);
k3 ~ exponential(bet);
k4 ~ exponential(bet);
k_1 ~ exponential(bet);
k_2 ~ exponential(bet);
k_3 ~ exponential(bet);
K4 ~ exponential(.025);
bet ~ gamma(1,1);

}
generated quantities{
  real dG1;
  real dG2;
  real dG3;
  real dG4;
  real R;
  real T;
  real k1;
  real k_4;
  real Hald;
  real K2;
  real K3;
  R = 1.987204e-3;//kcal/K.mol
  T = 298;//degrees K
  dG1 = -R*T*log(K1);
  dG2 = -R*T*log(k2/k_2);
  dG3 = -R*T*log(k3/k_3);
  dG4 = -R*T*log(K4);
  k1 = K1*k_1;
  k_4 = k4/K4;
  Hald = mukf*muKmr/(mukr*muKmf);
  K2 = k2/k_2;
  K3 = k3/k_3;
}
