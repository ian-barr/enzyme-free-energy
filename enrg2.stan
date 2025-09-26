
// The input data is a series of vectors of length 2.
data {
  vector[2] kl;
  vector[2] kr;
  vector[2] Kml;
  vector[2] Kmr;
  vector[2] Dkl;
  vector[2] DKmr;
  vector[2] DKml;
  vector[2] Dkr;
  vector[2] slope;
  vector[2] theta;
  vector[2] kez;
}
transformed data{
  vector[2] KIEl;
  vector[2] KIEr;
  KIEr[1] = (DKmr[1]-1.0)/(Dkr[1]-1.0);
  KIEl[1] = (DKml[1]-1.0)/(Dkl[1]-1.0);
  KIEr[2] = KIEr[1]*sqrt((DKmr[2]/DKmr[1])^2 + (Dkr[2]/Dkr[1])^2);
  KIEl[2] = KIEl[1]*sqrt((DKml[2]/DKml[1])^2 + (Dkl[2]/Dkl[1])^2);
}
// The parameters accepted by the model. 
// 
parameters {
  real<lower=0, upper=1e9> k1;
  real<lower=0, upper=1e12> k2;
  real<lower=0, upper=1e12> k3;
  real<lower=0, upper=1e12> k4;
  real<lower=0, upper=1e12> k_1; // reverse constants
  real<lower=0, upper=1e12> k_2;
  real<lower=0, upper=1e12> k_3;
  real<lower=0, upper=1e9> k_4;
  real<lower=0> Dk2;// intrinsic KIE For
  real<lower=0> Dk_3;// intrinsic KIE Rev
}
transformed parameters{
  real<lower=0> mukl;
  real<lower=0> mukr;
  real<lower=0> muKml;
  real<lower=0> muKmr;
  real<lower=0> muDkl;
  real<lower=0> muDkr;
  real<lower=0> muDKml;
  real<lower=0> muDKmr;
  real<lower=0> mutheta;
  real<lower=0> muslope;
  real<lower=0> mukez;
  real<lower=0> muKIEr;
  real<lower=0> muKIEl;
  mutheta = (k_2*(1+k_3/k4))/(k3*(1+k2/k_1));
  mukl = k2/(1 + k2/k4 + ((k2+k_2)/k3)*(1 + k_3/k4)) ;
  mukr = k_3/(1 + k_3/k_1 + ((k3+k_3)/k_2)*(1 + k2/k_1)) ;
  mukez = (k_2/k2 + k3/k_3 +1)^(-1);
  muslope = ( k_1*k_2*k_3 + k2*k3*k4)/(k_1*k_2*k_3 + k2*k3*k4 + k_1*k3*k4 + k_1*k_2*k4);
  muKml = (k_1/k1)*(1 + k2/k_1 +(k_2/k3)*(1+k_3/k4))/(1 + k2/k4 + ((k2 + k_2)/k3)*(1+k_3/k4));//
  muKmr = (k4/k_4)*(1 + k_3/k4 +(k3/k_2)*(1+k2/k_1))/(1 + k_3/k_1 + ((k3 + k_3)/k2)*(1+k2/k_1));//
  muDkl = ( Dk2 + (k2/k3)*(1+k_3/k4) + k2/k4 + Dk2*(k_2/k3)*(1+k_3/k4)) / ( 1 + (k2/k3)*(1+k_3/k4) + k2/k4 + (k_2/k3)*(1+k_3/k4));//
  muDkr = ( Dk_3 + (k_3/k_2)*(1+k2/k_1) + k_3/k_1 + Dk_3*(k3/k_2)*(1+k2/k_1))/( 1 + (k_3/k_2)*(1+k2/k_1) + k_3/k_1 + (k3/k_2)*(1+k2/k_1));//
  muDKml = ( Dk2 + (k2/k_1) + Dk2*(k_2/k3)*(1+k_3/k4)) / ( 1 + (k2/k_1) + (k_2/k3)*(1+k_3/k4));//
  muDKmr = ( Dk_3 + (k_3/k4) + Dk_3*(k3/k_2)*(1+k2/k_1))/(  1 + (k_3/k4) + (k3/k_2)*(1+k2/k_1));//
  muKIEl = ( 1 + (k2/k3)*(1+k_3/k4) + k2/k4 + (k_2/k3)*(1+k_3/k4))/(1 + k2/k_1+(k_2/k3)*(1+k_3/k4) );//
  muKIEr = ( 1 + (k_3/k_2)*(1+k2/k_1) + k_3/k_1 + (k3/k_2)*(1+k2/k_1))/(  1 + k_3/k4 +  (k3/k_2)*(1+k2/k_1) );//
}
// The model to be estimated. 
// 
// 
model {
  kl[1] ~ normal(mukl, kl[2]);
  kr[1] ~ normal(mukr, kr[2]);
  Kml[1] ~ normal(muKml, Kml[2]);
  Kmr[1] ~ normal(muKmr, Kmr[2]);
  slope[1] ~ normal(muslope, slope[2]);
  theta[1] ~ normal(mutheta, theta[2]);
  kez[1]  ~ normal(mukez, kez[2]);
  Dkl[1] ~  normal(muDkl, Dkl[2]);
  Dkr[1] ~ normal(muDkr, Dkr[2]);
  DKml[1] ~  normal(muDKml, DKml[2]);
  DKmr[1] ~ normal(muDKmr, DKmr[2]);
  KIEr[1] ~  normal(muKIEr, KIEr[2]);
  KIEl[1] ~ normal(muKIEl, KIEl[2]);
  // priors for ks
  k1 ~  lognormal(log(1e6),2);
  k2 ~ lognormal(log(1e6),2);
  k3 ~ lognormal(log(1e6),2);
  k4 ~ lognormal(log(1e6),2);
  k_1 ~ lognormal(log(1e6),2);
  k_2 ~ lognormal(log(1e6),2);
  k_3 ~ lognormal(log(1e6),2);
  k_4 ~ lognormal(log(1e6),2);
  Dk2 ~ lognormal(1,0.5);
  Dk_3 ~ lognormal(1,0.5);
  // k1 ~  uniform(0, 1e9);
  // k2 ~ uniform(0, 1e12);
  // k3 ~ uniform(0, 1e12);
  // k4 ~ uniform(0, 1e12);
  // k_1 ~ uniform(0, 1e12); 
  // k_2 ~ uniform(0, 1e12);
  // k_3 ~ uniform(0, 1e12);
  // k_4 ~ uniform(0, 1e9);
  // Dk2 ~ lognormal(1,0.5);
  // Dk_3 ~ lognormal(1,0.5);
}
generated quantities{
  
  real dG1;
  real dG2;
  real dG3;
  real dG4;
  real R;
  real T;
  R = 1.987204e-3;//kcal/K.mol
  T = 298;//degrees K
  dG1 = -R*T*log(0.001*k1/k_1);
  dG2 = -R*T*log(k2/k_2);
  dG3 = -R*T*log(k3/k_3);
  dG4 = -R*T*log(k4/(0.001*k_4));
  
}
