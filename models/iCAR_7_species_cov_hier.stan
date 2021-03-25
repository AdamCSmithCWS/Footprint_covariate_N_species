// This is a Stan implementation of a route-level slope model
// with an explicitly spatial prior structure on the 
// random effects for route-level intercepts and slopes
// AND it has route-level covariates on the slopes and intercepts that vary by species


//iCAR function
functions {
  real icar_normal_lpdf(vector bb, int nroutes, int[] node1, int[] node2) {
    return -0.5 * dot_self(bb[node1] - bb[node2])
      + normal_lpdf(sum(bb) | 0, 0.001 * nroutes); //soft sum to zero constraint on bb
 }
}


data {
   int<lower=1> nspecies;
   int<lower=1> ncounts; // sum of all ncounts variables used to define dimensions of log_lik calculations
   int<lower=1> nyears;
   int<lower=1> nobservers;
   int<lower=1> fixedyear; // centering value for years
   int<lower=0> nrts_cov; //total number of routes across all species
   vector[nrts_cov] rt_cov_dat; //single covariate value
   //matrix[nrts_cov,ncovs+1] rt_cov_mat_dat;//matrix of route indicators and covariate values
   // covariates start at column 2.


  //species specific components
  
//by species     
//species 1
  int<lower=1> nroutes1;
  int<lower=1> ncounts1;
  int<lower=0> count1[ncounts1];              // count observations
  int<lower=0> year1[ncounts1];              // year observations
  int<lower=0> route1[ncounts1];              // route observations - nested in species
  int<lower=0> observer1[ncounts1];              // observers 
  int<lower=0> route_cov1[ncounts1];              // covariate route observations
  int<lower=0> cov_routes1[nroutes1];              // covariate route indicators with species

 // spatial neighbourhood information - species 1
  int<lower=1> N_edges1;
  int<lower=1, upper=nroutes1> node11[N_edges1];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nroutes1> node21[N_edges1];  // and node1[i] < node2[i]


//species 2
  int<lower=1> nroutes2;
  int<lower=1> ncounts2;
  int<lower=0> count2[ncounts2];              // count observations
  int<lower=0> year2[ncounts2];              // year observations
  int<lower=0> route2[ncounts2];              // route observations - nested in species
  int<lower=0> route_cov2[ncounts2];              // rcovariate route observations
  int<lower=0> cov_routes2[nroutes2];              // covariate route indicators with species
  
  int<lower=0> observer2[ncounts2];              // observers 
 // spatial neighbourhood information - species 2
  int<lower=1> N_edges2;
  int<lower=1, upper=nroutes2> node12[N_edges2];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nroutes2> node22[N_edges2];  // and node1[i] < node2[i]


//species 3
  int<lower=1> nroutes3;
  int<lower=1> ncounts3;
  int<lower=0> count3[ncounts3];              // count observations
  int<lower=0> year3[ncounts3];              // year observations
  int<lower=0> route3[ncounts3];              // route observations - nested in species
  int<lower=0> route_cov3[ncounts3];              // covariate route observations
  int<lower=0> cov_routes3[nroutes3];              // covariate route indicators with species
  int<lower=0> observer3[ncounts3];              // observers 
 // spatial neighbourhood information - species 3
  int<lower=1> N_edges3;
  int<lower=1, upper=nroutes3> node13[N_edges3];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nroutes3> node23[N_edges3];  // and node1[i] < node2[i]


//species 4
  int<lower=1> nroutes4;
  int<lower=1> ncounts4;
  int<lower=0> count4[ncounts4];              // count observations
  int<lower=0> year4[ncounts4];              // year observations
  int<lower=0> route4[ncounts4];              // route observations - nested in species
  int<lower=0> route_cov4[ncounts4];              // covariate route observations
  int<lower=0> cov_routes4[nroutes4];              // covariate route indicators with species
  int<lower=0> observer4[ncounts4];              // observers 
 // spatial neighbourhood information - species 4
  int<lower=1> N_edges4;
  int<lower=1, upper=nroutes4> node14[N_edges4];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nroutes4> node24[N_edges4];  // and node1[i] < node2[i]


//species 5
  int<lower=1> nroutes5;
  int<lower=1> ncounts5;
  int<lower=0> count5[ncounts5];              // count observations
  int<lower=0> year5[ncounts5];              // year observations
  int<lower=0> route5[ncounts5];              // route observations - nested in species
  int<lower=0> route_cov5[ncounts5];              // covariate route observations
  int<lower=0> cov_routes5[nroutes5];              // covariate route indicators with species
  int<lower=0> observer5[ncounts5];              // observers 
 // spatial neighbourhood information - species 5
  int<lower=1> N_edges5;
  int<lower=1, upper=nroutes5> node15[N_edges5];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nroutes5> node25[N_edges5];  // and node1[i] < node2[i]

//species 6
  int<lower=1> nroutes6;
  int<lower=1> ncounts6;
  int<lower=0> count6[ncounts6];              // count observations
  int<lower=0> year6[ncounts6];              // year observations
  int<lower=0> route6[ncounts6];              // route observations - nested in species
  int<lower=0> route_cov6[ncounts6];              // covariate route observations
  int<lower=0> cov_routes6[nroutes6];              // covariate route indicators with species
  int<lower=0> observer6[ncounts6];              // observers 
 // spatial neighbourhood information - species 6
  int<lower=1> N_edges6;
  int<lower=1, upper=nroutes6> node16[N_edges6];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nroutes6> node26[N_edges6];  // and node1[i] < node2[i]


//species 7
  int<lower=1> nroutes7;
  int<lower=1> ncounts7;
  int<lower=0> count7[ncounts7];              // count observations
  int<lower=0> year7[ncounts7];              // year observations
  int<lower=0> route7[ncounts7];              // route observations - nested in species
  int<lower=0> route_cov7[ncounts7];              // covariate route observations
  int<lower=0> cov_routes7[nroutes7];              // covariate route indicators with species
  int<lower=0> observer7[ncounts7];              // observers 
 // spatial neighbourhood information - species 7
  int<lower=1> N_edges7;
  int<lower=1, upper=nroutes7> node17[N_edges7];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nroutes7> node27[N_edges7];  // and node1[i] < node2[i]



}

parameters {
//by species     
vector[nroutes1] alpha_raw1;
vector[nroutes1] beta_raw1;
vector[ncounts1] noise_raw1;

vector[nroutes2] alpha_raw2;
vector[nroutes2] beta_raw2;
vector[ncounts2] noise_raw2;

vector[nroutes3] alpha_raw3;
vector[nroutes3] beta_raw3;
vector[ncounts3] noise_raw3;

vector[nroutes4] alpha_raw4;
vector[nroutes4] beta_raw4;
vector[ncounts4] noise_raw4;

vector[nroutes5] alpha_raw5;
vector[nroutes5] beta_raw5;
vector[ncounts5] noise_raw5;

vector[nroutes6] alpha_raw6;
vector[nroutes6] beta_raw6;
vector[ncounts6] noise_raw6;


vector[nroutes7] alpha_raw7;
vector[nroutes7] beta_raw7;
vector[ncounts7] noise_raw7;




   vector[nobservers] obs_raw ;             // observer effects by species - consider pooling across all species
 

  real B_cov;         // covariate effects on slope
  real A_cov;         // covariate effects on intercept

  vector[nspecies] b_cov_raw; // species-level random effects on slopes
  vector[nspecies] a_cov_raw; // species-level random effects on slopes

  real<lower=0> sdb_cov;    // sd of over-dispersion 
  real<lower=0> sda_cov;    // sd of over-dispersion 

//real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  vector<lower=0> [nspecies]  sdnoise;    // sd of over-dispersion 
  vector<lower=0> [nspecies]  sdbeta;    // sd of slopes 
  vector<lower=0> [nspecies]  sdalpha;    // sd of intercepts
  vector[nspecies]  BETA;    // species mean slopes
  vector[nspecies]  ALPHA;    // species mean intercepts

  
}

// transformed parameters { 



//   }
  
model {

vector[nobservers] obs;             // observer effects by species - consider pooling across all species
matrix[nrts_cov,nspecies] cov_b_adj; // covariate effects on slope vary among species
matrix[nrts_cov,nspecies] cov_a_adj; // covariate effects on intercept vary among species

  vector[nspecies] b_cov; // species-level random effects on slopes
  vector[nspecies] a_cov; // species-level random effects on slopes


//by species     
real E1[ncounts1];           // log_scale additive likelihood
vector[nroutes1] alpha1;
vector[nroutes1] beta1;
vector[ncounts1] noise1;

real E2[ncounts2];           // log_scale additive likelihood
vector[nroutes2] alpha2;
vector[nroutes2] beta2;
vector[ncounts2] noise2;

real E3[ncounts3];           // log_scale additive likelihood
vector[nroutes3] alpha3;
vector[nroutes3] beta3;
vector[ncounts3] noise3;

real E4[ncounts4];           // log_scale additive likelihood
vector[nroutes4] alpha4;
vector[nroutes4] beta4;
vector[ncounts4] noise4;

real E5[ncounts5];           // log_scale additive likelihood
vector[nroutes5] alpha5;
vector[nroutes5] beta5;
vector[ncounts5] noise5;

real E6[ncounts6];           // log_scale additive likelihood
vector[nroutes6] alpha6;
vector[nroutes6] beta6;
vector[ncounts6] noise6;

real E7[ncounts7];           // log_scale additive likelihood
vector[nroutes7] alpha7;
vector[nroutes7] beta7;
vector[ncounts7] noise7;

  
  B_cov ~ normal(0,0.1);//covariate effect on slope parameters
  A_cov ~ normal(0,0.1);//covariate effect on intercept parameters
  
    sdb_cov ~ normal(0,0.1); //prior on sd of slope covariate parameters
  b_cov_raw ~ std_normal();//species level slope covariate effects
  sum(b_cov_raw) ~ normal(0,0.001*nspecies); //sum to zero constraint
  b_cov = (sdb_cov * b_cov_raw) + B_cov;
  
    sda_cov ~ std_normal(); //prior on sd of intercept covariate parameters
  a_cov_raw ~ std_normal();//species level intercept covariate effects
  sum(a_cov_raw) ~ normal(0,0.001*nspecies); //sum to zero constraint
  a_cov = (sda_cov * a_cov_raw) + A_cov;


for(s in 1:nspecies){
  cov_b_adj[,s] = b_cov[s]*rt_cov_dat;
  cov_a_adj[,s] = a_cov[s]*rt_cov_dat;

}


  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  
  sdobs ~ std_normal(); //prior on sd of gam hyperparameters
 
  obs_raw ~ std_normal();//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers); //sum to zero constraint

  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ normal(0,1);// prior on fixed effect mean intercept
  sdalpha ~ normal(0,1); //prior on sd of intercept variation
  sdbeta ~ normal(0,0.025); //prior on sd of slope variation - mild regularization

   obs = sdobs*obs_raw;

  //spatial iCAR intercepts and slopes by route and species


//by species     

  beta_raw1 ~ icar_normal_lpdf(nroutes1, node11, node21);
  alpha_raw1 ~ icar_normal_lpdf(nroutes1, node11, node21);
  noise_raw1 ~ std_normal();//~ student_t(4,0,1); //alternative heavy tailed extra Poisson log-normal variance
  alpha1 =  (sdalpha[1]*alpha_raw1) + ALPHA[1];
  beta1 =  (sdbeta[1]*beta_raw1) + BETA[1];
  noise1 =  (sdnoise[1]*noise_raw1);
   for(i in 1:ncounts1){
     E1[i] = (beta1[route1[i]] + cov_b_adj[route_cov1[i],1]) * (year1[i]-fixedyear) + alpha1[route1[i]] + obs[observer1[i]] + cov_a_adj[route_cov1[i],1] + noise1[i]; // 
  }
    count1 ~ poisson_log(E1); // count likelihood with log-transformation



  beta_raw2 ~ icar_normal_lpdf(nroutes2, node12, node22);
  alpha_raw2 ~ icar_normal_lpdf(nroutes2, node12, node22);
  noise_raw2 ~ std_normal();//~ student_t(4,0,1); //heavy tailed extra Poisson log-normal variance
  alpha2 =  (sdalpha[2]*alpha_raw2) + ALPHA[2];
  beta2 =  (sdbeta[2]*beta_raw2) + BETA[2];
  noise2 =  (sdnoise[2]*noise_raw2);
   for(i in 1:ncounts2){
     E2[i] = (beta2[route2[i]] + cov_b_adj[route_cov2[i],2])  * (year2[i]-fixedyear) + alpha2[route2[i]] + obs[observer2[i]] + cov_a_adj[route_cov2[i],2] + noise2[i]; // 
  }
    count2 ~ poisson_log(E2); // count likelihood with log-transformation



  beta_raw3 ~ icar_normal_lpdf(nroutes3, node13, node23);
  alpha_raw3 ~ icar_normal_lpdf(nroutes3, node13, node23);
  noise_raw3 ~ std_normal();//~ student_t(4,0,1); //heavy tailed extra Poisson log-normal variance
  alpha3 =  (sdalpha[3]*alpha_raw3) + ALPHA[3];
  beta3 =  (sdbeta[3]*beta_raw3) + BETA[3];
  noise3 =  (sdnoise[3]*noise_raw3);
   for(i in 1:ncounts3){
     E3[i] = (beta3[route3[i]] + cov_b_adj[route_cov3[i],3])  * (year3[i]-fixedyear) + alpha3[route3[i]] + obs[observer3[i]] + cov_a_adj[route_cov3[i],3] + noise3[i]; // 
  }
    count3 ~ poisson_log(E3); // count likelihood with log-transformation


  beta_raw4 ~ icar_normal_lpdf(nroutes4, node14, node24);
  alpha_raw4 ~ icar_normal_lpdf(nroutes4, node14, node24);
  noise_raw4 ~ std_normal();//~ student_t(4,0,1); //heavy tailed extra Poisson log-normal variance
  alpha4 =  (sdalpha[4]*alpha_raw4) + ALPHA[4];
  beta4 =  (sdbeta[4]*beta_raw4) + BETA[4];
  noise4 =  (sdnoise[4]*noise_raw4);
   for(i in 1:ncounts4){
     E4[i] = (beta4[route4[i]] + cov_b_adj[route_cov4[i],4])  * (year4[i]-fixedyear) + alpha4[route4[i]] + obs[observer4[i]] + cov_a_adj[route_cov4[i],4] + noise4[i]; // 
  }
    count4 ~ poisson_log(E4); // count likelihood with log-transformation


  beta_raw5 ~ icar_normal_lpdf(nroutes5, node15, node25);
  alpha_raw5 ~ icar_normal_lpdf(nroutes5, node15, node25);
  noise_raw5 ~ std_normal();//~ student_t(4,0,1); //heavy tailed extra Poisson log-normal variance
  alpha5 =  (sdalpha[5]*alpha_raw5) + ALPHA[5];
  beta5 =  (sdbeta[5]*beta_raw5) + BETA[5];
  noise5 =  (sdnoise[5]*noise_raw5);
   for(i in 1:ncounts5){
     E5[i] = (beta5[route5[i]] + cov_b_adj[route_cov5[i],5])  * (year5[i]-fixedyear) + alpha5[route5[i]] + obs[observer5[i]] + cov_a_adj[route_cov5[i],5] + noise5[i]; // 
  }
    count5 ~ poisson_log(E5); // count likelihood with log-transformation



  beta_raw6 ~ icar_normal_lpdf(nroutes6, node16, node26);
  alpha_raw6 ~ icar_normal_lpdf(nroutes6, node16, node26);
  noise_raw6 ~ std_normal();//~ student_t(4,0,1); //heavy tailed extra Poisson log-normal variance
  alpha6 =  (sdalpha[6]*alpha_raw6) + ALPHA[6];
  beta6 =  (sdbeta[6]*beta_raw6) + BETA[6];
  noise6 =  (sdnoise[6]*noise_raw6);
   for(i in 1:ncounts6){
     E6[i] = (beta6[route6[i]] + cov_b_adj[route_cov6[i],6])  * (year6[i]-fixedyear) + alpha6[route6[i]] + obs[observer6[i]] + cov_a_adj[route_cov6[i],6] + noise6[i]; // 
  }
    count6 ~ poisson_log(E6); // count likelihood with log-transformation


 
  beta_raw7 ~ icar_normal_lpdf(nroutes7, node17, node27);
  alpha_raw7 ~ icar_normal_lpdf(nroutes7, node17, node27);
  noise_raw7 ~ std_normal();//~ student_t(4,0,1); //heavy tailed extra Poisson log-normal variance
  alpha7 =  (sdalpha[7]*alpha_raw7) + ALPHA[7];
  beta7 =  (sdbeta[7]*beta_raw7) + BETA[7];
  noise7 =  (sdnoise[7]*noise_raw7);
   for(i in 1:ncounts7){
     E7[i] = (beta7[route7[i]] + cov_b_adj[route_cov7[i],7])  * (year7[i]-fixedyear) + alpha7[route7[i]] + obs[observer7[i]] + cov_a_adj[route_cov7[i],7] + noise7[i]; // 
  }
    count7 ~ poisson_log(E7); // count likelihood with log-transformation


}

  generated quantities {

 real log_lik[ncounts];

//by species     
real E1[ncounts1];           // log_scale additive likelihood
vector[nroutes1] alpha1;
vector[nroutes1] alpha1_adj;
vector[nroutes1] beta1;
vector[nroutes1] beta1_adj;
vector[ncounts1] noise1;

real E2[ncounts2];           // log_scale additive likelihood
vector[nroutes2] alpha2;
vector[nroutes2] alpha2_adj;
vector[nroutes2] beta2;
vector[nroutes2] beta2_adj;
vector[ncounts2] noise2;

real E3[ncounts3];           // log_scale additive likelihood
vector[nroutes3] alpha3;
vector[nroutes3] alpha3_adj;
vector[nroutes3] beta3;
vector[nroutes3] beta3_adj;
vector[ncounts3] noise3;

real E4[ncounts4];           // log_scale additive likelihood
vector[nroutes4] alpha4;
vector[nroutes4] alpha4_adj;
vector[nroutes4] beta4;
vector[nroutes4] beta4_adj;
vector[ncounts4] noise4;

real E5[ncounts5];           // log_scale additive likelihood
vector[nroutes5] alpha5;
vector[nroutes5] alpha5_adj;
vector[nroutes5] beta5;
vector[nroutes5] beta5_adj;
vector[ncounts5] noise5;

real E6[ncounts6];           // log_scale additive likelihood
vector[nroutes6] alpha6;
vector[nroutes6] alpha6_adj;
vector[nroutes6] beta6;
vector[nroutes6] beta6_adj;
vector[ncounts6] noise6;

real E7[ncounts7];           // log_scale additive likelihood
vector[nroutes7] alpha7;
vector[nroutes7] alpha7_adj;
vector[nroutes7] beta7;
vector[nroutes7] beta7_adj;
vector[ncounts7] noise7;

   vector[nobservers] obs ;             // observer effects by species - consider pooling across all species


matrix[nrts_cov,nspecies] cov_b_adj; // covariate effects on slope vary among species
matrix[nrts_cov,nspecies] cov_a_adj; // covariate effects on intercept vary among species

  vector[nspecies] b_cov; // species-level random effects on slopes
  vector[nspecies] a_cov; // species-level random effects on slopes


 
  obs = sdobs*obs_raw;


  b_cov = (sdb_cov * b_cov_raw) + B_cov;
  a_cov = (sda_cov * a_cov_raw) + A_cov;


for(s in 1:nspecies){
  cov_b_adj[,s] = b_cov[s]*rt_cov_dat;
  cov_a_adj[,s] = a_cov[s]*rt_cov_dat;

}



//by species     
  alpha1 =  (sdalpha[1]*alpha_raw1) + ALPHA[1];
  beta1 =  (sdbeta[1]*beta_raw1) + BETA[1];
  noise1 =  (sdnoise[1]*noise_raw1);
  for(i in 1:nroutes1){
    beta1_adj[i] = beta1[i]+cov_b_adj[cov_routes1[i],1];
    alpha1_adj[i] = alpha1[i]+cov_a_adj[cov_routes1[i],1];
  }
   for(i in 1:ncounts1){
     E1[i] =   (beta1[route1[i]] + cov_b_adj[route_cov1[i],1]) * (year1[i]-fixedyear) + alpha1[route1[i]] + obs[observer1[i]] + cov_a_adj[route_cov1[i],1] + noise1[i]; // 
     log_lik[i] = poisson_log_lpmf(count1[i] | E1[i]);
 }


  alpha2 =  (sdalpha[2]*alpha_raw2) + ALPHA[2];
  beta2 =  (sdbeta[2]*beta_raw2) + BETA[2];
  noise2 =  (sdnoise[2]*noise_raw2);
  for(i in 1:nroutes2){
    beta2_adj[i] = beta2[i]+cov_b_adj[cov_routes2[i],2];
    alpha2_adj[i] = alpha2[i]+cov_a_adj[cov_routes2[i],2];
  }
   for(i in 1:ncounts2){
     E2[i] =  (beta2[route2[i]] + cov_b_adj[route_cov2[i],2]) * (year2[i]-fixedyear) + alpha2[route2[i]] + obs[observer2[i]] + cov_a_adj[route_cov2[i],2] + noise2[i]; // 
     log_lik[i+ncounts1] = poisson_log_lpmf(count2[i] | E2[i]);
 }



  alpha3 =  (sdalpha[3]*alpha_raw3) + ALPHA[3];
  beta3 =  (sdbeta[3]*beta_raw3) + BETA[3];
  noise3 =  (sdnoise[3]*noise_raw3);
  for(i in 1:nroutes3){
    beta3_adj[i] = beta3[i]+cov_b_adj[cov_routes3[i],3];
    alpha3_adj[i] = alpha3[i]+cov_a_adj[cov_routes3[i],3];
  }
   for(i in 1:ncounts3){
     E3[i] = (beta3[route3[i]] + cov_b_adj[route_cov3[i],3]) * (year3[i]-fixedyear) + alpha3[route3[i]] + obs[observer3[i]] + cov_a_adj[route_cov3[i],3] + noise3[i]; // 
    log_lik[i+ncounts1+ncounts2] = poisson_log_lpmf(count3[i] | E3[i]);
}

  alpha4 =  (sdalpha[4]*alpha_raw4) + ALPHA[4];
  beta4 =  (sdbeta[4]*beta_raw4) + BETA[4];
  noise4 =  (sdnoise[4]*noise_raw4);
  for(i in 1:nroutes4){
    beta4_adj[i] = beta4[i]+cov_b_adj[cov_routes4[i],4];
    alpha4_adj[i] = alpha4[i]+cov_a_adj[cov_routes4[i],4];
  }
   for(i in 1:ncounts4){
     E4[i] = (beta4[route4[i]] + cov_b_adj[route_cov4[i],4]) * (year4[i]-fixedyear) + alpha4[route4[i]] + obs[observer4[i]] + cov_a_adj[route_cov4[i],4] + noise4[i]; // 
    log_lik[i+ncounts1+ncounts2+ncounts3] = poisson_log_lpmf(count4[i] | E4[i]);
}
 

  alpha5 =  (sdalpha[5]*alpha_raw5) + ALPHA[5];
  beta5 =  (sdbeta[5]*beta_raw5) + BETA[5];
  noise5 =  (sdnoise[5]*noise_raw5);
  for(i in 1:nroutes5){
    beta5_adj[i] = beta5[i]+cov_b_adj[cov_routes5[i],5];
    alpha5_adj[i] = alpha5[i]+cov_a_adj[cov_routes5[i],5];
  }
   for(i in 1:ncounts5){
     E5[i] = (beta5[route5[i]] + cov_b_adj[route_cov5[i],5]) * (year5[i]-fixedyear) + alpha5[route5[i]] + obs[observer5[i]] + cov_a_adj[route_cov5[i],5] + noise5[i]; // 
    log_lik[i+ncounts1+ncounts2+ncounts3+ncounts4] = poisson_log_lpmf(count5[i] | E5[i]);
 }

  alpha6 =  (sdalpha[6]*alpha_raw6) + ALPHA[6];
  beta6 =  (sdbeta[6]*beta_raw6) + BETA[6];
  noise6 =  (sdnoise[6]*noise_raw6);
  for(i in 1:nroutes6){
    beta6_adj[i] = beta6[i]+cov_b_adj[cov_routes6[i],6];
    alpha6_adj[i] = alpha6[i]+cov_a_adj[cov_routes6[i],6];
  }
   for(i in 1:ncounts6){
     E6[i] = (beta6[route6[i]] + cov_b_adj[route_cov6[i],6]) * (year6[i]-fixedyear) + alpha6[route6[i]] + obs[observer6[i]] + cov_a_adj[route_cov6[i],6] + noise6[i]; // 
    log_lik[i+ncounts1+ncounts2+ncounts3+ncounts4+ncounts5] = poisson_log_lpmf(count6[i] | E6[i]);
 }

  alpha7 =  (sdalpha[7]*alpha_raw7) + ALPHA[7];
  beta7 =  (sdbeta[7]*beta_raw7) + BETA[7];
  noise7 =  (sdnoise[7]*noise_raw7);
  for(i in 1:nroutes7){
    beta7_adj[i] = beta7[i]+cov_b_adj[cov_routes7[i],7];
    alpha7_adj[i] = alpha7[i]+cov_a_adj[cov_routes7[i],7];
  }
   for(i in 1:ncounts7){
     E7[i] = (beta7[route7[i]] + cov_b_adj[route_cov7[i],7]) * (year7[i]-fixedyear) + alpha7[route7[i]] + obs[observer7[i]] + cov_a_adj[route_cov7[i],7] + noise7[i]; // 
    log_lik[i+ncounts1+ncounts2+ncounts3+ncounts4+ncounts5+ncounts6] = poisson_log_lpmf(count7[i] | E7[i]);
 }


    }
    
    
  
  



