data {
  int S; // 41 (0-40dB) stimulus dynamic range
  int X; // 36 (0-35dB) true sensitivity range
  int K; // 42 (-1-40dB) measured sensitivity range
  int Y[X, K]; // counts of measued sensitivity K at locations X dB
  int N_response;
  int<lower=0, upper=1> seen[N_response];
  int<lower=0, upper=40> db[N_response];
  int N_sequence;
  int<lower=-1, upper=40> final[N_sequence];
  int start_row[N_sequence];
  int stop_row[N_sequence];
  real gamma;
}

parameters {
  vector<lower=-5, upper=40>[X] mu;
  vector<lower=-2, upper=3>[X] sigma_log;
  vector<lower=-5, upper=5>[X] lamda_logit;
}

transformed parameters {
  vector[X] sigma = exp(sigma_log);
  vector[X] lamda = inv_logit(lamda_logit);
  vector[K] theta[X];

  for (x in 1:X){
    vector[S] psi;
    theta[x] = rep_vector(0, K);

    for (s in 1:S)
      psi[s]=gamma+(1-lamda[x]-gamma)*(1-normal_cdf(s-1, mu[x], sigma[x]));
  
    for (n in 1:N_sequence){
      real P=1;
      
      for (i in start_row[n]:stop_row[n])
        P *= seen[i] ? psi[db[i]+1] : 1-psi[db[i]+1];
        
      theta[x][final[n]+2] += P;
    }
    
    theta[x] /= sum(theta[x]);
  }
}

model {
  for (x in 1:X)
    target += multinomial_lpmf(Y[x] | theta[x]);

  target += normal_lpdf(mu[2:X] | mu[1:(X-1)], 1);
  target += normal_lpdf(sigma_log[2:X] | sigma_log[1:(X-1)], 0.2);
  target += normal_lpdf(lamda_logit[2:X] | lamda_logit[1:(X-1)], 0.3);
}

generated quantities {
  vector[X] asym = 1-lamda;
}
