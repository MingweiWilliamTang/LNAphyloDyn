functions{
  real[] SIR_BD_period_ODE_onestep(real t,
    real [] X,
    real [] theta,
    real[] x_r, int [] x_i){
	/**
	one step differentiation of ODE
  t: time
  X: real {S,I,m[1],m[2], S[1,1], S[1,2], S[2,2]}

  theta = c(R0,gamma, mu, A,period,N)
  R0: R0 for infectious disease
  beta: infection rate
  gamma: recover rate
  mu: birth and death rate
  A: scale in the periodicity, i.e betat = beta * (1 + A * sin(2*pi*t/ period))
  period: period of infectious disease, one period denote one year
  N: the population size

*/

    real beta;
    real betat;
    real gamma;
    real mu;
    real N;
    real dX[7];
    matrix[4,2] A;
    matrix[2,2] S;
    matrix[2,2] F;
    vector[4] h;
    vector[2] l;

    A = to_matrix({{-1,1},{0,-1},{1,0},{0,-1}});
    S = to_matrix({{X[5],X[6]},{X[6],X[7]}});

    N = theta[6];
    beta = theta[1] * theta[2] / N;
    betat = beta * (1 + theta[4] * sin(2 * pi()* t / theta[5]));
    gamma = theta[2];
    mu = theta[3];

    h[1] = betat * X[1] * X[2];
	  h[2] = gamma * X[2];
	  h[3] = mu * (N-X[1]);
	  h[4] = mu * X[2];

    F[1,1] = - betat * X[2] - mu;
	  F[1,2] = - betat * X[1];
	  F[2,1] = betat * X[2];
	  F[2,2] = betat * X[1] - gamma - mu;

    l = F * to_vector({X[3],X[4]});

    S = F * S + S * F' + A' * diag_matrix(h) * A;
    dX[1] =  - betat * X[1] * X[2] + mu * N - mu  * X[1];
    dX[2] = betat * X[1] * X[2] - mu * X[2] - gamma * X[2];

    dX[3] = l[1];
    dX[4] = l[2];

    dX[5] = S[1,1];
    dX[6] = S[1,2];
    dX[7] = S[2,2];

    return dX;
	  }


 real coal_log_stan(vector w, vector C, vector y, vector tt, int[] tids,
      matrix SI,real [] theta,int nc){
  /*
  w: vector of mixed observed time interval of coalescent, sampling, grid
  C:  number of active lineages choose 2
  y: indicator for coalescent time: 1 coalesecent, 0 not a coalescent time
  tt: the time grid for each interval
  SI: counts of infected and recovered
  nc: number of terms for sum
  */
    real log_like;
    real Ne_inv;
    real t;
    real betat;
    real beta;
    beta = theta[1] * theta[2] / theta[6];
  log_like = 0;

  for(i in 1:nc){
    t = tt[tids[i]];
    betat = beta * (1 + theta[4] * sin(2 * pi()* t / theta[5]));
    Ne_inv = SI[tids[i]][1] / SI[tids[i]][2] * betat * theta[6];
    log_like = log_like + log(Ne_inv) * y[i] - 2 * w[i] * C[i] * Ne_inv;
  }

  return log_like;
  }

  matrix LNA_ns_onestep(row_vector SI, row_vector X0, real t0, real t1, real [] theta,real[] x_r, int [] x_i){
    real Init[7];
    vector [7] X1;
    matrix [4,2] res;
    Init = {X0[1], X0[2], SI[1] - X0[1], SI[2] - X0[2], 0, 0, 0};

    X1 = to_vector(integrate_ode_bdf(SIR_BD_period_ODE_onestep,
     Init,t0,{t1},theta,x_r,x_i)[1]);
/*
    // m + eta
    res[1,1] = X1[1] + X1[3];
    res[1,2] = X1[2] + X1[4];
    // eta
    res[2,1] = X1[1];
    res[2,2] = X1[2];
    //phi
    res[3,1] = X1[5];
    res[3,2] = X1[6];
    res[4,1] = X1[6];
    res[4,2] = X1[7];
*/
    return res;
  }


}


data {
  /**
  grid and time configuration

  */
  int<lower = 1> n;
  int<lower = 1> ngrid;
  int<lower = 1> gridsize;
  real period;
  real N;
  real ts[n];
  real tt[ngrid + 1];
  real mu;
  /**
  coalescent related data
  */

  int nc;
  int tids[nc];
  vector[nc] w;
  vector[nc] C;
  vector[nc] y;

  /**
  parmaters for prior distribution
  */

}

transformed data {
  real x_r[0];                 // no real data for ODE system
  int x_i[0];                  // no integer data for ODE system
}


parameters{
  /**
  R0: R0 in infectious disease
  gamma: recover rate
  A: scale of periodicity
  Alpha: the ratio of initial states, i.e S0/I0
  SI: vector array of S and I
  */

  real<lower=0,upper=10> R0;
  real<lower=0> gamma;
  real<lower=0,upper=1> A;
  real<lower=0> Alpha;
  //vector[2] SI[ngrid+1];
  matrix<lower=1> [ngrid,2] SI;
}



transformed parameters{
/**
  x_r, x_i unused varible for ode
  transform the parameters into
  eta: Ode trajectory for S and I
  B: linear transform matrices in LNA
  S: Covariance matrices in LNA
*/
 real Init_State[2];
  real theta[6];
 theta = {R0,gamma,mu,A,period,N};
 Init_State[1] = N * Alpha / (1 + Alpha);
 Init_State[2] = N - Init_State[1];

}


model{
  matrix [4,2] FT[ngrid];
  R0 ~ uniform(1,10);
 // gamma ~ lognormal(-3,0.4);
  //Alpha ~ lognormal(1,10);
  gamma ~ lognormal(-2.5,0.1);
  Alpha ~ lognormal(2,2);
  A ~ uniform(0,1);
  for(i in 1:1){
    if(i == 1){
      FT[i] = LNA_ns_onestep(to_row_vector(Init_State), to_row_vector(Init_State), tt[i], tt[i+1], theta, x_r, x_i);
       SI[i] ~ multi_normal(to_vector({1,1}),to_matrix({{1,0},{0,1}}));
     // SI[i] ~ multi_normal(to_vector(FT[i][1]),block(FT[i],3,3,4,4));
    }
    /*
    else{
      //print(to_vector(eta[1 + (i-1) * gridsize]) + FT[i-1] * (SI[i-1] - eta[1 + (i-2) * gridsize])');
      FT[i] = LNA_ns_onestep(SI[i-1], FT[i][2], tt[i], tt[i+1], theta, x_r, x_i);
      SI[i] ~ multi_normal(to_vector(FT[i][1]),block(FT[i],3,3,4,4));
    }
    */
  }
  increment_log_prob(coal_log_stan(w, C, y, to_vector(tt), tids, append_row(to_row_vector(Init_State),SI), theta, nc));
}


