functions{
  real[] SIR_BD_period_ODE_onestep(real t,
    real [] X,
    real [] theta,
    real[] x_r, int [] x_i){
	/**
	one step differentiation of ODE
  t: time
  X: real {S,I}

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
    real dX[2];

    N = theta[6];
    beta = theta[1] * theta[2] / N;
    betat = beta * (1 + theta[4] * sin(2 * pi()* t / theta[5]));
    gamma = theta[2];
    mu = theta[3];

    dX[1] =  - betat * X[1] * X[2] + mu * N - mu  * X[1];
    dX[2] = betat * X[1] * X[2] - mu * X[2] - gamma * X[2];
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

}


data {
  /**
  grid and time configuration

  */
 // int<lower = 1> n;
  int<lower = 1> ngrid;
//  int<lower = 1> gridsize;
  real period;
  real N;
  //real ts[ngrid + 1];
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
  //real<lower=1> SI[ngrid,2];
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
 matrix [ngrid ,2] SI;
 matrix [ngrid +1 ,2] SIs;
 real theta[6];

 theta = {R0,gamma,mu,A,period,N};
 Init_State[1] = N * Alpha / (1 + Alpha);
 Init_State[2] = N - Init_State[1];

 SI = to_matrix(integrate_ode_rk45(SIR_BD_period_ODE_onestep,Init_State,
  tt[1],tt[2:(ngrid+1)],theta,x_r,x_i));
SIs = append_row(to_row_vector(Init_State),SI);
}

model{
  R0 ~ uniform(3,8);
 // gamma ~ lognormal(-3,0.4);
  //Alpha ~ lognormal(1,10);
  gamma ~ lognormal(-2.5,0.2);
  Alpha ~ lognormal(2,2);
  A ~ uniform(0,1);
  increment_log_prob(coal_log_stan(w, C, y, to_vector(tt), tids, SIs, theta, nc));
}


