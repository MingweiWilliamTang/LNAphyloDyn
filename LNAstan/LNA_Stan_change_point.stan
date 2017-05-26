# phylodynamic with changing point
functions{

  real[] SIR_changepoint_onestep(real t,
    real [] X,
    real [] theta,
    real[] x_r, int [] x_i){
	/**
	one step differentiation of ODE
  t: time
  X: real {S,I,m[1],m[2], S[1,1], S[1,2], S[2,2]}

  theta = c(R0,gamma,changscale,lambda)
  R0: R0 for infectious disease
  beta: infection rate
  gamma: recover rate
  mu: birth and death rate
  A: scale in the periodicity, i.e betat = beta * (1 + A * sin(2*pi*t/ period))
  period: period of infectious disease, one period denote one year
  N: the population size

*/

    real betat;
    real gamma;
    real N;
    real Rt;
    real dX[5];
    matrix[2,2] A;
    matrix[2,2] S;
    matrix[2,2] F;
    vector [2] m;
    vector[2] h;
    int i;

    A = to_matrix({{-1,1},{0,-1}});
    S = to_matrix({{X[3],X[4]},{X[4],X[5]}});
    m = to_vector({X[6],X[7]});
    i = 1;
    N = x_r[1];
    Rt = theta[1];
    while(x_r[i+1]<= t){
      Rt = Rt * theta[i+3];
      if(i == x_i[1]) break;
      i = i + 1;
    }

    betat = Rt * theta[2] / N;
    gamma = theta[2];

    h[1] = betat * X[1] * X[2];
	  h[2] = gamma * X[2];

    F[1,1] = - betat * X[2];
	  F[1,2] = - betat * X[1];
	  F[2,1] = betat * X[2];
	  F[2,2] = betat * X[1] - gamma;


    S = F * S + S * F' + A' * diag_matrix(h) * A;
    dX[1] =  - betat * X[1] * X[2] ;
    dX[2] = betat * X[1] * X[2] - gamma * X[2];

    dX[3] = S[1,1];
    dX[4] = S[1,2];
    dX[5] = S[2,2];

    m = F * m;
    dX[6] = m[1];
    dX[7] = m[2];

    return dX;
	  }
  matrix LNA_rs_onestep(row_vector SI, real t0, real t1, real [] theta, real[] x_r, int [] x_i){
    real X1[1,5];
    real Init[5];
    matrix [3,2] res;
    Init = {SI[1], SI[2], 0.0, 0.0, 0.0};

    X1 = integrate_ode_rk45(SIR_changepoint_onestep,Init,t0,{t1},theta,x_r,x_i);

    // m + eta

    // eta
    res[1,1] = X1[1,1];
    res[1,2] = X1[1,2];
    //phi
    res[2,1] = X1[1,3];
    res[2,2] = X1[1,4];
    res[3,1] = X1[1,4];
    res[3,2] = X1[1,5];
    return res;
  }


 real coal_log_stan(vector w, vector C, vector y, vector tt, int[] tids,
      matrix SI,real [] theta,int nc,real[] x_r){
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
  log_like = 0;

  for(i in 1:nc){
    t = tt[tids[i]];
    Ne_inv = theta[4] / SI[tids[i]][2];
    log_like = log_like + log(Ne_inv) * y[i] - 2 * w[i] * C[i] * Ne_inv;
  }

  return log_like;
  }

}


data {
  /**
  grid and time configuration

  */
  int<lower = 1> ngrid;
  real x_r[2];   //real N; real ct;  time of changepoint
  real tt[ngrid + 1];
  int<lower = 0> nch;
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

transformed data {                  // no real data for ODE system
  int x_i[1];
  x_i[1] = nch; // no integer data for ODE system
}


parameters{
  /**
  R0: R0 in infectious disease
  gamma: recover rate
  ch: scale of changepoint
  ALpha: ratio of S0/I0
  SI: vector array of S and I
  */

  real<lower=0,upper=10> R0;
  real<lower=0> gamma;
  real<lower=0> ch;
   real<lower=0> Alpha;
      real<lower=1> lambda;
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
  real theta[4];
 theta = {R0,gamma,lambda,ch};
 Init_State[1] = x_r[1] * Alpha / (1 + Alpha);
 Init_State[2] = x_r[1] / (1 + Alpha) ;

}

model{
   matrix [3,2] FT[ngrid];
   R0 ~ uniform(1,10);
 // gamma ~ lognormal(-3,0.4);
  //Alpha ~ lognormal(1,10);
  gamma ~ lognormal(3,0.5);
  Alpha ~ lognormal(8,2);
  lambda ~ lognormal(10,3);
  ch ~ double_exponential(0,1);

   for(i in 1:ngrid){
    if(i == 1){
    FT[i] = LNA_rs_onestep(to_row_vector(Init_State), tt[i], tt[i+1], theta, x_r, x_i);
       //FT = integrate_ode_rk45(SIR_changepoint_onestep,{Init_State[1],Init_State[2],0.0,0.0,0.0},0,{0.3},theta,x_r,x_i);
       SI[i] ~ multi_normal(to_vector({1,1}),to_matrix({{1,0},{0,1}}));
     #print(FT);
      #SI[i] ~ multi_normal(to_vector(FT[i][1]),block(FT[i],2,1,2,2));
    }
    else{
      //print(to_vector(eta[1 + (i-1) * gridsize]) + FT[i-1] * (SI[i-1] - eta[1 + (i-2) * gridsize])');
      #FT[i] = LNA_rs_onestep(SI[i-1], tt[i], tt[i+1], theta, x_r, x_i);
      #SI[i] ~ multi_normal(to_vector(FT[i][1]),block(FT[i],2,1,3,2));
      SI[i] ~ multi_normal(to_vector({1,1}),to_matrix({{1,0},{0,1}}));
    }
  }
  increment_log_prob(coal_log_stan(w, C, y, to_vector(tt), tids, append_row(to_row_vector(Init_State),SI), theta, nc, x_r));
}
