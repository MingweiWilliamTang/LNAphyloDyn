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


  matrix submatrix(matrix X, int [] index ,int col_indexing){
    matrix[size(index),cols(X)] M1;
    matrix[rows(X),size(index)] M2;
    if(col_indexing == 1){
      for(i in 1:size(index)){
        M2[,i] = X[,index[i]];
      }
      return M2;
    }else{
      for(i in 1:size(index)){
        M1[i] = X[index[i]];
      }
      return M1;
    }
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


	matrix SIR_BD_F_onestep(real t, row_vector X, real [] theta){
		  /**
		  One step differentiation of residual process

		  t: time t
		  X: {S,I}
		  theta: parameter vector
		  */

		  real beta;
	     	real betat;
	     real gamma;
	     real mu;
        real N;
        matrix[2,2] F;
      N = theta[6];
	     beta = theta[1] * theta[2] / N;
	     betat = beta * (1 + theta[4] * sin(2 * pi()* t / theta[5]));
	     gamma = theta[2];
	     mu = theta[3];

	     F[1,1] = - betat * X[2] - mu;
	     F[1,2] = - betat * X[1];
	     F[2,1] = betat * X[2];
	     F[2,2] = betat * X[1] - gamma - mu;
	     return F;
	}

	matrix diff_Var(real t, row_vector X, real [] theta){
	  /**
	  calcuate the the variance matrix in SDE at each time t
	  Var = diag(h)

	  t: time t
		  X: {S,I}
		  theta: parameter vector
	  */
	  real beta;
	  real betat;
	  real gamma;
	  real mu;
    real N;
	  vector[4] h;

       N = theta[6];
	     beta = theta[1] * theta[2] / theta[6];
	     betat = beta * (1 + theta[4] * sin(2 * pi()* t / theta[5]));
	     gamma = theta[2];
	     mu = theta[3];
	     h[1] = betat * X[1] * X[2];
	     h[2] = gamma * X[2];
	     h[3] = mu * (N-X[1]);
	     h[4] = mu * X[2];
	     return diag_matrix(h);
	}

matrix [] SIR_FT(real [] ts, real [] Init_State,int gridsize, int ngrid,
      real R0, real gamma, real mu, real a, real period,real N,
      real [] x_r, int [] x_i){
  /**
  ts: the time grid for integration
  eta: Ode path
  gridsize: the gridsize for LNA, i.e LNA grid = dt * gridsize
  ngrid: how many LNA grids
  R0
  gamma: recover rate
  mu: birth and death rate
  a: scale in the periodicity, i.e betat = beta * (1 + A * sin(2*pi*t/ period))
  period: period for infectious disease
  N: population size
  */

  real dt;
  real theta[6];
  real t;
  int id;
  matrix[4,2] A;    // reaction matrix A
  matrix [2,2] FF;
  matrix[2,2] B[3 * ngrid]; // linear transform matrix in LNA   // variance matrix in LNA
  matrix[gridsize*ngrid + 1,2] eta;
  theta = {R0, gamma, mu, a, period, N};
  dt = ts[2] - ts[1];
/*
  A[1,1] = -1;
  A[1,2] = 1;
  A[2,1] = 0;
  A[2,2] = -1;
  A[3,1] = 1;
  A[3,2] = 0;
  A[4,1] = 0;
  A[4,2] = -1;
*/
  A = to_matrix({{-1,1},{0,-1},{1,0},{0,-1}});
  t = 0;
  eta = append_row(to_row_vector(Init_State),to_matrix(integrate_ode_rk45(SIR_BD_period_ODE_onestep,Init_State,
  ts[1],ts[2:(gridsize*ngrid+1)],theta,x_r,x_i)));
  for(i in 1:ngrid){
  /**
  B[1:ngrid] stores the 2*2 matrices for exp(integrate(F))
  B[(ngrid+1)] stores the 2*2 matrices for the gaussian noise

  Need firstly need to initialize all element to zero


    B[i][1,1] = 0;
    B[i][2,1] = 0;
    B[i][1,2] = 0;
    B[i][2,2] = 0;
    B[ngrid + i][1,1] = 0;
    B[ngrid + i][2,1] = 0;
    B[ngrid + i][1,2] = 0;
    B[ngrid + i][2,2] = 0;
*/
    B[i] = to_matrix({{0,0},{0,0}});
    B[ngrid + i] = to_matrix({{0,0}, {0,0}});

   for(j in 1:gridsize){
      id = (i-1) * gridsize + j;
      FF = (SIR_BD_F_onestep(ts[id], eta[id], theta));

      B[i] = B[i] + dt * FF; // numerical integration of F by hand

      // numerical integration of variance
      B[ngrid + i] = B[ngrid + i] + dt * (FF * B[ngrid + i] + B[ngrid + i] * FF' + A' * diff_Var(ts[id],eta[id],theta) * A);
    }
      B[i] = matrix_exp(B[i]);
      B[2 * ngrid + i][1] = to_row_vector(eta[id + 1]);
    }
    //print(B[ngrid + 1]);
    return B;
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
 //matrix[ngrid,2] eta;
 real betatN[n];
 matrix [2,2] FT[3*ngrid];
  real theta[6];

 theta = {R0,gamma,mu,A,period,N};
 Init_State[1] = N * Alpha / (1 + Alpha);
 Init_State[2] = N - Init_State[1];
 //eta = to_matrix(integrate_ode_rk45(SIR_BD_period_ODE_onestep,Init_State,ts[1],tt[2:(ngrid+1)],theta,x_r,x_i));
 FT = SIR_FT(ts,Init_State,gridsize,ngrid, R0, gamma, mu, A, period,N,x_r,x_i);
}


model{
  R0 ~ uniform(1,10);
 // gamma ~ lognormal(-3,0.4);
  //Alpha ~ lognormal(1,10);
  gamma ~ lognormal(-2.5,0.1);
  Alpha ~ lognormal(2,2);
  A ~ uniform(0,1);
  for(i in 1:ngrid){
    if(i == 1){
      SI[i] ~ multi_normal(to_vector(FT[i + 2*ngrid][1]), FT[i + ngrid]);
    }else{
      //print(to_vector(eta[1 + (i-1) * gridsize]) + FT[i-1] * (SI[i-1] - eta[1 + (i-2) * gridsize])');
      SI[i] ~ multi_normal(to_vector(FT[i + 2*ngrid][1]) + FT[i] * (SI[i] - FT[i + 2*ngrid][1])', FT[i + ngrid]);
    }
  }
  increment_log_prob(coal_log_stan(w, C, y, to_vector(tt), tids, append_row(to_row_vector(Init_State),SI), theta, nc));
}


/*
data{
  int<lower=1> T;
  real y0[2];
  real t0;
  real ts[T];
  real theta[6];
}
generated quantities {
  matrix [T+1,2] y_hat;
  y_hat = append_row(to_row_vector(y0),to_matrix(integrate_ode_rk45(SIR_BD_period_ODE_onestep, y0, t0, ts, theta, x_r, x_i)));
      // add measurement error
      for (t in 1:T) {
        y_hat[t, 1] = y_hat[t, 1] + normal_rng(0, 0.0001);
        y_hat[t, 2] = y_hat[t, 2] + normal_rng(0, 0.0001);
      }
*/


