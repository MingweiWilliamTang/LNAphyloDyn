# restarting LNA
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



  matrix IntSigma(real [] t_k,matrix eta,int gridsize, real [] theta){
    matrix [2,2] FF;
    matrix [2,2] S;
    matrix[4,2] A;    // reaction matrix A
    real dt;
    dt = t_k[2] - t_k[1];

    A[1,1] = -1;
    A[1,2] = 1;
    A[2,1] = 0;
    A[2,2] = -1;
    A[3,1] = 1;
    A[3,2] = 0;
    A[4,1] = 0;
    A[4,2] = -1;

    S[1,1] = 0;
    S[2,1] = 0;
    S[1,2] = 0;
    S[2,2] = 0;


    for(i in 1:gridsize){

       FF = SIR_BD_F_onestep(t_k[i], eta[i], theta);
       S = S + dt * (FF * S + S * FF' + A' * diff_Var(t_k[i],eta[i],theta) * A);
    }
    return S;
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
  matrix [gridsize,2] etas[ngrid];
  matrix[2,2] S[ngrid];
  R0 ~ uniform(1,10);
 // gamma ~ lognormal(-3,0.4);
  //Alpha ~ lognormal(1,10);
  gamma ~ lognormal(-2.5,0.1);
  Alpha ~ lognormal(2,2);
  A ~ uniform(0,1);

  for(i in 1:ngrid){
    if(i == 1){
      etas[i] = to_matrix(integrate_ode_rk45(SIR_BD_period_ODE_onestep,Init_State,ts[1],ts[2:(gridsize+1)],theta,x_r,x_i));
      S[i] = IntSigma(ts[(1+(i-1) * gridsize) : i * gridsize], etas[i], gridsize, theta);
      SI[i] ~ multi_normal(to_vector(etas[i][gridsize]),S[i]);
    }else{
      //print(to_vector(eta[1 + (i-1) * gridsize]) + FT[i-1] * (SI[i-1] - eta[1 + (i-2) * gridsize])');

      etas[i] = to_matrix(integrate_ode_rk45(SIR_BD_period_ODE_onestep,to_array_1d(SI[i-1]),
      ts[(i-1)*gridsize+ 1],ts[((i-1)*gridsize + 2):(i * gridsize+1)],theta,x_r,x_i));
      S[i] = IntSigma(ts[(1+(i-1) * gridsize) : i * gridsize], etas[i], gridsize, theta);
      //SI[i] ~ multi_normal(to_vector(etas[i][gridsize]),S[i]);

       SI[i] ~ multi_normal(to_vector(etas[i]),Si[i]);
    }
  }
  increment_log_prob(coal_log_stan(w, C, y, to_vector(tt), tids, append_row(to_row_vector(Init_State),SI), theta, nc));
}


