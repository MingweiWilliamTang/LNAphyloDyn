functions{

row_vector SIR_BD_period_ODE_onestep(real t,
vector X,
real [] theta)
{
	/**

	one step differentiation of ODE
  t: time
  X: vector {S,I}

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
row_vector[2] dX;

N = theta[6];
beta = theta[1] * theta[2] / N;
betat = beta * (1 + theta[4] * sin(2 * pi()* t / theta[5]));
gamma = theta[2];
mu = theta[3];

dX[1] =  - betat * X[1] * X[2] + mu * N - mu  * X[1];
dX[2] = betat * X[1] * X[2] - mu * X[2] - gamma * X[2];
return dX;
	}



matrix SIR_BD_ODE_stan(real[] ts, row_vector Init_State, real[] theta, int n_ode){
/**
Integrate the ode path

ts: the time grid for integration
Init_state: {S0,I0}
theta: parameter vector
n_ode: number of grids used the ode integration, i.e n_ode = length(ts)

return a ts * 2 matrix of ODE path, each row vector is {S,I} at time t
*/

   real dt;
   matrix [n_ode,2] traj;
  row_vector[2] k1;
  row_vector[2] k2;
  row_vector[2] k3;
  row_vector[2] k4;
  dt = ts[2] - ts[1];
  traj[1] = Init_State;
  for(i in 2:n_ode){
     // using 4th order ruge-kutta method for ode integration

     k1 = SIR_BD_period_ODE_onestep(ts[i-1],to_vector(traj[i-1]), theta);
     k2 = SIR_BD_period_ODE_onestep(ts[i-1],to_vector(traj[i-1] + k1 * dt / 2), theta);
     k3 = SIR_BD_period_ODE_onestep(ts[i-1],to_vector(traj[i-1] + k2 * dt / 2), theta);
     k4 = SIR_BD_period_ODE_onestep(ts[i-1],to_vector(traj[i-1] + k3 * dt / 2), theta);
     traj[i] = traj[i-1] + (k1/6 + k2/3 + k3/3 + k4/6) * dt;

    }
  return traj;
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




matrix [] SIR_FT(real [] ts, matrix eta, int gridsize, int ngrid, real R0, real gamma, real mu, real a, real period,real N){
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
  matrix[2,2] B[2 * ngrid]; // linear transform matrix in LNA   // variance matrix in LNA
  theta = {R0, gamma, mu, a, period, N};
  dt = ts[2] - ts[1];

  A[1,1] = -1;
  A[1,2] = 1;
  A[2,1] = 0;
  A[2,2] = -1;
  A[3,1] = 1;
  A[3,2] = 0;
  A[4,1] = 0;
  A[4,2] = -1;

t = 0;
for(i in 1:ngrid){
/**
B[1:ngrid] stores the 2*2 matrices for exp(integrate(F))
B[(ngrid+1)] stores the 2*2 matrices for the gaussian noise

Need firstly need to initialize all element to zero
*/

  B[i][1,1] = 0;
  B[i][2,1] = 0;
  B[i][1,2] = 0;
  B[i][2,2] = 0;
  B[ngrid + i][1,1] = 0;
  B[ngrid + i][2,1] = 0;
  B[ngrid + i][1,2] = 0;
  B[ngrid + i][2,2] = 0;

for(j in 1:gridsize){
  id = (i-1) * gridsize + j;
  FF = (SIR_BD_F_onestep(ts[id], eta[id], theta));

  B[i] = B[i] + dt * FF; // numerical integration of F

  // numerical integration of variance
  B[ngrid + i] = B[ngrid + i] + dt * (FF * B[ngrid + i] + B[ngrid + i] * FF' + A' * diff_Var(ts[id],eta[id],theta) * A);
        }
        B[i] = matrix_exp(B[i]);
      }
 //   print(B[ngrid + 1]);
  return B;
  }



real coal_log_stan(vector w, vector C, vector y, vector tt, int[] tids, matrix SI,real [] theta,int nc){
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
  int<lower = 1> n;
  int<lower = 1> ngrid;
  int<lower = 1> gridsize;
  real period;
  real N;
  real ts[n];
  vector[ngrid + 1] tt;
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
  matrix<lower=1> [ngrid + 1,2] SI;
}



transformed parameters{
/**
  x_r, x_i unused varible for ode
  transform the parameters into
  eta: Ode trajectory for S and I
  B: linear transform matrices in LNA
  S: Covariance matrices in LNA
*/
 vector[2] Init_State;
 matrix[n,2] eta;
 real betatN[n];
 matrix [2,2] FT[2*ngrid];
  real theta[6];

 theta = {R0,gamma,mu,A,period,N};
 Init_State[1] = N * Alpha / (1 + Alpha);
 Init_State[2] = N - Init_State[1];
 eta = SIR_BD_ODE_stan(ts, to_row_vector(Init_State), theta, n);
 FT = SIR_FT(ts,eta,gridsize,ngrid, R0, gamma, mu, A, period,N);
}







model{

  // specify the prior probability for the parameters

  R0 ~ uniform(3,8);
 // gamma ~ lognormal(-3,0.4);
  //Alpha ~ lognormal(1,10);
  gamma ~ normal(0.06,0.01);
  Alpha ~ normal(4,1);
  A ~ uniform(0,1);
  // construct the latent LNA trajector based on Ode trajectory, B and S
  for(i in 1:(ngrid+1)){
    if(i == 1){
      // this should be SI[1] = eta[1], however stan does not allows that, SI must be generated
      //from some distriubiton, so I choose a distribution approximation equals to
      // eta[1]
      SI[1] ~ multi_normal(to_vector(eta[1]), diag_matrix(to_vector({0.001,0.001})));
    }else{
      //print(to_vector(eta[1 + (i-1) * gridsize]) + FT[i-1] * (SI[i-1] - eta[1 + (i-2) * gridsize])');
      SI[i] ~ multi_normal(to_vector(eta[1 + (i-1) * gridsize]) + FT[i-1] * (SI[i-1] - eta[1 + (i-2) * gridsize])', FT[i + ngrid - 1]);
    }


  }

  //print(eta[n]' + FT[ngrid] * (SI[ngrid]-eta[n-gridsize])');
  // the likelihood model between observed data and latent trajectory

  // coalescent model
  increment_log_prob(coal_log_stan(w, C, y, to_vector(ts), tids, SI, theta, nc));

  // can also plugin other models
  // eg. incidence data, prevalence data

}


