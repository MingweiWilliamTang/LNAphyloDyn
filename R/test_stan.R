model_code <- 'functions{

row_vector SIR_BD_period_ODE_onestep(real t,
vector X,
real [] theta)
{
	/**
theta = c(R0,gamma, mu, A,period,N)
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

matrix SIR_BD_ODE_stan(real[] ts, row_vector Init_State, real[] theta, int N){
/**
Integrate the ode path
*/

   real dt;
   matrix [N,2] traj;
  row_vector[2] k1;
  row_vector[2] k2;
  row_vector[2] k3;
  row_vector[2] k4;
  dt = ts[2] - ts[1];
  traj[1] = Init_State;
  for(i in 2:N){
     k1 = SIR_BD_period_ODE_onestep(ts[i-1],to_vector(traj[i-1]), theta);
     k2 = SIR_BD_period_ODE_onestep(ts[i-1],to_vector(traj[i-1] + k1 * dt / 2), theta);
     k3 = SIR_BD_period_ODE_onestep(ts[i-1],to_vector(traj[i-1] + k2 * dt / 2), theta);
     k4 = SIR_BD_period_ODE_onestep(ts[i-1],to_vector(traj[i-1] + k3 * dt / 2), theta);
     traj[i] = traj[i-1] + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
    }
  return traj;
}

	matrix SIR_BD_F_onestep(real t, row_vector X, real [] theta){
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
B[i] = B[i] + dt * FF;
B[ngrid + i] = B[ngrid + i] + dt * (FF * B[ngrid + i] + B[ngrid + i] * FF\' + A\' * diff_Var(ts[id],eta[id],theta) * A) ;
        }
        B[i] = matrix_exp(B[i]);
      }
return B;
}

real coal_log_stan(vector w, vector C, vector y, vector tt, int[] tids, vector[] SI,real [] theta,int nc){
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
    print(Ne_inv);
  }

  return log_like;
}

}


model{}

'

#expose_stan_functions(stanc(model_code = model_code))
