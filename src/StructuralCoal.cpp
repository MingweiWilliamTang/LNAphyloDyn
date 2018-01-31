// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>

//#include <RcppArmadillo.h>
//#include "basic_util.h"
//#include "SIR_phylodyn.h"
#include "LNA_functional.h"

using namespace arma;
using namespace Rcpp;

static const int SAMPLE  = 0;
static const int CO = 1;

typedef std::vector<double> state_type;

arma::mat finite_size_correction2(const vec& p_a, const vec& A, const std::vector<bool> extant, arma::mat P)
{
  // NOTE mstates m X n
  int u;
  vec rho;
  vec rterm;
  //~ vec lterm;
  double lterm;
  vec p_u;
  vec Amin_pu;
  //~ for (int iu = 0; iu < extantLines.size(); iu++){
  for (int iu = 0; iu < extant.size(); iu++){
    if (extant.at(iu)){
      u = iu + 1;
      p_u = P.col(u-1);
      Amin_pu = clamp(( A - p_u), 1., INFINITY );
      rterm = p_a / Amin_pu ;
      rho = A / Amin_pu;
      lterm = dot( rho, p_a); //
      p_u = p_u % (lterm - rterm) ;
      p_u = p_u / sum(p_u ) ;
      P.col(u - 1) = p_u;
    }
  }
  return P;
}

void finite_size_correction3(const int ia, const vec& A, const std::vector<bool> extant, mat& P)
{
  // NOTE mstates m X n
  int u;
  vec rho;
  vec rterm;
  //~ vec lterm;
  double lterm;
  vec p_u;
  vec Amin_pu;
  vec p_a = P.col( ia );
  for (int iu = 0; iu < extant.size(); iu++){
    if (extant.at(iu) && iu!=ia){
      u = iu + 1;
      p_u = P.col(u-1);
      Amin_pu = clamp(( A - p_u), 1., INFINITY );
      rterm = p_a / Amin_pu ;
      rho = A / Amin_pu;
      lterm = dot( rho, p_a); //
      p_u = p_u % clamp((lterm - rterm), 0., INFINITY ) ; // l > r
      //~ if (sum(p_u)<=0.){
      //~ cout << " fsc3 messed up " << endl;
      //~ cout << p_u << endl;
      //~ cout << p_a << endl;
      //~ cout << lterm << endl;
      //~ cout << rterm << endl;
      //~ cout << Amin_pu << endl;
      //~ cout << rho << endl;
      //~ cout  << endl;
      //~ cout  << endl;
      //~ } //TODO check eval of A before calling this func
      p_u = p_u / sum(p_u ) ;
      P.col(u - 1) = p_u;
    }
  }
}


class DQAL{
  List Fs, Gs, Ys;
  int m;
  double hres;
  double treeT;
public:
  DQAL( List Fs, List  Gs, List Ys, int m, double hres, double treeT ) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT) {};
  void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
  {
    // time index
    //~ int i =  (int)min( 1+(int)( hres * (*t) / treeT ), hres);
    //~ NOTE hres = length(times)
    //~ NOTE treeT = max(times)
    int i =  (int)std::max(0., std::min( hres * t / treeT , (double)(hres-1.))); //TODO redefine max min? will use arma versin?
    mat F = as<mat>(Fs[i]);
    mat G = as<mat>(Gs[i]);
    vec Y = as<vec>(Ys[i]);

    int k,l,z,w;

    double a[m]; //normalized nlft
    //~ double sumA = 0.;
    //~ for (k = 0; k < m; k++) sumA += A(k);
    double r = 1. ; // Atotal / sumA; // TODO
    for (k = 0; k < m; k++) {
      //dA(k) = 0.;
      dxdt[Aind(k)] = 0.;
      if (Y(k) > 0) {
        //~ a[k] = r *  std::min(1., A(x, k)/ Y(k) );
        a[k] = r *  A(x, k)/ Y(k) ;
      } else{
        a[k] = 1.; //
      }
    }

    //dA
    // TODO shouldn't sum(A) be conserved over each interval?? then these are not the right eqns...
    for (k = 0; k < m; k++){
      for (l = 0; l < m; l++){
        if (k==l){
          //dA(k) -= a[l] * (F(i,l,k)) * a[k];
          dxdt[ Aind(k) ] -= a[l] * F(l,k) * a[k];
        } else{
          dxdt[ Aind(k) ]  += ( std::max(0., (1 - a[k])) * F(k,l) + G(k,l)) * a[l] ;
          dxdt[ Aind(k) ]  -= (F(l,k) + G(l,k)) * a[k];
        }
      }
    }
    //dQ
    for (z = 0; z < m; z++){ // col of Q
      for (k = 0; k < m; k++){ //row of Q
        dxdt[ Qind(k, z ) ] = 0. ;
        for (l = 0. ; l < m; l++){
          if (k!=l){
            if ( Q(x, l,z) > 0)
            {
              dxdt[ Qind(k, z ) ] += (F(k,l) + G(k,l)) *  Q(x, l,z)/  std::max(Q(x,l,z), Y(l));
            }
            if (Q(x, k,z) > 0)
            {
              dxdt[ Qind(k, z ) ] -= (F(l,k) + G(l,k)) *  Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
            }
          }
          // coalescent:
          //~ dQ(k,z) -= (F(i,k,l)+F(i,l,k)) * a[l] * Q(k,z)/Y(i,k);
          if (Q(x, k,z) > 0){
            dxdt[ Qind(k, z ) ] -= F(k,l) * a[l] * Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
          }
        }
      }
    }
    //dL
    double dL = 0.;
    double Ydenom;
    for (k = 0; k < m; k++){
      for (l =0 ; l < m; l++){
        if (k==l){
          if (Y(k) > 1 && A(x,k) > 1){
            Ydenom = std::max( (Y(k)-1.),(r*A(x,k)-1.) );//
            if (Ydenom > 0)
            {
              dL += std::max(0.,std::min(1.,a[k])) * (  (r*A(x, k)-1)/Ydenom) * F(k,l);
            }
          } else{
            dL += std::max(0.,std::min(1.,a[k])) * std::max(0.,std::min(1.,a[l])) * F(k,l);
          }
        } else{
          dL += std::max(0.,std::min(1.,a[k])) * std::max(0.,std::min(1.,a[l])) * F(k,l);
        }
      }
    }
    dL = std::max(dL, 0.);
    dxdt[Lind()] = dL;
  }

private:
  //~ double F(int ih, int k, int l){
  //~ return as<mat>(Fs[ih]).at(k,l);
  //~ }
  //~ double G(int ih, int k, int l){
  //~ return as<mat>(Gs[ih]).at(k,l);
  //~ }
  //~ double Y(int ih, int k){
  //~ return as<vec>(Ys[ih]).at(k);
  //~ }

  double Q( const state_type &x, int k, int l) {
    return x[l*m+k];
  }
  double A( const state_type &x, int k) {
    return x[(int)pow(m,2) + k];
  }
  double L( const state_type &x ) {
    return x[ (int)pow(m,2) + m] ;
  }

  int Qind( int k, int l ){
    return l*m + k ;
  }
  int Aind( int k ){
    return (int)pow(m,2) + k;
  }
  int Lind(){
    return (int)pow(m,2) + m ;
  }
};


state_type generate_initial_conditions( vec A){
  int m = A.size() ;
  state_type x( (int)pow(m,2) + m + 1, 0. );
  int k = 0;
  for (int i =0; i < m; i++){
    for (int j =0 ; j < m; j++){
      if (i == j){
        x[k] =1.;
      }
      k++;
    }
  }
  for (int i = 0; i < m; i++){
    x[k] = A.at(i) ;
    k++;
  }
  return x;
}

void Q_from_state( arma::mat &Q, state_type xfin ){
  int m = Q.n_rows;
  int k =0 ;
  for (int i = 0; i < m; i++){
    for (int j = 0; j < m; j++){
      Q.at(i,j) = xfin[k];
      k++;
    }
  }
  for (int i = 0; i < m; i++){
    Q.row(i) = Q.row(i) / sum(Q.row(i));
  }
  //~ Q = Q / arma::sum(Q, 1) ;
}

void A_from_state( vec &A, state_type xfin){
  int m = A.size();
  int k = 0;
  for ( int i = (int)pow(m,2); i < ((int)pow(m,2)+m); i++){
    A.at(k) = xfin[i];
    k++;
  }
}

double L_from_state( state_type xfin){
  return (double)xfin[xfin.size()-1];
}

//[[Rcpp::export()]]
double colik2cpp(const NumericVector heights, const List Fs, const List Gs, const List Ys
                   , const IntegerVector eventIndicator // sample or co
                   , const IntegerVector eventIndicatorNode // node involved at each event
                   , const NumericVector eventHeights
                   , const arma::mat sortedSampleStates
                   , const IntegerMatrix daughters // daughters of each node
                   , const int n
                   , const int Nnode
                   , const int m
                   , double AgtYboundaryCondition )
{
  double loglik = 0.;
  mat P(m, n + Nnode, fill::zeros);
  int nextant = 0;
  int u,v,w,z,a;
  int samplesAdded = 0;
  mat Q  = zeros(m, m );
  vec A_Y;
  vec Y ;
  vec A = zeros(m);
  mat F = zeros(m,m);
  mat G = zeros(m,m);
  vec puY, pvY, pa;
  vec m_rs_R;
  std::vector<bool> extant(n + Nnode, false);

  // instantiate solver
  double hres =  heights.size() ;
  double treeT = heights[heights.size()-1]; // note increasing order
  DQAL dqal(Fs, Gs, Ys, m, hres, treeT );
  state_type x0;//, x;
  //~ std::vector<state_type> x_vec;
  //~ std::vector<double> x_vec_times;

  // iterate events
  double nextEventHeight = eventHeights(0);
  int ievent = 0;

  double sterm0, sterm1;
  double h, h0, h1, hstar, dh;
  int ih; // index of nextEventHeight
  double default_step_size = std::abs(heights[1] - heights[0]);
  double L = 0.;

  h = 0.;
  while( nextEventHeight != INFINITY ){

    if (nextEventHeight > h ){
      A = sum(P, 1);
      A = arma::normalise(A,1.)  * ((double)nextant);
      x0 = generate_initial_conditions(A) ;
      //~ cout << "solving..." << endl;
      //~ cout << h << endl;
      //~ cout << nextEventHeight << endl;
      //~ cout << A << endl;
      //~ Rf_PrintValue( wrap(x0));
      //~ size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight , default_step_size, simple_observer(x)  );
      size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight
                                                          , std::min( default_step_size,  (nextEventHeight-h)/10.) );
      //~ size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight , default_step_size
      //~ , push_back_state_and_time( x_vec, x_vec_times)  );
      //~ size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight , (nextEventHeight-h)/10.
      //~ , push_back_state_and_time( x_vec, x_vec_times)  );
      //~ Rf_PrintValue( wrap(x0));
      //~ cout << steps << endl;
      //~ throw 1;
      //
      Q_from_state(Q, x0);
      A_from_state(A, x0);
      L = L_from_state(x0);
      A = arma::normalise(A,1.) * ((double)nextant);

      P  = abs(Q.t() * P);
      //~ cout << P.n_rows << " " << P.n_cols << endl;
      //~ cout << P.t() ;
    } else{
      L = 0.;
    }

    loglik -= std::max(0.,L);

    if (eventIndicator(ievent)==SAMPLE)
    {
      u = eventIndicatorNode(ievent);
      nextant++;
      //~ P.col(u-1) = arma::normalise(sortedSampleStates.col(samplesAdded)); //WRONG! 2-norm is default
      P.col(u-1) = arma::normalise(sortedSampleStates.col(samplesAdded), 1.);
      extant.at(u-1) = true;
      samplesAdded++;
    } else{
      // coalescent
      ih = (int)std::min( hres * nextEventHeight / treeT, (double)(heights.size()-1));
      F = as<mat>(Fs[ih]);
      G = as<mat>(Gs[ih]);
      Y = as<vec>(Ys[ih]);
      //~ if (sum(A) > sum(Y)) {
      //~ cout << Y ;
      //~ cout << nextEventHeight << " " << ih << endl;
      //~ cout << " sum A > sum Y  " << endl ;
      //~ }
      //~ if ( AgtYboundaryCondition && (sum(A) > sum(Y)) ) return -INFINITY;
      if (nextant > sum(Y)){
        loglik -= L * (AgtYboundaryCondition * (nextant-sum(Y)));
      }

      a = eventIndicatorNode(ievent);
      u = daughters(a -1 , 0);
      v = daughters( a - 1, 1 );
      puY = arma::normalise(  arma::min(Y,P.col(u- 1 )) ,1.) / arma::clamp(Y, 1e-6, INFINITY ) ;
      pvY = arma::normalise(  arma::min(Y,P.col(v- 1 )) ,1.) / arma::clamp( Y, 1e-6, INFINITY ) ;
      //~ loglik += log( (( puY * F) * pvY).at(0,0) ) ;
      //~ double llterm = log(  as_scalar( puY.t() * (F * pvY) )  + as_scalar( pvY.t() * (F * puY) ) ) ;
      //~ if (llterm==-INFINITY)
      //~ {
      //~ cout << " F " << endl;
      //~ cout << F;
      //~ cout << " G " << endl;
      //~ cout << G;
      //~ cout << " Y " << endl;
      //~ cout << Y;
      //~ cout << " puy " << endl;
      //~ cout << puY;
      //~ cout << " pvy " << endl;
      //~ cout << pvY;
      //~ cout << nextEventHeight << " " << ih << endl;
      //~ cout << heights.size() << endl;
      //~ cout << hres  << endl;
      //~ cout << treeT  << endl;
      //~ }
      //~ cout << " co event " << endl;
      //~ cout << puY;
      //~ cout << pvY;
      //~ cout << F;
      loglik += log(  as_scalar( puY.t() * (F * pvY) )  + as_scalar( pvY.t() * (F * puY) ) ) ;
      // state of ancestor
      pa =  arma::normalise( (F * puY) % pvY + (F * pvY) % puY ,1.) ;
      //~ pa = pa / sum(pa ) ;
      P.col(a - 1 ) = pa;
      //~ cout << loglik << endl;
      //~ cout << pa;
      //~ cout << endl << endl ;
      //~ if (any(is_nan(as<NumericVector>(wrap(pa))))) throw 1;
      P = finite_size_correction2(pa, A, extant, P);//TODO test
      // TODO make FSC optional?

      //bookkeeping
      nextant--;
      extant.at(u-1) = false;
      extant.at(v-1) = false;
      extant.at(a-1) = true;
      P.col(u-1) = zeros<colvec>(P.n_rows);
      P.col(v-1) = zeros<colvec>(P.n_rows);
    }
    //~ cout << endl;
    //~ cout << endl;
    //~ cout << endl;
    //~ cout << h << endl;
    //~ cout << nextEventHeight << endl;
    //~ cout << L << endl;
    //~ cout << loglik << endl;
    //~ cout << A << endl;
    //~ cout << Q << endl;

    //~ if (h > 1) throw 1;
    // prep next iter
    h = nextEventHeight;
    ievent++;
    if (ievent<eventHeights.size()){
      nextEventHeight = eventHeights(ievent);
    } else{
      nextEventHeight = INFINITY;
    }
  }

  return loglik;
}


////////////////////////////////////////////////////////////////////////////////
//~ TODO the following should be moved to phylo.source.attribution.cpp and should ref a
// dynamic library made out of colik (see colik.hpp).


// transition probabilities and survivor functions, conditioning on lineage not changing host
class DQpsiA{
  List Fs, Gs, Ys;
  int m;
  double hres;
  double treeT;
public:
  DQpsiA( List Fs, List  Gs, List Ys, int m, double hres, double treeT ) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT) {};
  void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
  {
    // time index
    //~ NOTE hres = length(times)
    //~ NOTE treeT = max(times)
    int i =  (int)std::max(0., std::min( hres * t / treeT , (double)(hres-1.)));
    mat F = as<mat>(Fs[i]);
    mat G = as<mat>(Gs[i]);
    vec Y = as<vec>(Ys[i]);

    int k,l,z,w;

    double a[m]; //normalized nlft
    for (k = 0; k < m; k++) {
      dxdt[Aind(k)] = 0.;
      if (Y(k) > 0) {
        a[k] = A(x, k)/ Y(k) ;
      } else{
        a[k] = 1.; //
      }
    }

    //dA
    for (k = 0; k < m; k++){
      for (l = 0; l < m; l++){
        if (k==l){
          //dA(k) -= a[l] * (F(i,l,k)) * a[k];
          dxdt[ Aind(k) ] -= a[l] * F(l,k) * a[k];
        } else{
          dxdt[ Aind(k) ]  += ( std::max(0., (1 - a[k])) * F(k,l) + G(k,l)) * a[l] ;
          dxdt[ Aind(k) ]  -= (F(l,k) + G(l,k)) * a[k];
        }
      }
    }
    //dQ
    // NOTE this is filling role of "rho" in the paper
    for (z = 0; z < m; z++){ // col of Q
      for (k = 0; k < m; k++){ //row of Q
        dxdt[ Qind(k, z ) ] = 0. ;
        for (l = 0. ; l < m; l++){
          if (k!=l){
            if ( Q(x, l,z) > 0)
            {
              dxdt[ Qind(k, z ) ] += ( G(k,l)) *  Q(x, l,z)/  std::max(Q(x,l,z), Y(l));
            }
            if (Q(x, k,z) > 0)
            {
              dxdt[ Qind(k, z ) ] -= ( G(l,k)) *  Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
            }
          }
        }
      }
    }
    // d psi
    for (z = 0; z < m; z++){ // col of Q, elem of psi- line starts in deme z
      dxdt[ psi_ind( z ) ] = 0. ;
      for (k = 0; k < m; k++){ //row of Q - current deme
        for (l = 0. ; l < m; l++){ //source of transm
          // TODO possibly should have (Y_l - A_l)/Y_l terms as in paper
          //~ dxdt[psi_ind(k,z)] -= F(i,l,k) * psi(x, k, z) / std::max(psi(x, k, z), Y(k));
          // NOTE Q is filling the role of rho_{ik} in the paper
          dxdt[psi_ind(z)] -= std::max(0., (1 - a[l])) * F(l,k) * std::max(0., psi(x, z)) * Q(x,k,z) / std::max(Q(x, k, z), Y(k));
        }
      }
    }
    // NO
    //~ for (k = 0; k < m; k++){
    //~ for (l = 0; l < m; l++){
    //~ dxdt[psi_ind( l )] -= F(i,l,k) * psi(x, k) / std::max(psi(x,k), Y(k));
    //~ }
    //~ }
  }
private:
  double Q( const state_type &x, int k, int l) {
    return x[l*m+k];
  }
  double psi( const state_type &x, int k) {
    return x[(int)pow(m,2) +  k];
  }
  double A( const state_type &x, int k) {
    return x[(int)pow(m,2) + m + k];
  }
  int Qind( int k, int l ){
    return l*m + k ;
  }
  int psi_ind( int k ){
    return (int)pow(m,2) + k;
  }
  int Aind( int k ){
    return (int)pow(m,2) + m + k;
  }
};



state_type generate_initial_conditions_dqpsia(vec A){
  int m = A.size() ;
  state_type x( (int)pow(m,2) + m + m, 0. );
  int k = 0;
  for (int i =0; i < m; i++){
    for (int j =0 ; j < m; j++){
      if (i == j){
        x[k] =1.; // Q or 'rho'
      }
      k++;
    }
  }
  for (int i = 0; i  < m ; i++){
    x[k] = 1.; //psi
    k++;
  }
  for (int i = 0; i < m; i++){
    x[k] = A.at(i) ; //A
    k++;
  }
  return x;
}


void Qrho_from_state( arma::mat &Qrho, state_type xfin ){
  int m = Qrho.n_rows;
  int k =0 ;
  for (int i = 0; i < m; i++){
    for (int j = 0; j < m; j++){
      Qrho.at(i,j) = std::min(1.,std::max(0., xfin[k]));
      k++;
    }
  }
  for (int i = 0; i < m; i++){
    Qrho.row(i) = Qrho.row(i) / sum(Qrho.row(i));
  }
}


void psi_from_state( vec &psi, state_type xfin ){
  int m = psi.size();
  int k = 0;
  for ( int i = (int)pow(m,2); i < ((int)pow(m,2)+m); i++){
    psi.at(k) = std::min(1.,std::max(0., xfin[i]));
    k++;
  }
}

//[[Rcpp::export()]]
List sourceAttribMultiDemeCpp( const NumericVector heights, const List Fs, const List Gs, const List Ys
                                 , const IntegerVector eventIndicator // sample or co
                                 , const IntegerVector eventIndicatorNode // node involved at each event
                                 , const NumericVector eventHeights
                                 , const arma::mat sortedSampleStates
                                 , const IntegerMatrix daughters // daughters of each node
                                 , const int n
                                 , const int Nnode
                                 , const int m
                                 , double AgtYboundaryCondition
                                 , const double maxHeight // terminate computation at this height
){
  double loglik = 0.;
  mat P(m, n + Nnode, fill::zeros);
  mat rho(m, n , fill::zeros);
  int nextant = 0;
  int u,v,w,z,a, k, l;
  int samplesAdded = 0;
  mat Q  = zeros(m, m );
  mat Qrho  = zeros(m, m );
  //~ vec Psi = zeros(n, n + Nnode); //records probability i'th sample is host at j'th node
  vec psi = zeros(m) ; // corresponding to each state over each interval
  vec psi_time = zeros(n); // for each tip, varies over time

  vec A_Y;
  vec Y ;
  vec A = zeros(m);
  mat F = zeros(m,m);
  mat G = zeros(m,m);
  vec puY, pvY, pa;
  vec m_rs_R;
  std::vector<bool> extant(n + Nnode, false);
  umat tipsDescendedFrom(n+Nnode, n, fill::zeros);
  for (u = 0; u < n; u++){
    tipsDescendedFrom.at(u,u) = 1;
  }

  // container for output
  //~ mat W = zeros(n,n);
  std::vector<int> donorW;
  donorW.reserve( (int)(n * std::sqrt((double)n) ) );
  std::vector<int> recipW;
  recipW.reserve( (int)(n * std::sqrt((double)n) ) );
  std::vector<double> W;
  W.reserve( (int)(n * std::sqrt((double)n) ) );

  // instantiate solver
  double hres =  heights.size() ;
  double treeT = heights[heights.size()-1]; // note increasing order
  DQAL dqal(Fs, Gs, Ys, m, hres, treeT );
  state_type x0;//, x;

  DQpsiA dqpsia( Fs, Gs, Ys, m, hres, treeT );
  state_type x0sa;

  // iterate events
  double nextEventHeight = eventHeights(0);
  int ievent = 0;

  double sterm0, sterm1;
  double h, h0, h1, hstar, dh;
  int ih; // index of nextEventHeight
  double default_step_size = std::abs(heights[1] - heights[0]);
  double L = 0.;

  //
  //~ boost::numeric::odeint::euler< state_type > eu_stepper;

  h = 0.;
  while( nextEventHeight != INFINITY && nextEventHeight < maxHeight ){
    //~ std::cout << h << " " << nextEventHeight << " " << maxHeight << " "  << std::endl;
    if (nextEventHeight > h )
    {
      A = sum(P, 1);
      A = arma::normalise(A,1.)  * ((double)nextant);
      x0 = generate_initial_conditions(A) ;
      //~ cout << "solving..." << endl;
      //~ cout << h << endl;
      //~ cout << nextEventHeight << endl;
      //~ cout << A << endl;
      //~ Rf_PrintValue( wrap(x0));
      size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight
                                                          , std::min( default_step_size,  (nextEventHeight-h)/10.) );
      //~ size_t steps = boost::numeric::odeint::integrate_const(eu_stepper, dqal ,  x0 , h , nextEventHeight
      //~ , std::min( default_step_size,  (nextEventHeight-h)/10.) );
      //~ Rf_PrintValue( wrap(x0));
      //~ cout << steps << endl;
      //~ throw 1;
      //

      //  same for sa ...
      x0sa = generate_initial_conditions_dqpsia( A) ;
      //~ Rf_PrintValue( wrap(x0sa) );
      steps = boost::numeric::odeint::integrate( dqpsia ,  x0sa , h , nextEventHeight
                                                   , std::min( default_step_size,  (nextEventHeight-h)/10.) );
      //~ steps = boost::numeric::odeint::integrate_const(eu_stepper, dqpsia ,  x0sa , h , nextEventHeight
      //~ , std::min( default_step_size,  (nextEventHeight-h)/10.) );
      //~ Rf_PrintValue( wrap( x0sa));
      //~ throw 1;
      Q_from_state(Q, x0);
      A_from_state(A, x0);
      L = L_from_state(x0);
      A = arma::normalise(A,1.) * ((double)nextant);

      Qrho_from_state( Qrho, x0sa );
      psi_from_state( psi, x0sa );

      //update psi_time
      // note need to run this before rho updated
      //~ double psi_factor;
      //~ for (u = 0; u < n ; u++){
      //~ psi_factor = 0.;
      //~ for (k = 0; k < m; k++){
      //~ psi_factor += rho.at(k,u) * psi.at(k);
      //~ }
      //~ //psi_factor /= (double)m;
      //~ psi_time.at(u) *= psi_factor;
      //~ }
      //vec psi_factor =  rho.t() * psi  ;
      psi_time = psi_time % (rho.t() * psi );

      P  = abs(Q.t() * P);
      P = arma::normalise( P, 1., 0 );
      rho = abs( Qrho.t() * rho );
      rho = arma::normalise( rho, 1., 0); // p=1, dim=0
      //~ std::cout << h <<std::endl;
      //~ std::cout << P ;
      //~ std::cout << rho ;
      //~ std::cout << psi;
      //~ std::cout << Qrho;
      //~ std::cout << psi_time ;
      //~ std::cout << std::endl ;
      //~ std::cout << std::endl ;
      //~ std::cout << std::endl ;
    } else{
      L = 0.;
    }

    if (eventIndicator(ievent)==SAMPLE)
    {
      u = eventIndicatorNode(ievent);
      nextant++;
      P.col(u-1) = arma::normalise(sortedSampleStates.col(samplesAdded) ,1.);
      rho.col(u-1) = P.col(u-1);
      psi_time.at(u-1) = 1.;
      extant.at(u-1) = true;
      samplesAdded++;
    } else {
      // coalescent
      ih = (int)std::min( hres * nextEventHeight / treeT, (double)(heights.size()-1));
      F = as<mat>(Fs[ih]);
      G = as<mat>(Gs[ih]);
      Y = as<vec>(Ys[ih]);
      Y = arma::clamp(Y, 1e-6, INFINITY );
      //~ if (sum(A) > sum(Y)) {
      //~ cout << Y ;
      //~ cout << nextEventHeight << " " << ih << endl;
      //~ cout << " sum A > sum Y  " << endl ;
      //~ }

      a = eventIndicatorNode(ievent);
      u = daughters(a -1 , 0);
      v = daughters( a - 1, 1 );
      puY = arma::normalise(  arma::min(Y,P.col(u- 1 )) ,1.) / arma::clamp( Y, 1e-6, INFINITY ) ;
      pvY = arma::normalise(  arma::min(Y,P.col(v- 1 )) ,1.) / arma::clamp( Y, 1e-6, INFINITY ) ;
      //~ loglik += log( (( puY * F) * pvY).at(0,0) ) ;

      // state of ancestor
      pa =  arma::normalise( (F * puY) % pvY + (F * pvY) % puY ,1.) ;
      //~ pa = pa / sum(pa ) ;
      P.col(a - 1 ) = pa;
      //P = finite_size_correction2(pa, A, extant, P);

      //bookkeeping
      nextant--;
      extant.at(u-1) = false;
      extant.at(v-1) = false;
      extant.at(a-1) = true;
      P.col(u-1) = zeros<colvec>(P.n_rows);
      P.col(v-1) = zeros<colvec>(P.n_rows);

      // sa stuff ; upate psi , rho
      vec rho_w__Y = zeros(m);
      vec rho_z__Y = zeros(m);
      // update W(iw, iz)
      for (int iw = 0; iw < n; iw++){
        if (tipsDescendedFrom.at(u-1, iw)==1){
          rho_w__Y = arma::normalise(  arma::min(Y,rho.col(iw)) ,1.) / arma::clamp(Y, 1e-6, INFINITY ) ;
          for (int iz = iw+1; iz < n; iz++){
            if (tipsDescendedFrom.at(v-1, iz)==1){
              rho_z__Y = arma::normalise(  arma::min(Y,rho.col(iz)) ,1.) / arma::clamp(Y, 1e-6, INFINITY ) ;
              double pwz0 = sum( rho_w__Y % (F * rho_z__Y) );
              double pzw0 = sum( rho_z__Y % (F * rho_w__Y ));
              double pwz = psi_time.at(iw) * psi_time.at(iz) * pwz0 / (pwz0 + pzw0) ;
              double pzw = psi_time.at(iw) * psi_time.at(iz) *  pzw0 / (pwz0 + pzw0 ) ;

              donorW.push_back( iw + 1);
              recipW.push_back( iz+1 );
              W.push_back( pwz) ;

              donorW.push_back( iz + 1 );
              recipW.push_back( iw+ 1);
              W.push_back( pzw );

              //~ W.at( iw, iz) = pwz;
              //~ W.at(iz, iw) = pzw;
            }
          }
        }
      }
      tipsDescendedFrom.row(a-1) = tipsDescendedFrom.row(u-1) + tipsDescendedFrom.row(v-1);
      //update rho and psi
      for (int iw = 0; iw < n; iw++){// w & z tips
        if (tipsDescendedFrom.at(u-1, iw)==1){
          //tip descended from a and u
          rho_w__Y = arma::normalise( arma::clamp(  arma::min(Y,rho.col(iw))  / Y, 1e-6, 1. ) ,1.) ;
          //update psi(iw)
          double puv = sum( rho_w__Y % (F * pvY) );
          double pvu = sum( pvY % (F * rho_w__Y ));
          puv = puv / (puv + pvu ) ;
          psi_time.at(iw) *= puv ;
          //update rho(iw)
          rho.col(iw) = arma::normalise( arma::clamp( rho_w__Y % (F * pvY) , 1e-6, 1.) ,1.);
        } else if (tipsDescendedFrom.at(v-1, iw)==1){
          //tip descended from a and v
          rho_w__Y = arma::normalise(arma::clamp(  arma::min(Y,rho.col(iw))  / Y, 1e-6, 1. ) ,1.) ;
          //update psi(iw)
          double pvu = sum( rho_w__Y % (F * puY) );
          double puv = sum( puY % (F * rho_w__Y ));
          pvu = pvu / (puv + pvu ) ;
          psi_time.at(iw) *= pvu ;
          //update rho(iw)
          rho.col(iw) = arma::normalise( arma::clamp( rho_w__Y % (F * puY), 1e-6, 1.) ,1.);
        }
      }

    }
    h = nextEventHeight;
    ievent++;
    if (ievent<eventHeights.size()){
      nextEventHeight = eventHeights(ievent);
    } else{
      nextEventHeight = INFINITY;
    }
    if ( h > maxHeight ) {
      break;
    }
  }
  List wlist;
  wlist["donor"] = wrap(donorW);
  wlist["recip"] = wrap(recipW);
  wlist["infectorProbability"] = wrap(W);
  return wlist;
}



//[[Rcpp::export()]]
arma::mat birthMx(arma::vec states, arma::vec thetas, arma::vec x_r, std::string model = "SEIR"){
  arma::mat M;
  if(model == "SEIR"){
    M << 0 << 0 << arma::endr
      << thetas[0] * states[0] * states[2] << 0 << arma::endr;
  }else{
  M << 0 << 0 << arma::endr
    << thetas[0] * x_r[0] * states[1] << 0 << arma::endr;
  }
  return M;
}


//[[Rcpp::export()]]
arma::mat MigrationMx(arma::vec states, arma::vec thetas, arma::vec x_r,std::string model = "SEIR"){
  arma::mat M;
  if(model == "SEIR"){
    M << 0 << thetas[1] * states[1] << arma::endr << 0 << 0 << arma::endr;
  }else{
    M << 0 << thetas[1] * states[0] << arma::endr << 0 << 0 << arma::endr;
  }
  return M;
}

//[[Rcpp::export()]]
Rcpp::NumericVector DeathMX(arma::vec states, arma::vec thetas, arma::vec x_r, std::string model = "SEIR"){
  NumericVector D(2);
  D(0) = 0;
  if(model == "SEIR"){
    D(1) = states(2) * thetas[2];
  }else{
    D(1) = states(1) * thetas[2];
  }
  return D;
}

//[[Rcpp::export()]]
List TFGY_list(arma::mat LatentTraj, arma::vec param,
               arma::vec x_r, arma::ivec x_i,
               std::string transP = "changepoint", std::string model = "SEIR",std::string transX = "standard"){
  int n = LatentTraj.n_rows;
  int p = LatentTraj.n_cols - 1;
  double t;
  List Births;
  List Migrations;
  List Deaths;
  List Ys;
  arma::mat trajm(n,p+1);
  trajm.col(0) = LatentTraj.col(0);
  trajm.col(1) = LatentTraj.col(2);
  trajm.col(2) = LatentTraj.col(3);
  trajm.col(3) = LatentTraj.col(1);
  arma::vec times(n);
  XPtr<parat> param_trans = transformPtr(transP);
  for(int i = n - 1; i >=0; i --){
    t = LatentTraj(i,0);
    times(n - 1 - i) = t;
    arma::vec states = LatentTraj.submat(i,1,i,p).t();
    arma::vec state2 = LatentTraj.submat(i,2,i,p).t();
    arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
    Births.push_back(birthMx(states, thetas, x_r));
    Migrations.push_back(MigrationMx(states, thetas, x_r));
    Ys.push_back(Rcpp::NumericVector(state2.begin(), state2.end()));
    Deaths.push_back(DeathMX(states, thetas, x_r));
  }

  List tfgy;
  //tfgy["Fs"] = Births;
  //tfgy["Gs"] = Migrations;
  //tfgy["Ds"] = Deaths;
  //tfgy["Ys"] = Ys;

  tfgy["times"] = Rcpp::NumericVector(times.begin(),times.end());
  tfgy["births"] = Births;
  tfgy["migrations"] = Migrations;
  tfgy["sizes"] = Ys;
  tfgy["traj"] = trajm;
  tfgy["deaths"] = Deaths;

  return tfgy;
}





//[[Rcpp::export()]]
List TFGY_list2(arma::mat LatentTraj, std::string model = "SEIR"){
  int n = LatentTraj.n_rows;
  int p = LatentTraj.n_cols - 1;
  double t;
  List Births;
  List Migrations;
  List Deaths;
  List Ys;
  List SZ;
  arma::vec times(n);
  for(int i = n - 1; i >=0; i --){
    t = LatentTraj(i,0);
    times(i) = t;
    arma::vec states = LatentTraj.submat(i,1,i,p).t();
    arma::mat B(2,2);
    arma::vec D(2);
    arma::mat M(2,2);
    if(i < n - 1){
      B(1,0) = LatentTraj(i,1) - LatentTraj(i+1, 1);
      M(0,1) = LatentTraj(i,2) - LatentTraj(i+1, 2) + B(1,0);
      D(1) = M(0,1) + LatentTraj(i,3) - LatentTraj(i+1, 3);
    }
    Births.push_back(B);
    Migrations.push_back(M);
    Ys.push_back(Rcpp::NumericVector(states.begin(), states.end()));
    Deaths.push_back(D);
  }
  List tfgy;
  tfgy["times"] = times;
  tfgy["births"] = Births;
  tfgy["migrations"] = Migrations;
  tfgy["sizes"] = Ys;
  tfgy["traj"] = LatentTraj;
  tfgy["deaths"] = Deaths;
  return tfgy;
}


//[[Rcpp::export]]
double Structural_Coal_lik(List Init_Details, arma::mat LatentTraj,
                           arma::vec param,
                           arma::vec x_r, arma::ivec x_i,
                           std::string transP = "changepoint",
                           std::string model = "SEIR2",
                           std::string transX = "standard",
                           double AgtYboundCondition = 1){
  /*
  * Init_Details: contains the informations related to the tree structure and the
  *
  */
  List tfgy = TFGY_list(LatentTraj, param, x_r, x_i, transP, model, transX);
  List Fs = tfgy["births"];
  List Gs = tfgy["migrations"];
  List Ys = tfgy["sizes"];
  NumericVector heights = Init_Details[0];
  IntegerVector eventIndicator = Init_Details[1];
  IntegerVector eventIndicatorNode = Init_Details[2];
  NumericVector eventHeights = Init_Details[3];
  arma::mat sortedSampleStates = Init_Details[4];
  IntegerMatrix daughters = Init_Details[5];
  int n = Init_Details[6];
  int Nnode = Init_Details[7];
  int m = Init_Details[8];

  return colik2cpp(heights, Fs, Gs, Ys
                     ,   eventIndicator // sample or co
                     , eventIndicatorNode // node involved at each event
                     , eventHeights
                     , sortedSampleStates
                     , daughters // daughters of each node
                     , n
                     , Nnode
                     ,  m
                     , AgtYboundCondition);
}




//[[Rcpp::export()]]
List ESlice_general_NC_Structural(arma::mat f_cur, arma::mat OdeTraj, List FTs,
                                  List init, arma::vec param, arma::vec x_r, arma::ivec x_i,
                                  double coal_log = 0,
                                  std::string model = "SEIR2",
                                  std::string transX = "standard"){
  // OdeTraj is the one with low resolution

  int p = f_cur.n_cols;
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;

  arma::mat v_traj = arma::randn(f_cur.n_rows, f_cur.n_cols);
  double u = R::runif(0,1);
  //  if(funname == "standard"){
  // logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
  if(coal_log != 0){
    logy = coal_log + log(u);
  }else{
    logy = Structural_Coal_lik(init, TransformTraj(OdeTraj,f_cur, FTs), param, x_r, x_i); + log(u);
  }

  double theta = R::runif(0,2 * pi);

  double theta_min = theta - 2*pi;
  double theta_max = theta;

  arma::mat f_prime = f_cur * cos(theta) + v_traj * sin(theta);
  newTraj = TransformTraj(OdeTraj,f_prime, FTs);
  int i = 0;
  double loglike;

  loglike = Structural_Coal_lik(init, newTraj, param, x_r, x_i);

  while(newTraj.cols(1,p).min() <0 || loglike <= logy){
    // shrink the bracket
    i += 1;
   // Rcout << -i << endl;
    if(i>20){
      newTraj = TransformTraj(OdeTraj,f_cur, FTs);
      loglike = Structural_Coal_lik(init, newTraj, param, x_r, x_i);
      f_prime = f_cur;
      Rcout<<"theta = "<<theta<<endl;
      break;
    }
    if(theta < 0){
      theta_min = theta;
    }else{
      theta_max = theta;
    }
    theta = R::runif(theta_min,theta_max);
    f_prime = f_cur * cos(theta) + v_traj * sin(theta);
    newTraj = TransformTraj(OdeTraj,f_prime, FTs);
    //Rcout << -i * 10 << endl;
    loglike = Structural_Coal_lik(init, newTraj, param, x_r, x_i);

  }
  List Res;
  double logOrigin = 0;
  for(int j = 0; j < f_cur.n_rows - 1; j ++){
    for(int k = 0; k < p; k ++){
      logOrigin -= 0.5 * f_prime(j,k) * f_prime(j,k);
    }
  }
  Res["LatentTraj"] = newTraj;
  Res["OriginTraj"] = f_prime;
  Res["logOrigin"] = logOrigin;
  Res["CoalLog"] = loglike;
  return Res;
}


