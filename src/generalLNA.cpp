#include "basic_util.h"
#include "SIR_phylodyn.h"
using namespace Rcpp;
using namespace arma;

//' @export Foo
class Foo{
  public:
    Foo(arma::mat x_, double y_, double z_ ):x(x_), y(y_), z(z_) {}
    arma::mat x;
    double y;
    double get_z() {return z;}
    arma::mat get_chol_x() {return chols(x);}
    void set_z( double z_ ) { z = z_; }private:double z;
      };

RCPP_MODULE(mod_Foo){
  using namespace Rcpp;
  class_<Foo>( "Foo" )
  .constructor<arma::mat,double,double>()
  .field( "x", &Foo::x )
  .field_readonly( "y", &Foo::y )
  .property( "z", &Foo::get_z, &Foo::set_z )
  .method( "get_chol_x", &Foo::get_chol_x)
  .method( "get_z", &Foo::get_z );
  }
