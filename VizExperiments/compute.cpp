#include <iostream>
#include <cmath>

const long double PI = 3.14159265359;

void compute_val(long double & xval, long double & yval, long long omega, long long n, long long y){
  long double xtotal = 0;
  long double ytotal = 0;
  xtotal = std::cos(2*PI*omega*y/n);
  ytotal = std::sin(2*PI*omega*y/n);
  long long omegapower = (omega*omega) % n;
  while (omegapower != 0 && omegapower != 1 && omegapower != omega){
    xtotal += std::cos(2*PI*omegapower*y/n);
    ytotal += std::sin(2*PI*omegapower*y/n);
    omegapower = (omega * omegapower) % n;
  }
  xval = xtotal;
  yval = ytotal;
}


void compute_all_vals(long long omega, long long n){
  long double x;
  long double y;
  for (long long i = 1; i < n; ++i){
    compute_val(x, y, omega, n, i);
    std::cout << x << "   " << y << "\n";
  }
}

int main(){
  //compute_all_vals(20719, 551905);
  compute_all_vals(3124, 478125);
  return 0;
}
