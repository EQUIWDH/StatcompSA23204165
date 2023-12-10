#include <Rcpp.h>
#include <random>

using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix gibbs(int a, int b, int n, int N) {
  NumericMatrix output(N,2);
  std::random_device rd;  // obtain a random number from hardware
  std::mt19937 gen(rd());
  std::gamma_distribution<> gamma1(0, 1);
  std::gamma_distribution<> gamma2(0, 1);
  std::binomial_distribution<> binomial(n, 0);
  for (int i = 1; i < N ; i++){
    binomial.param(std::binomial_distribution<>::param_type(n,output(i-1,1)));
    output(i,0) = binomial(gen);
    gamma1.param(std::gamma_distribution<>::param_type(a + output(i,1),1));
    gamma2.param(std::gamma_distribution<>::param_type(n-output(i,1)+b,1));
    output(i,1) =  gamma1(gen) / (gamma1(gen)+gamma2(gen));
  }
  return output;
  
}


