#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sumRcpp(NumericVector x) {
    int n = x.size();
    double res;
    for (int i = 0; i < n; i++) res += x[i];
    return res;
}

// [[Rcpp::export]]
double sumRcppBis(NumericVector x) {
    double res = std::accumulate(x.begin(), x.end(), 0.0);
    return res;
}

// [[Rcpp::export]]
NumericVector sumRcppBisUpAndDown(NumericVector x) {
    NumericVector res(2);
    double up = 0, down = 0;
    int n = x.size();
    for (int i = 0; i < n; i++) up += x[i];
    for (i = n - 1; i > -1; i--) down += x[i];
    res[0] = up; res[1] = down;
    return res;
}


