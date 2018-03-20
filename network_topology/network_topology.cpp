#include<RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

mat a31(int n) {
	int n2 = pow(n, 2), n0 = 0;
	mat a1 = zeros(n2, n2), a2(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a2(i, j) = n0;
			n0++;
		}
	}
	for (int i = 0; i < n; i++) {
		a1(a2(i, 0), a2(i, n - 1)) = 1;
		a1(a2(i, n - 1), a2(i, 0)) = 1;
		a1(a2(0, i), a2(n - 1, i)) = 1;
		a1(a2(n - 1, i), a2(0, i)) = 1;
		for (int j = 0; j < n - 1; j++) {
			a1(a2(i, j), a2(i, j + 1)) = 1;
			a1(a2(i, j + 1), a2(i, j)) = 1;
			a1(a2(j, i), a2(j + 1, i)) = 1;
			a1(a2(j + 1, i), a2(j, i)) = 1;
		}
	}
	return a1;
}

//[[Rcpp::export]]

mat a32(int n1, int n2)
{
	int n22 = n1*n1*n2, n0 = 0;
	mat a1 = zeros(n22, n22);
	cube a2(n1, n1, n2);
	for (int k = 0; k < n2; k++) {
		for (int i = 0; i < n1; i++) {
			for (int j = 0; j < n1; j++) {
				a2(i, j, k) = n0;
				n0++;
			}
		}
	}
	for (int i = 0; i < n1; i++) {
		for (int k = 0; k < n2; k++) {
			a1(a2(i, 0, k), a2(i, n1 - 1, k)) = 1;
			a1(a2(i, n1 - 1, k), a2(i, 0, k)) = 1;
			a1(a2(0, i, k), a2(n1 - 1, i, k)) = 1;
			a1(a2(n1 - 1, i, k), a2(0, i, k)) = 1;
		}
		for (int j = 0; j < n1; j++) {
			a1(a2(i, j, 0), a2(i, j, n2 - 1)) = 1;
			a1(a2(i, j, n2 - 1), a2(i, j, 0)) = 1;
		}
		for (int k = 0; k < n2; k++) {
			for (int j = 0; j < n1 - 1; j++) {
				a1(a2(i, j, k), a2(i, j + 1, k)) = 1;
				a1(a2(i, j + 1, k), a2(i, j, k)) = 1;
				a1(a2(j, i, k), a2(j + 1, i, k)) = 1;
				a1(a2(j + 1, i, k), a2(j, i, k)) = 1;
			}
		}
		for (int j = 0; j < n1; j++) {
			for (int k = 0; k < n2 - 1; k++) {
				a1(a2(i, j, k), a2(i, j, k + 1)) = 1;
				a1(a2(i, j, k + 1), a2(i, j, k)) = 1;
			}
		}
	}
	return a1;
}

//[[Rcpp::export]]

mat a33(int n) {
	int n2 = pow(2, n), n00 = 1, n0 = 2;
	mat a(n2, n2), temp;
	temp = zeros(1, 1);
	for (int i = 1; i <= n; i++) {
		a(span(n00, n0 - 1), span(n00, n0 - 1)) = temp;
		a(span(0, n00 - 1), span(n00, n0 - 1)) = eye<mat>(n00, n00);
		a(span(n00, n0 - 1), span(0, n00 - 1)) = eye<mat>(n00, n00);
		temp = a(span(0, n0 - 1), span(0, n0 - 1));
		n00 *= 2;
		n0 *= 2;
	}
	return a;
}