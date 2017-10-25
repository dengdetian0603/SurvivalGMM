#include <Rcpp.h>
#include <iostream>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
int binSearch(NumericVector& A, double x) {
  // return i s.t. A[i] is the largest number that < x
	int start = 0;
	int end = A.size() - 1;
	int mid;
	if (x < A[0]) {
		return -1;
	}
	while (start + 1 < end) {
		mid = start + (end - start)/2;
		if (x < A[mid]) {
			end = mid;
		} else {
			start = mid;
		}
	}
	if (A[end] < x) {
		return end;
	} else {
		return start;
	}
}

// [[Rcpp::export]]
int binSearch2(NumericVector& A, double x) {
  // return i s.t. A[i] is the smallest number that >= x
  int start = 0;
  int end = A.size() - 1;
  int mid;
  if (x > A[end]) {
    return end + 1;
  }
  while (start + 1 < end) {
    mid = start + (end - start)/2;
    if (x < A[mid]) {
      end = mid;
    } else {
      start = mid;
    }
  }
  if (A[start] >= x) {
    return start;
  } else {
    return end;
  }
}

// [[Rcpp::export]]
NumericMatrix hazardMat(NumericVector& cumSumVec,
                        NumericVector& eVec,
                        NumericMatrix& epsilonMat) {
	NumericMatrix Alpha(epsilonMat.nrow(), epsilonMat.ncol());
  int index;
  for (int i = 0; i < Alpha.nrow(); ++i) {
    for (int t = 0; t < Alpha.ncol(); ++t) {
      index = binSearch(eVec, epsilonMat(i, t));
// cout << "i: " << i << " t: " << t << " index: " << index << endl;
      Alpha(i, t) = (index < 0) ? 0 : cumSumVec[index];
    }
  }
  return Alpha;
}




// [[Rcpp::export]]
NumericMatrix psiMat(NumericVector& cumSumVec2,
                     NumericVector& eVec,
                     NumericVector& dVec,
                     NumericMatrix& epsilonMat,
                     NumericMatrix& alphaMat,
                     IntegerVector& grpID,
                     NumericMatrix& survProbs) {
  int N = epsilonMat.nrow();
  int nGrp = survProbs.nrow();
  int nT = survProbs.ncol();
  NumericMatrix psiMoments(N, nGrp * nT);
  NumericMatrix vMat(N, N);
  NumericMatrix hMat(N, N);
  int index, psiCol;
  double Ehi, Esumhij;

  for (int t = 0; t < nT; ++t) { 
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        index = binSearch(eVec, min(epsilonMat(i, t), eVec[j]));
        vMat(i, j) = alphaMat(i, t);
        vMat(i, j) -= (index < 0) ? 0: cumSumVec2[index];
        vMat(i, j) += (eVec[j] > epsilonMat(i, t)) ? 0 : (N * dVec[j]/(N - j)); 
      }
    }
    for (int k = 1; k <= nGrp; ++k) {
      psiCol = t * nGrp + k - 1;
      Esumhij = 0.0;
      for (int i = 0; i < N; ++i) {
        psiMoments(i, psiCol) = 0;
        if (grpID[i] == k) {
          psiMoments(i, psiCol) += 
            (1 + alphaMat(i, t)) * exp(-alphaMat(i, t)) -
            survProbs(k - 1, t) - 1/N * exp(-alphaMat(i, t)) * vMat(i, i);
        }
        Ehi = 0.0;
        for (int j = 0; j < N; ++j) {       
          hMat(i, j) = 
            ((grpID[i] == k) ? (exp(-alphaMat(i, t)) * vMat(i, j)) : 0) +
            ((grpID[j] == k) ? (exp(-alphaMat(j, t)) * vMat(j, i)) : 0);
          Ehi += hMat(i, j);
          if (j > i) {
            Esumhij += 
              ((grpID[i] == k) ? (exp(-alphaMat(i, t)) * alphaMat(i, t)) : 0) +
              ((grpID[j] == k) ? (exp(-alphaMat(j, t)) * alphaMat(j, t)) : 0);
          }
        }
        Ehi = Ehi/N;
        psiMoments(i, psiCol) -= Ehi;
      }
      Esumhij = Esumhij/(N * N - N);
      for (int i = 0; i < N; ++i) {
        psiMoments(i, psiCol) += Esumhij;
      }
    }
  }
  return psiMoments;
}

