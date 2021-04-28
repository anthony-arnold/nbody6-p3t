#include <algorithm>

namespace KNN {
   double norm(double xi[3], double xj[3]) {
      double d = 0;
      for (int k = 0; k < 3; k++) {
	 double r = xi[k] - xj[k];
	 d += r*r;
      }
      return d;
   }

   void knn(int i, int n, double x[][3], int k, double* rs) {
      double* xi = x[i];
      double* ri = rs + k*i;
      int kk = 0;

      for (int j = 0; j < n; j++) {
	 if (j == i) continue;
	 double r = norm(xi, x[j]);
	 double* q = ri + kk;
	 double* p = std::lower_bound(ri, q, r);
	 if (p < ri+k) {
	    kk = std::min(k, kk+1);
	    std::copy_backward(p, q, ri+kk);
	    *p = r;
	 }
      }
   }
}

extern "C" void knn_(int* l, int* n, double x[][3], int* k, double* rs) {
#pragma omp parallel for
   for (int i = 0; i < *l; i++) {
      KNN::knn(i, *n, x, *k, rs);
   }
}
