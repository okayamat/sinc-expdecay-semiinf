#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double t)
{
  return sqrt( (1 + (1 - 2*exp(-t))*(1 - 2*exp(-t))) ) * t * exp(-t) / (1 + t);
}

double p(double x)
{
  return asinh(exp(x));
}

double pinv(double t)
{
  return log(sinh(t));
}

double sinc(double x)
{
  double val = 1.0;

  if (x != 0) {
    val = sin(M_PI * x) / (M_PI * x);
  }

  return val;
}

double S(int k, double h, double x)
{
  return sinc((x/h) - k);
}

double fapp(int n, double d, double alpha, double beta, double t)
{
  int k, M, N;
  double mu= fmin(alpha, beta);
  double h = sqrt(M_PI * d / (mu*n));
  if (alpha <= beta) {
    M = n;
    N = (int)ceil(alpha * n / beta);
  } else {
    M = (int)ceil(beta * n / alpha);
    N = n;
  }
  double x = pinv(t);
  double val1 = 0;
  double val2 = 0;

  for (k = -M; k < 0; k++) {
    val1 += f(p(k*h)) * S(k,h,x);
  }
  for (k = N; k >= 0; k--) {
    val2 += f(p(k*h)) * S(k,h,x);
  }

  return val1 + val2;
}

double err_bound(int n, double K, double d, double alpha, double beta)
{
  double mu  = fmin(alpha, beta);
  double pdm = sqrt(M_PI * d * mu);
  double sqn = sqrt(n);
  double val = 1.0;
  val += pow(2, 1+(alpha+beta)*0.5) / (pdm * (1 - exp(-2 * pdm)) * pow(cos(0.5*d),alpha+beta));
  val *= 2 * K / pdm;

  return val * sqn * exp(- pdm * sqn);
}

int main()
{
  int i, n;
  double d = atan(3.0);
  double alpha = 1.0;
  double beta  = 1.0;
  double K = M_SQRT2;

  double t, err, maxerr;

  for (n = 2; n <= 200; n += 5) {
    maxerr = 0;

    for (i = -100; i <= 100; i++) {
      t = pow(2, 0.5*i);
      err = fabs(f(t) - fapp(n, d, alpha, beta, t));

      if (maxerr < err) {
        maxerr = err;
      }
    }

    printf("%d\t%e\t%e\n", n, maxerr, err_bound(n, K, d, alpha, beta));
  }

  return EXIT_SUCCESS;
}
