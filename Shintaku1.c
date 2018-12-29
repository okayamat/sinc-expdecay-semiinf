#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double t)
{
  return pow(t, M_PI_4)*exp(-t);
}

double p(double x)
{
  double val;

  if (x >= 0) {
    val = x + log1p(exp(- x));
  } else {
    val = log1p(exp(x));
  }

  return val;
}

double pinv(double t)
{
  return log(expm1(t));
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
  val += 2 * pow(M_E / (M_E - 1), 0.5*mu) / (pdm * (1 - exp(-2 * pdm)) * pow(cos(0.5*d),alpha+beta));
  val *= 2 * K / pdm;

  return val * sqn * exp(- pdm * sqn);
}

int main()
{
  int i, n;
  double d = 3;
  double alpha = M_PI_4;
  double beta  = 1.0 - 0.5*alpha*M_1_PI;
  double gamma = - log(cos(0.5*d));
  double K = pow(((1 - gamma)*(1 - gamma) + M_PI*M_PI)*exp(gamma*M_1_PI),0.5*alpha);

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
