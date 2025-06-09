#include <R.h>

void sum_C(double *x, unsigned int *n, double *sum)
{
    *sum = 0; /* in case not initialized in R */
    for (int i = 0; i < *n; i++) *sum += x[i];
}
