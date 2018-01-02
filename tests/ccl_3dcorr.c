#include "ccl.h"
#include "ccl_correlation.h"
#include "ccl_utils.h"
#include "ctest.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

// Define cosmological parameter values
#define OC 0.25
#define OB 0.05
#define OK 0.00
#define ON 0.00
#define HH 0.70
#define W0 -1.0
#define WA 0.00
#define NS 0.96
#define NORMPS 0.80
#define ZD 0.5
#define NZ 128
#define Z0_GC 0.50
#define SZ_GC 0.05
#define NL 512
#define PS 0.1 
#define NREL 3.046
#define NMAS 0
#define MNU 0.0
#define K_MIN 0.001
#define K_MAX 100
#define N_MAX_PK 20000


int main(int argc,char **argv)
{
  // Use the default configuration, plus the cosmological parameters that were 
  // defined above
  int status=0;
  ccl_configuration config = default_config;
  ccl_parameters params = ccl_parameters_create(OC, OB, OK, NREL, NMAS, MNU, 
                                                W0, WA, HH, NORMPS, NS,
						14.079181246047625, 0.5, 55, 0, NULL, NULL, &status);
  ccl_cosmology *cosmo = ccl_cosmology_create(params,config);

  int i;
  double *pkarr = malloc(N_MAX_PK*sizeof(double));
  double *karr;

  karr=ccl_log_spacing(K_MIN,K_MAX,N_MAX_PK);
  for(i=0; i < N_MAX_PK; i++){
    // Calculate power spectrum
    pkarr[i] = ccl_nonlin_matter_power(cosmo, karr[i], 1.0, &status); 
  }
  
  // Define cosine tapering, to reduce ringing. The first two numbers are 
  // [lmin, lmax] for the low-ell taper, and the last two are [lmin, lmax] for 
  // the high-ell taper.
  double taper_pk_limits[4] = {0.001, 0.002, 50, 75};
  
  double *xi, *r;
  int n_r = 20;
  r = ccl_log_spacing(0.01, 500., n_r); // New array with log spacing
  xi = malloc(n_r*sizeof(double));

  // Calculate 3dcorrelation function from power spectrum
  ccl_correlation3d(cosmo, N_MAX_PK, karr, pkarr, 
                  n_r, r, xi,
                  1, taper_pk_limits, 
                  &status);
  
  // Print results
  for(int it=0; it < n_r; it++)
    printf("%le %le\n", r[it], xi[it]);
  
  // Free allocated memory
  ccl_cosmology_free(cosmo);
  free(xi);
  free(karr);
  free(pkarr);

  return 0;
}
