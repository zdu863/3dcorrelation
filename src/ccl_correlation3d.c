#include "gsl/gsl_integration.h"
#include "ccl_pk.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_spline.h"
#include "ccl_error.h"
#include "ccl_utils.h"
#include "ccl_correlation.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "ccl_power.h"
#include "ccl.h"
#include "fftlog.h"


/*--------ROUTINE: taper_pk ------
TASK:n Apply cosine tapering to Cls to reduce aliasing
INPUT: number of k bins for pk, k vector, pk vector, limits for tapering
       e.g., k_limits=[low_k_limit_lower,low_k_limit_upper,high_k_limit_lower,high_k_limit_upper]
*/
static int taper_pk(int n_k,double *k,double *pk, double *k_limits)
{

  for(int i=0;i<n_k;i++) {
    if(k[i]<k_limits[0] || k[i]>k_limits[3]) {
      k[i]=0;//k outside desirable range
      continue;
    }
    if(k[i]>=k_limits[1] && k[i]<=k_limits[2])
      continue;//k within good k range
    
    if(k[i]<k_limits[1])//tapering low k
      pk[i]*=cos((k[i]-k_limits[1])/(k_limits[1]-k_limits[0])*M_PI/2.);
    
    if(k[i]>k_limits[2])//tapering high k
      pk[i]*=cos((k[i]-k_limits[2])/(k_limits[3]-k_limits[2])*M_PI/2.);
  }

  return 0;
}

/*--------ROUTINE: ccl_tracer_corr_fftlog ------
TASK: For a given tracer, get the correlation function
      Following function takes a function to calculate angular pk as well.
      By default above function will call it using ccl_angular_pk
INPUT: type of tracer, number of r values to evaluate = NL, r vector
 */
#define K_MIN 0.0001
#define K_MAX 500
#define N_K 5000

typedef struct {
  int nell;
  double ell0;
  double ellf;
  double cl0;
  double clf;
  int extrapol_0;
  int extrapol_f;
  double tilt0;
  double tiltf;
  SplPar *pk_spl;
  int i_bessel;
  double th;
} corr_int_par;

void ccl_correlation3d(ccl_cosmology *cosmo,
		     int n_k,double *k,double *pk,
		     int n_r,double *r,double *xi,
		     int do_taper_pk,double *taper_pk_limits,
		     int *status)
{
  int i;
  double *k_arr,*pk_arr,*r_arr,*xi_arr;

  k_arr=ccl_log_spacing(K_MIN,K_MAX,N_K);
  if(k_arr==NULL) {
    *status=CCL_ERROR_LINSPACE;
    strcpy(cosmo->status_message,"ccl_3dcorrelation.c: ran out of memory\n");
    return;
  }
  pk_arr=malloc(N_K*sizeof(double));
  if(pk_arr==NULL) {
    free(k_arr);
    *status=CCL_ERROR_MEMORY;
    strcpy(cosmo->status_message,"ccl_3dcorrelation.c: ran out of memory\n");
    return;
  }

  //Interpolate input pk into array needed for FFTLog
  SplPar *pk_spl=ccl_spline_init(n_k,k,pk,pk[0],0);
  if(pk_spl==NULL) {
    free(k_arr);
    free(pk_arr);
    *status=CCL_ERROR_MEMORY;
    strcpy(cosmo->status_message,"ccl_3dcorrelation.c: ran out of memory\n");
    return;
  }

  double pk_tilt,k_edge,pk_edge;
  k_edge=k[n_k-1];
  if((pk[n_k-1]*pk[n_k-2]<0) || (pk[n_k-2]==0)) {
    pk_tilt=0;
    pk_edge=0;
  }
  else {
    pk_tilt=log(pk[n_k-1]/pk[n_k-2])/log(k[n_k-1]/k[n_k-2]);
    pk_edge=pk[n_k-1];
  }
  for(i=0;i<N_K;i++) {
    if(k_arr[i]>=k_edge)
      pk_arr[i]=pk_edge*pow(k_arr[i]/k_edge,pk_tilt);
    else
      pk_arr[i]=ccl_spline_eval(k_arr[i],pk_spl);
  }
  ccl_spline_free(pk_spl);

  if (do_taper_pk)
    taper_pk(N_K,k_arr,pk_arr,taper_pk_limits);

  r_arr=malloc(sizeof(double)*N_K);
  if(r_arr==NULL) {
    free(k_arr);
    free(pk_arr);
    *status=CCL_ERROR_MEMORY;
    strcpy(cosmo->status_message,"ccl_correlation.c: ccl_tracer_corr_fftlog ran out of memory\n");
    return;
  }
  xi_arr=malloc(sizeof(double)*N_K);
  if(xi_arr==NULL) {
    free(k_arr); free(pk_arr); free(r_arr);
    *status=CCL_ERROR_MEMORY;
    strcpy(cosmo->status_message,"ccl_correlation.c: ccl_tracer_corr_fftlog ran out of memory\n");
    return;
  }

  for(i=0;i<N_K;i++)
    r_arr[i]=0;
 
  pk2xi(N_K,k_arr,pk_arr,r_arr,xi_arr);

  // Interpolate to output values of r
  SplPar *xi_spl=ccl_spline_init(N_K,r_arr,xi_arr,xi_arr[0],0);
  for(i=0;i<n_r;i++)
    xi[i]=ccl_spline_eval(r[i],xi_spl);
  ccl_spline_free(xi_spl);

  free(k_arr); free(pk_arr);
  free(r_arr); free(xi_arr);

  ccl_check_status(cosmo,status);

  return;
}
