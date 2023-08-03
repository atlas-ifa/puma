/* Header file for sortlib */

#ifdef __cplusplus
extern "C" {
#endif

/* Quicksort, possibly carry an integer array: 16 bit int data */
   void tsort_s(int n, short int *x, int *idx);
/* Quicksort, possibly carry an integer array: 16 bit unsigned int data */
   void tsort_u(int n, unsigned short *x, int *idx);
/* Quicksort, possibly carry an integer array: 32 bit integer data */
   void tsort_i(int n, int *x, int *idx);
/* Quicksort, possibly carry an integer array: 64 bit integer data */
   void tsort_I(int n, long long int *x, int *idx);
/* Quicksort, possibly carry an integer array: 32 bit float data */
   void tsort_f(int n, float *x, int *idx);
/* Quicksort, possibly carry an integer array: 64 bit double data */
   void tsort_d(int n, double *x, int *idx);
/* Quicksort, possibly carry an integer array: 128 bit double data */
   void tsort_D(int n, long double *x, int *idx);
/* Quicksort, possibly carry an integer array: C string pointers */
   void tsort_str(int n, char **x, int *idx);


/* Quicksort, possibly carry a pointer array: 16 bit int data */
   void tsort_sp(int n, short int *x, void* *idx);
/* Quicksort, possibly carry a pointer array: 16 bit unsigned int data */
   void tsort_up(int n, unsigned short *x, void* *idx);
/* Quicksort, possibly carry a pointer array: integer data */
   void tsort_ip(int n, int *x, void* *idx);
/* Quicksort, possibly carry a pointer array: 64 bit integer data */
   void tsort_Ip(int n, long long int *x, void* *idx);
/* Quicksort, possibly carry a pointer array: 32 bit float data */
   void tsort_fp(int n, float *x, void* *idx);
/* Quicksort, possibly carry a pointer array: 64 bit double data */
   void tsort_dp(int n, double *x, void* *idx);
/* Quicksort, possibly carry a pointer array: 128 bit double data */
   void tsort_Dp(int n, long double *x, void* *idx);
/* Quicksort, possibly carry a pointer array: 128 bit double data */
   void tsort_strp(int n, char **x, void* *idx);


/* quantile (e.g. median): 16 bit int data */
   short int quantile_s(int n, short int *x, double *wgt, double q);
/* quantile (e.g. median): unsigned short data */
   unsigned short quantile_u(int n, unsigned short *x, double *wgt, double q);
/* quantile (e.g. median): int data */
   int quantile_i(int n, int *x, double *wgt, double q);
/* quantile (e.g. median): 64 bit int data */
   long long int quantile_I(int n, long long int *x, double *wgt, double q);
/* quantile (e.g. median): 32 bit float data */
   float quantile_f(int n, float *x, double *wgt, double q);
/* quantile (e.g. median): 64 bit double data */
   double quantile_d(int n, double *x, double *wgt, double q);
/* quantile (e.g. median): 128 bit long double data */
   long double quantile_D(int n, long double *x, double *wgt, double q);

#if 1
/* insertion sort, possibly carry an integer array, 16 bit int data */
   void insort2_s(int n, short *x, int *idx);
/* insertion sort, possibly carry an integer array, 16 bit unsigned int data */
   void insort2_u(int n, unsigned short *x, int *idx);
/* insertion sort, possibly carry an integer array, int data */
   void insort2_i(int n, int *x, int *idx);
/* insertion sort, possibly carry an integer array, 64 bit int data */
   void insort2_I(int n, long long int *x, int *idx);
/* insertion sort, possibly carry an integer array, 32 bit float data */
   void insort2_f(int n, float *x, int *idx);
/* insertion sort, possibly carry an integer array, 64 bit double data */
   void insort2_d(int n, double *x, int *idx);
/* insertion sort, possibly carry an integer array, 128 bit double data */
   void insort2_D(int n, long double *x, int *idx);
#endif

/* Simple, fast median, short integer data */
short median_s(int n, short *x);
/* Simple, fast median, unsigned short integer data */
unsigned short median_u(int n, unsigned short *x);
/* Simple, fast median, integer data */
int median_i(int n, int *x);
/* Simple, fast median, 64 bit integer data */
int median_I(int n, long long int *x);
/* Simple, fast median, float data */
float median_f(int n, float *x);
/* Simple, fast median, double data */
double median_d(int n, double *x);

/* Simple, fast median, long double data */
long double median_D(int n, long double *x);


// FORTRAN wrappers for sort, quantile, and median
/* Quantiles */
   short int quantile_s_(int *n, short int *x, double *q);
   int quantile_i_(int *n, int *x, double *q);
   float quantile_f_(int *n, float *x, double *q);
   double quantile_d_(int *n, double *x, double *q);
/* weighted quantiles */
   short int wquantile_s_(int *n, short int *x, double *wgt, double *q);
   int wquantile_i_(int *n, int *x, double *wgt, double *q);
   float wquantile_f_(int *n, float *x, double *wgt, double *q);
   double wquantile_d_(int *n, double *x, double *wgt, double *q);

/* sort */
   void tsort_s_(int *n, short int *x);
   void tsort_i_(int *n, int *x);
   void tsort_f_(int *n, float *x);
   void tsort_d_(int *n, double *x);
/* sort and carry an auxiliary array */
   void tsort_aux_s_(int *n, short int *x, int *idx);
   void tsort_aux_i_(int *n, int *x, int *idx);
   void tsort_aux_f_(int *n, float *x, int *idx);
   void tsort_aux_d_(int *n, double *x, int *idx);

/* median */
   short int median_s_(int *n, short int *x);
   int median_i_(int *n, int *x);
   float median_f_(int *n, float *x);
   double median_d_(int *n, double *x);

/* Quicksort, possibly carry an integer array: 128 bit double data */
   void tsort_str_cpp(int n, char **x, int *idx);


#ifdef __cplusplus
}
#endif
