void covariance ( int m, int n, int x[], double *average, double *std, double *covc );
double d_epsilon ( void );
double d_huge ( void );
int d_nint ( double x );
double d_uniform_01 ( int *seed );
char digit_to_ch ( int i );
double dvec_average ( int n, double a[] );
double dvec_std ( int n, double a[] );
int get_seed ( void );
int i_log_10 ( int i );
int i_max ( int i1, int i2 );
int i_min ( int i1, int i2 );
char *i_to_s ( int i );
int i_uniform ( int b, int c, int *seed );
void ihs ( int m, int n, int d, int *seed, int x[] );
void ihs_write ( int m, int n, int d, int seed_init, int seed, 
  int r[], char *file_out_name );
void timestamp ( void );
char *timestring ( void );
