#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <queue>
#include <time.h>

using namespace std;

#define ABS(x) (((x)>0)?(x):(-(x)))
#define MAX(x,y) ((x)>(y))?(x):(y)
#define MIN(x,y) ((x)<(y))?(x):(y)
#define PI acos(-1.0)
#define RP(i,n) for(i=0; i<(n); i++)
#define FR(i,a,b) for(i=(a); i<=(b); i++)
#define SQR(x) (x)*(x)
#define EPS 1e-9

#ifdef WIN32
#pragma warning(disable: 4786 4267)
#endif

#ifdef WIN32
typedef __int64 int64;
#else
typedef long long int64;
#endif