#include <Python.h>
#include <math.h>
#include <float.h>
#include <numpy/arrayobject.h>


#define pi       3.141592654

long idum;


/*
 *  Routines from the Numerical Recipes
 */
double gammln(float xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}


double factln(int n)
{
  double gammln(float xx);
  static double a[101];
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
  else return gammln(n+1.0);
}

float factrl(int n)
{
  double gammln(float xx);
  static int ntop=4;
  static float a[33]={1.0,1.0,2.0,6.0,24.0};
  int j;

  if (n > 32) return exp(gammln(n+1.0));
  while (ntop<n) {
    j=ntop++;
    a[ntop]=a[j]*ntop;
  }
  return a[n];
}


/* Routine ran2 from Numerical Recipes */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/*Routine to generate Poisson deviates*/
float poidev(float xm, long *idum)
{
  double gammln(float xx);
  float ran2(long *idum);
  static float sq,alxm,g,oldm=(-1.0);
  float em,t,y;
  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em = -1;
    t=1.0;
    do {
      ++em;
      t *= ran2(idum);
    } while (t > g);
  } else {
    if (xm != oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
	y=tan(pi*ran2(idum));
	em=sq*y+xm;
      } while (em < 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (ran2(idum) > t);
  }
  return em;
}

/*
 *  Compute 1sigma error bars of signal via Neyman-Pearson confidence interval construction
 *  by Feldman & Cousins 1998,Phys.Rev.D.,57,3873:
 *  Features: proper treatment of limits independent of data (no flip-flopping), correct coverage
 *  Construction for single Poisson channel: probability P=((mu+b)^n)*exp(-(mu+b))/n!
 *  (1) For all possible counts n compute P(n|mu) for signal mu
 *  (2) For all possible counts n compute P(n|mu') with physically allowed signal mu'=max(0,n-b)
 *  (3) Compute probability ratio R=P(n|mu)/P(n|mu')
 *  (4) Sum P(n|mu) in decreasing order of R until desired confidence level is reached, 
 *      confidence interval for mu is given by minimum and maximum of included n
 *  (5) Repeat (1)-(4) for all possible signals mu
 *  (6) Confidence interval of gross counts is minimum and maximum mu with matching n
 *  (7) Compute error bar from confidence interval
 */
void feldmancousinserror(float grosstot, float backtot, float *errorup, float *errordown)
{
  double factln(int n);
  //float  *prob,*ratio,*mumin,*mumax,*mu;
  float p,rmax,probtest,ptot,mubest;
  long i,j,k,l,m,n;//*rsel,*mutest;
  char text[80];

  /*
   *  These vector lengths define the range of tested counts and mu,
   *  range is limited to reasonable values (speed)
   */

  m = (int) 3*grosstot+100;
  n = (int) (grosstot-backtot+5*sqrt(grosstot)-fmax(grosstot-backtot-5*sqrt(grosstot),0))*100;
  float prob[m];
  double ratio[m];
  long rsel[m];
  long mutest[m];
  double mumin[n];
  double mumax[n];
  double mu[n];
  for (l=1; l<=n; l++) mu[l] = fmax(grosstot-backtot-5*sqrt(grosstot),0)+0.01*(l-1);
  for (l=1; l<=n; l++)
    {
      ptot = 0;
      for (i=1; i<=m; i++)
        {
        rsel[i] = 0;
        prob[i] = exp((i-1)*log(mu[l]+backtot)-(mu[l]+backtot)-factln(i-1));
        ptot += prob[i];
        mubest = (i-1)-backtot;
        if (mubest<0) mubest = 0;
        mutest[i] = i-1;
        probtest = exp((i-1)*log(mubest+backtot)-(mubest+backtot)-factln(i-1));
        ratio[i] = prob[i]/probtest;
        }
      if (ptot<0.99)
        {
        printf("Warning: Not enough coverage of the probability distribution: %f",ptot);
        //SCTPUT(text);
        }
      p = 0;
      mumax[l] = 0;
      mumin[l] = 1e20;
      do
        {
        rmax = 0;
        for (k=1; k<=m; k++)
            {
            if (rsel[k]==0 && ratio[k]>rmax)
            {
            rmax = ratio[k];
            j = k;
            }
            }
        rsel[j] = 1;
        p += prob[j];
        if (mutest[j]>mumax[l]) mumax[l] = mutest[j];
        if (mutest[j]<mumin[l]) mumin[l] = mutest[j];
        }
      while (p<0.6826);
    }
  *errordown = 1e20;
  *errorup = 0;
  for (l=1; l<=n; l++)
    {
      if (mumax[l]>grosstot-1e-6 && mumax[l]<grosstot+1e-6)
        {
        if (mu[l]<*errordown) *errordown=mu[l];
        }
      if (mumin[l]>grosstot-1e-6 && mumin[l]<grosstot+1e-6)
        {
        if (mu[l]>*errorup) *errorup=mu[l];
        }
    }
  if (*errordown>1e19 || *errorup<1e-20)
    {
      *errorup = 1e20;
      *errordown = 0;
      for (l=1; l<=n; l++)
        {
        if (mumax[l]>grosstot-1e-6 && mumax[l]<grosstot+1e-6)
            {
            if (mu[l]>*errordown) *errordown=mu[l];
            }
        if (mumin[l]>grosstot-1e-6 && mumin[l]<grosstot+1e-6)
            {
            if (mu[l]<*errorup) *errorup=mu[l];
            }
        }
    }
  *errordown = grosstot-backtot-*errordown;
  *errorup = *errorup-(grosstot-backtot);
}

/*
 *  Compute bootstrap error of the background in case of negative signal
 *  Feldman-Cousins not defined here due to negative signal, observed is a Poisson fluctuation of the background
 */
void bootstrapbkgerror(float gross, float back, int m, float *errorup, float *errordown)
{
  float poidev(float xm, long *idum);
  float compute_percentile(float *x, int nselect, float percent);
  long j;
  //float *signal;
  //signal = vector(1,m);
  float signal[m];
  for (j=1; j<=m; j++) signal[j] = poidev(back,&idum)-back;
  *errordown = (gross-back)-compute_percentile(signal,m,15.87);
  *errorup = compute_percentile(signal,m,84.13)-(gross-back);
  if (*errordown<0) *errordown = 0;
  if (*errorup<0) *errorup = 0;
  //free_vector(signal,1,m);
}



/*
 *  Function to compute a given percentile of an array
 */
float compute_percentile(float *x, int nselect, float percent)
{
  float selectval(unsigned long k, unsigned long n, float arr[]);
  float percentile;
  unsigned long nom,denom;

  denom = 0;
  do denom++; while (fabs(percent/100*denom-floor(percent/100*denom))>1e-6);
  nom = floor(percent/100*denom); 
  if (div(nselect*nom,denom).rem!=0)
    {
      percentile = selectval(div(nselect*nom,denom).quot+1,nselect,x);
    }
  else
    {
      percentile = (selectval(div(nselect*nom,denom).quot,nselect,x)+selectval(div(nselect*nom,denom).quot+1,nselect,x))/2;
    }
  return percentile;
}


/*
 *  find the kth smallest value in array arr[1..n], input array rearranged to
 *  have this value in location arr[k] with smaller elements moved to arr[1..k-1]
 *  and larger elements moved to arr[k+1..n] in arbitrary order
 *
 *  Routine from Numerical Recipes p. 342
 */
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
float selectval(unsigned long k, unsigned long n, float arr[])
{
  unsigned long i,ir,j,l,mid;
  float a,temp;
  
  l = 1;
  ir = n;
  for (;;)
    {
      if (ir<=l+1)
	{
	  if (ir==l+1 && arr[ir]<arr[l])
	    {
	      SWAP(arr[l],arr[ir]);
	    }
	  return arr[k];
	}
      else
	{
	  mid=(l+ir)>>1;
	  SWAP(arr[mid],arr[l+1])
	  if (arr[l]>arr[ir])
	    {
	      SWAP(arr[l],arr[ir])
	    }
	  if (arr[l+1]>arr[ir])
	    {
	      SWAP(arr[l+1],arr[ir])
	    }
	  if (arr[l]>arr[l+1])
	    {
	      SWAP(arr[l],arr[l+1])
	    }
	  i = l+1;
	  j = ir;
	  a = arr[l+1];
	  for (;;)
	    {
	      do i++; while (arr[i]<a);
	      do j--; while (arr[j]>a);
	      if (j<i) break;
	      SWAP(arr[i],arr[j])
	    }
	  arr[l+1]=arr[j];
	  arr[j]=a;
	  if (j>=k) ir=j-1;
	  if (j<=k) l=i;
	}
    }
}


// confidence limits according to Feldman & Cousins 1998
static PyObject* feldman_cousins_conf_lim(PyObject* self, PyObject* args)
{
    float N, B, err_up, err_down;
    if(!PyArg_ParseTuple(args, "ff", &N, &B))
        return NULL;
    feldmancousinserror(N, B, &err_up, &err_down);
    return Py_BuildValue("(ff)", err_down, err_up);
}


static PyObject* bootstrap_bkg_conf_lim(PyObject* self, PyObject* args)
{
    float gross,back,errorup,errordown;
    int m;
    if(!PyArg_ParseTuple(args, "ffi", &gross, &back, &m))
        return NULL;
    bootstrapbkgerror(gross, back, m, &errorup, &errordown);
    return Py_BuildValue("(ff)", errordown, errorup);
}



// Our Module's Function Definition struct
// We require this `NULL` to signal the end of our method
// definition
static PyMethodDef myMethods[] = {
    { "feldman_cousins_conf_lim", feldman_cousins_conf_lim, 
        METH_VARARGS, "Calculates conf limits using Feldman & Cousins 1998" },
     { "bootstrap_bkg_conf_lim", bootstrap_bkg_conf_lim, 
        METH_VARARGS, "Bootstrap" },
    { NULL, NULL, 0, NULL }
    
};

// Our Module Definition struct
static struct PyModuleDef CustomConfLim = {
    PyModuleDef_HEAD_INIT,
    "CustomConfLim",
    "Test Module",
    -1,
    myMethods
};

// Initializes our module using our above struct
PyMODINIT_FUNC PyInit_CustomConfLim(void)
{
    import_array();
    return PyModule_Create(&CustomConfLim);
}
