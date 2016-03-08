#define SIGN(a,b) ( (b) >= 0.0 ? fabs(a) : -fabs(a) )
#define NR_END 1
#define FREE_ARG char*


void nrerror(char error_text[])
{
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}


double *vector(long nl, long nh)
     /* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) 
    nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}


void free_vector(double *v, long nl, long nh)
     /* free a double vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
		   (lminarg1) : (lminarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
		   (lmaxarg1) : (lmaxarg2))

static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))


static double dminarg1, dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))

#define SQR(x) ( (x) * (x) )




void rk4(double y[], double dydx[], long n, double x, double h,
	 double yout[], void (*derivs)(double, double [], double []))
{
  long i;
  double xh, hh, h6, *dym, *dyt, *yt;

  dym = (double *) malloc( (size_t)(n+1)*sizeof(double) );
  dyt = (double *) malloc( (size_t)(n+1)*sizeof(double) );
  yt = (double *) malloc( (size_t)(n+1)*sizeof(double) );
  hh = h*0.5;
  h6 = h/6.0;
  xh = x+hh;
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
  (*derivs)(xh,yt,dyt);
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
  (*derivs)(xh,yt,dym);
  for (i=1;i<=n;i++)
    {
      yt[i] = y[i]+h*dym[i];
      dym[i] += dyt[i];
    }
  (*derivs)(x+h,yt,dyt);
  for (i=1;i<=n;i++)
    yout[i] = y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  free(dym);
  free(dyt);
  free(yt);
}



double *xp, **yp;


void rkdumb(double vstart[], long nvar, double x1, double x2, long nstep,
	    void (*derivs)(double, double [], double []),
	    void (*print_result)(double, double []))
{
  void rk4(double y[], double dydx[], long n, double x, double h,
	   double yout[], void (*derivs)(double, double [], double []));
  long i, k;
  double x, h;
  double *v, *vout, *dv;
  
  v = (double *) malloc( (size_t)(nvar+1)*sizeof(double) );
  vout = (double *) malloc( (size_t)(nvar+1)*sizeof(double) );
  dv = (double *) malloc( (size_t)(nvar+1)*sizeof(double) );

  for (i=1;i<=nvar;i++)
    v[i] = vstart[i];
  /*
    for (i=1; i<=nvar; i++)
    {
    v[i] = vstart[i];
    yp[i][1] = v[i];
    }
    xp[1] = x1;
  */

  x = x1;
  h = (x2-x1)/nstep;
  for (k=1;k<nstep;k++)
    {
      (*derivs)(x,v,dv);
      rk4(v,dv,nvar,x,h,vout,derivs);
      if ((double)(x+h)==x) nrerror("Step size too small in routine rkdumb");
      x += h;
      for (i=1;i<=nvar;i++)
	v[i] = vout[i];

      (*print_result)(x, v);

      /*
	xp[k+1] = x;
	for (i=1;i<=nvar;i++)
	{
	v[i] = vout[i];
	yp[i][k+1] = v[i];
	}
	*/
    }
  free(dv);
  free(vout);
  free(v);
}  








#define KMAXX 8
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

double **d, *x;

void bsstep( double y[], double dydx[], long nv, double *xx, double htry,
	     double eps, double yscal[], double *hdid, double *hnext,
	     void (*derivs)(double, double [], double []) )
{
  void mmid(double y[], double dydx[], long nvar, double xs, double htot,
	    long nstep, double yout[], 
	    void (*derivs)(double, double [], double []));
  void pzextr(long iest, double xest, double yest[], double yz[], double dy[],
	      long nv);
  void pzextr1(long iest, double xest, double yest[], double yz[], double dy[],
	      long nv);

  long i, iq, k, kk, km;
  static long first=1, kmax, kopt;
  static double epsold=-1.0, xnew;
  double eps1, errmax, fact, h, red, scale, work, wrkmin, xest;
  double *err, *yerr, *ysav, *yseq;
  static double a[IMAXX+1];
  static double alf[KMAXX+1][KMAXX+1];
  static long nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
  long reduct, exitflag=0;

  d = (double **) malloc( (size_t)(nv+1)*sizeof(double*));
  d[0] = (double *) malloc( (size_t)((KMAXX+1)*(nv+1))*sizeof(double) );
  for (i=1; i<=nv; i++)
    d[i] = d[i-1]+KMAXX+1;
  err = (double *) malloc( (size_t)(KMAXX+1)*sizeof(double) );
  x = (double *) malloc( (size_t)(KMAXX+1)*sizeof(double) );
  yerr = (double *) malloc( (size_t)(nv+1)*sizeof(double) );
  ysav = (double *) malloc( (size_t)(nv+1)*sizeof(double) );
  yseq = (double *) malloc( (size_t)(nv+1)*sizeof(double) );

  if (eps!=epsold)
    {
      *hnext = xnew = -1.0e29;
      eps1 = SAFE1*eps;
      a[1] = nseq[1]+1;
      for (k=1; k<=KMAXX; k++) a[k+1]=a[k]+nseq[k+1];
      for (iq=2; iq<=KMAXX; iq++)
	{
	  for (k=1; k<iq; k++)
	    alf[k][iq] = pow(eps1, (a[k+1]-a[iq+1])/
			     ((a[iq+1]-a[1]+1.0)*(2*k+1)));
	}
      epsold = eps;
      for (kopt=2; kopt<KMAXX; kopt++)
	if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
      kmax = kopt;
    }
  h = htry;
  for (i=1; i<=nv; i++) ysav[i] = y[i];
  if (*xx != xnew || h!=(*hnext))
    {
      first = 1;
      kopt = kmax;
    }
  reduct = 0;
  
  for(;;)
    {
      for (k=1; k<=kmax; k++)
	{
	  xnew = (*xx)+h;
	  if (xnew== (*xx)) nrerror("step size underflow in bsstep");
	  mmid(ysav, dydx, nv, *xx, h, nseq[k], yseq, derivs);
	  xest=SQR(h/nseq[k]);
	  /*	  pzextr(k, xest, yseq, y, yerr, nv); */
	  pzextr1(k, xest, yseq, y, yerr, nv);
	  if (k!=1)
	    {
	      errmax = TINY;
	      for (i=1; i<=nv; i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
	      errmax /= eps;
	      km = k-1;
	      err[km] = pow(errmax/SAFE1, 1.0/(2*km+1));
	    }
	  if (k!=1 && (k>=kopt-1 || first))
	    {
	      if (errmax<1.0)
		{
		  exitflag = 1;
		  break;
		}
	      if (k==kmax || k==kopt+1)
		{
		  red = SAFE2/err[km];
		  break;
		}
	      else if (k==kopt && alf[kopt-1][kopt]<err[km])
		{
		  red = 1.0/err[km];
		  break;
		}
	      else if (kopt==kmax && alf[km][kmax-1]<err[km])
		{
		  red = alf[km][kmax-1]*SAFE2/err[km];
		  break;
		}
	      else if (alf[km][kopt] < err[km])
		{
		  red = alf[km][kopt-1]/err[km];
		  break;
		}
	    }
	}
      if (exitflag) break;
      red=DMIN(red,REDMIN);
      red=DMAX(red,REDMAX);
      h *= red;
      reduct = 1;
    }
  *xx = xnew;
  *hdid = h;
  first = 0;
  wrkmin = 1.0e35;
  for (kk=1; kk<=km; kk++)
    {
      fact=DMAX(err[kk],SCALMX);
      work=fact*a[kk+1];
      if (work<wrkmin)
	{
	  scale=fact;
	  wrkmin=work;
	  kopt=kk+1;
	}
    }
  *hnext = h/scale;
  if (kopt>=k && kopt!=kmax && !reduct)
    {
      fact = DMAX(scale/alf[kopt-1][kopt], SCALMX);
      if (a[kopt+1]*fact <= wrkmin)
	{
	  *hnext=h/fact;
	  kopt++;
	}
    }
  free(yseq);
  free(ysav);
  free(yerr);
  free(x);
  free(err);
  free(d[0]);
  free(d);
}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX





void pzextr(long iest, double xest, double yest[], double yz[], double dy[],
	    long nv)
{
  long k1, j;
  double q, f2, f1, delta, *c;

  c = (double *) malloc( (size_t)(nv+1)*sizeof(double) );
  x[iest]=xest;
  for (j=1; j<=nv; j++) dy[j]=yz[j]=yest[j];
  if (iest==1)
    for (j=1; j<=nv; j++) d[j][1]=yest[j];
  else
    {
      for (j=1; j<=nv; j++) c[j]=yest[j];
      for (k1=1; k1<iest; k1++)
	{
	  delta=1.0/(x[iest-k1]-xest);
	  f1=xest*delta;
	  f2=x[iest-k1]*delta;
	  for (j=1; j<=nv; j++)
	    {
	      q=d[j][k1];
	      d[j][k1]=dy[j];
	      delta=c[j]-q;
	      dy[j]=f1*delta;
	      c[j]=f2*delta;
	      yz[j]+=dy[j];
	    }
	}
      for (j=1; j<=nv; j++) d[j][iest]=dy[j];
    }
  free(c);
}




void pzextr1(long iest, double xest, double yest[], double yz[], double dy[],
	    long nv)
{
  long k, j;
  double yy, v, ddy, c, b1, b, *fx;


  fx = (double *) malloc( (size_t)(iest+1)*sizeof(double) );
  x[iest]=xest;
  if (iest==1)
    for (j=1; j<=nv;j++)
      {
	yz[j]=yest[j];
	d[j][1]=yest[j];
	dy[j]=yest[j];
      }
  else
    {
      for (k=1; k<iest; k++)
	fx[k+1]=x[iest-k]/xest;
      for (j=1;j<=nv;j++)
	{
	  v=d[j][1];
	  d[j][1]=yy=c=yest[j];
	  for (k=2;k<=iest;k++)
	    {
	      b1=fx[k]*v;
	      b=b1-c;
	      if (b)
		{
		  b=(c-v)/b;
		  ddy=c*b;
		  c=b1*b;
		}
	      else
		ddy=v;
	      if (k!=iest) v=d[j][k];
	      d[j][k]=ddy;
	      yy+=ddy;
	    }
	  dy[j]=ddy;
	  yz[j]=yy;
	}
    }
  free(fx);
}
      






void mmid(double y[], double dydx[], long nvar, double xs, double htot,
	  long nstep, double yout[], 
	  void (*derivs)(double, double [], double []))
{
  long n, i;
  double x, swap, h2, h, *ym, *yn;
  
  ym = (double *) malloc( (size_t)(nvar+1)*sizeof(double) );
  yn = (double *) malloc( (size_t)(nvar+1)*sizeof(double) );
  h=htot/nstep;
  for (i=1;i<=nvar;i++)
    {
      ym[i] = y[i];
      yn[i] = y[i]+h*dydx[i];
    }
  x=xs+h;
  (*derivs)(x,yn,yout);
  h2=2.0*h;
  for (n=2; n<=nstep; n++)
    {
      for (i=1; i<=nvar; i++)
	{
	  swap=ym[i]+h2*yout[i];
	  ym[i]=yn[i];
	  yn[i]=swap;
	}
      x+=h;
      (*derivs)(x,yn,yout);
    }
  for (i=1; i<=nvar; i++)
    yout[i] = 0.5*(ym[i]+yn[i]+h*yout[i]);
  free(yn);
  free(ym);
}
      









#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqs(double y[], double dydx[], long n, double *x, double htry, 
	  double eps, double yscal[], double *hdid, double *hnext, 
	     void (*derivs)(double, double [], double []))
{
  void rkck(double y[], double dydx[], long n, double x, double h,
	    double yout[], double yerr[], 
	    void (*derivs)(double, double [], double []));
  long i;
  double errmax,h,htemp,xnew,*yerr,*ytemp;

  yerr = (double *) malloc( (size_t)(n+2)*sizeof(double) );
  ytemp = (double *) malloc( (size_t)(n+2)*sizeof(double) );
  h = htry;
  for(;;)
    {
      rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
      errmax = 0.0;
      for (i=1;i<=n;i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
      errmax /= eps;
      if (errmax <= 1.0) break;
      htemp = SAFETY*h*pow(errmax, PSHRNK);
      h = (h >= 0.0 ? DMAX(htemp,0.1*h) : DMIN(htemp, 0.1*h));
      xnew = (*x)+h;
      if (xnew == *x) nrerror("Stepsize underflow in rkqs");
    }
  if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax, PGROW);
  else *hnext=5.0*h;
  *x += (*hdid=h);
  for (i=1; i<=n; i++) y[i]=ytemp[i];
  free(ytemp);
  free(yerr);
}


#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON




void rkck(double y[], double dydx[], long n, double x, double h,
	  double yout[], double yerr[], 
	  void (*derivs)(double, double [], double []))
{
  long i;
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,b31=3.0/40.0
    ,b32=9.0/40.0,b41=0.3,b42=-0.9,b43=1.2,b51=-11.0/54.0,b52=2.5,
    b53=-70.0/27.0,b54=35.0/27.0,b61=1631.0/55296.0,b62=175.0/512.0,
    b63=575.0/13824.0,b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,dc5=-277.0/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,dc4=c4-13525.0/55296.0
    ,dc6=c6-0.25;
  double *ak2, *ak3, *ak4, *ak5, *ak6, *ytemp;

  ak2 = (double *) malloc( (size_t)(n+2)*sizeof(double) );
  ak3 = (double *) malloc( (size_t)(n+2)*sizeof(double) );
  ak4 = (double *) malloc( (size_t)(n+2)*sizeof(double) );
  ak5 = (double *) malloc( (size_t)(n+2)*sizeof(double) );
  ak6 = (double *) malloc( (size_t)(n+2)*sizeof(double) );
  ytemp = (double *) malloc( (size_t)(n+2)*sizeof(double) );

  for (i=1;i<=n;i++)                 /* First step  */
    ytemp[i] = y[i]+b21*h*dydx[i];
  (*derivs)(x+a2*h,ytemp,ak2);       /* Second step  */
  for (i=1;i<=n;i++)                 
    ytemp[i] = y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3);       /* Third step  */
  for (i=1;i<=n;i++)          
    ytemp[i] = y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4);       /* Fourth step  */
  for (i=1;i<=n;i++)          
    ytemp[i] = y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5);       /* Fifth step  */
  for (i=1;i<=n;i++)          
    ytemp[i] = y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]
		       +b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6);       /* Sixth step  */
  for (i=1;i<=n;i++)          
    yout[i] = y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=1;i<=n;i++)          
    yerr[i] = h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

  free(ytemp);
  free(ak6);
  free(ak5);
  free(ak4);
  free(ak3);
  free(ak2);
}


#define MAXSTP 50000000
#define TINY 1.0e-30;

long kmax, kout=0;
double dxsav;

void odeint( double ystart[], long nvar, double x1, double x2, double eps,
	     double h1, double hmin, long *nok, long *nbad, 
	     void (*derivs)(double, double [], double []),
	     void (*rkqs)(double [], double [], long, double *, double, 
			  double, double [], double *, double *, 
			  void (*)(double, double [], double [])),
	     void (*print_result)(double, double []) )
{
  long nstp, i;
  double xsav, x, hnext, hdid, h;
  double *yscal, *y, *dydx;

  yscal = (double *) malloc( (size_t)(nvar+2)*sizeof(double) );
  y = (double *) malloc( (size_t)(nvar+2)*sizeof(double) );
  dydx = (double *) malloc( (size_t)(nvar+2)*sizeof(double) );
  x = x1;
  h = SIGN(h1, x2-x1);
  *nok = (*nbad) = kout = 0;
  for (i=1; i<=nvar; i++) y[i]=ystart[i];
  if (kmax > 0) xsav=x-dxsav*2.0;
  for (nstp=1;nstp<=MAXSTP;nstp++)
    {
      (*derivs)(x,y,dydx);
      for (i=1;i<=nvar;i++)
	yscal[i] = fabs(y[i])+fabs(dydx[i]*h)+TINY;
	/* yscal[i] = fabs(eps*dydx[i]*h)+TINY;  */

      if (kmax > 0 && kout < kmax-1 && fabs(x-xsav) > fabs(dxsav))
	{
	  /*
	    xp[++kout] = x;
	    for (i=1; i<=nvar; i++) 
	    {
	    Its possible to print here the results 
	    yp[i][kout]=y[i];
	    }
	    */
	  (*print_result)(x, y);
	  xsav = x;
	}
      if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
      (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
      if (hdid == h) ++(*nok); else ++(*nbad);
      if ((x-x2)*(x2-x1) >= 0.0)
	{
	  for (i=1;i<=nvar;i++) ystart[i]=y[i];
	  if (kmax)
	    {
	      (*print_result)(x, y);	      
	      /*
		xp[++kout]=x;
		for (i=1;i<=nvar;i++) yp[i][kout]=y[i];
		Its possible to print here the results 
		*/
	    }
	  free(dydx);
	  free(y);
	  free(yscal);
	  return;
	}
      if (fabs(hnext) <= hmin) 
	{
	  printf("\n\nStep size too small in odeint\n\n");
	  /*	  nrerror("Step size too small in odeint");  */
	  free(dydx);
	  free(y);
	  free(yscal);
	  return;
	}
      h=hnext;
    }

  printf( "\n\nToo many steps in routine odeint\n\n" );
/*
  nrerror("Too many steps in routine odeint");
*/
}



#undef MAXSTP
#undef TINY

#define MAXSTP 50000000
#define TINY 1.0e-30;

void my_odeint( double ystart[], long nvar, double x1, double x2, double eps,
	     double h1, double hmin, long *nok, long *nbad, 
	     void (*derivs)(double, double [], double []),
	     void (*rkqs)(double [], double [], long, double *, double, 
			  double, double [], double *, double *, 
			  void (*)(double, double [], double [])))
{
  long nstp, i;
  double xsav, x, hnext, hdid, h;
  double *yscal, *y, *dydx;

  yscal = (double *) malloc( (size_t)(nvar+2)*sizeof(double) );
  y = (double *) malloc( (size_t)(nvar+2)*sizeof(double) );
  dydx = (double *) malloc( (size_t)(nvar+2)*sizeof(double) );
  x = x1;
  h = SIGN(h1, x2-x1);
  *nok = (*nbad) = kout = 0;
  for (i=1; i<=nvar; i++) y[i]=ystart[i];
  if (kmax > 0) xsav=x-dxsav*2.0;
  for (nstp=1;nstp<=MAXSTP;nstp++)
    {
      (*derivs)(x,y,dydx);
      for (i=1;i<=nvar;i++)
	yscal[i] = fabs(y[i])+fabs(dydx[i]*h)+TINY;
	/* yscal[i] = fabs(eps*dydx[i]*h)+TINY;  */

      if (kmax > 0 && kout < kmax-1 && fabs(x-xsav) > fabs(dxsav))
	{
	  xp[++kout] = x;
	  for (i=1; i<=nvar; i++) 
	    {
	      /*	    Its possible to print here the results  */
	      yp[i][kout]=y[i];
	    }
	  /*	  (*print_result)(x, y); */
	  xsav = x;
	}
      if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
      (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
      if (hdid == h) ++(*nok); else ++(*nbad);
      if ((x-x2)*(x2-x1) >= 0.0)
	{
	  for (i=1;i<=nvar;i++) ystart[i]=y[i];
	  if (kmax)
	    {
/*	      (*print_result)(x, y);  */
	      xp[++kout]=x;
	      for (i=1;i<=nvar;i++) yp[i][kout]=y[i];
	      /*		Its possible to print here the results   */
	    }
	  free(dydx);
	  free(y);
	  free(yscal);
	  return;
	}
      if (fabs(hnext) <= hmin) 
	{
	  printf("\n\nStep size too small in odeint\n\n");
	  /*	  nrerror("Step size too small in odeint");  */
	  free(dydx);
	  free(y);
	  free(yscal);
	  return;
	}
      h=hnext;
    }

  printf( "\n\nToo many steps in routine odeint\n\n" );
/*
  nrerror("Too many steps in routine odeint");
*/
}



#undef MAXSTP
#undef TINY





/********************************************************/
/*                                                      */
/*                    Wave lets                         */
/*                                                      */
/********************************************************/



void wt1(double a[], unsigned long n, int isign,
	 void (*wtstep)(double [], unsigned long, int))
{
  unsigned long nn;

  if (n < 4) return;
  if (isign >= 0)
    {
      for (nn=n;nn>=4;nn>>=1)
	(*wtstep)(a,nn,isign);
    }
  else
    {
      for (nn=4;nn<=n;nn<<=1)
	(*wtstep)(a,nn,isign);
    }
}



void my_wt1(double a[], unsigned long n, unsigned long m, int isign,
	 void (*wtstep)(double [], unsigned long, int))
{
  unsigned long nn, i;

  if (m==0) return;

  i=(n>>(m-1));

  if (isign >= 0)
    {
      for (nn=n;nn>=i;nn>>=1)
	(*wtstep)(a,nn,isign);
    }
  else
    {
      for (nn=i;nn<=n;nn<<=1)
	(*wtstep)(a,nn,isign);
    }
}




typedef struct
{
  int ncof,ioff,joff;
  double *cc, *cr;
} wavefilt;

wavefilt wfilt;


void pwtset(int n)
{
  void nrerror(char error_text[]);
  int k;
  double sig = -1.0;
  static double c2[3]={0.0,0.7071067811865475,0.7071067811865475};
  static double c4[5]={0.0,0.4829629131445341,0.8365163037378079,
		       0.2241438680420134,-0.1294095225512604};
  static double c6[7]={0.0,-.3326705529500826, -.8068915093110594, 
		       -.459877502118546, .1350110200101011, 
		       .8544127388180240e-1, -.3522629188572628e-1};
  static double c8[9]={0.0,-.2303778133088965, -.7148465705529875, 
		       -.6308807679308086, .2798376941690793e-1, 
		       .1870348117194144, -.3084138183561371e-1, 
		       -.3288301166694168e-1, .1059740178508723e-1};
  static double c10[11]={0.0,0.160102397974,.603829269797,.724308528438,
			 .138428145901,-.242294887066,-.032244869585,
			 0.077571493840,-.006241490213,-.01258075199,
			 .003335725285};
  static double c12[13]={0.0,0.111540743350,0.494623890398,
			 0.751133908021,0.315250351709,
			 -0.226264693965,-0.129766867567,
			 0.097501605587,0.027522865530,
			 -0.031582039318,0.000553842201,
			 0.004777257511,-0.001077301085};

  static double c14[15]={0.0,0.077852054085,0.396539319482,
			 0.729132090846,0.469782287405,
			 -0.143906003929,-0.224036184994,
			 0.071309219267,0.080612609151,
			 -0.038029936935,-0.016574541631,
			 0.012550998556,0.000429577973,
			 -0.001801640704,0.000353713800};

  static double c16[17]={0.0,0.054415842243,0.312871590914,
			 0.675630736297,0.585354683654,
			 -0.015829105256,-0.284015542962,
			 0.000472484574,0.128747426620,
			 -0.017369301002,-0.044088253931,
			 0.013981027917,0.008746094047,
			 -0.004870352993,-0.000391740373,
			 0.000675449406,-0.000117476784};

  static double c18[19]={0.0,0.038077947364,0.243834674613,
			 0.604823123690,0.657288078051,
			 0.133197385825,-0.293273783279,
			 -0.096840783223,0.148540749338,
			 0.030725681479,-0.067632829061,
			 0.000250947115,0.022361662124,
			 -0.004723204758,-0.004281503682,
			 0.001847646883,0.000230385764,
			 -0.000251963189,0.000039347320};

  static double c20[21]={0.0,0.026670057901,0.188176800078,
			 0.527201188932,0.688459039454,
			 0.281172343661,-0.249846424327,
			 -0.195946274377,0.127369340336,
			 0.093057364604,-0.071394147166,
			 -0.029457536822,0.033212674059,
			 0.003606553567,-0.010733175483,
			 0.001395351747,0.001992405295,
			 -0.000685856695,-0.000116466855,
			 0.000093588670,-0.000013264203};
  

  static double c2r[3],c4r[5],c6r[7],c8r[9],c10r[11],c12r[13],c14r[15],
    c16r[17],c18r[19],c20r[21];

  wfilt.ncof=n;
  switch(n)
    {
    case 2:
      wfilt.cc=c2;
      wfilt.cr=c2r;
      break;
    case 4:
      wfilt.cc=c4;
      wfilt.cr=c4r;
      break;
    case 6:
      wfilt.cc=c6;
      wfilt.cr=c6r;
      break;
    case 8:
      wfilt.cc=c8;
      wfilt.cr=c8r;
      break;
    case 10:
      wfilt.cc=c10;
      wfilt.cr=c10r;
      break;
    case 12:
      wfilt.cc=c12;
      wfilt.cr=c12r;
      break;
    case 14:
      wfilt.cc=c14;
      wfilt.cr=c14r;
      break;
    case 16:
      wfilt.cc=c16;
      wfilt.cr=c16r;
      break;
    case 18:
      wfilt.cc=c18;
      wfilt.cr=c18r;
      break;
    case 20:
      wfilt.cc=c20;
      wfilt.cr=c20r;
      break;
    default :
      nrerror("unimplemented value n in pwtset");
      break;
    }      
  
  for (k=1;k<=n;k++)
    {
      wfilt.cr[wfilt.ncof+1-k]=sig*wfilt.cc[k];
      sig = -sig;
    }
  wfilt.ioff = wfilt.joff = -(n >> 1);
}


extern wavefilt wfilt;

void pwt(double a[], unsigned long n, int isign)
{
  double ai, ai1, *wksp;
  unsigned long i, ii, j, jf, jr, k, n1, ni, nj, nh, nmod;

  if (n<4) return;
  wksp = (double *)malloc( (size_t)(n+1)*sizeof(double) );
  nmod = wfilt.ncof*n;
  n1 = n-1;
  nh = n >> 1;
  for (j=1;j<=n;j++) wksp[j] = 0.0;
  if (isign>=0)
    {
      for (ii=1,i=1;i<=n;i+=2,ii++)
	{
	  ni = i+nmod+wfilt.ioff;
	  nj = i+nmod+wfilt.joff;
	  for (k=1;k<=wfilt.ncof;k++)
	    {
	      jf=n1 & (ni+k);
	      jr=n1 & (nj+k);
	      wksp[ii] += wfilt.cc[k]*a[jf+1];
	      wksp[ii+nh] += wfilt.cr[k]*a[jr+1];
	    }
	}
    }
  else
    {
      for (ii=1,i=1;i<=n;i+=2,ii++)
	{
	  ai = a[ii];
	  ai1 = a[ii+nh];
	  ni = i+nmod+wfilt.ioff;
	  nj = i+nmod+wfilt.joff;
	  for (k=1;k<=wfilt.ncof;k++)
	    {
	      jf=(n1 & (ni+k))+1;
	      jr=(n1 & (nj+k))+1;
	      wksp[jf] += wfilt.cc[k]*ai;
	      wksp[jr] += wfilt.cr[k]*ai1;
	    }
	}
    }
  for (j=1;j<=n;j++) a[j]=wksp[j];
  free(wksp);
}


#define MOD(a,b) ((a)-((unsigned long)floor((double)(a)/(double)(b))*(b)))

void my_pwt(double a[], unsigned long n, int isign)
{
  double ai, ai1, *wksp;
  unsigned long i, ii, j, jf, jr, k, n1, ni, nj, nh, nmod;

  if (n<4) return;
  wksp = (double *)malloc( (size_t)(n+1)*sizeof(double) );
  nmod = wfilt.ncof*n;
  n1 = n-1;
  nh = n >> 1;
  for (j=1;j<=n;j++) wksp[j] = 0.0;
  if (isign>=0)
    {
      for (ii=1,i=1;i<=n;i+=2,ii++)
	{
	  ni = i+nmod+wfilt.ioff;
	  nj = i+nmod+wfilt.joff;
	  for (k=1;k<=wfilt.ncof;k++)
	    {
	      jf=MOD(ni+k,n);
	      jr=MOD(nj+k,n);
	      wksp[ii] += wfilt.cc[k]*a[jf+1];
	      wksp[ii+nh] += wfilt.cr[k]*a[jr+1];
	    }
	}
    }
  else
    {
      for (ii=1,i=1;i<=n;i+=2,ii++)
	{
	  ai = a[ii];
	  ai1 = a[ii+nh];
	  ni = i+nmod+wfilt.ioff;
	  nj = i+nmod+wfilt.joff;
	  for (k=1;k<=wfilt.ncof;k++)
	    {
	      jf=MOD(ni+k,n)+1;
	      jr=MOD(nj+k,n)+1;
	      wksp[jf] += wfilt.cc[k]*ai;
	      wksp[jr] += wfilt.cr[k]*ai1;
	    }
	}
    }
  for (j=1;j<=n;j++) a[j]=wksp[j];
  free(wksp);
}

#undef MOD




void wtn(double a[], unsigned long nn[], int ndim, int isign, 
	 void (*wtstep)(double [], unsigned long, int)) 
/*     Replaces a by its ndim-dimensional discrete wavelet transform, if 
       isign is input as 1. Here nn[1..ndim] is an integer array containing 
       the lengths of each dimension (number of real values), which MUST all 
       be powers of 2. a is a real array of length equal to the product of 
       these lengths, in which the data are stored as in a multidimensional 
       real array. If isign is input as -1, a is replaced by its inverse 
       wavelet transform. The routine wtstep, whose actual name must be 
       supplied in calling this routine, is the underlying wavelet filter. 
       Examples of wtstep are daub4 and (preceded by pwtset) pwt. 
*/
{ 
  unsigned long i1,i2,i3,k,n,nnew,nprev=1,nt,ntot=1,max_n=1; 
  int idim; 
  double *wksp; 

  for (idim=1;idim<=ndim;idim++) 
    {
      ntot *= nn[idim]; 
      if (nn[idim]>max_n)
        max_n=nn[idim];
    }
  wksp=vector(1,max_n); 

  for (idim=1;idim<=ndim;idim++) 
    { 
/*      Main loop over the dimensions. */
      n=nn[idim]; 
      nnew=n*nprev; 
      if (n > 4) 
	{ 
	  for (i2=0;i2<ntot;i2+=nnew) 
	    { 
	      for (i1=1;i1<=nprev;i1++) 
		{ 
		  for (i3=i1+i2,k=1;k<=n;k++,i3+=nprev) 
		    wksp[k]=a[i3]; 
		  /* Copy the relevant row or column or etc. into workspace.*/
		  if (isign >= 0) 
		    { 
		      /* Do one-dimensional wavelet transform. */
		      for(nt=n;nt>=4;nt >>= 1) 
			(*wtstep)(wksp,nt,isign); 
		    } 
		  else 
		    { 
		      /* Or inverse transform. */
		      for(nt=4;nt<=n;nt <<= 1) 
			(*wtstep)(wksp,nt,isign); 
		    } 
		  for (i3=i1+i2,k=1;k<=n;k++,i3+=nprev) 
		    a[i3]=wksp[k]; /* Copy back from workspace. */
		} 
	    }
	} 
      nprev=nnew; 
    } 
  free_vector(wksp,1,ntot); 
}





/******************************************************************/
/*                                                                */
/*   Statistics                                                   */
/*                                                                */
/******************************************************************/




double betai(double a, double b, double x) 
     /* Returns the incomplete beta function I x (a; b). */
{ 
  double betacf(double a, double b, double x); 
  double gammln(double xx); 
  void nrerror(char error_text[]); 
  double bt; 
  if (x < 0.0 || x > 1.0) 
    nrerror("Bad x in routine betai"); 
  if (x == 0.0 || x == 1.0) 
    bt=0.0; 
  else /* Factors in front of the continued fraction. */
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x)); 
  if (x < (a+1.0)/(a+b+2.0)) /* Use continued fraction directly. */
    return bt*betacf(a,b,x)/a; 
  else /* Use continued fraction after making the symmetry transformation. */
    return 1.0-bt*betacf(b,a,1.0-x)/b; 
}


double gammln(double xx) 
     /* Returns the value ln[gamma(xx)] for xx>0  */
{ 
  /* Internal arithmetic will be done in double precision, a nicety that you 
     can omit if  ve- gure accuracy is good enough. */
  double x,y,tmp,ser; 
  static double cof[6]={76.18009172947146,-86.50532032941677, 
			  24.01409824083091,-1.231739572450155, 
			  0.1208650973866179e-2,-0.5395239384953e-5}; 
  long j; 

  y=x=xx; 
  tmp=x+5.5; 
  tmp -= (x+0.5)*log(tmp); 
  ser=1.000000000190015; 
  for (j=0;j<=5;j++) 
    ser += cof[j]/++y; 
  return -tmp+log(2.5066282746310005*ser/x); 
}


#define MAXIT 100 
#define EPS 3.0e-7 
#define FPMIN 1.0e-30 

double betacf(double a, double b, double x) 
/* Used by betai: Evaluates continued fraction for incomplete beta function 
   by modi ed Lentz's method ( x 5.2). */
{ 
  void nrerror(char error_text[]);
  long m,m2; 
  double aa,c,d,del,h,qab,qam,qap; 

  qab=a+b; /* These q's will be used in factors that occur in the coe cients 
	      (6.4.6). */
  qap=a+1.0; 
  qam=a-1.0; 
  c=1.0; /* First step of Lentz's method. */
  d=1.0-qab*x/qap; 
  if (fabs(d) < FPMIN) 
    d=FPMIN; 
  d=1.0/d; 
  h=d; 
  for (m=1;m<=MAXIT;m++) 
    { 
      m2=2*m; 
      aa=m*(b-m)*x/((qam+m2)*(a+m2)); 
      d=1.0+aa*d; /* One step (the even one) of the recurrence. */
      if (fabs(d) < FPMIN) 
	d=FPMIN; 
      c=1.0+aa/c; 
      if (fabs(c) < FPMIN) 
	c=FPMIN; 
      d=1.0/d; 
      h *= d*c; 
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2)); 
      d=1.0+aa*d; /* Next step of the recurrence (the odd one). */
      if (fabs(d) < FPMIN) 
	d=FPMIN; 
      c=1.0+aa/c; 
      if (fabs(c) < FPMIN) 
	c=FPMIN; 
      d=1.0/d; 
      del=d*c; 
      h *= del; 
      if (fabs(del-1.0) < EPS) 
	break; /* Are we done? */
    } 
  if (m > MAXIT) 
    nrerror("a or b too big, or MAXIT too small in betacf"); 
  return h; 
}

#undef MAXIT
#undef EPS
#undef FPMIN




void avevar(double data[], unsigned long n, double *ave, double *var) 
     /* Given array data[1..n], returns its mean as ave and its variance as 
	var. */
{ 
  unsigned long j; 
  double s,ep; 
  for (*ave=0.0,j=1;j<=n;j++) 
    *ave += data[j]; 
  *ave /= n; 
  *var=ep=0.0; 
  for (j=1;j<=n;j++) 
    { 
      s=data[j]-(*ave); 
      ep += s; 
      *var += s*s; 
    } 
  *var=(*var-ep*ep/n)/(n-1); /* Corrected two-pass formula (14.1.8). */
}


void moment(double data[], long n, double *ave, double *adev, double *sdev, 
	    double *var, double *skew, double *curt) 
     /* Given an array of data[1..n], this routine returns its mean ave, 
	average deviation adev, standard deviation sdev, variance var, 
	skewness skew, and kurtosis curt. */
{ 
  void nrerror(char error_text[]); 
  long j; 
  double ep=0.0,s,p; 
  if (n <= 1) 
    nrerror("n must be at least 2 in moment"); 
  s=0.0; /* First pass to get the mean. */
  for (j=1;j<=n;j++) 
    s += data[j]; 
  *ave=s/n; 
  *adev=(*var)=(*skew)=(*curt)=0.0; /* Second pass to get the  rst (absolute), 
				       second, third, and fourth moments of 
				       the deviation from the mean. */
  for (j=1;j<=n;j++) 
    { 
      *adev += fabs(s=data[j]-(*ave)); 
      ep += s; 
      *var += (p=s*s); 
      *skew += (p *= s); 
      *curt += (p *= s); 
    } 
  *adev /= n; 
  *var=(*var-ep*ep/n)/(n-1); /* Corrected two-pass formula. */
  *sdev=sqrt(*var); /* Put the pieces together according to the conventional 
		       definitions. */
  if (*var) 
    { 
      *skew /= (n*(*var)*(*sdev)); 
      *curt=(*curt)/(n*(*var)*(*var))-3.0; 
    } 
  else 
    nrerror("No skew/kurtosis when variance = 0 (in moment)"); 
}



void ttest(double data1[], unsigned long n1, double data2[], unsigned long n2, 
	   double *t, double *prob) 
     /* Given the arrays data1[1..n1] and data2[1..n2], this routine returns 
	Student's t as t, and its significance as prob, small values of prob 
	indicating that the arrays have significantly different means. 
	The data arrays are assumed to be drawn from populations with the
	same true variance. */
{ 
  void avevar(double data[], unsigned long n, double *ave, double *var); 
  double betai(double a, double b, double x); 
  double var1,var2,svar,df,ave1,ave2; 

  avevar(data1,n1,&ave1,&var1); 
  avevar(data2,n2,&ave2,&var2); 
  df=n1+n2-2; /* Degrees of freedom. */
  svar=((n1-1)*var1+(n2-1)*var2)/df; /* Pooled variance. */
  *t=(ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2)); 
  *prob=betai(0.5*df,0.5,df/(df+(*t)*(*t))); /* See equation (6.4.9). */
}




void tutest(double data1[], unsigned long n1, double data2[], 
	    unsigned long n2, double *t, double *prob) 
     /* Given the arrays data1[1..n1] and data2[1..n2], this routine returns 
	Student's t as t, and its signi cance as prob, small values of prob 
	indicating that the arrays have significantly differ-ent means. The 
	data arrays are allowed to be drawn from populations with unequal 
	variances.  */
{ 
  void avevar(double data[], unsigned long n, double *ave, double *var); 
  double betai(double a, double b, double x); 
  double var1,var2,df,ave1,ave2;

  avevar(data1,n1,&ave1,&var1); 
  avevar(data2,n2,&ave2,&var2); 
  *t=(ave1-ave2)/sqrt(var1/n1+var2/n2); 
  df=SQR(var1/n1+var2/n2)/(SQR(var1/n1)/(n1-1)+SQR(var2/n2)/(n2-1)); 
  *prob=betai(0.5*df,0.5,df/(df+SQR(*t))); 
}




void tptest(double data1[], double data2[], unsigned long n, double *t, 
	    double *prob) 
/* Given the paired arrays data1[1..n] and data2[1..n], this routine returns 
   Student's t for paired data as t, and its significance as prob, small 
   values of prob indicating a significant difference of means. */
{ 
  void avevar(double data[], unsigned long n, double *ave, double *var); 
  double betai(double a, double b, double x); 
  unsigned long j; 
  double var1,var2,ave1,ave2,sd,df,cov=0.0; 

  avevar(data1,n,&ave1,&var1); 
  avevar(data2,n,&ave2,&var2); 
  for (j=1;j<=n;j++) 
    cov += (data1[j]-ave1)*(data2[j]-ave2); 
  cov /= df=n-1; 
  sd=sqrt((var1+var2-2.0*cov)/n); 
  *t=(ave1-ave2)/sd; 
  *prob=betai(0.5*df,0.5,df/(df+(*t)*(*t))); 
}




void ftest(double data1[], unsigned long n1, double data2[], unsigned long n2, 
	   double *f, double *prob) 
/* Given the arrays data1[1..n1] and data2[1..n2], this routine returns the 
   value of f, and its significance as prob. Small values of prob indicate 
   that the two arrays have significantly different variances. */
{ 
  void avevar(double data[], unsigned long n, double *ave, double *var); 
  double betai(double a, double b, double x); 
  double var1,var2,ave1,ave2,df1,df2; 

  avevar(data1,n1,&ave1,&var1); 
  avevar(data2,n2,&ave2,&var2); 
  if (var1 > var2) 
    { 
      /* Make F the ratio of the larger variance to the smaller one. */
      *f=var1/var2; 
      df1=n1-1; 
      df2=n2-1; 
    } 
  else 
    { 
      *f=var2/var1; 
      df1=n2-1; 
      df2=n1-1; 
    } 
  *prob = 2.0*betai(0.5*df2,0.5*df1,df2/(df2+df1*(*f))); 
  if (*prob > 1.0) *prob=2.0-*prob; 
}





/******************************************************************/
/*                                                                */
/*   Fourier transform                                            */
/*                                                                */
/******************************************************************/





#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr 
void four1(double data[], unsigned long nn, long isign) 
/* Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input 
   as 1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier 
   transform, if isign is input as 1; or replaces data[1..2*nn] by nn times 
   its inverse discrete Fourier transform, if isign is input as -1. data is a 
   complex array of length nn or, equivalently, a real array of length 2*nn. 
   nn MUST be an longeger power of 2 (this is not checked for!). */
{ 
  unsigned long n,mmax,m,j,istep,i; 
  double wtemp,wr,wpr,wpi,wi,theta; /* Double precision for the trigonometric 
				       recurrences. */
  double tempr,tempi; 
  n=nn << 1; j=1; 
  for (i=1;i<n;i+=2) 
    { /*This is the bit-reversal section of the routine. */
      if (j > i) 
	{ 
	  SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
	  SWAP(data[j+1],data[i+1]); 
	} 
      m=n >> 1; 
      while (m >= 2 && j > m) 
	{ 
	  j -= m; 
	  m >>= 1; 
	} 
      j += m; 
    } /* Here begins the Danielson-Lanczos section of the routine. */
  mmax=2; 
  while (n > mmax) 
    { /* Outer loop executed log 2 nn times. */
      istep=mmax << 1; 
      theta=isign*(6.28318530717959/mmax); 
      /* Initialize the trigonometric recurrence. */
      wtemp=sin(0.5*theta); 
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta); 
      wr=1.0; 
      wi=0.0; 
      for (m=1;m<mmax;m+=2) 
	{ /* Here are the two nested inner loops. */
	  for (i=m;i<=n;i+=istep) 
	    { 
	      j=i+mmax; /* This is the Danielson-Lanczos formula: */
	      tempr=wr*data[j]-wi*data[j+1]; 
	      tempi=wr*data[j+1]+wi*data[j]; 
	      data[j]=data[i]-tempr; 
	      data[j+1]=data[i+1]-tempi; 
	      data[i] += tempr; 
	      data[i+1] += tempi; 
	    } 
	  wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
	  wi=wi*wpr+wtemp*wpi+wi; 
	} 
      mmax=istep;
    }
}

#undef SWAP



#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr 
void fourn(double data[], unsigned long nn[], long ndim, long isign) 
/* Replaces data by its ndim-dimensional discrete Fourier transform, if isign 
   is input as 1. nn[1..ndim] is an longeger array containing the lengths of 
   each dimension (number of complex values), which MUST all be powers of 2. 
   data is a real array of length twice the product of these lengths, in 
   which the data are stored as in a multidimensional complex array: real and 
   imaginary parts of each element are in consecutive locations, and the 
   rightmost index of the array increases most rapidly as one proceeds along 
   data. For a two-dimensional array, this is equivalent to storing the array 
   by rows. If isign is input as -1, data is replaced by its inverse transform
   times the product of the lengths of all dimensions.*/
{ 
  long idim; 
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2; 
  unsigned long ibit,k1,k2,n,nprev,nrem,ntot; 
  double tempi,tempr; 
  double theta,wi,wpi,wpr,wr,wtemp; /* Double precision for trigonometric 
				       recurrences. */
  for (ntot=1,idim=1;idim<=ndim;idim++) /* Compute total number of complex 
					   values. */
    ntot *= nn[idim]; 
  nprev=1; 
  for (idim=ndim;idim>=1;idim--) 
    { /* Main loop over the dimensions. */
      n=nn[idim]; 
      nrem=ntot/(n*nprev); 
      ip1=nprev << 1; 
      ip2=ip1*n; 
      ip3=ip2*nrem; 
      i2rev=1; 
      for (i2=1;i2<=ip2;i2+=ip1) 
	{ /* This is the bit-reversal section of the routine. */
	  if (i2 < i2rev) 
	    { 
	      for (i1=i2;i1<=i2+ip1-2;i1+=2) 
		{ 
		  for (i3=i1;i3<=ip3;i3+=ip2) 
		    { 
		      i3rev=i2rev+i3-i2; 
		      SWAP(data[i3],data[i3rev]); 
		      SWAP(data[i3+1],data[i3rev+1]); 
		    } 
		} 
	    } 
	  ibit=ip2 >> 1; 
	  while (ibit >= ip1 && i2rev > ibit) 
	    { 
	      i2rev -= ibit; 
	      ibit >>= 1;
	    } 
	  i2rev += ibit; 
	} 
      ifp1=ip1; /* Here begins the Danielson-Lanczos section of the routine.*/
      while (ifp1 < ip2) 
	{ 
	  ifp2=ifp1 << 1; 
	  theta=isign*6.28318530717959/(ifp2/ip1); 
	  /* Initialize for the trig. recurrence. */
	  wtemp=sin(0.5*theta); 
	  wpr = -2.0*wtemp*wtemp; 
	  wpi=sin(theta); 
	  wr=1.0; 
	  wi=0.0; 
	  for (i3=1;i3<=ifp1;i3+=ip1) 
	    { 
	      for (i1=i3;i1<=i3+ip1-2;i1+=2) 
		{ 
		  for (i2=i1;i2<=ip3;i2+=ifp2) 
		    { 
		      k1=i2; /* Danielson-Lanczos formula: */
		      k2=k1+ifp1; 
		      tempr=(double)wr*data[k2]-(double)wi*data[k2+1]; 
		      tempi=(double)wr*data[k2+1]+(double)wi*data[k2]; 
		      data[k2]=data[k1]-tempr; 
		      data[k2+1]=data[k1+1]-tempi; 
		      data[k1] += tempr; 
		      data[k1+1] += tempi; 
		    } 
		} 
	      wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence.*/
	      wi=wi*wpr+wtemp*wpi+wi; 
	    } 
	  ifp1=ifp2; 
	} 
      nprev *= n; 
    } 
}

#undef SWAP






void realft(double data[], unsigned long n, long isign) 
     /* Calculates the Fourier transform of a set of n real-valued data 
	polongs. Replaces this data (which is stored in array data[1..n]) by 
	the positive frequency half of its complex Fourier transform. The 
	real-valued  rst and last components of the complex transform are 
	returned as elements data[1] and data[2], respectively. n must be a 
	power of 2. This routine also calculates the inverse transform of a 
	complex data array if it is the transform of real data. (Result in 
	this case must be multiplied by 2/n.) */
{ 
  void four1(double data[], unsigned long nn, long isign); 
  unsigned long i,i1,i2,i3,i4,np3; 
  double c1=0.5,c2,h1r,h1i,h2r,h2i; 
  double wr,wi,wpr,wpi,wtemp,theta; 
  /* Double precision for the trigonometric recurrences. */
  theta=3.141592653589793/(double) (n>>1); /* Initialize the recurrence. */
  if (isign == 1) 
    { 
      c2 = -0.5; 
      four1(data,n>>1,1); /* The forward transform is here. */
    } 
  else 
    { 
      c2=0.5; /* Otherwise set up for an inverse transform. */
      theta = -theta; 
    } 
  wtemp=sin(0.5*theta); 
  wpr = -2.0*wtemp*wtemp; 
  wpi=sin(theta); 
  wr=1.0+wpr; 
  wi=wpi; 
  np3=n+3; 
  for (i=2;i<=(n>>2);i++) 
    { /* Case i=1 done separately below. */
      i4=1+(i3=np3-(i2=1+(i1=i+i-1))); 
      h1r=c1*(data[i1]+data[i3]); /* The two separate transforms are 
				     separated out of data. */
      h1i=c1*(data[i2]-data[i4]); 
      h2r = -c2*(data[i2]+data[i4]); 
      h2i=c2*(data[i1]-data[i3]); 
      data[i1]=h1r+wr*h2r-wi*h2i; /* Here they are recombined to form the 
				    true transform of the original real data.*/
      data[i2]=h1i+wr*h2i+wi*h2r; 
      data[i3]=h1r-wr*h2r+wi*h2i; 
      data[i4] = -h1i+wr*h2i+wi*h2r; 
      wr=(wtemp=wr)*wpr-wi*wpi+wr; /* The recurrence. */
      wi=wi*wpr+wtemp*wpi+wi; 
    } 
  if (isign == 1) 
    {
      data[1] = (h1r=data[1])+data[2]; /*Squeeze the  rst and last data 
					 together to get them all within the 
					 original array. */
      data[2] = h1r-data[2]; 
  } 
  else 
    { 
      data[1]=c1*((h1r=data[1])+data[2]); 
      data[2]=c1*(h1r-data[2]); 
      four1(data,n>>1,-1); /* This is the inverse transform for the case 
			      isign=-1. */
    } 
}








#define TWOPID 6.2831853071795865 
void period(double x[], double y[], long n, double ofac, double hifac, 
	    double px[], double py[], long np, long *nout, long *jmax, 
	    double *prob) 
/* Given n data polongs with abscissas x[1..n] (which need not be equally 
spaced) and ordinates y[1..n], and given a desired oversampling factor ofac 
(a typical value being 4 or larger), this routine  lls array px[1..np] with 
an increasing sequence of frequencies (not angular frequencies) up to hifac 
times the \average" Nyquist frequency, and  lls array py[1..np] with the 
values of the Lomb normalized periodogram at those frequencies. The arrays x 
and y are not altered. np, the dimension of px and py, must be large enough 
to contain the output, or an error results. The routine also returns jmax 
such that py[jmax] is the maximum element in py, and prob, an estimate of the 
signi cance of that maximum against the hypothesis of random noise. A small 
value of prob indicates that a signi cant periodic signal is present. */
{ 
  void avevar(double data[], unsigned long n, double *ave, double *var); 
  long i,j; 
  double ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh, 
  sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy; 
  double arg,wtemp,*wi,*wpi,*wpr,*wr; 
  wi=vector(1,n); 
  wpi=vector(1,n); 
  wpr=vector(1,n); 
  wr=vector(1,n); 
  *nout=0.5*ofac*hifac*n; 
  if (*nout > np) 
    nrerror("output arrays too short in period"); 
  avevar(y,n,&ave,&var); /* Get mean and variance of the input data. */
  xmax=xmin=x[1]; /* Go through data to get the range of abscis- sas. */
  for (j=1;j<=n;j++) 
    { 
      if (x[j] > xmax) 
	xmax=x[j]; 
      if (x[j] < xmin) 
	xmin=x[j]; 
    } 
  xdif=xmax-xmin; 
  xave=0.5*(xmax+xmin); 
  pymax=0.0; 
  pnow=1.0/(xdif*ofac); /* Starting frequency. */
  for (j=1;j<=n;j++) 
    { 
      /* Initialize values for the trigonometric recurrences at each data 
	 polong. The recur- rences are done in double precision. */
      arg=TWOPID*((x[j]-xave)*pnow); 
      wpr[j] = -2.0*SQR(sin(0.5*arg));/*rences are done in double precision.*/
      wpi[j]=sin(arg); 
      wr[j]=cos(arg); 
      wi[j]=wpi[j]; 
    } 
  for (i=1;i<=(*nout);i++) 
    { 
      /* Main loop over the frequencies to be evaluated. */
      px[i]=pnow; 
      sumsh=sumc=0.0;/*First,loop over the data to get and related quantities*/
      for (j=1;j<=n;j++) 
	{ 
	  c=wr[j]; 
	  s=wi[j]; 
	  sumsh += s*c; 
	  sumc += (c-s)*(c+s); 
	} 
      wtau=0.5*atan2(2.0*sumsh,sumc); 
      swtau=sin(wtau); 
      cwtau=cos(wtau); 
      sums=sumc=sumsy=sumcy=0.0; 
      /* Then, loop over the data again to get the periodogram value. */
      for (j=1;j<=n;j++) 
	{ 
	  s=wi[j]; 
	  c=wr[j]; 
	  ss=s*cwtau-c*swtau; 
	  cc=c*cwtau+s*swtau; 
	  sums += ss*ss; 
	  sumc += cc*cc; 
	  yy=y[j]-ave; 
	  sumsy += yy*ss; 
	  sumcy += yy*cc; 
	  wr[j]=((wtemp=wr[j])*wpr[j]-wi[j]*wpi[j])+wr[j]; 
	  /* Update the trigonometric recurrences. */
	  wi[j]=(wi[j]*wpr[j]+wtemp*wpi[j])+wi[j]; 
	} 
      py[i]=0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var; 
      if (py[i] >= pymax) 
	pymax=py[(*jmax=i)];
      pnow += 1.0/(ofac*xdif); /* The next frequency. */
    } 
  expy=exp(-pymax); /* Evaluate statistical signi cance of the maximum. */
  effm=2.0*(*nout)/ofac; 
  *prob=effm*expy; 
  if (*prob > 0.01) 
    *prob=1.0-pow(1.0-expy,effm); 
  free_vector(wr,1,n); 
  free_vector(wpr,1,n); 
  free_vector(wpi,1,n); 
  free_vector(wi,1,n); 
} 

#undef TWOPID 








#define MOD(a,b) while(a >= b) a -= b; /* Positive numbers only. */
#define MACC 4 
/* Number of longerpolation polongs per 1/4 cycle of highest frequency. */
void fasper(double x[], double y[], unsigned long n, double ofac, 
	    double hifac, double wk1[], double wk2[], unsigned long nwk, 
	    unsigned long *nout, unsigned long *jmax, double *prob)
     /* Given n data polongs with abscissas x[1..n] (which need not be equally 
	spaced) and ordinates y[1..n], and given a desired oversampling 
	factor ofac (a typical value being 4 or larger), this routine  lls 
	array wk1[1..nwk] with a sequence of nout increasing frequencies (not 
	angular frequencies) up to hifac times the \average" Nyquist 
	frequency, and  lls array wk2[1..nwk] with the values of the Lomb 
	normalized periodogram at those frequencies. The arrays x and y are 
	not altered. nwk, the dimension of wk1 and wk2, must be large enough 
	for longermediate work space, or an error results. The routine also 
	returns jmax such that wk2[jmax] is the maximum element in wk2, and 
	prob, an estimate of the signi cance of that maximum against the 
	hypothesis of random noise. A small value of prob indicates that a 
	significant periodic signal is present.  */
{ 
  void avevar(double data[], unsigned long n, double *ave, double *var); 
  void realft(double data[], unsigned long n, long isign); 
  void spread(double y, double yy[], unsigned long n, double x, long m); 
  unsigned long j,k,ndim,nfreq,nfreqt; 
  double ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt; 
  double hs2wt,hypo,pmax,sterm,swt,var,xdif,xmax,xmin;
  *nout=0.5*ofac*hifac*n; 
  nfreqt=ofac*hifac*n*MACC; /* Size the FFT as next power of 2 above nfreqt. */
  nfreq=64; 
  while (nfreq < nfreqt) 
    nfreq <<= 1; 
  ndim=nfreq << 1; 
  if (ndim > nwk) 
    { 
      printf("nwk should be larger than %d\n",ndim/n);
      nrerror("workspaces too small in fasper"); 
    }
  avevar(y,n,&ave,&var); /* Compute the mean, variance,and range of the data.*/
  xmin=x[1]; 
  xmax=xmin; 
  for (j=2;j<=n;j++) 
    { 
      if (x[j] < xmin) 
	xmin=x[j]; 
      if (x[j] > xmax) 
	xmax=x[j]; 
    } 
  xdif=xmax-xmin; 
  for (j=1;j<=ndim;j++) 
    wk1[j]=wk2[j]=0.0; /* Zero the workspaces. */
  fac=ndim/(xdif*ofac); 
  fndim=ndim; 
  for (j=1;j<=n;j++) 
    { 
      /* Extirpolate the data longo the workspaces. */
      ck=(x[j]-xmin)*fac; 
      MOD(ck,fndim) 
	ckk=2.0*(ck++); 
      MOD(ckk,fndim) 
	++ckk; 
      spread(y[j]-ave,wk1,ndim,ck,MACC); 
      spread(1.0,wk2,ndim,ckk,MACC); 
    } 
  realft(wk1,ndim,1); /* Take the Fast Fourier Transforms. */
  realft(wk2,ndim,1); 
  df=1.0/(xdif*ofac); 
  pmax = -1.0; 
  for (k=3,j=1;j<=(*nout);j++,k+=2) 
    { 
      /* Compute the Lomb value for each frequency. */
      hypo=sqrt(wk2[k]*wk2[k]+wk2[k+1]*wk2[k+1]); 
      hc2wt=0.5*wk2[k]/hypo; 
      hs2wt=0.5*wk2[k+1]/hypo; 
      cwt=sqrt(0.5+hc2wt); 
      swt=SIGN(sqrt(0.5-hc2wt),hs2wt); 
      den=0.5*n+hc2wt*wk2[k]+hs2wt*wk2[k+1]; 
      cterm=SQR(cwt*wk1[k]+swt*wk1[k+1])/den; 
      sterm=SQR(cwt*wk1[k+1]-swt*wk1[k])/(n-den); 
      wk1[j]=j*df; 
      wk2[j]=(cterm+sterm)/(2.0*var); 
      if (wk2[j] > pmax) pmax=wk2[(*jmax=j)]; 
    } 
  expy=exp(-pmax); /* Estimate signi cance of largest peak value. */
  effm=2.0*(*nout)/ofac; 
  *prob=effm*expy; 
  if (*prob > 0.01) 
    *prob=1.0-pow(1.0-expy,effm); 
}

#undef MOD 
#undef MACC 


void spread(double y, double yy[], unsigned long n, double x, long m) 
     /* Given an array yy[1..n], extirpolate (spread) a value y longo m actual 
	array elements that best approximate the \ ctional" (i.e., possibly 
	nonlongeger) array element number x. The weights used are coe cients 
	of the Lagrange longerpolating polynomial. */
{ 
  long ihi,ilo,ix,j,nden; 
  static long nfac[11]={0,1,1,2,6,24,120,720,5040,40320,362880}; 
  double fac; 
  if (m > 10) 
    nrerror("factorial table too small in spread");
  ix=(long)x; 
  if (x == (double)ix) 
    yy[ix] += y; 
  else 
    { 
      ilo=LMIN(LMAX((long)(x-0.5*m+1.0),1),n-m+1); 
      ihi=ilo+m-1; 
      nden=nfac[m]; 
      fac=x-ilo; 
      for (j=ilo+1;j<=ihi;j++) 
	fac *= (x-j); 
      yy[ihi] += y*fac/(nden*(x-ihi)); 
      for (j=ihi-1;j>=ilo;j--) 
	{ 
	  nden=(nden/(j+1-ilo))*(j-ihi); 
	  yy[j] += y*fac/(nden*(x-j)); 
	} 
    } 
}




/******************************************************************/
/*                                                                */
/*   linear algebra                                                 */
/*                                                                */
/******************************************************************/

#define TINY 1.0e-20; /*  all number. */

long ludcmp(double **a, long n, long *indx, double *d) 
/* 
Given a matrix a[1..n][1..n], this routine replaces it by the LU 
decomposition of a rowwise permutation of itself. a and n are input. a is 
output, arranged as in equation (2.3.14) above; indx[1..n] is an output 
vector that records the row permutation e ected by the partial pivoting; 
d is output as   1 depending on whether the number of row longerchanges was 
even or odd, respectively. This routine is used in combination with lubksb to 
solve linear equations or invert a matrix.  */
{ 
  long i,imax,j,k; 
  double big,dum,sum,temp; 
  double *vv; /* vv stores the implicit scaling of each row. */
  
  vv=(double *) malloc( (size_t)(n+1)*sizeof(double) );
  
  *d=1.0;          /* No row longerchanges yet. */
  for (i=1;i<=n;i++) 
    { 
      /*  Loop over rows to get the implicit scaling informa- tion. */
      big=0.0; 
      for (j=1;j<=n;j++) 
	if ((temp=fabs(a[i][j])) > big) 
	  big=temp; 
      if (big == 0.0) 
	{
/* 	nrerror("Singular matrix in routine ludcmp"); */
	  return 0;
	}
      /* No nonzero largest element. */
      vv[i]=1.0/big;      /* Save the scaling. */
    } 
  for (j=1;j<=n;j++) 
    { 
      /*This is the loop over columns of Crout's method. */
      for (i=1;i<j;i++) 
	{ 
	  /* This is equation (2.3.12) except for i = j. */
	  sum=a[i][j]; 
	  for (k=1;k<i;k++) 
	    sum -= a[i][k]*a[k][j]; 
	  a[i][j]=sum; 
	} 
      big=0.0; 
      /* Initialize for the search for largest pivot element. */
      for (i=j;i<=n;i++) 
	{ 
	  /* This is i = j of equation (2.3.12) and i = j +1 : : : N of 
	     equation (2.3.13). */
	  sum=a[i][j]; 
	  for (k=1;k<j;k++)
	    sum -= a[i][k]*a[k][j]; 
	  a[i][j]=sum; 
	  if ( (dum=vv[i]*fabs(sum)) >= big) 
	    { 
	      /* Is the  gure of merit for the pivot better than the best 
		 so far? */
	      big=dum; 
	      imax=i; 
	    } 
	} 
      if (j != imax) 
	{ 
	  /* Do we need to longerchange rows? */
	  for (k=1;k<=n;k++) 
	    { 
	      /* Yes, do so... */
	      dum=a[imax][k]; 
	      a[imax][k]=a[j][k]; 
	      a[j][k]=dum; 
	    } 
	  *d = -(*d); /* ...and change the parity of d. */
	  vv[imax]=vv[j]; 
	  /* Also longerchange the scale factor. */
	} 
      indx[j]=imax; 
      if (a[j][j] == 0.0) 
	a[j][j]=TINY; 
      /* If the pivot element is zero the matrix is singular (at least to 
	 the precision of the algorithm). For some applications on 
	 singular matrices, it is desirable to substitute TINY for zero. */
      if (j != n) 
	{ 
	  /*Now,  nally, divide by the pivot element.  */
	  dum=1.0/(a[j][j]); 
	  for (i=j+1;i<=n;i++) 
	    a[i][j] *= dum; 
	} 
    } 
  /* Go back for the next column in the reduction. */
  free(vv); 
  return 1;
}

#undef TINY

void lubksb(double **a, long n, long *indx, double b[]) 
/* Solves the set of n linear equations A   X = B. Here a[1..n][1..n] is 
   input, not as the matrix A but rather as its LU decomposition, determined 
   by the routine ludcmp. indx[1..n] is input as the permutation vector 
   returned by ludcmp. b[1..n] is input as the right-hand side vector B, and 
   returns with the solution vector X. a, n, and indx are not modi ed by this 
   routine and can be left in place for successive calls with di erent 
   right-hand sides b. This routine takes longo account the possibility that b 
   will begin with many zero elements, so it is e cient for use in matrix 
   inversion.  */
{ 
  long i,ii=0,ip,j; 
  double sum; 

  for (i=1;i<=n;i++) 
    { 
      /* When ii is set to a positive value, it will become the index of the  
	 rst nonvanishing element of b. We now do the forward substitution, 
	 equation (2.3.6). The only new wrinkle is to unscramble the 
	 permutation as we go. */
      ip=indx[i]; 
      sum=b[ip]; 
      b[ip]=b[i]; 
      if (ii) 
	for (j=ii;j<=i-1;j++) 
	  sum -= a[i][j]*b[j]; 
      else 
	if (sum) 
	  ii=i; 
      /* A nonzero element was encountered, so from now on we will have to do 
	 the sums in the loop above. */
      b[i]=sum; 
    } 
  for (i=n;i>=1;i--) 
    { 
      /* Now we do the backsubstitution, equation (2.3.7). */
      sum=b[i]; 
      for (j=i+1;j<=n;j++) 
	sum -= a[i][j]*b[j]; 
      b[i]=sum/a[i][i]; 
      /* Store a component of the solution vector X. */
    } 
  /* All done! */
}


/*
   To summarize, this is the preferred way to solve the linear set of 
   equations A   x = b: 
   double **a,*b,d; 
   long n,*indx; 
   ... 
   ludcmp(a,n,indx,&d); 
   lubksb(a,n,indx,b); 
   The answer x will be given back in b. Your original matrix A will have 
   been destroyed. If you subsequently want to solve a set of equations with 
   the same A but a different right-hand side b, you repeat only 
   lubksb(a,n,indx,b); 
   not, of course, with the original matrix A, but with a and indx as were 
   already set by ludcmp.
*/



/* Inverse of a Matrix  */
/* Using the above LU decomposition and backsubstitution routines, it is 
   com-pletely straightforward to find the inverse of a matrix column by 
   column.  */

/*
#define N  ... 
double **a,**y,d,*col; 
long i,j,*indx; 
... 
ludcmp(a,N,indx,&d); 
Decompose the matrix just once. 
for(j=1;j<=N;j++) 
{ 
  Find inverse by columns. 
    for(i=1;i<=N;i++) 
      col[i]=0.0; 
  col[j]=1.0; 
  lubksb(a,N,indx,col); 
  for(i=1;i<=N;i++) 
    y[i][j]=col[i]; 
}
*/
/*
 The matrix y will now contain the inverse of the original matrix a, which 
 will have been destroyed. Alternatively, there is nothing wrong with using a 
 Gauss-Jordan routine like gaussj ( x 2.1) to invert a matrix in place, again 
 destroying the original. Both methods have practically the same operations 
 count.
*/



/* Calculation of a determinant thus requires one call to ludcmp, with no 
   subse-quent backsubstitutions by lubksb.   */
/* #define N ...
double **a,d; 
long j,*indx; 
... 
ludcmp(a,N,indx,&d); 
This returns d as   1.
for(j=1;j<=N;j++) 
     d *= a[j][j]; 
*/
/* The variable d now contains the determinant of the original matrix a, 
   which will have been destroyed.  */




void svbksb(double **u, double w[], double **v, long m, long n, double b[], 
	    double x[])
     /* Solves A  X = B for a vector X, where A is specied by the arrays 
	u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by svdcmp. m and n 
	are the dimensions of a, and will be equal for square matrices. 
	b[1..m] is the input right-hand side. x[1..n] is the output solution 
	vector. No input quantities are destroyed, so the routine may be 
	called sequentially with diferent b's. */
{
  long jj,j,i;
  double s,*tmp;
  tmp=vector(1,n);
  for (j=1;j<=n;j++) 
    { 
      /* Calculate U T B. */
      s=0.0;
      if (w[j]) 
	{ 
	  /* Nonzero result only if wj is nonzero. */
	  for (i=1;i<=m;i++) 
	    s += u[i][j]*b[i];
	  s /= w[j];        /*  This is the divide by wj . */
	}
      tmp[j]=s;
    }
  for (j=1;j<=n;j++) 
    { 
      /* Matrix multiply by V to get answer. */
      s=0.0;
      for (jj=1;jj<=n;jj++) 
	s += v[j][jj]*tmp[jj];
      x[j]=s;
    }
  free_vector(tmp,1,n);
}



double pythag(double a, double b)
/* Computes (a 2 + b 2 ) 1=2 without destructive under ow or over ow. */
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) 
    return absa*sqrt(1.0+SQR(absb/absa));
  else 
    return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


void svdcmp(double **a, long m, long n, double w[], double **v)
     /*
	Given a matrix a[1..m][1..n], this routine computes its singular 
	value decomposition, A = U  W  V T . The matrix U replaces a on 
	output. The diagonal matrix of singular values W is output as a 
	vector w[1..n]. The matrix V (not the transpose V T ) is output as 
	v[1..n][1..n].
	*/
{
  double pythag(double a, double b);
  long flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  rv1=vector(1,n);
  g=scale=anorm=0.0; /* Householder reduction to bidiagonal form. */
  for (i=1;i<=n;i++) 
    {
      l=i+1;
      rv1[i]=scale*g;
      g=s=scale=0.0;
      if (i <= m) 
	{
	  for (k=i;k<=m;k++) 
	    scale += fabs(a[k][i]);
	  if (scale) 
	    {
	      for (k=i;k<=m;k++) 
		{
		  a[k][i] /= scale;
		  s += a[k][i]*a[k][i];
		}
	      f=a[i][i];
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      a[i][i]=f-g;
	      for (j=l;j<=n;j++) 
		{
		  for (s=0.0,k=i;k<=m;k++) 
		    s += a[k][i]*a[k][j];
		  f=s/h;
		  for (k=i;k<=m;k++) 
		    a[k][j] += f*a[k][i];
		}
	      for (k=i;k<=m;k++) 
		a[k][i] *= scale;
	    }
	}
      w[i]=scale *g;
      g=s=scale=0.0;
      if (i <= m && i != n) 
	{
	  for (k=l;k<=n;k++) 
	    scale += fabs(a[i][k]);
	  if (scale) 
	    {
	      for (k=l;k<=n;k++) 
		{
		  a[i][k] /= scale;
		  s += a[i][k]*a[i][k];
		}
	      f=a[i][l];
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      a[i][l]=f-g;
	      for (k=l;k<=n;k++) 
		rv1[k]=a[i][k]/h;
	      for (j=l;j<=m;j++) 
		{
		  for (s=0.0,k=l;k<=n;k++) 
		    s += a[j][k]*a[i][k];
		  for (k=l;k<=n;k++) 
		    a[j][k] += s*rv1[k];
		}
	      for (k=l;k<=n;k++) 
		a[i][k] *= scale;
	    }
	}
      anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
  for (i=n;i>=1;i--) 
    { 
      /* Accumulation of right-hand transformations.  */
      if (i < n) 
	{
	  if (g) 
	    {
	      /* Double division to avoid possible under ow. */
	      for (j=l;j<=n;j++) 
		v[j][i]=(a[i][j]/a[i][l])/g;
	      for (j=l;j<=n;j++) 
		{
		  for (s=0.0,k=l;k<=n;k++) 
		    s += a[i][k]*v[k][j];
		  for (k=l;k<=n;k++) 
		    v[k][j] += s*v[k][i];
		}
	    }
	  for (j=l;j<=n;j++) 
	    v[i][j]=v[j][i]=0.0;
	}
      v[i][i]=1.0;
      g=rv1[i];
      l=i;
    }
  for (i=LMIN(m,n);i>=1;i--) 
    {
      /* Accumulation of left-hand transformations. */
      l=i+1;
      g=w[i];
      for (j=l;j<=n;j++) 
	a[i][j]=0.0;
      if (g) 
	{
	  g=1.0/g;
	  for (j=l;j<=n;j++) 
	    {
	      for (s=0.0,k=l;k<=m;k++) 
		s += a[k][i]*a[k][j];
	      f=(s/a[i][i])*g;
	      for (k=i;k<=m;k++) 
		a[k][j] += f*a[k][i];
	    }
	  for (j=i;j<=m;j++) 
	    a[j][i] *= g;
      } 
      else 
	for (j=i;j<=m;j++) 
	  a[j][i]=0.0;
      ++a[i][i];
    }
  for (k=n;k>=1;k--) 
    { 
      /* Diagonalization of the bidiagonal form: Loop over
	 singular values, and over allowed iterations. */
      for (its=1;its<=30;its++) 
	{
	  flag=1;
	  for (l=k;l>=1;l--) 
	    { 
	      /* Test for splitting. */
	      nm=l-1;     /* Note that rv1[1] is always zero. */
	      if ((double)(fabs(rv1[l])+anorm) == anorm) 
		{
		  flag=0;
		  break;
		}
	      if ((double)(fabs(w[nm])+anorm) == anorm) 
		break;
	    }
	  if (flag) 
	    {
	      c=0.0; /* Cancellation of rv1[l], if l > 1. */
	      s=1.0;
	      for (i=l;i<=k;i++) 
		{
		  f=s*rv1[i];
		  rv1[i]=c*rv1[i];
		  if ((double)(fabs(f)+anorm) == anorm) 
		    break;
		  g=w[i];
		  h=pythag(f,g);
		  w[i]=h;
		  h=1.0/h;
		  c=g*h;
		  s = -f*h;
		  for (j=1;j<=m;j++) 
		    {
		      y=a[j][nm];
		      z=a[j][i];
		      a[j][nm]=y*c+z*s;
		      a[j][i]=z*c-y*s;
		    }
		}
	    }
	  z=w[k];
	  if (l == k) 
	    { 
	      /*  Convergence. */
	      if (z < 0.0) 
		{ 
		  /* Singular value is made nonnegative.  */
		  w[k] = -z;
		  for (j=1;j<=n;j++) 
		    v[j][k] = -v[j][k];
		}
	      break;
	    }
	  if (its == 30) 
	    nrerror("no convergence in 30 svdcmp iterations");
	  x=w[l];  /* Shift from bottom 2-by-2 minor.  */
	  nm=k-1;
	  y=w[nm];
	  g=rv1[nm];
	  h=rv1[k];
	  f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	  g=pythag(f,1.0);
	  f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	  c=s=1.0;  /* Next QR transformation: */
	  for (j=l;j<=nm;j++) 
	    {
	      i=j+1;
	      g=rv1[i];
	      y=w[i];
	      h=s*g;
	      g=c*g;
	      z=pythag(f,h);
	      rv1[j]=z;
	      c=f/z;
	      s=h/z;
	      f=x*c+g*s;
	      g = g*c-x*s;
	      h=y*s;
	      y *= c;
	      for (jj=1;jj<=n;jj++) 
		{
		  x=v[jj][j];
		  z=v[jj][i];
		  v[jj][j]=x*c+z*s;
		  v[jj][i]=z*c-x*s;
		}
	      z=pythag(f,h);
	      w[j]=z;     /*  Rotation can be arbitrary if z = 0. */
	      if (z) 
		{
		  z=1.0/z;
		  c=f*z;
		  s=h*z;
		}
	      f=c*g+s*y;
	      x=c*y-s*g;
	      for (jj=1;jj<=m;jj++) 
		{
		  y=a[jj][j];
		  z=a[jj][i];
		  a[jj][j]=y*c+z*s;
		  a[jj][i]=z*c-y*s;
		}
	    }
	  rv1[l]=0.0;
	  rv1[k]=f;
	  w[k]=x;
	}
    }
  free_vector(rv1,1,n);
}







#define TOL 1.0e-5

void svdfit(double x[], double y[], double sig[], long ndata, double a[], 
	    double ma, double **u, double **v, double w[], double *chisq,
	    void (*funcs)(double, double [], long))
     /*
	Given a set of data polongs x[1..ndata],y[1..ndata] with individual 
	standard deviations sig[1..ndata], use  2 minimization to determine 
	the coeficients a[1..ma] of the fitting function y = P i ai  
	afunci(x). Here we solve the fitting equations using singular value 
	decomposition of the ndata by ma matrix, as in x 2.6. Arrays 
	u[1..ndata][1..ma], v[1..ma][1..ma], and w[1..ma] provide workspace 
	on input; on output they define the singular value decomposition, and 
	can be used to obtain the covariance matrix. The pro-gram returns 
	values for the ma fit parameters a, and  2 , chisq. The user 
	supplies a routine funcs(x,afunc,ma) that returns the ma basis 
	functions evaluated at x = x in the array afunc[1..ma].
	*/
{
  void svbksb(double **u,double w[],double **v,long m,long n,double b[],
	      double x[]);
  void svdcmp(double **a, long m, long n, double w[], double **v);
  long j,i;
  double wmax,tmp,thresh,sum,*b,*afunc;

  b=vector(1,ndata);
  afunc=vector(1,ma);

  for (i=1;i<=ndata;i++) 
    { 
      /* Accumulate coeficients of the fitting matrix.*/
      (*funcs)(x[i],afunc,ma);
      tmp=1.0/sig[i];
      for (j=1;j<=ma;j++) 
	u[i][j]=afunc[j]*tmp;
      b[i]=y[i]*tmp;
    }
  svdcmp(u,ndata,ma,w,v); /* Singular value decomposition. */
  wmax=0.0; 
  /* Edit the singular values, given TOL from the #define statement, 
     between here ... */
  for (j=1;j<=ma;j++)
    if (w[j] > wmax) 
      wmax=w[j];
  thresh=TOL*wmax;
  for (j=1;j<=ma;j++)
    if (w[j] < thresh) 
      w[j]=0.0;         /*  ...and here. */
  svbksb(u,w,v,ndata,ma,b,a);
  *chisq=0.0;             /*  Evaluate chi-square. */
  for (i=1;i<=ndata;i++) 
    {
      (*funcs)(x[i],afunc,ma);
      for (sum=0.0,j=1;j<=ma;j++) 
	sum += a[j]*afunc[j];
      *chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
    }
  free_vector(afunc,1,ma);
  free_vector(b,1,ndata);
}


#undef TOL



#define TOL 1.0e-10

void my_svdfit(double x[], double y[], long ndata, double a[], 
	    double ma, double **u, double **v, double w[],
	    void (*funcs)(double, double [], long))
     /*
	Given a set of data polongs x[1..ndata],y[1..ndata] with individual 
	standard deviations sig[1..ndata], use chi^2 minimization to determine 
	the coeficients a[1..ma] of the fitting function y = P i ai  
	afunci(x). Here we solve the fitting equations using singular value 
	decomposition of the ndata by ma matrix, as in x 2.6. Arrays 
	u[1..ndata][1..ma], v[1..ma][1..ma], and w[1..ma] provide workspace 
	on input; on output they define the singular value decomposition, and 
	can be used to obtain the covariance matrix. The program returns 
	values for the ma fit parameters a, and chi^2 , chisq. The user 
	supplies a routine funcs(x,afunc,ma) that returns the ma basis 
	functions evaluated at x = x in the array afunc[1..ma].
	*/
{
  void svbksb(double **u,double w[],double **v,long m,long n,double b[],
	      double x[]);
  void svdcmp(double **a, long m, long n, double w[], double **v);
  long j,i;
  double wmax,tmp,thresh,sum,*b,*afunc;

  b=vector(1,ndata);
  afunc=vector(1,ma);

  for (i=1;i<=ndata;i++) 
    { 
      /* Accumulate coeficients of the fitting matrix.*/
      (*funcs)(x[i],afunc,ma);
      for (j=1;j<=ma;j++) 
	u[i][j]=afunc[j];
      b[i]=y[i];
    }

  svdcmp(u,ndata,ma,w,v); /* Singular value decomposition. */
  wmax=0.0; 
  /* Edit the singular values, given TOL from the #define statement, 
     between here ... */
  for (j=1;j<=ma;j++)
    if (w[j] > wmax) 
      wmax=w[j];
  thresh=TOL*wmax;
  for (j=1;j<=ma;j++)
    if (w[j] < thresh) 
      w[j]=0.0;         /*  ...and here. */


  svbksb(u,w,v,ndata,ma,b,a);

  free_vector(afunc,1,ma);
  free_vector(b,1,ndata);
}


#undef TOL



void svdvar(double **v, long ma, double w[], double **cvm)
/*
   To evaluate the covariance matrix cvm[1..ma][1..ma] of the fit for ma 
   parameters obtained by svdfit, call this routine with matrices 
   v[1..ma][1..ma], w[1..ma] as returned from svdfit. */
{
  long k,j,i;
  double sum,*wti;
  wti=vector(1,ma);
  for (i=1;i<=ma;i++) 
    {
      wti[i]=0.0;
      if (w[i]) 
	wti[i]=1.0/(w[i]*w[i]);
    }
  for (i=1;i<=ma;i++) 
    { 
      /* Sum contributions to covariance matrix (15.4.20). */
      for (j=1;j<=i;j++) 
	{
	  for (sum=0.0,k=1;k<=ma;k++) 
	    sum += v[i][k]*v[j][k]*wti[k];
	  cvm[j][i]=cvm[i][j]=sum;
	}
    }
  free_vector(wti,1,ma);
}






void tred2(double **a, long n, double d[], double e[])
/*
Householder reduction of a real, symmetric matrix a[1..n][1..n]. On
output, a is replaced by the orthogonal matrix Q e ecting the
transformation. d[1..n] returns the diagonal elements of the
tridiagonal matrix, and e[1..n] the off-diagonal elements, with
e[1]=0. Several statements, as noted in comments, can be omitted if only
eigenvalues are to be found, in which case a contains no useful
information on output. Otherwise they are to be included. 
*/
{ 
  long l,k,j,i;
  double scale,hh,h,g,f; 
  for (i=n;i>=2;i--) 
    { 
      l=i-1; 
      h=scale=0.0; 
      if (l > 1) 
	{ 
	  for (k=1;k<=l;k++) 
	    scale += fabs(a[i][k]); 
	  if (scale == 0.0) 
	    /* Skip transformation. */
	    e[i]=a[i][l]; 
	  else 
	    { 
	      for (k=1;k<=l;k++) 
		{ 
		  a[i][k] /= scale; /* Use scaled a's for transformation. */
		  h += a[i][k]*a[i][k]; /* Form   in h. */
		} 
	      f=a[i][l]; 
	      g=(f >= 0.0 ? -sqrt(h) : sqrt(h)); 
	      e[i]=scale*g; 
	      h -= f*g; /* Now h is equation (11.2.4). */
	      a[i][l]=f-g; /* Store u in the ith row of a. */
	      f=0.0; 
	      for (j=1;j<=l;j++) 
		{ 
		/* Next statement can be omitted if eigenvectors not wanted */ 
		  a[j][i]=a[i][j]/h; 
		  /* Store u=H in ith column of a. */
		  g=0.0; /* Form an element of A   u in g. */
		  for (k=1;k<=j;k++) 
		    g += a[j][k]*a[i][k]; 
		  for (k=j+1;k<=l;k++) 
		    g += a[k][j]*a[i][k]; 
		  e[j]=g/h; 
		  /* Form element of p in temporarily unused element of e.  */
		  f += e[j]*a[i][j]; 
		} 
	      hh=f/(h+h); /* Form K, equation (11.2.11). */
	      for (j=1;j<=l;j++) 
		{ 
		/*  Form q and store in e overwriting p. */
		  f=a[i][j]; 
		  e[j]=g=e[j]-hh*f; 
		  for (k=1;k<=j;k++) /* Reduce a, equation (11.2.13). */
		    a[j][k] -= (f*e[k]+g*a[i][k]); 
		} 
	    } 
	} 
      else 
	e[i]=a[i][l]; 
      d[i]=h; 
    } 
  /* Next statement can be omitted if eigenvectors not wanted */ 
  d[1]=0.0; 
  e[1]=0.0; 
  /* Contents of this loop can be omitted if eigenvectors not wanted except 
     for statement d[i]=a[i][i]; */ 
  for (i=1;i<=n;i++) 
    { 
      /* Begin accumulation of transformation matrices. */
      l=i-1; 
      if (d[i]) 
	{ 
	  /* This block skipped when i=1. */
	  for (j=1;j<=l;j++) 
	    { 
	      g=0.0; 
	      for (k=1;k<=l;k++) /* Use u and u=H stored in a to form P Q.*/
		g += a[i][k]*a[k][j]; 
	      for (k=1;k<=l;k++) 
		a[k][j] -= g*a[k][i]; 
	    } 
	} 
      d[i]=a[i][i]; /* This statement remains. */
      a[i][i]=1.0; 
      /* Reset row and column of a to identitymatrix for next iteration. */
      for (j=1;j<=l;j++) 
	a[j][i]=a[i][j]=0.0; 
    } 
}



void tqli(double d[], double e[], long n, double **z)
/* QL algorithm with implicit shifts, to determine the eigenvalues and 
   eigenvectors of a real, symmetric, tridiagonal matrix, or of a real, 
   symmetric matrix previously reduced by tred2 x 11.2. On input, d[1..n] 
   contains the diagonal elements of the tridiagonal matrix. 
   On output, it returns the eigenvalues. The vector e[1..n] inputs the 
   subdiagonal elements of the tridiagonal matrix, with e[1] arbitrary. 
   On output e is destroyed. When ending only the eigenvalues, several lines 
   may be omitted, as noted in the comments. If the eigenvectors of a 
   tridiagonal matrix are desired, the matrix z[1..n][1..n] is input as the 
   identity matrix. If the eigenvectors of a matrix that has been reduced by 
   tred2 are required, then z is input as the matrix output by tred2. In 
   either case, the kth column of z returns the normalized eigenvector 
   corresponding to d[k]. */
{ 
  double pythag(double a, double b); 
  long m,l,iter,i,k; 
  double s,r,p,g,f,dd,c,b; 

  for (i=2;i<=n;i++) 
    e[i-1]=e[i]; /* Convenient to renumber the elements of e. */
  e[n]=0.0; 
  for (l=1;l<=n;l++) 
    { 
      iter=0; 
      do 
	{ 
	  for (m=l;m<=n-1;m++) 
	    { 
      /* Look for a single small subdiagonal element to split the matrix. */
	      dd=fabs(d[m])+fabs(d[m+1]); 
	      if ((double)(fabs(e[m])+dd) == dd) 
		break; 
	    } 
	  if (m != l) 
	    { 
	      if ((iter++) == 30) 
		nrerror("Too many iterations in tqli"); 
	      g=(d[l+1]-d[l])/(2.0*e[l]); /* Form shift. */
	      r=pythag(g,1.0); 
	      g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); /* This is dm */
	      s=c=1.0;
	      p=0.0;
	      for (i=m-1;i>=l;i--) 
		{ 
		  /* A plane rotation as in the original QL, followed by 
		     Givens rotations to restore tridiagonal form. */
		  f=s*e[i]; 
		  b=c*e[i]; 
		  e[i+1]=(r=pythag(f,g)); 
		  if (r == 0.0) 
		    { 
		      /* Recover from under ow. */
		      d[i+1] -= p; 
		      e[m]=0.0; 
		      break; 
		    } 
		  s=f/r; 
		  c=g/r; 
		  g=d[i+1]-p; 
		  r=(d[i]-g)*s+2.0*c*b; 
		  d[i+1]=g+(p=s*r); 
		  g=c*r-b; 
		  /* Next loop can be omitted if eigenvectors not wanted*/ 
		  for (k=1;k<=n;k++) 
		    { /* Form eigenvectors. */
		      f=z[k][i+1]; 
		      z[k][i+1]=s*z[k][i]+c*f; 
		      z[k][i]=c*z[k][i]-s*f; 
		    } 
		} 
	      if (r == 0.0 && i >= l) 
		continue; 
	      d[l] -= p; 
	      e[l]=g; 
	      e[m]=0.0; 
	    } 
	} 
      while (m != l); 
    } 
}







/*************************************************************************/
/****                                                                 ****/
/****                 RANDOM NUMBERS                                  ****/
/****                                                                 ****/
/*************************************************************************/


#define IA 16807 
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define NTAB 32 
#define NDIV (1+(IM-1)/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS) 

double ran1(long *idum) 
/*
"Minimal" random number generator of Park and Miller with Bays-Durham shufle 
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 
(exclusive of the endpoint values). Call with idum a negative integer to 
initialize; thereafter, do not alter idum between successive deviates in a 
sequence. RNMX should approximate the largest floating value that is less 
than 1. */
{ 
  long j; 
  long k; 
  static long iy=0; 
  static long iv[NTAB]; double temp; 
  if (*idum <= 0 || !iy) 
    { /* Initialize. */
      if (-(*idum) < 1) 
	*idum=1; /* Be sure to prevent idum = 0.*/ 
      else *idum = -(*idum); 
      for (j=NTAB+7;j>=0;j--) 
	{ /* Load the shuffle table (after 8 warm-ups). */
	  k=(*idum)/IQ; 
	  *idum=IA*(*idum-k*IQ)-IR*k; 
	  if (*idum < 0) 
	    *idum += IM; 
	  if (j < NTAB) 
	    iv[j] = *idum; 
	} 
      iy=iv[0]; 
    } 
  k=(*idum)/IQ; /* Start here when not initializing. */
  *idum=IA*(*idum-k*IQ)-IR*k;/* Compute idum=(IA*idum) % IM without overflows 
				by Schrage's method. */
  if (*idum < 0) 
    *idum += IM; 
  j=iy/NDIV; /* Will be in the range 0..NTAB-1. */
  iy=iv[j]; /* Output previously stored value and re ll the shuffle table. */
  iv[j] = *idum; 
  if ((temp=AM*iy) > RNMX) 
    return RNMX; /* Because users don't expect endpoint values. */
  else return temp; 
}





double gasdev(long *idum) 
/* Returns a normally distributed deviate with zero mean and unit variance, 
using ran1(idum) as the source of uniform deviates. */
{ 
  double ran1(long *idum); 
  static int iset=0; 
  static double gset; 
  double fac,rsq,v1,v2; 
  if (*idum < 0) 
    iset=0; /* Reinitialize. */
  if (iset == 0) 
    { /* We don't have an extra deviate handy, so */
      do 
	{ 
	  v1=2.0*ran1(idum)-1.0; /* pick two uniform numbers in the square 
				  extending from -1 to +1 in each direction, */
	  v2=2.0*ran1(idum)-1.0; 
	  rsq=v1*v1+v2*v2; /* see if they are in the unit circle, */
	} 
      while (rsq >= 1.0 || rsq == 0.0); /* and if they are not, try again. */
      fac=sqrt(-2.0*log(rsq)/rsq); /* Now make the Box-Muller transformation 
				      to get two normal deviates. Return one 
				      and save the other for next time. */
      gset=v1*fac; 
      iset=1; /* Set flag. */
      return v2*fac; 
  } 
  else 
    { /* We have an extra deviate handy, */
      iset=0; /* so unset the flag, */
      return gset; /* and return it. */
    } 
}
