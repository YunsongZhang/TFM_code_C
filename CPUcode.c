#include <math.h>               /* standard include files */
#include <gsl/gsl_rng.h>
#include "./mypsov2.h"
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include "nrutil.h"         /* Numerical Recipies */
#include "nrutil.c"         /* Numerical Recipies */

#define TOL 2.0e-4
#define TOLF 1.0e-4        /* function tolerance for conjugate gradient */
#define ITMAX 10000
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define GOLD 1.6180339887
#define GLIMIT 100.0
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define EPS 1.0e-10
#define KAPPA0 0.05

#define NMAX 10000
#define GET_PSUM \
					for (j=0;j<ndim;j++) {\
					for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];\
					psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

#define _ADIABATIC_


void                            /* function prototypes */
   SetUpProblem(),	
   initial_cell(double,double,double),
   setVisible(double Xcenter,double Ycenter,double VisibleRadius),
   deleteBond(int,int,int,int,int**,int**),
   calc_all_energy(),
   CALC_FORCE(),
   stretch(double,double,double,double),
   CellContraction(double,double,int*,int,double*,double*),
   PlotNet(char*),
   conjugate_gradient(),
   recover_net();   // for the sake of adiabatic optimization, recover the net to unstretched state
int
   frprmn(double);

void linmin(double p[],double xi[]);
void dlinmin(double p[],double xi[]);
double brent(double ax, double bx, double cx,double (*f)(double), double tol, double *xmin);
double dbrent(double ax, double bx, double cx,double (*f)(double), double (*df)(double), double tol,double *xmin);
double f1dim(double x);
double df1dim(double x);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,double *fc, double (*func)(double));
double amotry(double **p, double y[], double psum[], int ndim, double (*funk)(double [],int, void*), int ihi, double fac);
void amoeba(double **p, double y[], int ndim, double ftol, double (*funk)(double [],int, void*), int *nfunk); 


void SaveVector(double *X,double *X1,int N);
double Fourier(double theta,double *Apara,double *Bpara,int num_orders);
double deviation_func(double*,int,void*);
double deviation_fourier(double*,int,void*);
double cart2pol(double dX,double dY);

int compare_double(const void *a, const void *b)
{
	const double *da = (const double*)a;
	const double *db = (const double*)b;
	return( (*da>*db) - (*db-*da) );
}

int num_realizations,num_atoms,num_bonds,numBoundaryPts;
int Nx,Ny;
int Fourier_Orders;
double kappa,stiffness,A_lat,*A_lat0,CellRadius0,percolation,fret;
double *X,*Y;
double energy[3];
double *X_observations,*Y_observations;
double *X_initial,*Y_initial;
double Xcenter,Ycenter;
double *AX,*AY;
double *angle_known;
int **Bonds,**Junctions,*LinkedPts;
int *isfixed,*isVisible;


int main()
{
  FILE *pFileDisSet,*pFileDisPredicted,*pFileLinks,*pFileLinks0,*pFileLinks1;
  int numPt,indPt,i,j;
  double dis,xold,yold;
  double noise = 0.02; // noise of the magnitude on the netowrk deformation
  
  percolation=0.66;
  printf("percolation p=%.2f\n",percolation);
  unsigned int randseed = (unsigned)time(NULL);
  printf("Random number seed: %d\n",randseed);
  printf(" %.2f %% noise will be added to the network deformation\n",100.0*noise);
  srand48(randseed);
//  srand48(1467);
  stiffness=1.0;
  kappa=0.05;
  printf("bending coefficient kappa=%f\n",kappa);
  A_lat=1.0;
  Nx=59;
  Ny=59;
  CellRadius0=8.0*A_lat;
  double shrink_percentage; 
  Xcenter=(double)(Nx)*0.5*A_lat;
  Ycenter=(double)(Ny)*0.5*A_lat;
  shrink_percentage=0.3;
  double shrink_dis=CellRadius0*shrink_percentage;

  SetUpProblem();
#ifdef CLOSEST_BLIND
  for(i=0;i<num_atoms;i++)
  {
	  isVisible[i] = 1;
  }
#endif
  PlotNet("./data/BlankNet.txt");
  initial_cell(Xcenter,Ycenter,CellRadius0);
  setVisible(Xcenter,Ycenter,100.0); // all visible
#ifdef CLOSE_NEIGHBOR_ONLY
  setVisible(Xcenter,Ycenter,15.0);
#endif
  PlotNet("./data/unrelaxedNet.txt");
  X_initial=(double*)malloc(num_atoms*sizeof(double));
  Y_initial=(double*)malloc(num_atoms*sizeof(double));
  SaveVector(X,X_initial,num_atoms);
  SaveVector(Y,Y_initial,num_atoms);



  pFileLinks0=fopen("./data/linkedPtsInitial.txt","w");
  for(numPt = 0; numPt<numBoundaryPts; numPt++)
  {
   fprintf(pFileLinks0,"%.3f %.3f\n",X[LinkedPts[numPt]],Y[LinkedPts[numPt]]);
  }
  fclose(pFileLinks0);


  printf("%d points are attached to the cell\n",numBoundaryPts);
/* This module sets the shift for all attached points*/
  double *contract_dis = (double*)malloc(numBoundaryPts*sizeof(double));
  double *contract_angle = (double*)malloc(numBoundaryPts*sizeof(double));
  
  Fourier_Orders = 4;

  double *Acoef = (double*)malloc((Fourier_Orders+1)*sizeof(double));
  double *Bcoef = (double*)malloc((Fourier_Orders+1)*sizeof(double));
     Acoef[0] = 2.0*drand48(); //shrink_percentage*CellRadius0;
     Bcoef[0] = 0.0;

     for(i = 1; i<= Fourier_Orders; i++)
     {
	     Acoef[i] = 1.0 * drand48();
	     Bcoef[i] = 1.0 * drand48();
     }

		  

  for(numPt=0;numPt<numBoundaryPts;numPt++)
  {
	  indPt=LinkedPts[numPt];
	  xold=X[indPt]; yold=Y[indPt];
	  //dis=sqrt((xold-Xcenter)*(xold-Xcenter)+(yold-Ycenter)*(yold-Ycenter)+DBL_EPSILON);

//	  contract_dis[numPt] = shrink_percentage*CellRadius0;
	  contract_angle[numPt] = cart2pol(xold-Xcenter,yold-Ycenter);
	  contract_dis[numPt] = Fourier(contract_angle[numPt],Acoef,Bcoef,Fourier_Orders);
  }

  pFileDisSet = fopen("./data/contraction_set.txt","w");
  if (pFileDisSet )
  {
	  for( numPt = 0; numPt<numBoundaryPts; numPt++)
	  {
		  fprintf(pFileDisSet,"%.5f %.5f\n",contract_dis[numPt],contract_angle[numPt]);
	  }
	  fclose(pFileDisSet);
  }
  else
  {
	  printf("file contraction_set open failed!\n");
	  exit(-1);
  }

CellContraction(Xcenter,Ycenter,LinkedPts,numBoundaryPts,
		  contract_dis,contract_angle);


  pFileLinks=fopen("./data/linkedPts.txt","w");
  for(numPt = 0; numPt<numBoundaryPts; numPt++)
  {
   fprintf(pFileLinks,"%.3f %.3f\n",X[LinkedPts[numPt]],Y[LinkedPts[numPt]]);
  }
  fclose(pFileLinks);


 // free(contract_dis);
 // stretch(Xcenter,Ycenter,CellRadius0,shrink_percentage);
  conjugate_gradient();
  PlotNet("./data/RelaxedNet.txt");

  X_observations=(double*)malloc(num_atoms*sizeof(double));
  Y_observations=(double*)malloc(num_atoms*sizeof(double));
  for( i = 0 ; i < num_atoms ; i++)
  {
	  X_observations[i] *= (1.0+noise*2.0*(drand48()-0.5));
	  Y_observations[i] *= (1.0+noise*2.0*(drand48()-0.5));
  }
  SaveVector(X,X_observations,num_atoms);
  SaveVector(Y,Y_observations,num_atoms);
  angle_known=(double*)malloc(numBoundaryPts*sizeof(double));
  SaveVector(contract_angle,angle_known,numBoundaryPts);


  /* Now we start the PSO portion of our CPU code! */
  
  pso_obj_fun_t obj_fun = deviation_fourier;
  pso_settings_t settings;
  pso_set_default_settings(&settings);
  
  // set PSO settings manually
      settings.dim=2*Fourier_Orders+1; //numBoundaryPts;
      settings.size=100;// pso_calc_swarm_size(settings.dim);
      printf(" %d particles in PSO algorithm!\n",settings.size);
      settings.print_every=1;
      settings.nhood_strategy = PSO_NHOOD_RANDOM;
      settings.nhood_size = 10;
     // settings.w_strategy = PSO_W_LIN_DEC;
      settings.x_lo=-2.0;
      settings.x_hi=2.0;
      settings.goal=2.0e-2;
      settings.steps= 1000;
      settings.w_min=1e-4;
      settings.x_min=(double*)malloc(settings.dim*sizeof(double));
      settings.x_max=(double*)malloc(settings.dim*sizeof(double));
      
      int num;
      for(num = 0; num<settings.dim; num++)
      {
	settings.x_min[num]=-3.5;
	settings.x_max[num]=3.5;
      }

   
  pso_result_t solution;
  
  solution.gbest=malloc(settings.dim*sizeof(double));
  
    struct timespec start,finish;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC,&start);



    pso_solve(obj_fun,NULL,&solution,&settings);

    clock_gettime(CLOCK_MONOTONIC,&finish);
    elapsed = (finish.tv_sec-start.tv_sec);
    elapsed+= (finish.tv_nsec-start.tv_nsec)*1.0e-9;
    printf(" Particle Swarm Optimization time=%lf sec\n",elapsed);


    printf("A=%f,Apredicted=%f\n",Acoef[0],*solution.gbest);
    for(i = 1 ;i<(settings.dim+1)/2 ; i++)
    {
	    printf("A=%f,Apredicted=%f ",Acoef[i],*(solution.gbest+i)); 
	    printf("B=%f,Bpredicted=%f\n",Bcoef[i],*(solution.gbest+i+(settings.dim-1)/2));
    }
    
    
    double *contract_angle_p = (double*)malloc(numBoundaryPts*sizeof(double));
    double *contract_dis_p   = (double*)malloc(numBoundaryPts*sizeof(double));
    double *A_latest = (double*)malloc((settings.dim+1)/2*sizeof(double));
    double *B_latest = (double*)malloc((settings.dim+1)/2*sizeof(double));
    A_latest[0] = *(solution.gbest);  
    B_latest[0] = 0.0;
    for ( i = 1; i<(settings.dim+1)/2 ; i++)
    {
	    A_latest[i] = *(solution.gbest+i);
	    B_latest[i] = *(solution.gbest+i+(settings.dim-1)/2);
    }

    for( numPt =0; numPt<numBoundaryPts; numPt++)
    {
	    contract_angle_p[numPt] = contract_angle[numPt];
	    contract_dis_p[numPt]   = Fourier(contract_angle[numPt],A_latest,B_latest,(settings.dim-1)/2);
    }
	    
  pFileDisPredicted = fopen("./data/contraction_predicted.txt","w");
  if (pFileDisPredicted )
  {
	  for( numPt = 0; numPt<numBoundaryPts; numPt++)
	  {
		  fprintf(pFileDisPredicted,"%.5f %.5f\n",contract_dis_p[numPt],contract_angle_p[numPt]);
	  }
	  fclose(pFileDisPredicted);
  }
  else
  {
	  printf("file contraction_set open failed!\n");
	  exit(-1);
  }

  /* This module outputs the predicted network*/
  CellContraction(Xcenter,Ycenter,LinkedPts,numBoundaryPts,
		  contract_dis_p,contract_angle_p);
  conjugate_gradient();
  PlotNet("./data/predictedNet.txt");

  double tmp_x,tmp_y;
  pFileLinks1=fopen("./data/linkedPtsPredicted.txt","w");
  for(numPt = 0; numPt<numBoundaryPts; numPt++)
  {
    tmp_x = X_initial[LinkedPts[numPt]];
    tmp_y = Y_initial[LinkedPts[numPt]];
    tmp_x-= contract_dis_p[numPt]*cos(contract_angle_p[numPt]);
    tmp_y-= contract_dis_p[numPt]*sin(contract_angle_p[numPt]);
   fprintf(pFileLinks1,"%.3f %.3f\n",tmp_x,tmp_y);
  }
  fclose(pFileLinks1);
  
#ifdef _ADIABATIC_
  if ( solution.error > settings.goal*25.0)
  {
	  printf(" PSO fails, no need to go further!\n");
	  exit(-1);
  }
  printf("\n \n \n \n \n ");
  printf("From now on, I adiabatically decrease the bending modulus of the material!\n");
  int num_AdiabaticSteps = 49,numStep;
  double delta_kappa = 0.001;
  FILE *pFileRelax,*pFilePredicct;
  char filename[50];
  int namelength;
  double *startpoint = (double*)malloc((2*Fourier_Orders+1)*sizeof(double));
  double *values_simplex = (double*)malloc((2*Fourier_Orders+2)*sizeof(double));
  double **pts_simplex = (double**)malloc((2*Fourier_Orders+2)*sizeof(double*));
  int *nfunk = (int*)malloc(sizeof(int));
  double *predicted_dis = (double*)malloc(numBoundaryPts*sizeof(double));
  double *deviation_percentage = (double*)malloc(numBoundaryPts*sizeof(double));
  double *abs_deviation_percentage = (double*)malloc(numBoundaryPts*sizeof(double));

  for( i = 0 ; i < 2*Fourier_Orders+2 ; i++)
  {
	  pts_simplex[i] = (double*)malloc((2*Fourier_Orders+1)*sizeof(double));
  }



  for( numStep = 0 ; numStep <= num_AdiabaticSteps ; numStep++)
  {
	  kappa -= delta_kappa * numStep; 
	  printf("#Adiabatic Decrease = %d , bending modulus = %.4f\n",numStep, kappa);


	  // recover the net to the unstretched state
	  SaveVector(X_initial,X,num_atoms);
	  SaveVector(Y_initial,Y,num_atoms);
	  // redo the polar to cartesian transform to double check
           for(numPt=0;numPt<numBoundaryPts;numPt++)
           {
                   indPt=LinkedPts[numPt];
                   xold=X[indPt]; yold=Y[indPt];
                   contract_angle[numPt] = cart2pol(xold-Xcenter,yold-Ycenter);
                   contract_dis[numPt] = Fourier(contract_angle[numPt],Acoef,Bcoef,Fourier_Orders);
           }


           CellContraction(Xcenter,Ycenter,LinkedPts,numBoundaryPts,
            		  contract_dis,contract_angle);
	   conjugate_gradient();
	   namelength = sprintf(filename,"./data/relax_%d.txt",numPt);
	   PlotNet(filename);
           for( i = 0 ; i < num_atoms ; i++)
           {
                   X_observations[i] *= (1.0+noise*2.0*(drand48()-0.5));
                   Y_observations[i] *= (1.0+noise*2.0*(drand48()-0.5));
           }
	   SaveVector(X,X_observations,num_atoms);
	   SaveVector(Y,Y_observations,num_atoms);

	   // Use the pso predicted values A_latest and B_latest as initial guess to do the outer loop optimization
	   // Here due to a lack of derivatives 
	   // (also one shall not numerically calculate the derivatives 
	   // due to the potential highly discontinuous nature of the problem ) 
	   // Here a Nelder-Mead simplex algorithm is applied

	   startpoint[0] = A_latest[0];
	   for( i = 1 ; i <= Fourier_Orders ; i++)
	   {
		   startpoint[i] = A_latest[i];
		   startpoint[i+Fourier_Orders] = B_latest[i];
	   }

	   for( i = 0 ; i < 2*Fourier_Orders+2 ; i++)
	   {
		   for( j = 0 ; j < 2*Fourier_Orders+1 ; j++)
		   {
			   pts_simplex[i][j] = startpoint[j];
		   }
	   }

	   for( i = 0 ; i < 2*Fourier_Orders+1 ; i++)
	   {
		   pts_simplex[i][i] +=0.05;
	   }

	   for( i = 0 ; i<2*Fourier_Orders+2; i++)
	   {
		   values_simplex[i] = deviation_fourier(pts_simplex[i],2*Fourier_Orders+1,NULL); 
	   }
	   
            amoeba(pts_simplex,values_simplex,2*Fourier_Orders+1,TOLF,deviation_fourier,nfunk);


	    printf(" adiabatic change # %d finished!\n",numStep);
	    printf(" optimization value = %e \n",values_simplex[0]);
	    printf(" Object functions are evaluated %d times\n",*nfunk);
	    printf(" Compare results:\n");
	    printf(" A0 = %f  , A0_pred=%f\n",Acoef[0],pts_simplex[0][0]);
	    for( i = 1 ; i <= Fourier_Orders; i++)
	    {
		    printf(" A%d = %f  , A%d_pred=%f   ",i,Acoef[i],i,pts_simplex[0][i]);
		    printf(" B%d = %f  , B%d_pred=%f  \n ",i,Bcoef[i],i,pts_simplex[0][i+Fourier_Orders]);
	    }
	    printf("\n \n \n \n");

	    if( values_simplex[0] < 1.0e-6)
	    {
		    printf("Victory belongs to us!\n");
		    exit(0);
	    }


	    // update A_latest and B_latest
		    printf("Update A_latest and B_latest with downhill simplex optimization\n");
		    A_latest[0] = pts_simplex[0][0];
		    B_latest[0] = 0.0;
		    for( i = 1; i<= Fourier_Orders ; i++)
		    {
			    A_latest[i] = pts_simplex[0][i];
			    B_latest[i] = pts_simplex[0][i+Fourier_Orders];
		    }
             // update A_latest and B_latest

	     for( i = 0; i < numBoundaryPts; i++)
	     {
		     predicted_dis[i] = Fourier(contract_angle[i],A_latest,B_latest,Fourier_Orders);
		     deviation_percentage[i] = (predicted_dis[i]-contract_dis[i])/predicted_dis[i];
		     printf(" Boundary Point # %d : deviation=%.2f%% \n",i,deviation_percentage[i]*100.0);
		     abs_deviation_percentage[i] = fabs(deviation_percentage[i]);

	     }

	     qsort(abs_deviation_percentage,numBoundaryPts,sizeof(double),compare_double);

	     printf("\n\n prediction results evaluation:\n");
	     printf(" smallest deviation %.2f %%, largest deviation %.2f %%, median deviation %.2f %%\n",
			     abs_deviation_percentage[0]*100.0,abs_deviation_percentage[numBoundaryPts-1]*100.0,
			     abs_deviation_percentage[numBoundaryPts/2-1]*100.0 );
	     printf("\n\n");





	   


  }
		    


     // free memory occupation for ADIABATIC situations
         free(startpoint);
	 free(values_simplex);

	 for( i =0 ; i<2*Fourier_Orders+2 ; i++) { free(pts_simplex[i]); }
	 free(pts_simplex);
	 free(nfunk);
	 free(deviation_percentage);
	 free(abs_deviation_percentage);
	 free(predicted_dis);


#endif



  // Before finishing the main function
  // Clear all memory occupation

   free(A_latest);
   free(B_latest);
   free(contract_angle_p);
   free(contract_dis_p);


   free(Acoef);
   free(Bcoef);

   free(solution.gbest);
   free(settings.x_min);
   free(settings.x_max);
   
 
   free(X_initial);
   free(Y_initial);
   free(X_observations);
   free(Y_observations);
   free(A_lat0);
   free(isfixed);
   free(X);
   free(Y);
   free(AX);
   free(AY);
   free(contract_dis);
   free(contract_angle);
   free(LinkedPts);
 
   int pnum;
   for(pnum=0;pnum<num_bonds;pnum++)
   {
 	 free(Bonds[pnum]);
   }
   for(pnum=0;pnum<num_atoms;pnum++)
   {
 	 free(Junctions[pnum]);
   }
   free(Bonds);
   free(Junctions);
}


/* All functions definition */
void SetUpProblem()
{
	num_atoms=(Nx+1)*(Ny+1);
	num_bonds=(Ny+1)*Nx+(Nx+1)*Ny+Nx*Ny;

        AX=(double*)malloc(num_atoms*sizeof(double));
        AY=(double*)malloc(num_atoms*sizeof(double));
        X=(double*)malloc(num_atoms*sizeof(double));
        Y=(double*)malloc(num_atoms*sizeof(double));
	LinkedPts=(int*)malloc(num_atoms*sizeof(int));
	isfixed=(int*)malloc(num_atoms*sizeof(int));
	isVisible=(int*)malloc(num_atoms*sizeof(int));

	int i;
	Bonds=(int**)malloc(num_bonds*sizeof(int*));
	   for(i=0;i<num_bonds;i++)  { Bonds[i] = (int*) malloc(2*sizeof(int)); }

	Junctions=(int**)malloc(num_atoms*sizeof(int*));
	   for(i=0;i<num_atoms;i++)  { Junctions[i]=(int*) malloc(6*sizeof(int)); }

	/* Start constructing networks */
	/*******5*****0**********/
	/*****4**Node**1********/
	/******3*****2*********/
	int row,column,index;
	for(row=0;row<=Nx;row+=2)
	{
		for(column=0;column<=Ny;column++)
		{
			index=row*(Nx+1)+column;
			Y[index]=sqrt(3.0)*A_lat/2.0*row;
			X[index]=column*A_lat;
			  Junctions[index][0]=index+Nx+1;
			  Junctions[index][1]=index+1;
			  Junctions[index][2]=index-Nx-1;
			  Junctions[index][3]=index-Nx-2;
			  Junctions[index][4]=index-1;
			  Junctions[index][5]=index+Nx;
			  
			  if (row ==0) 
			  { Junctions[index][2]=-1;
		            Junctions[index][3]=-1; }
			  if (row ==Ny)
			  { Junctions[index][0]=-1;
		            Junctions[index][5]=-1; }
			  if (column==0)
			  { Junctions[index][3]=-1;
		            Junctions[index][4]=-1; 
			    Junctions[index][5]=-1; }
			  if (column==Nx)
			  { Junctions[index][1]=-1; }
		}
	}


	for(row=1;row<=Nx;row+=2)
	{
		for(column=0;column<=Ny;column++)
		{
			index=row*(Nx+1)+column;
			Y[index]=sqrt(3.0)*A_lat/2.0*row;
			X[index]=column*A_lat+0.5*A_lat;
			  Junctions[index][0]=index+Nx+2;
			  Junctions[index][1]=index+1;
			  Junctions[index][2]=index-Nx;
			  Junctions[index][3]=index-Nx-1;
			  Junctions[index][4]=index-1;
			  Junctions[index][5]=index+Nx+1;

			  if (row ==0) 
			  { Junctions[index][2]=-1;
		            Junctions[index][3]=-1; }
			  if (row ==Ny)
			  { Junctions[index][0]=-1;
		            Junctions[index][5]=-1; }
			  if (column==Nx)
			  { Junctions[index][0]=-1;
		            Junctions[index][2]=-1;
		            Junctions[index][1]=-1; }
			  if (column==0)
			  { Junctions[index][4]=-1; }
		}
	} 
	/*End of construction networks */
	

	int count_bond=0;
	int J,temp;

	A_lat0=(double*)malloc(num_bonds*sizeof(double));
	int ind[3]={5,0,1};

	for(row=0;row<=Nx;row++)
	{
		for(column=0;column<=Ny;column++)
		{
			index=row*(Nx+1)+column;
			if(row==0||row==Ny||column==0 || column==Nx)
			{
				isfixed[index]=1;
			}
			else
			{
				isfixed[index]=0;
			}

			for(J=0;J<3;J++)
			{
				if(Junctions[index][ind[J]]>0)
				{
					
					if( count_bond>=num_bonds ) 
					{
						printf("Trying to get nonexistant bond!\n");
						exit(-1);
					}
					temp=Junctions[index][ind[J]];
					Bonds[count_bond][0]=(temp>index)?index:temp;
					Bonds[count_bond][1]=(temp<index)?index:temp;

					A_lat0[count_bond]=A_lat;
					count_bond++;
				}
			}
		}
	} /* end constructing all bonds*/

	for(J=0;J<3;J++)
	{
		energy[J]=0.0;
	}

	calc_all_energy();
					
	printf("Problem Initialized!\n");
	printf("%d atoms,%d bonds in total\n",num_atoms,num_bonds);
	printf("Total Energy:%.2f, Stretching Energy:%.2f, Bending Energy:%.2f \n",energy[2],energy[0],energy[1]);




	/* cut a certain portion of bonds in the network */
	int num_cuts,num_removed,p,a,b;

	num_cuts=(int)((double)(num_bonds)*(1-percolation));
        num_removed=0;	

      while(num_removed < num_cuts)
      {
/*        p=min((int)(drand48()*(double)(num_bonds)+1.0),num_bonds-1);*/
	p=(int)(drand48()*(double)num_bonds);
        a=Bonds[p][0];
        b=Bonds[p][1];
          if((a>0)&&(b>0))
          {
	    num_removed++;
	    Bonds[p][0]=-1;
	    Bonds[p][1]=-1;
               for(J=0; J<6; J++)
               {
                 if(Junctions[a][J]==b) { Junctions[a][J]=-1; }
                 if(Junctions[b][J]==a) { Junctions[b][J]=-1; }
               } /* end of junction adjustment */
          }  /* end of one bond elimination */
       }  /* End of all bond elimination */

}  /* End of SetUpProblem function */

void calc_all_energy()
{
	int i,a,b,k;
	double r,R1,R2,x0,x1,x2,y0,y1,y2,cos_theta,theta;
	int vec1[3]={0,1,2};
	int vec2[3]={3,4,5};

	for(i=0;i<=3;i++) { energy[i]=0.0; }

	/* stretching energy */
	for(i=0;i<num_bonds;i++)
	{
		if( (Bonds[i][0]>-1)&&(Bonds[i][1]>-1) )
		{
			a=Bonds[i][0];
			b=Bonds[i][1];
			x1=X[a]; y1=Y[a];
			x2=X[b]; y2=Y[b];
			r=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
			energy[0]+=0.5*stiffness*(r-A_lat0[i])*(r-A_lat0[i]);

			if ( isnan(energy[0])) 
			{
				printf("At this moment,stretching energy is NAN!\n");
				printf("%d (%.2f,%.2f) and %d (%.2f,%.2f)\n",a,x1,y1,b,x2,y2);
				exit(-1);
			}
		}
	} /*end of calculating stretching energy */

	/* bending energy*/
	if ( kappa>1.0e-5)
	{
		for(i=0;i<num_atoms;i++)
		{
			for(k=0;k<3;k++)
			{
				if ( (Junctions[i][vec1[k]]>0) && (Junctions[i][vec2[k]]>0) )
				{
					a=Junctions[i][vec1[k]];
					b=Junctions[i][vec2[k]];
					x0=X[i]; y0=Y[i];
					x1=X[a]; y1=Y[a];
					x2=X[b]; y2=Y[b];
					R1=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
					R2=sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
					cos_theta=((x1-x0)*(x2-x0)+(y1-y0)*(y2-y0))/(R1*R2);
                                        if( isnan(cos_theta) )
                                        {
                                            printf("cos_theta is NAN!\n");
						printf("positions:(%.2f,%.2f),(%.2f,%.2f),(%.2f,%.2f)\n",x1,y1,x0,y0,x2,y2);
						printf("R1=%f,R2=%f,cos=%f\n",R1,R2,cos_theta);
						printf("theta:  %f \n",acos((cos_theta)));
						exit(-1);
					}

					if (cos_theta<-1.0) {cos_theta=-1.0;}
					if (cos_theta>1.0)  {cos_theta=1.0; }
					theta=fabs(M_PI-acos((cos_theta)));
					energy[1]+=0.5*kappa*theta*theta;

					if( isnan(energy[1]) )
					{
						printf("At this moment,bending energy is NAN!\n");
						printf("positions:(%.2f,%.2f),(%.2f,%.2f),(%.2f,%.2f)\n",x1,y1,x0,y0,x2,y2);
						printf("newly added bending energy %f, %f\n",0.5*kappa*theta*theta,theta);
						printf("R1=%f,R2=%f,cos=%f\n",R1,R2,cos_theta);
						printf("theta:  %f \n",acos((cos_theta)));

						exit(-1);
					}

				}
			}
		}
	} /* end of calculating bending energy */

	energy[2]=energy[0]+energy[1];

} /* end of function calc_all_energy() */

void initial_cell(double Xcenter,double Ycenter,
		  double Radius0)
{
	int a,b,i,j;
	double r1,r2,a1,a2,a3,s1,s2,s,x1,x2,y1,y2;
	
	numBoundaryPts=0;

	for(i=0;i<num_bonds;i++)
	{ 
		a=Bonds[i][0];
		b=Bonds[i][1];
		if ( (a>-1) && (b>-1) )
		{
                     r1=sqrt((X[a]-Xcenter)*(X[a]-Xcenter)+(Y[a]-Ycenter)*(Y[a]-Ycenter));
                     r2=sqrt((X[b]-Xcenter)*(X[b]-Xcenter)+(Y[b]-Ycenter)*(Y[b]-Ycenter));

		     if( (r1<=Radius0) && (r2<=Radius0) )
		     {
		   	  Bonds[i][0]=-1; Bonds[i][1]=-1; 
		   	  for(j=0;j<6;j++)
		   	  {
		   		  if( Junctions[a][j]==b) { Junctions[a][j]=-1; }
		   		  if( Junctions[b][j]==a) { Junctions[b][j]=-1; }
		   	  }
		     }

		     else if( (r1<Radius0 ) && (r2>Radius0) )
		     {
			  LinkedPts[numBoundaryPts]=a;
			  numBoundaryPts++;
		   	  x1=X[a]; y1=Y[a];
		   	  x2=X[b]; y2=Y[b];
		   	  a1=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
		   	  a2=2*((x1-Xcenter)*(x2-x1)+(y1-Ycenter)*(y2-y1));
		   	  a3=r1*r1-Radius0*Radius0;
		   	  s1=(-a2+sqrt(a2*a2-4*a1*a3))/(2*a1);
		   	  s2=(-a2-sqrt(a2*a2-4*a1*a3))/(2*a1);
		   	  s=(s1<=1 && s1>=0)?s1:s2;
		   	  X[a]=x1+s*(x2-x1);
		   	  Y[a]=y1+s*(y2-y1);
		   	  A_lat0[i]=(1-s)*sqrt(a1);
			  isfixed[a]=1;
#ifdef CLOSEST_BLIND
			  isVisible[b]=0;
#endif

			  for(j=0;j<6;j++)
			  {
				  if ( (Junctions[a][j]!=b) && ( Junctions[a][j]>-1) )
				  {
				 	  deleteBond(a,Junctions[a][j],num_bonds,num_atoms,Bonds,Junctions);
				  }
			  }

		   	  if( (s<0) || (s>1))  { printf("Error1, no crossing,s=%.10f\n",s); } 
/*				  printf("Node moved from (%.3f,%.3f) to (%.3f,%.3f), distance (go up)from %.3f to %.3f s=%.5f\n",x1,y1,X[a],Y[a],r1,sqrt((X[a]-Xcenter)*(X[a]-Xcenter)+(Y[a]-Ycenter)*(Y[a]-Ycenter)),s);*/
		     }
		     else if ( (r1>Radius0) && (r2<Radius0) )
		     {
			     LinkedPts[numBoundaryPts]=b;
			     numBoundaryPts++;
		             x1=X[a]; y1=Y[a];
		             x2=X[b]; y2=Y[b];
		             a1=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
		             a2=2*((x1-Xcenter)*(x2-x1)+(y1-Ycenter)*(y2-y1));
		             a3=r1*r1-Radius0*Radius0;
		             s1=(-a2+sqrt(a2*a2-4*a1*a3))/(2*a1);
		             s2=(-a2-sqrt(a2*a2-4*a1*a3))/(2*a1);
		             s=(s1<=1 && s1>=0)?s1:s2;
		             X[b]=x1+s*(x2-x1);
		             Y[b]=y1+s*(y2-y1);
		             A_lat0[i]=s*sqrt(a1);
			     isfixed[b]=1;
#ifdef CLOSEST_BLIND
			  isVisible[a]=0;
#endif
			     

			  for(j=0;j<6;j++)
			  {
				  if ( (Junctions[b][j]!=a) && (Junctions[b][j]>-1) )
				  {
				 	  deleteBond(b,Junctions[b][j],num_bonds,num_atoms,Bonds,Junctions);
				  }
			  }
		             if( (s<0) | (s>1) ) { printf("Error2, no crossing,s=%.10f\n Notice:",s); }
/*				  printf("Node moved from (%.3f,%.3f) to (%.3f,%.3f), distance (go down )from %.3f to %.3f s=%.10f \n",x1,y1,X[b],Y[b],r1,sqrt((X[b]-Xcenter)*(X[b]-Xcenter)+(Y[b]-Ycenter)*(Y[b]-Ycenter)),s);*/
		     }
		  }
              }


	/* This module will delete all dangling bonds*/
  	int numLinks,link_num;
  	for(i=0;i<num_atoms;i++)
  	{
  		numLinks=0;
  	   if( isfixed[i]==0)
  	   {
  		for(j=0;j<6;j++)
  		{
  			if( Junctions[i][j]>-1)
  			{
  				numLinks++;
  				link_num=Junctions[i][j];
  			}
  		}
  		if (numLinks==1 )
  		{
  			deleteBond(i,link_num,num_bonds,num_atoms,Bonds,Junctions);
  	//		printf("dangling bond deleted: %d,(%f,%f) and %d, (%f,%f)\n",i,X[i],Y[i],link_num,X[link_num],Y[link_num]);
  		}
  	   } // end isfixed....
  	}
	/* end of dangling bond remove */

        calc_all_energy();
	printf("\n\nCell Inserted!\n");
	printf("Center (%f,%f),Radius=%f\n",Xcenter,Ycenter,Radius0);
	printf("Total Energy:%.2f, Stretching Energy:%.2f, Bending Energy:%.2f \n\n",energy[2],energy[0],energy[1]);
		  
} /* end of function void initial_cell() */


void stretch(double Xcenter,double Ycenter,double Radius0,double shrink_percentage)
{
	double r,xnew,ynew,xold,yold;
	int count;

	for(count=0;count<num_atoms;count++)
	{
		r=sqrt((X[count]-Xcenter)*(X[count]-Xcenter)+(Y[count]-Ycenter)*(Y[count]-Ycenter));
		if (fabs(r-Radius0)<=1e-10)
		{
			xold=X[count];yold=Y[count];
			xnew=Xcenter+(1-shrink_percentage)*(X[count]-Xcenter);
			ynew=Ycenter+(1-shrink_percentage)*(Y[count]-Ycenter);
			X[count]=xnew;
			Y[count]=ynew;
			isfixed[count]=1;
			printf("stretch: %d (%.3f,%.3f) to (%.3f,%.3f)\n",count,xold,yold,xnew,ynew);
		}
	}


	calc_all_energy();
	printf("Network Stretched by %.0f percent!\n Stretching Energy=%f, Bending Energy=%f,Total Energy=%f\n\n",100*shrink_percentage,energy[0],energy[1],energy[2]);

} /* end of function stretch(Xcenter,Ycenter,R0,shrink_percentage) */

void CellContraction(double Xcenter,double Ycenter,int* LinkedPts,int numBoundaryPts,double* dis_shift,double* angle_shift)
{
	int i,a;
	double xold,yold;
	for(i=0;i<num_atoms;i++)
	{
		X[i]=X_initial[i];
		Y[i]=Y_initial[i];
	}
	for(i=0;i<numBoundaryPts;i++)
	{
		a=LinkedPts[i];
		xold=X[a]; yold=Y[a];
		X[a]-=dis_shift[i]*cos(angle_shift[i]);
		Y[a]-=dis_shift[i]*sin(angle_shift[i]);
		isfixed[a]=1;
		//	printf("contraction: %d (%.3f,%.3f) to (%.3f,%.3f),angle=%f\n",a,xold,yold,X[a],Y[a],angle_shift[i]);
	}
	//calc_all_energy();
	//printf("Network Stretched  \n Stretching Energy=%f, Bending Energy=%f,Total Energy=%f\n\n",energy[0],energy[1],energy[2]);
}/* end of function CellContraction */



void deleteBond(int a,int b,int num_bonds,int num_atoms,int** Bonds,int** Junctions)
{
	int count;

	if ( a<0 || b<0)
	{
		printf("Error: One of the bond nodes does NOT exist!\n");
		return;
	}
	else
	{
		for(count=0;count<num_bonds;count++)
		{
			if (  (Bonds[count][0]==a && Bonds[count][1]==b) || ( Bonds[count][0]==b && Bonds[count][1]==a) )
			{
				Bonds[count][0]=-1;
				Bonds[count][1]=-1;
			}
		}
		for(count=0;count<6;count++)
		{
			if( Junctions[a][count]==b ) {Junctions[a][count]=-1; }
			if( Junctions[b][count]==a ) {Junctions[b][count]=-1; }
		}
	}

}
		

void conjugate_gradient()
{
	int its;
	double ftolf;
	double elapsed;
	struct timespec start,finish;

	ftolf=TOLF;

		clock_gettime(CLOCK_MONOTONIC,&start);
		its=frprmn(ftolf);
		//printf("its=%d\n",its);
		clock_gettime(CLOCK_MONOTONIC,&finish);
		elapsed = (finish.tv_sec-start.tv_sec);
		elapsed+= (finish.tv_nsec-start.tv_nsec)*1.0e-9;
                //printf("conjugate-gradient time = %lf sec\n",elapsed);
                //printf("Using conjugate-gradient algorithm, its=%d, energy=%lf\n",its,energy[2]);

} /* end of function conjugate_gradient() */
	
int frprmn(double ftol)
{
	int j,its,i;
	double gg,gam,fp,dgg;
	double *g,*h,*xi,*p;

	g=dvector(1,2*num_atoms);
	h=dvector(1,2*num_atoms);
	xi=dvector(1,2*num_atoms);
	p=dvector(1,2*num_atoms);
	calc_all_energy();
	fp=energy[2];
	CALC_FORCE();
	for(i=1;i<=num_atoms;i++)
	{
	  if(isfixed[i-1]==0)
 	  {
	    p[2*i-1]=X[i-1];
	    p[2*i]=Y[i-1];
	    xi[2*i-1]=-AX[i-1];
	    xi[2*i]=-AY[i-1];
	  }
	}
	for (j=1;j<=2*num_atoms;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {  /* ITMAX*/

		  linmin(p,xi);

		if (2.0*fabs(fret-fp) <= ftol*(fabs(fret)+fabs(fp)+EPS)) {
			free_dvector(xi,1,2*num_atoms);
			free_dvector(h,1,2*num_atoms);
			free_dvector(g,1,2*num_atoms);
			free_dvector(p,1,2*num_atoms);
			return its;
		}     		       
		calc_all_energy();
		fp=energy[2];
		CALC_FORCE();
	        for(i=1;i<=num_atoms;i++)
	        {
	          xi[2*i-1]=-AX[i-1];
	          xi[2*i]=-AY[i-1];
	        }
		dgg=gg=0.0;
		for (j=1;j<=num_atoms;j++) {
			if(isfixed[j-1]==0)
			{
			  gg += g[2*j-1]*g[2*j-1];
			  dgg += (xi[2*j-1]+g[2*j-1])*xi[2*j-1];
			  gg += g[2*j]*g[2*j];
			  dgg += (xi[2*j]+g[2*j])*xi[2*j];
			}
		}
		if (gg == 0.0) {
			free_dvector(xi,1,2*num_atoms);
			free_dvector(h,1,2*num_atoms);
			free_dvector(g,1,2*num_atoms);
			free_dvector(p,1,2*num_atoms);
			return its;
		}
		gam=dgg/gg;
		for (j=1;j<=num_atoms;j++) {
			if(isfixed[j-1]==0)
			{
			  g[2*j-1] = -xi[2*j-1];
			  xi[2*j-1]=h[2*j-1]=g[2*j-1]+gam*h[2*j-1];
			  g[2*j] = -xi[2*j];
			  xi[2*j]=h[2*j]=g[2*j]+gam*h[2*j];
			}
		}
	//	printf("iterNo=%d,Es=%.3f,Eb=%.3f,Et=%.3f\n",its,energy[0],energy[1],energy[2]);
	}

			free_dvector(xi,1,2*num_atoms);
			free_dvector(h,1,2*num_atoms);
			free_dvector(g,1,2*num_atoms);
			free_dvector(p,1,2*num_atoms);
	nrerror("Too many iterations in frprmn");
  return its;
}



int ncom;
double *pcom,*xicom;

void linmin(double p[], double xi[])
{
	int j,i;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=2*num_atoms;
	pcom=dvector(1,num_atoms*2);
	xicom=dvector(1,num_atoms*2);
	for (j=1;j<=num_atoms*2;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=2*num_atoms;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	for(i=1;i<=num_atoms;i++)
	{
	  if(isfixed[i-1]==0)
	  {
	    X[i-1]=p[2*i-1];
	    Y[i-1]=p[2*i];
	  }
	}	
	free_dvector(xicom,1,2*num_atoms);
	free_dvector(pcom,1,2*num_atoms);
}


double brent(double ax, double bx, double cx, double (*f)(double), double tol,double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	return fx;
}





extern int ncom;
extern double *pcom,*xicom;

double f1dim(double x)
{
	int j,i;
	double f,*xt;

	xt=dvector(1,ncom);
	for (j=1;j<=ncom/2;j++)
	{  
	  if(isfixed[j-1]==0)
	  {
	    xt[2*j-1]=pcom[2*j-1]+x*xicom[2*j-1];
	    xt[2*j]=pcom[2*j]+x*xicom[2*j];
	  }
	}
	for(i=1;i<=num_atoms;i++)
	{
	  if(isfixed[i-1]==0)
	  {
	    X[i-1]=xt[2*i-1];
	    Y[i-1]=xt[2*i];
	  }
	}
	calc_all_energy();
	f=energy[2];
	
	for(i=1;i<=num_atoms;i++)
	{
	  if(isfixed[i-1]==0)
	  {
	    X[i-1]=pcom[2*i-1];
	    Y[i-1]=pcom[2*i];
	  }
	}
	
	free_dvector(xt,1,ncom);
	return f;
}



void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,double (*func)(double))
{
	double ulim,u,r,q,fu,dum;


	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),1.0e-20),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}



void CALC_FORCE()
{
	int count;
	int i,pi,pk;
	double R12,x_i,y_i,x_j,y_j,x_k,y_k,x1,y1,x2,y2;
	double xij,yij,xkj,ykj,yik;
	double Force;
	double temp,delta_theta,common_factor;
	double c1,c2,c3,c4,c5,c6,c7,c8;
	double *aX,*aY;

	aX=(double*)malloc(num_atoms*sizeof(double));
	aY=(double*)malloc(num_atoms*sizeof(double));
	for(count=0;count<num_atoms;count++) { aX[count]=0; aY[count]=0;}
/*
	printf("check node positions before force calculation!\n");
	for(count=0;count<num_atoms;count++)
	{
		printf("Node #%d,(%.1f,%.1f)\n",count,X[count],Y[count]);
	}*/
/*
	printf("Check all bonds!\n");
	for(count=0;count<num_bonds;count++)
	{
		printf("bond No%d,%d,%d,(%.1f,%.1f),(%.1f,%.1f)\n",count,Bonds[count][0],Bonds[count][1],\
				X[Bonds[count][0]],Y[Bonds[count][0]],X[Bonds[count][1]],Y[Bonds[count][1]]);
	}
	*/
	/* calculate stretching energy */
	for(count=0;count<num_bonds;count++)
	{
		
	    if( Bonds[count][0]>-1 && Bonds[count][1]>-1)
	    {
      		x1=X[Bonds[count][0]];
		y1=Y[Bonds[count][0]];
		x2=X[Bonds[count][1]];
		y2=Y[Bonds[count][1]];
		R12=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

		Force=(R12-A_lat0[count])*stiffness;
		aX[Bonds[count][0]]+=(x2-x1)/R12*Force;
		aY[Bonds[count][0]]+=(y2-y1)/R12*Force;
		aX[Bonds[count][1]]+=(x1-x2)/R12*Force;
		aY[Bonds[count][1]]+=(y1-y2)/R12*Force;
		if ( isnan(Force/R12) )
		{
			printf("singularity found for stretching force!,");
			printf("Points:%d,(%.2f,%.2f) and %d,(%.2f,%.2f),A_lat0=%.2f\n",Bonds[count][0],x1,y1,Bonds[count][1],x2,y2,A_lat0[count]);
			exit(-1);
		}
		if ( isnan(aX[Bonds[count][0]]) || isnan(aY[Bonds[count][0]]) || isnan(aX[Bonds[count][1]]) || isnan(aY[Bonds[count][1]]) )
		{
			printf("Oops!aX or aY goes to NAN!\n");
			printf("aX0+=%.2f,aY0+=%.2f,aX1+=%.2f,aY1+=%.2f\n",(x2-x1)/R12*Force,(y2-y1)/R12*Force,(x1-x2)/R12*Force,(y1-y2)/R12*Force);
			printf("x1=%.2f,y1=%.2f,x2=%.2f,y2=%.2f  bondNo %d, node:%d %d\n",x1,y1,x2,y2,count,Bonds[count][0],Bonds[count][1]);
			exit(-1);
		}
	    }
	}

/*	printf("At this point, all stretching forces are well calculated!\n"); */

	double Fix,Fiy,Fkx,Fky,Fjx,Fjy;
	if(kappa > 1.0e-5)
	{
		int vec1[3]={0,1,2};
		int vec2[3]={3,4,5};
		for(count=0;count<num_atoms;count++)
		{
			for(i=0;i<3;i++)
			{
				if(  (Junctions[count][vec1[i]]>=0) && (Junctions[count][vec2[i]]>=0) )
				{
					pi=Junctions[count][vec1[i]];
					pk=Junctions[count][vec2[i]];
					x_i=X[pi];  y_i=Y[pi];
					x_j=X[count]; y_j=Y[count];
					x_k=X[pk];  y_k=Y[pk];
					
					if( isnan(x_i)||isnan(y_i)||isnan(x_j)||isnan(y_j)||isnan(x_k)||isnan(y_k) )
					{
						printf("Nan found for coordinations: %d:(%.1f,%.1f) %d:(%.1f,%.1f) %d:(%.1f,%.1f)\n",pi,x_i,y_i,count,x_j,y_j,pk,x_k,y_k);
						exit(-1);
					}

					temp=((x_i-x_j)*(x_k-x_j)+(y_i-y_j)*(y_k-y_j));
					temp=temp/(sqrt((x_i-x_j)*(x_i-x_j)+(y_i-y_j)*(y_i-y_j)+DBL_EPSILON));
					temp=temp/(sqrt((x_k-x_j)*(x_k-x_j)+(y_k-y_j)*(y_k-y_j)+DBL_EPSILON));

					if (temp<-1.0) { temp=-1.0; }
					else if ( temp>1.0 ) {temp =1.0;}

					delta_theta=acos(temp)-M_PI;
					common_factor=-kappa*delta_theta/(sqrt(1-temp*temp)+DBL_EPSILON);

					xij=x_i-x_j;
					yij=y_i-y_j;
					xkj=x_k-x_j;
					ykj=y_k-y_j;
					yik=y_i-y_k;

					c1=x_k*(-yij)+x_j*(yik)+x_i*(ykj);   
					c2=xij*xij+yij*yij+LDBL_EPSILON;
					c3=xkj*xkj+ykj*ykj+LDBL_EPSILON;
   			                c4=(-c1)/(c2*sqrt(c2)*sqrt(c3));
				        c5=x_k*x_k*yij+x_j*x_j*yik+(x_i*x_i+yij*yik)*(-ykj)+2*x_j*(x_k*(-yij)+x_i*ykj);
				        c6=x_i*x_i*(-xkj)+x_j*x_j*x_k-x_k*yij*yij+x_i*(-x_j*x_j+x_k*x_k+ykj*ykj)-x_j*(x_k*x_k-yik*(yij+ykj));
				        c7=c2*sqrt(c2)*c3*sqrt(c3);
				        c8=(-c1)/(sqrt(c2)*c3*sqrt(c3));


					aX[pi]-=common_factor*(yij)*c4;
					aY[pi]-=-common_factor*(xij)*c4;
					aX[pk]-=-common_factor*(ykj)*c8;
					aY[pk]-=common_factor*(xkj)*c8;
					aX[count]-=common_factor*c1*c5/c7;
					aY[count]-=-common_factor*c1*c6/c7;
					
					Fix=-common_factor*(yij)*c4;
				        Fiy=common_factor*(xij)*c4;
					Fkx=common_factor*(ykj)*c8;
					Fky=-common_factor*(xkj)*c8;
					Fjx=-common_factor*c1*c5/c7;
					Fjy=common_factor*c1*c6/c7;
					if( isnan(Fix)||isnan(Fiy)||isnan(Fjx)||isnan(Fjy)||isnan(Fkx)||isnan(Fky) )
					{
						printf("numerical errors in calculating bending force!\n");
						printf("Junction:(%.2f,%.2f),(%.2f,%.2f),(%.2f,%.2f)\n",x_i,y_i,x_j,y_j,x_k,y_k);
						printf("Fix=%.2f,Fiy=%.2f,Fjx=%.2f,Fjy=%.2f,Fkx=%.2f,Fky=%.2f\n",Fix,Fiy,Fjx,Fjy,Fkx,Fky);
						exit(-1);
					}

				}
			}
		}
	} /* end of bending energy */
/*
	printf("At this point! All bending forces have been calculated without any singularities!\n");
	printf("Let us check the forces one by one!\n");
	for(count=0;count<num_atoms;count++)
	{
		printf("     node %d ax=%.2f,ay=%.2f\n",count,aX[count],aY[count]);
	}
	*/


	for(count=0;count<num_atoms;count++)
	{
		AX[count]=aX[count];  AY[count]=aY[count]; 
	}
	
	 free(aX);  free(aY);
}/* end of function CALC_FORCE() */

void SaveVector(double *X,double *X1,int N)
{
	int i;
	for(i=0;i<N;i++)
	{
		X1[i]=X[i];
	}
}

double deviation_func(double *solution_guess,int dim,void *params)
{
   int i;
   double *shift_guess,*angle_guess;
   shift_guess=(double*)malloc(numBoundaryPts*sizeof(double));
   angle_guess=(double*)malloc(numBoundaryPts*sizeof(double));
   for(i = 0; i<numBoundaryPts; i++)
   {
	   shift_guess[i] = solution_guess[i];
//	   angle_guess[i] = solution_guess[i+numBoundaryPts];
           angle_guess[i] = angle_known[i];
   }

   CellContraction(Xcenter,Ycenter,LinkedPts,numBoundaryPts,shift_guess,angle_guess);	
   conjugate_gradient();
   int num;
   double sum=0.0;
   for (num=0;num<num_atoms;num++)
   {
	   if (isfixed[num]==0 && isVisible[num] == 1)
	   {
		   sum+=(X[num]-X_observations[num])*(X[num]-X_observations[num]);
		   sum+=(Y[num]-Y_observations[num])*(Y[num]-Y_observations[num]);
	   }
   }
   return sum;
}

void PlotNet(char *filename)
{
	FILE *pFile;
	int k,a,b;
	pFile=fopen(filename,"w");
	if ( pFile==NULL)
	{
		//printf("File %s open failed!\n",filename);
	}
	else
	{
//		printf("Writing Net data into file %s ...\n",filename);
		for(k=0;k<num_bonds;k++)
		{
			a=Bonds[k][0];
			b=Bonds[k][1];
			if( a>-1 && b>-1) 
			{
			   fprintf(pFile,"%.3f %.3f %.3f %.3f\n",X[a],Y[a],X[b],Y[b]);
			}
		}
		fclose(pFile);
	}
}

double cart2pol(double dX,double dY)
{
	double theta;
	double tmpx,tmpy;
	if( fabs(dX)<DBL_EPSILON && fabs(dY)<DBL_EPSILON )
	{
		printf("0 vector!\n");
		return(0.);
	}
	
	tmpx=dX/sqrt(dX*dX+dY*dY);
	tmpy=dY/sqrt(dX*dX+dY*dY);
	theta=acos(tmpx);
	theta*=(tmpy<0)?(-1):1;
	theta=(theta>-M_PI+DBL_EPSILON)?theta:(-M_PI+DBL_EPSILON);
	
	return(theta);
}


	
double Fourier(double theta,double *Apara,double *Bpara,int num_orders)
{
	int k;
	double value = 0.0;
	Bpara[0] = 0.0;
	for(k = 0; k<=num_orders ; k++)
	{
		value+= Apara[k] * cos((double)k*theta);
		value+= Bpara[k] * sin((double)k*theta);
	}
	return(value);
}

double deviation_fourier(double *guess,int dim,void *params)
{
	int i,num_orders;
	double deviation_value;

	num_orders = (dim-1)/2;
	double *sol_guess = (double*)malloc(numBoundaryPts*sizeof(double));
	double *Ap = (double*)malloc((num_orders+1)*sizeof(double));
	double *Bp = (double*)malloc((num_orders+1)*sizeof(double));
	Ap[0] = guess[0];
	Bp[0] = 0.0;
	for(i = 1; i<=num_orders; i++)
	{
		Ap[i] = guess[i];
		Bp[i] = guess[i+num_orders];
	}


	for(i = 0 ; i<numBoundaryPts; i++)
	{
		sol_guess[i] = Fourier(angle_known[i],Ap,Bp,num_orders);
	}
//	printf("guess %.2e",sol_guess[i]);
//	printf("\n");

	deviation_value = deviation_func(sol_guess,numBoundaryPts,NULL);



	free(sol_guess);
	free(Ap);
	free(Bp);

	return deviation_value;
}
	

void setVisible(double Xcenter,double Ycenter,double VisibleRadius)
{
	int i;
	double dis;
	for(i = 0; i < num_atoms; i++)
	{
		dis = sqrt( (X[i]-Xcenter)*(X[i]-Xcenter)+(Y[i]-Ycenter)*(Y[i]-Ycenter) );
		isVisible[i] = (dis<VisibleRadius)?1:0;
	}
}




void recover_net()
{
	int i;
	for( i = 0 ; i < num_atoms ; i++)
	{
		X[i] = X_initial[i];
		Y[i] = Y_initial[i];
	}
}


#ifdef _ADIABATIC_ 


/* Attention: amoeba and amotry functions are taken from numerical recipes (C)
 * Here I modify it a little bit due to the fact a Numerical Recipes array starts from 1 
 * and a regular C array starts from 0 */


void amoeba(double **p, double y[], int ndim, double ftol,
	double (*funk)(double [],int, void*), int *nfunk)
{
	int i,ihi,ilo,inhi,j,mpts=ndim+1;
	double rtol,sum,swap,ysave,ytry,*psum;

	psum = (double*)malloc(ndim*sizeof(double));
	*nfunk=0;
	GET_PSUM
	for (;;) {
		ilo=0;
		ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
		for (i=0;i<mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
		if (rtol < ftol) {
			SWAP(y[0],y[ilo])
			for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i])
			break;
		}
		if (*nfunk >= NMAX) 
		{
			printf("NMAX exceeded!\n");
			SWAP(y[0],y[ilo])
			for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i])
			break;
		}
		*nfunk += 2;
		ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
		if (ytry <= y[ilo])
			ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
			if (ytry >= ysave) {
				for (i=0;i<mpts;i++) {
					if (i != ilo) {
						for (j=0;j<ndim;j++)
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						y[i]=(*funk)(psum,2*Fourier_Orders+1,NULL);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		} else --(*nfunk);
	}
	free(psum);
}



double amotry(double **p, double y[], double psum[], int ndim,
	double (*funk)(double [],int, void*), int ihi, double fac)
{
	int j;
	double fac1,fac2,ytry,*ptry;

	ptry = (double*)malloc(ndim*sizeof(double));
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry,2*Fourier_Orders+1,NULL);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=0;j<ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free(ptry);
	return ytry;
}

#endif
