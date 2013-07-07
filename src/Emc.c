# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <time.h>
# include <R.h>
# include <Rinternals.h>
# include <Rdefines.h>

// The procedure performs a Monte Carlo resampling under a null hypothesis to calculate
// the E measure of interindividual variation and the value of the \eqn{C_{ws}} measure
// of modularity following Araujo et al. (2008).
//
// Author: Dr. Nicola ZACCARELLI (nicola.zaccarelli@gmail.com)
//
// Version 1
// Date: 10/11/2012


SEXP Emc(SEXP Proportions, SEXP IndType, SEXP PopDietType, SEXP DietToTInd, SEXP nreplicates)
{

time_t t;
int NInds, NPrey, i, j, k, x, R, Ki, TypeInd, nreps;
float item;
double Wmedio, Wmax, CWt, Wpot3, F, Si, CwS, E;
double tmp, tmp2, cumulativep, lowerbound;
double *totaldieti, *populationdiet, *DataTmp, *RisBoot;
double **MatrixP, **MatrixW;
SEXP Rris, Rdim;

// Get dimension of data
Rdim = getAttrib(Proportions, R_DimSymbol);
NInds = INTEGER(Rdim)[0];
NPrey = INTEGER(Rdim)[1];

// Coerce data matrix into vector for data access form C
PROTECT(Proportions = coerceVector(Proportions, REALSXP));
DataTmp = REAL(Proportions);

// Coerce other paramters
IndType = coerceVector(IndType, INTSXP);
TypeInd = INTEGER(IndType)[0];
PROTECT(PopDietType = coerceVector(PopDietType, REALSXP));
populationdiet = REAL(PopDietType);
PROTECT(DietToTInd = coerceVector(DietToTInd, REALSXP));
totaldieti = REAL(DietToTInd);
nreplicates = coerceVector(nreplicates, INTSXP);
nreps = INTEGER(nreplicates)[0];

// Create the return vector
PROTECT(Rris = allocMatrix(REALSXP, (nreps + 1), 4));
RisBoot = REAL(Rris);
// Initialize the vector
for (i = 0; i < (4 *(nreps + 1)); i++) { RisBoot[i] = 0; }

// create other matrices needed
MatrixP = calloc(NInds, sizeof(double *));
MatrixW = calloc(NInds, sizeof(double *));
for (i=0; i < NInds; i++)
    { MatrixP[i] = calloc(NPrey, sizeof(double));
      MatrixW[i] = calloc(NInds, sizeof(double)); }

// Initialize some matrices
for (j=0; j<NPrey; j++) {
  for (i=0; i<NInds; i++){
    MatrixP[i][j] = DataTmp[i + NInds*j];
}}


srand((unsigned) time(&t));

// Calculation and Monte Carlo resampling
for (R = 0; R < (nreps +1); R++)
    {
     if (R > 0)
	 {
	  for (i=0; i<NInds; i++) { for (j=0; j<NPrey; j++) {
				 MatrixP[i][j] = 0; }}
      for (i=0; i<NInds; i++)
          {for (x=0; x< totaldieti[i]; x++)
              {item = (double)rand()/(double)RAND_MAX;
               cumulativep = 0;
               for (j=0; j<NPrey; j++)
                   {lowerbound = cumulativep;
                    cumulativep = cumulativep + populationdiet[j];
                    if (item>=lowerbound && item<cumulativep) { MatrixP[i][j]= MatrixP[i][j] + 1; }
                    }
              } // for x
           }; // for i
     // Transform into proportions
     for (i=0; i<NInds; i++) for (j=0; j< NPrey; j++) MatrixP[i][j] = MatrixP[i][j] / totaldieti[i];
     } // end if (R>0)
    
     Wmedio = 0;
     Wmax=0;
     
// calculate MatrixW
for (i=0; i< NInds; i++)
   { 
     MatrixW[i][i] = 0;
	  for (k = 0; k<NInds; k++)
	     {
		 if (k != i)
			 {
			  tmp = 0;
              tmp2 = 0;
			 for (j=0; j<NPrey; j++)
				 {
				  if (MatrixP[i][j] > MatrixP[k][j])
					 tmp += (MatrixP[i][j] - MatrixP[k][j]);
				  else
					 tmp += (MatrixP[k][j] - MatrixP[i][j]);
				 }
			 if (TypeInd == 1) 
                            { tmp2 = (double)1.0 - (double)0.5*tmp;
                              MatrixW[i][k] = pow(tmp2, (double)1.0/ (double)3.0);
                             } else { 
                              tmp2 = (double)1.0 - (double)0.5*tmp;
                              MatrixW[i][k] = tmp2; }
			 Wmedio = Wmedio + tmp2;
			 if (tmp2 > Wmax) Wmax = tmp2;
			 } //else
	     } 
    }
    
// Calculate some indices
Wmedio = Wmedio /(double)(NInds*(NInds - 1));
E = 0;
E = 1- Wmedio;

// Initialise some variables
CWt = 0;
tmp2 = 0;
if (TypeInd == 1)
   { for (i=0; i<NInds; i++)
	 { Ki=0;
	   F=0;
	   tmp=0;
	   for (j=0; j<NInds; j++) {if (MatrixW[i][j] > 0) Ki = Ki + 1;}
	   for (j=0; j<NInds; j++)
		 {for (k=0; k<NInds; k++)
			{if ((i != j) && (i != k) && (j != k))
			   {if ((MatrixW[i][j] > 0) && (MatrixW[i][k] > 0) && (MatrixW[j][k] > 0))
				  { Wpot3 = (MatrixW[i][j]*MatrixW[j][k]*MatrixW[k][i]) / Wmax;
					  F = F + Wpot3;
				  }
			   }
			}
		 }
		 tmp = (Ki - 1.0);
		 if (tmp == 0) tmp = 1.0;
         CWt = CWt + F /(Ki*tmp);
	  }
   } else
   { for (i=0; i<NInds; i++)
	  {  Ki=0;
		 Si=0;
		 F = 0;
		 for (j=0; j<NInds; j++) { if (MatrixW[i][j] > 0) 
                                      { Ki = Ki + 1;
                                        Si = Si + MatrixW[i][j]; }
		 }
		 for (j=0; j<NInds; j++)
		 {for (k=0; k<NInds; k++)
			{if ((i != j) && (i != k) && (j != k))
			   {if ((MatrixW[i][j] > 0) && (MatrixW[i][k] > 0) && (MatrixW[j][k] > 0))
				  {Wpot3 = (MatrixW[i][j] + MatrixW[i][k])/((double)2);
				   F = (F + Wpot3);
				  }
			   }
			}
		 }
		 tmp = (Ki - 1.0);
		 if (tmp == 0) tmp = 1.0;
		 tmp2=(1.0/(Si*tmp));
         CWt = CWt + (tmp2 * F);
      }
   }
CWt = CWt /(double)NInds;
CwS = 0;
if ((CWt + Wmedio) != 0) {CwS = (CWt - Wmedio)/(CWt + Wmedio);} else {CwS= 0;}

  //Write results
  RisBoot[R + (nreps + 1) * 0] = Wmedio;
  RisBoot[R + (nreps + 1) * 1] = E;
  RisBoot[R + (nreps + 1) * 2] = CWt;
  RisBoot[R + (nreps + 1) * 3] = CwS;
} // end R cycle
UNPROTECT(4);

return Rris;}
