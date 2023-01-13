/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  
  int ii, jj, kk;
 
  for (jj=0; jj<(*la) ;jj++){
    
    kk = jj*(*lab);
    
    if (*kv>=0){
      for (ii=0; ii < *kv;ii++){
	    AB[kk+ii]=0.0;
      }
    }
    
    AB[kk+ *kv]   =  -1.0;
    AB[kk+ *kv+1] =   2.0;
    AB[kk+ *kv+2] =  -1.0;
  
  }

  AB[0]=0.0;
  if (*kv == 1) {
    AB[1] = 0;
    }
  
  AB[(*lab)*(*la)-1] = 0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]      = *BC0;
  RHS[(*la)-1]= *BC1;
  
  for (jj=1; jj < (*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

/*void eig_poisson1D(double* eigval, int *la){
  //dgeev
  double* WR   = (double*) malloc(sizeof(double) * la) ; //we have eigen values as many as our dimension maximum
  double* WI   = (double*) malloc(sizeof(double) * la) ;
  double* VL   = (double*) malloc(sizeof(double) * la) ;
  double* VR   = (double*) malloc(sizeof(double) * la) ;
  double* WORK = (double*) malloc(sizeof(double) * la) ;
  int LDLV  = 1;
  int LDVR  = 1; 
  int LWORK = 1;
  int info = 0;
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  //comme on va pas faire un changement sur AB, on ne le passe pas comme un parametre 
  dgeev("N","N", &la, AB, &lab, WR, WI, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &info);

}*/

double eigmax_poisson1D(int *la){
//max

  return 0;
}

double eigmin_poisson1D(int *la){

  //min 

  return 0;
}

double richardson_alpha_opt(int *la){
  
 // b = 2 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
  return 0; //b;
}
/*REPOSANT SUR LES APPELS BLAS*/

void richardson_alpha(double *AB,                 //matrix stocké format compacte
                      double *RHS,                //right hand side - B
                      double *X,                  //SOL est passé comme argument 
                      double *alpha_rich,         //alpha
                      int *lab,                   //leading dim of AB
                      int *la, 
                      int *ku, int *kl,           //info about the matrices 
                      double *tol,                //precision
                      int *maxit,                 //maximum number of iterations
                      double *resvec,             //residu_vecteur
                      int *nbite){ //??
//x_k+1 = x_k + alpha *( b - A*x_k)
//X_k+1 = X_k + (RHS - AB*x_k)
  int ii; 
  double* rhs_copy = (double*) malloc(sizeof(RHS));
  for (ii = 0; ii < *maxit; ++ii){
    memcpy(rhs_copy, RHS, sizeof(RHS));
    //resvec[ii] = (iteration_avant - iteration apres);
    if (/*resvec[ii] == tol*/0){
      break;
    }else{
      //dgemv va ecraser les valeurs donc il faut faire une copie

      cblas_dgemv("N",la, la, alpha_rich, AB, la, X, 1, 1, rhs_copy, 1);
      //resvec[ii] = 
    }
  }


}

void extract_MB_jacobi_tridiag(double *AB, 
                               double *MB, //La matrice M de Jacobi stocké format GB  M**-1 = D**-1
                               int *lab,  //kv + ku +kl +1
                               int *la,
                               int *ku, int*kl, int *kv){
  int ii; 
  for (ii = 0; ii < *la; ++ii){
    MB[(*lab)*ii+1] = AB[(*lab)*ii+1];
    //printf("%f\n", MB[(*lab)*ii+1]);
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, 
                                     double *MB, //matrice de gauss seidel = (D-E)**-1 = M**-1
                                     int *lab, 
                                     int *la,
                                     int *ku, int*kl, int *kv){

  int ii, jj; 
   extract_MB_jacobi_tridiag(AB, MB, lab, la, ku, kl, kv);
   for (ii = 0; ii< *la; ++ii){
    MB[(*lab)*ii+2] = AB[(*lab)*ii+2];
    //printf("%f\n", MB[(*lab)*ii+2]);
   }
  
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i; //exam
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
