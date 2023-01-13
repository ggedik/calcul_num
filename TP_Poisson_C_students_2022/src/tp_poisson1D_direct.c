/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;
  double alpha = 1; 
  double temp, relres;

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2; //leading dimension of RHS ?? 
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS   =(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X     =(double *) malloc(sizeof(double)*la);

  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;           //decalage
  ku=1;           //upper
  kl=1;           //lower
  lab=kv+kl+ku+1; //leading dimension

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

  /*
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  */

  printf("Solution with LAPACK\n");
  /* LU Factorization */
  info=0;
  ipiv = (int *) calloc(la, sizeof(int));
  //ierr = dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  /*
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  */
  
  /* Solution (Triangular) */
  if (info==0){
    //ierr = 
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info, 1);
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }

  /* It can also be solved with dgbsv (dgbtrf+dgbtrs) */
  /*  ierr = dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);*/

//EXO 4.
  //dgbmv("N", &la, &la, &kl, &ku, alpha, AB, &lab, EX_SOL, 1,0, RHS, 1 ); 


  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;

/*resoudre le systeme lineaire*/
/*DGBTRF: compute LU, partial pivot 
  DGBTRS: resoudre system in, fait appel à DGBTRF => so it should be slower than it 
  DGBSV : resoudre system */
  /*system : poisson * x = b => AB * x = rhs */
  /*prendre des mesures*/
  /*copying from above so that i can see evrything clearly*/
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);  //== AB 
  set_grid_points_1D(X, &la);                              //== X
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);                   //== RHS
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);//== analytical solution to our problem 

/**exo5**/
/*get time*/
  dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info, 1 ); //without 1: too few arguments error why? it's not in the doc
/*get time*/

/*get time*/
  dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
/*get time*/
  
/*A = LU => LUx = b => first calculate Ux then resoudre L (Ux) = b */ /*dgbtrf?? */
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

/*exo6*/


  printf("\nThe relative forward error is relres = %e\n",relres);

  //free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
