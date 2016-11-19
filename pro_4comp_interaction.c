#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "nrutil.c"
//#include <lapack.h>
/*#include "f2c.h"
#include "cblas.h"
#include "clapack.h"*/
#include <time.h>
//#include <omp.h>

#define MIN 10e-5
#define MAX 10e5
#define MAXLEN 2048
#define NTIMES 1 // change later
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);



int ngroup;
double stat;
double **GG_tem;
double **GG;
clock_t begin, end;
double time_spent;
//double **combY;
//double ***mu;
//double ***mu_read;
double sigma_t;
double sigma_t_initial;
double *YG;
double **GG_inv;
double **MsigmacombY;
double **GsigmacombY_acc;
double ***GsigmacombX;
double ***MsigmacombX;
double ****MsigmacombX_bac;
double ****MsigmacombX_acc;
double ***MsigmacombX_total;
double ***MsigmacombM_acc;
double ***MsigmacombM_bac;
double ****MsigmacombM;
double ****MsigmacombM_total;
double *beta;
double **Fisher_inv;
double **Fisher;
double ***Fisher_bac;
double ***Fisher_acc;
double **combY_initial;
double **sigmacombY_initial;
double *beta_add;
double **beta_add_bac;
double **beta_add_acc;
double log_det;
double det_Ka;
double det_Kb;
double det_Na;
double det_Nb;
int snp_acc;
int Nbeta;
int Npeo;
//int bac_snp_total=25900;
int bac_snp_total=1;

int acc_snp_total=1;
//int acc_snp_total=5000;

//int bac_snp_total=10;
//int acc_snp_total=1;

//int acc_snp_total=3712;
//int acc_snp_total=3612;

//int bac_snp_total=1;
//int acc_snp_total=1;

double ****combX_acc;
double ****combX_bac;
double ***combX;
int ngroup_start=2;
const int Nparams=4;
const int Nbeta_start=6;
int nrep=3;
double **matrix1;
double **matrix2;
double **phi1;
double **phi2;
double *w1;
double **u1;
double **u1_tem;
double *w2;
double **u2;
double **u2_tem;
int *nPlant;
int **nComb;
int *nBac;
int ***nTriple;
double ***Y;
double ***Y_read;
double ***Y_sq_read;
double ***Cov_read[MAXLEN];
double ***Cov[MAXLEN];
double *params;
double **new_nComb_nPair;
double **new_nComb;
double ***new_nTriple;
//double ***sigmacombX_initial;
//double ****sigmacombX_initial_acc;
//double ****sigmacombX_initial_bac;
double **extend_Matrix1;
double **extend_Matrix2;
double **extend_matrix1;
double **extend_matrix1_inv;
double **extend_matrix2;
double **extend_matrix2_inv;
double **cross;
double **doublekernal;
double **quad1;
double **quad2;
double *tem_x;
double *tem_m;
double *new_nPlant;
double *new_nBac;
/*integer info;
integer info2;*/
int info;
int info2;
double **temcombX;
double **temcombX2;
double ***temcombX_acc;
double ***temcombX2_acc;
double ***temcombX_bac;
double ***temcombX2_bac;

int n1;
int n2;
FILE *matrixf1,*matrixf2,*covf,*outf_acc,*outf_bac,*genof,*countf,*memberf,*genof_bac,*mfile;

char matrixfile[MAXLEN]="matrix",covfile[MAXLEN]="covaraite",matrixfile2[MAXLEN]="matrix2",genofile[MAXLEN]="genofile",countfile[MAXLEN]="countfile",memberfile[MAXLEN]="memberfile",genofile_bac[MAXLEN]="genotype_bac.txt",majorfile[MAXLEN]="majoreffect",outfile_acc[MAXLEN]="stat_acc_min.txt",outfile_bac[MAXLEN]="stat_bac_min.txt";

char letter[MAXLEN];


int inverse(double *params);
//int solve(double **y, double **a,int n, double *m,double *x,integer *info, integer *info2);
//int solve(double **y, double **a,int n, double *m,double *x,int *info, int *info2);
int matinvchol(double **a, int n, double **inv, double *det);
void readin();
void initial(void);
void jacobi(double **a, int n, double d[], double **v, int *nrot);

int main(int argc, char *argv[])
{
    int i,j,p,k,n,m,t;
    
    if (argc>1){
        for (i=1;i<argc && argv[i][0]=='-';i++){
            switch(argv[i][1])
            {
                case 'n':
                    strncpy(countfile,argv[++i],MAXLEN);
                    printf("user specified count: %s\n",countfile);
                    //mfile=1;
                    break;
                case 'm':
                    strncpy(matrixfile,argv[++i],MAXLEN);
                    printf("user specified matrix: %s\n",matrixfile);
                    //mfile=1;
                    break;
                case 'p':
                    strncpy(matrixfile2,argv[++i],MAXLEN);
                    printf("user specified matrix: %s\n",matrixfile2);
                    //pfile=1;
                    break;
                case 'c':
                    strncpy(covfile,argv[++i],MAXLEN);
                   printf("user specified covaraite: %s\n",covfile);
                    //cfile=1;
                    break;
                case 'g':
                    strncpy(genofile,argv[++i],MAXLEN);
                 printf("user specified genofile: %s\n",genofile);
                  //  ofile=1;
                    break;
                case 's':
                    strncpy(genofile_bac,argv[++i],MAXLEN);
                 printf("user specified bacteria genofile: %s\n",genofile_bac);
                   // ofile=1;
                    break;
                case 'a':
                    strncpy(outfile_acc,argv[++i],MAXLEN);
                   printf("user specified output: %s\n",outfile_acc);
                    //ofile=1;
                    break;
                case 'b':
                    strncpy(outfile_bac,argv[++i],MAXLEN);
                    printf("user specified output: %s\n",outfile_bac);
                    //ofile=1;
                    break;
                default:
                  //  printf ("Unknown option");
                    exit(1);
            }
        }
    }
                    
                    
    
    if((countf=fopen(countfile,"r"))==NULL){
        printf("Cannot open count file.\n");
        exit(1);
    }
    
    fscanf(countf,"%d %d %d",&Npeo,&n1,&n2);
    
    if((matrixf1=fopen(matrixfile,"r"))==NULL){
        printf("Cannot open matrix file.\n");
        exit(1);
    }
    
    
    if((matrixf2=fopen(matrixfile2,"r"))==NULL){
        printf("Cannot open matrix file.\n");
        exit(1);
    }
    
    
    if((covf=fopen(covfile,"r"))==NULL){
        printf("Cannot open covariate file.\n");
        exit(1);
    }
    
    if((genof=fopen(genofile,"r"))==NULL){
        printf("Cannot open genotype file.\n");
        exit(1);
    }
    
    if((genof_bac=fopen(genofile_bac,"r"))==NULL){
        printf("Cannot open bacteria genotype file.\n");
        exit(1);
    }
    
    
    
    
    
    if((outf_bac=fopen(outfile_bac,"w"))==NULL){
        printf("Cannot open output file.\n");
        exit(1);
    }
    
    if((outf_acc=fopen(outfile_acc,"w"))==NULL){
        printf("Cannot open output file.\n");
        exit(1);
    }
    
    
    double **stat_acc=dmatrix(1,acc_snp_total,1,2);
    int **inter_acc=imatrix(1,acc_snp_total,1,2);
    
    double *stat_bac=dvector(1,bac_snp_total);
    int *inter_bac=ivector(1,bac_snp_total);

    
    GG_tem=dmatrix(0,ngroup_start,1,Nbeta_start);
    GG=dmatrix(1,ngroup_start,1,ngroup_start);
    //combY=dmatrix(1,n1,1,n2);
    //mu=d3tensor(1,n1,1,n2,1,3);
   // mu_read=d3tensor(1,n1,1,n2,1,3);
    YG=dvector(1,ngroup_start);
    GG_inv=dmatrix(1,ngroup_start,1,ngroup_start);
    GsigmacombY_acc=dmatrix(1,n2,1,acc_snp_total);//changed
    MsigmacombY=dmatrix(1,n1,1,n2);//changed
    GsigmacombX=d3tensor(1,n2,1,Nbeta_start,1,acc_snp_total);
    MsigmacombM_total=d4tensor(1,n1,0,ngroup_start,1,n1,0,ngroup_start);
    MsigmacombM_acc=d3tensor(1,n1,1,n2,1,n1);
    MsigmacombM_bac=d3tensor(1,n1,1,n2,1,n2);
    MsigmacombM=d4tensor(1,n1,1,n2,1,n1,1,n2);//changed
    MsigmacombX_acc=d4tensor(1,n1,1,n2,Nbeta_start,Nbeta_start,1,acc_snp_total);
    MsigmacombX_bac=d4tensor(1,n1,1,n2,Nbeta_start-2,Nbeta_start-1,1,bac_snp_total);
    MsigmacombX=d3tensor(1,n1,1,n2,1,Nbeta_start);
    Fisher_inv=dmatrix(1,Nbeta_start,1,Nbeta_start);
    Fisher=dmatrix(1,Nbeta_start,1,Nbeta_start);
    beta=dvector(1,Nbeta_start);
    Fisher_bac=d3tensor(1,Nbeta_start-1,1,Nbeta_start-1,1,bac_snp_total);
    Fisher_acc=d3tensor(1,Nbeta_start,1,Nbeta_start,1,acc_snp_total);
    combX_acc=d4tensor(1,n1,1,n2,Nbeta_start,Nbeta_start,1,acc_snp_total);
    combX_bac=d4tensor(1,n1,1,n2,Nbeta_start-2,Nbeta_start-1,1,bac_snp_total);
    beta_add_bac=dmatrix(Nbeta_start-2,Nbeta_start-1,1,bac_snp_total);
    beta_add_acc=dmatrix(Nbeta_start,Nbeta_start,1,acc_snp_total);
    beta_add=dvector(1,Nbeta_start);
    matrix1=dmatrix(1,n1,1,n1);
    matrix2=dmatrix(1,n2,1,n2);
    phi1=dmatrix(1,n1,1,n1);
    w1=dvector(1,n1);
    u1=dmatrix(1,n1,1,n1);
    u1_tem=dmatrix(1,n1,1,n1);
    phi2=dmatrix(1,n2,1,n2);
    w2=dvector(1,n2);
    u2=dmatrix(1,n2,1,n2);
    u2_tem=dmatrix(1,n2,1,n2);
    nPlant=ivector(1,n1);
    nBac=ivector(1,n2);
    nComb=imatrix(1,n1,1,n2);
    nTriple=i3tensor(1,n1,1,n2,1,3);
    Y=d3tensor(1,n1,1,n2,1,3);
    Y_read=d3tensor(1,n1,1,n2,1,3);
     Y_sq_read=d3tensor(1,n1,1,n2,1,3);
    for (i=1;i<=n1;i++){
        Cov[i]=d3tensor(1,n2,1,3,1,Nbeta_start);
        Cov_read[i]=d3tensor(1,n2,1,3,1,Nbeta_start);
    }
    params=dvector(1,Nparams);
    new_nComb_nPair=dmatrix(1,n1,1,n2);
    new_nComb=dmatrix(1,n1,1,n2);
    new_nTriple=d3tensor(1,n1,1,n2,1,3);
    //sigmacombX_initial=d3tensor(1,n1,1,n2,1,Nbeta_start);
   // sigmacombX_initial_acc=d4tensor(1,n1,1,n2,Nbeta_start,Nbeta_start,1,acc_snp_total);
    //sigmacombX_initial_bac=d4tensor(1,n1,1,n2,4,5,1,bac_snp_total);
    extend_Matrix1=dmatrix(1,n1,1,n1);
    extend_Matrix2=dmatrix(1,n2,1,n2);
    extend_matrix1=dmatrix(1,n1,1,n1);
    extend_matrix1_inv=dmatrix(1,n1,1,n1);
    extend_matrix2=dmatrix(1,n2,1,n2);
    extend_matrix2_inv=dmatrix(1,n2,1,n2);
    cross=dmatrix(1,n1,1,n2);
    doublekernal=dmatrix(1,n1,1,n2);
    quad1=dmatrix(1,n1,1,n1);
    quad2=dmatrix(1,n1,1,n2);
    tem_x=dvector(0,Npeo*Npeo-1);
    tem_m=dvector(0,Npeo*Npeo-1);
    new_nPlant=dvector(1,n1);
    new_nBac=dvector(1,n2);
    combY_initial=dmatrix(1,n1,1,n2);
    sigmacombY_initial=dmatrix(1,n1,1,n2);
    combX=d3tensor(1,n1,1,n2,1,Nbeta_start);
    temcombX=dmatrix(1,n1,1,Nbeta_start);
    temcombX2=dmatrix(1,n2,1,Nbeta_start);
    temcombX_acc=d3tensor(1,n1,Nbeta_start,Nbeta_start,1,acc_snp_total);
    temcombX2_acc=d3tensor(1,n2,Nbeta_start,Nbeta_start,1,acc_snp_total);
    temcombX_bac=d3tensor(1,n1,Nbeta_start-2,Nbeta_start-1,1,bac_snp_total);
    temcombX2_bac=d3tensor(1,n2,Nbeta_start-2,Nbeta_start-1,1,bac_snp_total);
    
    for (i=1;i<=n1;i++){
        nPlant[i]=0;
        for (j=1;j<=n2;j++){
            nComb[i][j]=0;
            for (k=1;k<=nrep;k++){
                Y[i][j][k]=0;
                Y_read[i][j][k]=0;
                Y_sq_read[i][j][k]=0;
                nTriple[i][j][k]=0;
            }
        }
    }
    
    
    for (i=1;i<=n2;i++){
        nBac[i]=0;
    }
    

    
    
    readin();
    
    initial();
    
    params[1]=0.0864937724;
    params[2]=1.8199200899;
    params[3]=0.2409041527;
    params[4]=0.0645489199;
    
    
    //new kinship for bac (more genotypes)
    params[1]=0.0865004127;
    params[2]=1.2614736635;
    params[3]=0.2409026570;
    params[4]=0.0645161743;
    
        
    //prepare for taking inverse
for (i=1;i<=n1;i++){
    new_nPlant[i]=0;
          for (j=1;j<=n2;j++){
              new_nComb[i][j]=0;
              for (k=1;k<=nrep;k++){
                   new_nTriple[i][j][k]=nTriple[i][j][k]/(1+params[3]*nTriple[i][j][k]);
                  new_nComb[i][j]+= new_nTriple[i][j][k];
              }
              new_nComb_nPair[i][j]=new_nComb[i][j]/(1+params[4]*new_nComb[i][j]);
              new_nPlant[i]+=new_nComb_nPair[i][j];
          }
}
 
    for (j=1;j<=n2;j++){
        new_nBac[j]=0;
        for (i=1;i<=n1;i++){
            new_nBac[j]+=new_nComb_nPair[i][j];
            
        }
    }

    
    
    //taking inverse
    double error;
    error=inverse(params);

    
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            combY_initial[i][j]=0;
            for (k=1;k<=nrep;k++){
                Y[i][j][k]=1/(1+params[3]*nTriple[i][j][k])*Y_read[i][j][k];
                combY_initial[i][j]+=1/(1+params[4]*new_nComb[i][j])*Y[i][j][k];
            }
        }
    }
    
    
        
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            sigmacombY_initial[i][j]=0;
            
            for(n=1;n<=n1;n++){
                for(m=1;m<=n2;m++){
                    
                    sigmacombY_initial[i][j]-=(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[i][m]-doublekernal[n][j])*combY_initial[n][m];
                    
                }
            }
        }
        
    }
    
    //prepare sigma_t_initial
    sigma_t_initial=0;
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            for (k=1;k<=nrep;k++){
                
                 sigma_t_initial+=Y_sq_read[i][j][k];
                
                 sigma_t_initial-=params[3]*(1+nTriple[i][j][k]*params[3])*Y[i][j][k]*Y[i][j][k];
                
                                
                
            }
            sigma_t_initial-=params[4]*combY_initial[i][j]*combY_initial[i][j]*(1+params[4]*new_nComb[i][j]);
            
             sigma_t_initial+=combY_initial[i][j]*sigmacombY_initial[i][j];
            
        }
    }
    //printf("sigma %lf\n",sigma_t);

    
//prepare Cov 1-3
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            for(p=1;p<=3;p++){
                combX[i][j][p]=0;
                for (k=1;k<=nrep;k++){
                    Cov[i][j][k][p]=1/(1+params[3]*nTriple[i][j][k])*Cov_read[i][j][k][p];
                    combX[i][j][p]+=1/(1+params[4]*new_nComb[i][j])*nTriple[i][j][k]*Cov[i][j][k][p];
                }

            }
        }
    }

    
//prepare coefficient beta_add[1.2.3]
    
    for (p=1;p<=3;p++){
        beta_add[p]=0;
        
        for (i=1;i<=n1;i++){
            temcombX[i][p]=0;
            for (j=1;j<=n2;j++){
                temcombX[i][p]+=combX[i][j][p];
                for (k=1;k<=nrep;k++){
                    beta_add[p]+=Cov_read[i][j][k][p]*Y[i][j][k];
                }
                
                beta_add[p]-=params[4]*combY_initial[i][j]*combX[i][j][p]*(1+params[4]*new_nComb[i][j]);
                
                beta_add[p]+=combX[i][j][p]*sigmacombY_initial[i][j];
                
            }
            
        }
        for(j=1;j<=n2;j++){
            temcombX2[j][p]=0;
            for(i=1;i<=n1;i++){
                temcombX2[j][p]+=combX[i][j][p];
            }
            
        }
        
        
    }
    

//prepare Fisher [1.2.3][1.2.3]
    for (p=1;p<=3;p++){
        for (t=p;t<=3;t++){
            
            Fisher[p][t]=0;
            
            
            for (j=1;j<=n2;j++){
                for(m=1;m<=n2;m++){
                    Fisher[p][t]-=extend_matrix2[j][m]*temcombX2[j][p]*temcombX2[m][t];
                    
                }
            }
            
            for (i=1;i<=n1;i++){
                
                for(n=1;n<=n1;n++){
                    Fisher[p][t]-=extend_matrix1[i][n]*temcombX[i][p]*temcombX[n][t];
                }
                
                
                for (j=1;j<=n2;j++){
                    for (k=1;k<=nrep;k++){
                        Fisher[p][t]+=new_nTriple[i][j][k]*Cov_read[i][j][k][p]*Cov_read[i][j][k][t];
                    }
                    
                    Fisher[p][t]-=params[4]*(1+params[4]*new_nComb[i][j])*combX[i][j][p]*combX[i][j][t];
                    
                    
                    Fisher[p][t]-=-doublekernal[i][j]*temcombX[i][p]*temcombX2[j][t]-doublekernal[i][j]*temcombX2[j][p]*temcombX[i][t];
                    
                }
                
            }
            Fisher[t][p]=Fisher[p][t];
        }
        
    }

    
    //prepare others...
    /*for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            for(p=1;p<=3;p++){
                sigmacombX_initial[i][j][p]=0;
                for(n=1;n<=n1;n++){
                    for(m=1;m<=n2;m++){
                        sigmacombX_initial[i][j][p]-=(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[i][m]-doublekernal[n][j])*combX[n][m][p];
                    }
                }
            }
        }
    }*/

    
    //MsigmacombY
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            MsigmacombY[i][j]= sigmacombY_initial[i][j]*new_nComb[i][j]/(1+params[4]*new_nComb[i][j])-params[4]*new_nComb[i][j]*combY_initial[i][j];
            ;
            for (k=1;k<=nrep;k++){
                MsigmacombY[i][j]+=Y[i][j][k];
            }

        }
        
    }
    
    
    //MsigmacombX[n][m][1.2.3]
    for(n=1;n<=n1;n++){
        for(m=1;m<=n2;m++){
            for(p=1;p<=3;p++){
                MsigmacombX[n][m][p]=-params[4]*new_nComb[n][m]*combX[n][m][p];
            }
                   for(k=1;k<=nrep;k++){
                            for(p=1;p<=3;p++){
                            MsigmacombX[n][m][p]+=new_nTriple[n][m][k]*Cov_read[n][m][k][p];
                            }
            }
        }
        }

    
    
    
    // MsigmacombM
    
    for(n=1;n<=n1;n++){
        for(j=1;j<=n2;j++){
            for(k=1;k<=n2;k++){
                MsigmacombM_bac[n][j][k]=0;
                for(i=1;i<=n1;i++){
                    MsigmacombM[n][j][i][k]=0;
                    MsigmacombM_acc[n][j][i]=0;
                }
            }
        }
        
    }
    for(n=1;n<=n1;n++){
        for(m=1;m<=n2;m++){
            MsigmacombM[n][m][n][m]=new_nComb[n][m]-params[4]/(1+params[4]*new_nComb[n][m])*new_nComb[n][m]*new_nComb[n][m];
        }
    }
    
    

    double tem_term;
    for(n=1;n<=n1;n++){
        
        for(m=1;m<=n2;m++){
            for(j=1;j<=n2;j++){
                
                for(i=1;i<=n1;i++){
                  tem_term=new_nComb_nPair[n][m]*(extend_matrix1[n][i]+extend_matrix2[m][j]-doublekernal[n][j]-doublekernal[i][m]);
                    
                    MsigmacombM[n][m][i][j]-=tem_term*new_nComb_nPair[i][j];
                    
                    for(p=1;p<=3;p++){
                        MsigmacombX[n][m][p]-=tem_term*combX[i][j][p];
                        
                    }
                    MsigmacombM_acc[n][m][i]+=MsigmacombM[n][m][i][j];
                    MsigmacombM_bac[n][m][j]+=MsigmacombM[n][m][i][j];

                }
            }
        }
    }
    
      printf("read bac\n");
    
//store bac info
//snp_bac=0;
fscanf(genof_bac,"%s ",letter);

for(i=1;i<=n2;i++){
    fscanf(genof_bac,"%s ", letter);
     //printf("%s\n",letter);
}
   
    int *count=ivector(1,bac_snp_total);
    int snp_bac;
    int chr1;
    double freq1,freq2;
    int gene;
    double ***G_bac=d3tensor(1,2,1,n2,1,bac_snp_total);
    int **groupmember_bac=imatrix(1,n2,1,bac_snp_total);
    int *ngroup_bac=ivector(1,bac_snp_total);
    
//while(fscanf(genof_bac,"%d ",&chr1)==1){
//snp_bac++;
for(snp_bac=1;snp_bac<=bac_snp_total;snp_bac++){
    
    fscanf(genof_bac,"%d ",&chr1);
    count[snp_bac]=chr1;
    freq1=freq2=0;
    
    for(i=1;i<=n2;i++){
        //below is for letter coding of genotype (original pattern)
        fscanf(genof_bac,"%d ", &gene);
        //printf("%d\n",gene);
        
        G_bac[1][i][snp_bac]=(double)(gene==1);//mutation
        
        G_bac[2][i][snp_bac]=(double)(gene==2);//deletion
        
        
        freq1+=G_bac[1][i][snp_bac];//mutation
        freq2+=G_bac[2][i][snp_bac];//deletion
        
        if(gene==1)
            groupmember_bac[i][snp_bac]=1;
        else if (gene==2)
            groupmember_bac[i][snp_bac]=2;
        else
            groupmember_bac[i][snp_bac]=0;
        //0 means original allele only, 1 means mutation as well, 2 mean deletion;
        
    }
    
   // printf("df is %d\n",(freq2!=0)+(freq1!=0)+((n2-freq1-freq2)!=0)-1);
    
    ngroup_bac[snp_bac]=(freq2!=0)+(freq1!=0)+((n2-freq1-freq2)!=0)-1;//degree of freedom
    
    /*if(freq2==0){
        ngroup_bac[snp_bac]=1;//degree of freedom
    }
    else
        ngroup_bac[snp_bac]=2;*/
    
    
    
        //add main effect
      for(i=1;i<=n1;i++){
            for(j=1;j<=n2;j++){
               combX_bac[i][j][4][snp_bac]= combX_bac[i][j][5][snp_bac]=0;
                
                for(k=1;k<=nrep;k++){
                    Cov_read[i][j][k][4]=G_bac[1][j][snp_bac];
                    Cov[i][j][k][4]=1/(1+params[3]*nTriple[i][j][k])*Cov_read[i][j][k][4];
                    
                    Nbeta=4;
                    if(ngroup_bac[snp_bac]==2){
                        Nbeta=5;
                        Cov_read[i][j][k][5]=G_bac[2][j][snp_bac];
                        Cov[i][j][k][5]=1/(1+params[3]*nTriple[i][j][k])*Cov_read[i][j][k][5];
                        
                        
                    }
                }
            }
      }

 
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            for(p=4;p<=Nbeta;p++){
                combX_bac[i][j][p][snp_bac]=new_nComb_nPair[i][j]*G_bac[p-3][j][snp_bac];
                
            }
            
        }
    }
    
    //prepare regression coefficient beta_add[4.5]

    for (p=4;p<=Nbeta;p++){
        beta_add_bac[p][snp_bac]=0;
        for (i=1;i<=n1;i++){
            temcombX_bac[i][p][snp_bac]=0;
            for (j=1;j<=n2;j++){
                temcombX_bac[i][p][snp_bac]+=combX_bac[i][j][p][snp_bac];
                for (k=1;k<=nrep;k++){
                    beta_add_bac[p][snp_bac]+=Cov_read[i][j][k][p]*Y[i][j][k];
                }
                
                beta_add_bac[p][snp_bac]-=params[4]*combY_initial[i][j]*combX_bac[i][j][p][snp_bac]*(1+params[4]*new_nComb[i][j]);
                
                beta_add_bac[p][snp_bac]+=combX_bac[i][j][p][snp_bac]*sigmacombY_initial[i][j];
                
            }
        }
        
        for(j=1;j<=n2;j++){
            temcombX2_bac[j][p][snp_bac]=0;
            for(i=1;i<=n1;i++){
                temcombX2_bac[j][p][snp_bac]+=combX_bac[i][j][p][snp_bac];
            }
        }
        
    }
    
     //prepare Fisher_bacsnp[4.5][4.5] and Fisher_bacsnp[4.5][4.5]
    
    for (p=4;p<=Nbeta;p++){
        for (t=p;t<=Nbeta;t++){
            Fisher_bac[p][t][snp_bac]=0;
            
            
            for (j=1;j<=n2;j++){
                for(m=1;m<=n2;m++){
                    Fisher_bac[p][t][snp_bac]-=extend_matrix2[j][m]*temcombX2_bac[j][p][snp_bac]*temcombX2_bac[m][t][snp_bac];
                    
                }
            }
            
           // printf("okay%d\t%d\n",p,t);
            for (i=1;i<=n1;i++){
                
                for(n=1;n<=n1;n++){
                    Fisher_bac[p][t][snp_bac]-=extend_matrix1[i][n]*temcombX_bac[i][p][snp_bac]*temcombX_bac[n][t][snp_bac];
                }
                
                
                for (j=1;j<=n2;j++){
                    for (k=1;k<=nrep;k++){
                        Fisher_bac[p][t][snp_bac]+=new_nTriple[i][j][k]*Cov_read[i][j][k][p]*Cov_read[i][j][k][t];
                    }
                    
                    Fisher_bac[p][t][snp_bac]-=params[4]*(1+params[4]*new_nComb[i][j])*combX_bac[i][j][p][snp_bac]*combX_bac[i][j][t][snp_bac];
                    
                    
                    Fisher_bac[p][t][snp_bac]-=-doublekernal[i][j]*temcombX_bac[i][p][snp_bac]*temcombX2_bac[j][t][snp_bac]-doublekernal[i][j]*temcombX2_bac[j][p][snp_bac]*temcombX_bac[i][t][snp_bac];
                    
                }
                
            }
            Fisher_bac[t][p][snp_bac]=Fisher_bac[p][t][snp_bac];
        }
    }
        // printf("\n");
   
    
    //Fisher[1.2.3][4.5] and Fisher[4.5][1.2.3]
    for (p=4;p<=Nbeta;p++){
        for (t=1;t<=3;t++){
            Fisher_bac[p][t][snp_bac]=0;
            
            
            for (j=1;j<=n2;j++){
                for(m=1;m<=n2;m++){
                    Fisher_bac[p][t][snp_bac]-=extend_matrix2[j][m]*temcombX2_bac[j][p][snp_bac]*temcombX2[m][t];
                    
                }
            }
 
            for (i=1;i<=n1;i++){
                
                for(n=1;n<=n1;n++){
                    Fisher_bac[p][t][snp_bac]-=extend_matrix1[i][n]*temcombX_bac[i][p][snp_bac]*temcombX[n][t];
                }
                
                
                for (j=1;j<=n2;j++){
                    for (k=1;k<=nrep;k++){
                        Fisher_bac[p][t][snp_bac]+=new_nTriple[i][j][k]*Cov_read[i][j][k][p]*Cov_read[i][j][k][t];
                    }
                    
                    Fisher_bac[p][t][snp_bac]-=params[4]*(1+params[4]*new_nComb[i][j])*combX_bac[i][j][p][snp_bac]*combX[i][j][t];
                    
                    
                    Fisher_bac[p][t][snp_bac]-=-doublekernal[i][j]*temcombX_bac[i][p][snp_bac]*temcombX2[j][t]-doublekernal[i][j]*temcombX2_bac[j][p][snp_bac]*temcombX[i][t];
                    
                }
                
            }
            Fisher_bac[t][p][snp_bac]=Fisher_bac[p][t][snp_bac];
        }
    }
    
    
    
    //prepare other ....
   /* for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            for(p=4;p<=Nbeta;p++){
                sigmacombX_initial_bac[i][j][p][snp_bac]=0;
                for(n=1;n<=n1;n++){
                    for(m=1;m<=n2;m++){
                        sigmacombX_initial_bac[i][j][p][snp_bac]-=(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[i][m]-doublekernal[n][j])*combX_bac[n][m][p][snp_bac];
                    }
                }
            }
        }
    }*/
    
    
    
    //MsigmacombX_bac[n1][n2][4.5]
    //new_nComb_nPair[i][j]
    
    for(n=1;n<=n1;n++){
        for(m=1;m<=n2;m++){
             for(p=4;p<=Nbeta;p++){
                 MsigmacombX_bac[n][m][p][snp_bac]=0;
                  for(j=1;j<=n2;j++){
            MsigmacombX_bac[n][m][p][snp_bac]+=MsigmacombM_bac[n][m][j]*G_bac[p-3][j][snp_bac];
        }
             }
    }
    }
    
       
    
}
    
    printf("read acc\n");
    
    // iterate by snp_acc
   // snp_acc=1;
    
    //while(fscanf(genof,"%d %d ",&chr[snp_acc],&pos[snp_acc])==2){
   
    int *chr=ivector(1,acc_snp_total);
    int *pos=ivector(1,acc_snp_total);
    int *nis1;
    int **is1;
    nis1=ivector(1,acc_snp_total);
    is1=imatrix(1,n1,1,acc_snp_total);
    double value;
    double **G_acc=dmatrix(1,n1,1,acc_snp_total);
    
    for(snp_acc=1;snp_acc<=acc_snp_total;snp_acc++){
          //begin=clock();
        
         /*if((snp_acc%1000)==1){
             printf("%d\n",snp_acc);
        }*/
        
        fscanf(genof,"%d %d ",&chr[snp_acc],&pos[snp_acc]);
        nis1[snp_acc]=0;
        
        for(i=1;i<=n1;i++){
            fscanf(genof,"%lf",&value);
            G_acc[i][snp_acc]=value/2;
            
            if((int)G_acc[i][snp_acc]==1){
                nis1[snp_acc]++;
                is1[nis1[snp_acc]][snp_acc]=i;
            }
                        
        }
        
        for (i=1;i<=n1;i++){
            for (j=1;j<=n2;j++){
                
                combX_acc[i][j][6][snp_acc]=new_nComb_nPair[i][j]*G_acc[i][snp_acc];
                
            }
        }

        //prepare beta_add beta[6]
        beta_add_acc[6][snp_acc]=0;
        
        for (i=1;i<=n1;i++){
            temcombX_acc[i][6][snp_acc]=0;
            for (j=1;j<=n2;j++){
                temcombX_acc[i][6][snp_acc]+=combX_acc[i][j][6][snp_acc];
                for (k=1;k<=nrep;k++){
                    beta_add_acc[6][snp_acc]+=G_acc[i][snp_acc]*Y[i][j][k];
                }
                
                beta_add_acc[6][snp_acc]-=params[4]*combY_initial[i][j]*combX_acc[i][j][6][snp_acc]*(1+params[4]*new_nComb[i][j]);
                
                beta_add_acc[6][snp_acc]+=combX_acc[i][j][6][snp_acc]*sigmacombY_initial[i][j];
                
            }
        }
        
        
        for(j=1;j<=n2;j++){
            temcombX2_acc[j][6][snp_acc]=0;
            for(i=1;i<=n1;i++){
                temcombX2_acc[j][6][snp_acc]+=combX_acc[i][j][6][snp_acc];
            }
        }
                
        
        //prepare Fisher[6][6]
        Fisher_acc[6][6][snp_acc]=0;
        
        
        for (j=1;j<=n2;j++){
            for(m=1;m<=n2;m++){
                Fisher_acc[6][6][snp_acc]-=extend_matrix2[j][m]*temcombX2_acc[j][6][snp_acc]*temcombX2_acc[m][6][snp_acc];
                
            }
        }
        
        
        for (i=1;i<=n1;i++){
            
            for(n=1;n<=n1;n++){
                Fisher_acc[6][6][snp_acc]-=extend_matrix1[i][n]*temcombX_acc[i][6][snp_acc]*temcombX_acc[n][6][snp_acc];
            }
            
            
            for (j=1;j<=n2;j++){
                    Fisher_acc[6][6][snp_acc]+=new_nComb[i][j]*G_acc[i][snp_acc]*G_acc[i][snp_acc];
            
                
                
                //printf("okay%d\t%d\t%d\t?\n",snp_acc,i,j);
                Fisher_acc[6][6][snp_acc]-=params[4]*(1+params[4]*new_nComb[i][j])*combX_acc[i][j][6][snp_acc]*combX_acc[i][j][6][snp_acc];
                
              //  printf("okay%d\t%d\t%d\t?\n",snp_acc,i,j);
                
                Fisher_acc[6][6][snp_acc]-=-doublekernal[i][j]*temcombX_acc[i][6][snp_acc]*temcombX2_acc[j][6][snp_acc]-doublekernal[i][j]*temcombX2_acc[j][6][snp_acc]*temcombX_acc[i][6][snp_acc];
                
               // printf("okay%d\t%d\t%d\t?\n",snp_acc,i,j);
                
            }
            
        }
        
        
        
        
        //prepare Fisher_accsnp[snp_acc][6][1.2.3] and Fisher_accsnp[snp_acc][1.2.3][6]
        
         for (t=1;t<=3;t++){
            
            Fisher_acc[6][t][snp_acc]=0;
            
            
            for (j=1;j<=n2;j++){
                for(m=1;m<=n2;m++){
                    Fisher_acc[6][t][snp_acc]-=extend_matrix2[j][m]*temcombX2_acc[j][6][snp_acc]*temcombX2[m][t];
                    
                }
            }
            
            for (i=1;i<=n1;i++){
                
                for(n=1;n<=n1;n++){
                    Fisher_acc[6][t][snp_acc]-=extend_matrix1[i][n]*temcombX_acc[i][6][snp_acc]*temcombX[n][t];
                }
                
                
                for (j=1;j<=n2;j++){
                    
                    for (k=1;k<=nrep;k++){
                        Fisher_acc[6][t][snp_acc]+=new_nTriple[i][j][k]*G_acc[i][snp_acc]*Cov_read[i][j][k][t];
                    }
                    
                    Fisher_acc[6][t][snp_acc]-=params[4]*(1+params[4]*new_nComb[i][j])*combX_acc[i][j][6][snp_acc]*combX[i][j][t];
                    
                    
                    Fisher_acc[6][t][snp_acc]-=-doublekernal[i][j]*temcombX_acc[i][6][snp_acc]*temcombX2[j][t]-doublekernal[i][j]*temcombX2_acc[j][6][snp_acc]*temcombX[i][t];
                    
                }
                
            }
            Fisher_acc[t][6][snp_acc]=Fisher_acc[6][t][snp_acc];
        }
    

        
             
        //prepare others
       /* for (i=1;i<=n1;i++){
            for (j=1;j<=n2;j++){
                sigmacombX_initial_acc[i][j][6][snp_acc]=0;
                for(n=1;n<=n1;n++){
                    for(m=1;m<=n2;m++){
                        sigmacombX_initial_acc[i][j][6][snp_acc]-=(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[i][m]-doublekernal[n][j])*combX_acc[n][m][6][snp_acc];
                    }
                }
            }
        }*/
        
        
        

        
        for(n=1;n<=n1;n++){
            for(m=1;m<=n2;m++){
                MsigmacombX_acc[n][m][6][snp_acc]=0;
                for(j=1;j<=n1;j++){
                    MsigmacombX_acc[n][m][6][snp_acc]+=MsigmacombM_acc[n][m][j]*G_acc[j][snp_acc];
                }
            }
        }
        
  
    
    //GsigmacombX= G^T M Sigma^{-1} M X [1.2.3]
      for(m=1;m<=n2;m++){
        for(p=1;p<=3;p++){
            GsigmacombX[m][p][snp_acc]=0;
            for(n=1;n<=n1;n++){
        GsigmacombX[m][p][snp_acc]+=G_acc[n][snp_acc]*MsigmacombX[n][m][p];
    }
        }
    }
        
        
        
        
         //GsigmacombX= G^T M Sigma^{-1} M X [6]
        for(m=1;m<=n2;m++){
                GsigmacombX[m][6][snp_acc]=0;
                for(n=1;n<=n1;n++){
                    GsigmacombX[m][6][snp_acc]+=G_acc[n][snp_acc]*MsigmacombX_acc[n][m][6][snp_acc];
                }
        }
        
        
    
            for(m=1;m<=n2;m++){
                GsigmacombY_acc[m][snp_acc]=0;
                for(n=1;n<=n1;n++){
                    GsigmacombY_acc[m][snp_acc]+=MsigmacombY[n][m]*G_acc[n][snp_acc];
            }
        }
    
               

        
    }
    
    /*fprintf(outf4,"%d\t%d\t",chr[1],1);
    for(snp_acc=1;snp_acc<=acc_snp_total;snp_acc++){
    
        fprintf(outf4,"%d\t",pos[snp_acc]);
    };
    fprintf(outf4,"\n");*/
    
    //test for pair-wise interaction
    
     fprintf(outf_bac,"index\tcount\tdf\tstat\tacc_chr\tacc_pos\tbp\n");
    
    printf("start testing\n");
    for(snp_bac=1;snp_bac<=bac_snp_total;snp_bac++){
        
        
        Nbeta=3+ngroup_bac[snp_bac];
        ngroup=ngroup_bac[snp_bac];
        
        //fprintf(outf4,"%d\t%d\t",snp_bac,ngroup);
        
        //prepare Fisher [1.2.3.4.5][4.5] and Fisher[4.5][1.2.3.4.5]
        for(p=4;p<=Nbeta;p++){
            for(t=1;t<=Nbeta;t++){
                Fisher[t][p]=Fisher[p][t]=Fisher_bac[p][t][snp_bac];
            }
        }

        
 
        for(n=1;n<=n1;n++){
            for(j=1;j<=ngroup;j++){
                for(i=1;i<=n1;i++){
                for(m=1;m<=ngroup;m++){
                    MsigmacombM_total[n][j][i][m]=0;
            }
        }
            }
        }
        
        for(n=1;n<=n1;n++){
            for(m=1;m<=n2;m++){
                for(i=1;i<=n1;i++){
                    for(j=1;j<=n2;j++){
                        MsigmacombM_total[n][groupmember_bac[m][snp_bac]][i][groupmember_bac[j][snp_bac]]+=MsigmacombM[n][m][i][j];
                    }
                }

            }
        }
        
        
        Nbeta=Nbeta+1;
        
        
        begin=clock();
        
        //for acc
        
      /* omp_set_num_threads(4);
        #pragma omp parallel 
        #pragma omp threadprivate(snp_acc,snp_bac,Fisher,Fisher_acc,Nbeta,t,i,j,p,m,n,k,beta_add,Fisher_inv,tem_m,tem_x,beta,GsigmacombX,stat,GG_tem,GG,GG_inv,YG,sigma_t,combY,mu,mu_read,ngroup,chr,pos,groupmember_bac,info,info2)
          {
        
        #pragma omp for*/
        
        for(snp_acc=1;snp_acc<=acc_snp_total;snp_acc++){
            
            //#pragma omp critical
            //{
            //prepare Fisher [1.2.3][6] and Fisher [6][1.2.3] and
            for(t=1;t<=3;t++){
                Fisher[t][Nbeta]=Fisher[Nbeta][t]=Fisher_acc[6][t][snp_acc];
            }
            
            
           /* #pragma omp critical
            printf("%d\t",snp_acc);
            {
           for(i=1;i<=3;i++){
               // for(n=1;n<=Nbeta;n++){
                   // printf("%d\t%d\t%lf\n",snp_acc,i,beta_add[i]);
                   printf("%lf\t",Fisher[i][Nbeta]);
               // }
               // printf("\n");
            }
printf("\n");
        }*/
            
            // prepare Fisher [6][6]
            Fisher[Nbeta][Nbeta]=Fisher_acc[6][6][snp_acc];

            //prepare Fisher[4.5.][6] and Fisher[6][4.5.] or prepare Fisher[4][5] and Fisher [5][4]
            for (p=4;p<=Nbeta-1;p++){
                
                Fisher[p][Nbeta]=0;
                
                
                for (j=1;j<=n2;j++){
                    for(m=1;m<=n2;m++){
                        Fisher[p][Nbeta]-=extend_matrix2[j][m]*temcombX2_bac[j][p][snp_bac]*temcombX2_acc[m][6][snp_acc];
                        
                    }
                }
                
                
                for (i=1;i<=n1;i++){
                    
                    for(n=1;n<=n1;n++){
                        Fisher[p][Nbeta]-=extend_matrix1[i][n]*temcombX_bac[i][p][snp_bac]*temcombX_acc[n][6][snp_acc];
                    }
                    
                    for (j=1;j<=n2;j++){
                        for (k=1;k<=nrep;k++){
                            Fisher[p][Nbeta]+=new_nTriple[i][j][k]*G_bac[p-3][j][snp_bac]*G_acc[i][snp_acc];
                        }
                        
                        Fisher[p][Nbeta]-=params[4]*(1+params[4]*new_nComb[i][j])*combX_bac[i][j][p][snp_bac]*combX_acc[i][j][6][snp_acc];
                        
                        
                        Fisher[p][Nbeta]-=-doublekernal[i][j]*temcombX_bac[i][p][snp_bac]*temcombX2_acc[j][6][snp_acc]-doublekernal[i][j]*temcombX2_bac[j][p][snp_bac]*temcombX_acc[i][6][snp_acc];
                        
                    }
                    
                }
                
                Fisher[Nbeta][p]=Fisher[p][Nbeta];
            }
          
            //prepare beta_add
           for (i=4;i<=Nbeta-1;i++){
                beta_add[i]=beta_add_bac[i][snp_bac];
            }
            beta_add[Nbeta]=beta_add_acc[6][snp_acc];
            
/*#pragma omp critical
            {
                printf("%d\n",snp_acc);
          for(i=1;i<=Nbeta;i++){
                for(n=1;n<=Nbeta;n++){
              //printf("%d\t%d\t%lf\n",snp_acc,i,beta_add[i]);
                    printf("%lf\t",Fisher[i][n]);
                }
                printf("\n");
            }
        }*/

            //error=solve(Fisher_inv,Fisher,Nbeta,tem_m,tem_x,&info,&info2);
            matinvchol(Fisher, Nbeta, Fisher_inv, &det_Na);

            
            //regression coefficient
            for (i=1;i<=Nbeta;i++){
                beta[i]=0;
                for (j=1;j<=Nbeta;j++){
                    beta[i]+=Fisher_inv[i][j]*beta_add[j];
                }
            }
              
            //GsigmacombX =G^T M Sigma^{-1} M X [4.5];
            for(p=4;p<=Nbeta-1;p++){
                for(m=1;m<=n2;m++){
                    GsigmacombX[m][p][snp_acc]=0;
                    for(n=1;n<=n1;n++){
                      GsigmacombX[m][p][snp_acc]+= G_acc[n][snp_acc]*MsigmacombX_bac[n][m][p][snp_bac];
                    }
                }
            }
            
          
            
            stat=0;
            

            for(p=1;p<=Nbeta;p++){
                for(m=1;m<=ngroup;m++){
                    GG_tem[m][p]=0;
                }
            }
            
             for(m=1;m<=n2;m++){
             for(p=1;p<=Nbeta-1;p++){
                     GG_tem[groupmember_bac[m][snp_bac]][p]+=GsigmacombX[m][p][snp_acc];
             }
                  GG_tem[groupmember_bac[m][snp_bac]][Nbeta]+=GsigmacombX[m][6][snp_acc];
            }
            
                                   
            for(p=1;p<=ngroup;p++){
                for(k=1;k<=p;k++){
                    GG[p][k]=GG_inv[p][k]=GG_inv[k][p]=0;
                     for(n=1;n<=nis1[snp_acc];n++){
                        for(i=1;i<=nis1[snp_acc];i++){
                            GG[p][k]+=MsigmacombM_total[is1[n][snp_acc]][p][is1[i][snp_acc]][k];
                        }
                    }
                    GG[k][p]=GG[p][k];
               }
            }

            
                
            for(p=1;p<=ngroup;p++){
                    for(k=1;k<=p;k++){
                        for(m=1;m<=Nbeta;m++){
                            for(t=1;t<=Nbeta;t++){
                                GG[p][k]-=GG_tem[p][m]*Fisher_inv[m][t]*GG_tem[k][t];
                            }
                        }
                        GG[k][p]=GG[p][k];
                    }
                }
                
                


            
           error=matinvchol(GG, ngroup, GG_inv, &det_Ka);
            
            
            
            for(k=1;k<=ngroup;k++){
                YG[k]=0;
                for(p=1;p<=Nbeta;p++){
                    YG[k]-=beta[p]*GG_tem[k][p];
                }
            }
            
            
            for(m=1;m<=n2;m++){
                    YG[groupmember_bac[m][snp_bac]]+=GsigmacombY_acc[m][snp_acc];
                }
      
           /* end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
            printf("time spent in total 1 %lf\n", time_spent);*/
            

            
            //sigma_t
            
             // sigma_t=0;
            
           // begin=clock();
           /* for (i=1;i<=n1;i++){
                for (j=1;j<=n2;j++){
                   // combY[i][j]=0;
                    for (k=1;k<=nrep;k++){
                        mu[i][j][k]=0;
                        mu_read[i][j][k]=0;
                        for (p=1;p<=3;p++){
                            mu_read[i][j][k]+=Cov_read[i][j][k][p]*beta[p];
                        }
                        for(p=4;p<=Nbeta-1;p++){
                            mu_read[i][j][k]+=G_bac[p-3][j][snp_bac]*beta[p];
                        }
                       // mu_read[i][j][k]+=G_acc[i][snp_acc]*beta[Nbeta];
                        
                       // mu[i][j][k]=1/(1+params[3]*nTriple[i][j][k])*mu_read[i][j][k];
                        
                        //combY[i][j]+=1/(1+params[4]*new_nComb[i][j])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k]);
                    }
                }
            }*/
            
          
        

            
           // printf("check %lf\n", sigma_t);
            sigma_t=sigma_t_initial;
            //sigma_t=0;
            for (p=1;p<=Nbeta;p++){
                // sigma_t-=combY_initial[i][j]*beta[6]*sigmacombX_initial_acc[i][j][6][snp_acc];
                sigma_t-=beta_add[p]*beta[p];
            }
            
           // printf("this is sigma_t third check %.10lf\n",sigma_t);
            

     
            
            
            sigma_t=sigma_t*(1+params[1]+params[2]+params[3]+params[4])/Npeo;

            
            
            for(p=1;p<=ngroup;p++){
                for(k=1;k<=p;k++){
                    stat+=(2-(p==k))*YG[p]*GG_inv[p][k]*YG[k]*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
                }
            }
         // printf("this is stat %d\t%d\t%lf\n",snp_acc,snp_bac,stat);
            
            //omit below
          /*
           if(ngroup==1){
               if(stat>15.14){ // p<10^(-4)
                    fprintf(outf4,"%d\t%d\t%d\t",chr[snp_acc],pos[snp_acc],snp_bac);
                    fprintf(outf4,"%lf\t%d\n",stat,ngroup);
                   // printf("%d\t%d\t%d\t",chr[snp_acc],pos[snp_acc],snp_bac);
                   // printf("%lf\t%d\n",stat,ngroup);
               }
            }
            if(ngroup==2){
                if(stat>18.42){
                    fprintf(outf4,"%d\t%d\t%d\t",chr[snp_acc],pos[snp_acc],snp_bac);
                    fprintf(outf4,"%lf\t%d\n",stat,ngroup);
                   // printf("%d\t%d\t%d\t",chr[snp_acc],pos[snp_acc],snp_bac);
                    //printf("%lf\t%d\n",stat,ngroup);
              }
            }*/
             //fflush(outf4);
        
            
        //}
            
           // printf("%d\t%d\t%lf\n",pos[snp_acc],snp_bac,stat);
            
            
            if(stat>stat_acc[snp_acc][ngroup]){
                inter_acc[snp_acc][ngroup]=snp_bac;
                stat_acc[snp_acc][ngroup]=stat;
            }
            
            if(stat>stat_bac[snp_bac]){
                stat_bac[snp_bac]=stat;
                inter_bac[snp_bac]=snp_acc;
            }
            
           

        } //finish loop on acc
        //  }
        // fflush(outf4);
        //fprintf(outf4,"\n");
        //fflush(outf4);
        
        /*end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("time spent in total %lf\n", time_spent);*/
    
        fprintf(outf_bac,"%d\t",snp_bac);
        fprintf(outf_bac,"%d\t",count[snp_bac]);
        fprintf(outf_bac,"%d\t",ngroup_bac[snp_bac]);
        fprintf(outf_bac,"%lf\t",stat_bac[snp_bac]);
        fprintf(outf_bac,"%d\t%d\t%d\n",chr[1],inter_bac[snp_bac],pos[inter_bac[snp_bac]]);
        
         fflush(outf_bac);
        
    }//finish loop on bac
    
    
    
    
    fprintf(outf_acc,"chr\tpos\tdf1\tbac_pos\tdf2\tbac_pos\n");
    for(snp_acc=1;snp_acc<=acc_snp_total;snp_acc++){
        fprintf(outf_acc,"%d\t",chr[snp_acc]);
        fprintf(outf_acc,"%d\t",pos[snp_acc]);
        fprintf(outf_acc,"%lf\t",stat_acc[snp_acc][1]);
        fprintf(outf_acc,"%d\t",inter_acc[snp_acc][1]);
        fprintf(outf_acc,"%lf\t",stat_acc[snp_acc][2]);
        fprintf(outf_acc,"%d\n",inter_acc[snp_acc][2]);
    }
     fflush(outf_acc);
    
    //fprintf(outf_bac,"index\tcount\tdf\tstat\tacc_chr\tacc_pos\tbp\n");
   /* for(snp_bac=1;snp_bac<=bac_snp_total;snp_bac++){
        fprintf(outf_bac,"%d\t",snp_bac);
        fprintf(outf_bac,"%d\t",count[snp_bac]);
        fprintf(outf_bac,"%d\t",ngroup_bac[snp_bac]);
        fprintf(outf_bac,"%lf\t",stat_bac[snp_bac]);
        fprintf(outf_bac,"%d\t%d\t%d\n",chr[1],inter_bac[snp_bac],pos[inter_bac[snp_bac]]);
    }

    fflush(outf_bac);*/
    
    
    
   // free_dmatrix(combY,1,n1,1,n2);
    //free_d3tensor(mu,1,n1,1,n2,1,3);
   // free_d3tensor(mu_read,1,n1,1,n2,1,3);
    free_dvector(YG,1,ngroup_start);
    free_dmatrix(GG_inv,1,ngroup_start,1,ngroup_start);
    free_dmatrix(GsigmacombY_acc,1,n2,1,acc_snp_total);//changed
    free_dmatrix(MsigmacombY,1,n1,1,n2);//changed
    free_d3tensor(GsigmacombX,1,n2,1,Nbeta_start,1,acc_snp_total);
    //free_d4tensor(MsigmacombM_total,1,n1,0,ngroup_start,1,n1,0,ngroup_start);
    free_d3tensor(MsigmacombM_acc,1,n1,1,n2,1,n1);
    free_d3tensor(MsigmacombM_bac,1,n1,1,n2,1,n2);
    //free_d4tensor(MsigmacombM,1,n1,1,n2,1,n1,1,n2);//changed
   // free_d4tensor(MsigmacombX_acc,1,n1,1,n2,Nbeta_start,Nbeta_start,1,bac_snp_total);
   // free_d4tensor(MsigmacombX_bac,1,n1,1,n2,Nbeta_start-2,Nbeta_start-1,1,acc_snp_total);
    free_d3tensor(MsigmacombX,1,n1,1,n2,1,Nbeta_start);
    free_dmatrix(Fisher_inv,1,Nbeta_start,1,Nbeta_start);
    free_dmatrix(Fisher,1,Nbeta_start,1,Nbeta_start);
    free_dvector(beta,1,Nbeta_start);
    free_d3tensor(Fisher_bac,1,Nbeta_start-1,1,Nbeta_start-1,1,bac_snp_total);
    free_d3tensor(Fisher_acc,1,Nbeta_start,1,Nbeta_start,1,acc_snp_total);
   // free_d4tensor(combX_acc,1,n1,1,n2,Nbeta_start,Nbeta_start,1,acc_snp_total);
    //free_d4tensor(combX_bac,1,n1,1,n2,Nbeta_start-2,Nbeta_start-1,1,bac_snp_total);
    free_dmatrix(beta_add_bac,Nbeta_start-2,Nbeta_start-1,1,bac_snp_total);
    free_dmatrix(beta_add_acc,Nbeta_start,Nbeta_start,1,acc_snp_total);
    free_dvector(beta_add,1,Nbeta_start);
    free_dmatrix(matrix1,1,n1,1,n1);
    free_dmatrix(matrix2,1,n2,1,n2);
    free_dmatrix(phi1,1,n1,1,n1);
    free_dvector(w1,1,n1);
    free_dmatrix(u1,1,n1,1,n1);
    free_dmatrix(u1_tem,1,n1,1,n1);
    free_dmatrix(phi2,1,n2,1,n2);
    free_dvector(w2,1,n2);
    free_dmatrix(u2,1,n2,1,n2);
    free_dmatrix(u2_tem,1,n2,1,n2);
    free_ivector(nPlant,1,n1);
    free_ivector(nBac,1,n2);
    free_imatrix(nComb,1,n1,1,n2);
    free_i3tensor(nTriple,1,n1,1,n2,1,3);
    free_d3tensor(Y,1,n1,1,n2,1,3);
    free_d3tensor(Y_read,1,n1,1,n2,1,3);
    free_d3tensor(Y_sq_read,1,n1,1,n2,1,3);
    for (i=1;i<=n1;i++){
        free_d3tensor(Cov[i],1,n2,1,3,1,Nbeta_start);
        free_d3tensor(Cov_read[i],1,n2,1,3,1,Nbeta_start);
    }
    free_dvector(params,1,Nparams);
    free_dmatrix(new_nComb_nPair,1,n1,1,n2);
    free_dmatrix(new_nComb,1,n1,1,n2);
    free_d3tensor(new_nTriple,1,n1,1,n2,1,3);
    //free_d3tensor(sigmacombX_initial,1,n1,1,n2,1,Nbeta_start);
    //free_d4tensor(sigmacombX_initial_acc,1,n1,1,n2,Nbeta_start,Nbeta_start,1,acc_snp_total);
    //free_d4tensor(sigmacombX_initial_bac,1,n1,1,n2,4,5,1,bac_snp_total);
    free_dmatrix(extend_Matrix1,1,n1,1,n1);
    free_dmatrix(extend_Matrix2,1,n2,1,n2);
    free_dmatrix(extend_matrix1,1,n1,1,n1);
    free_dmatrix(extend_matrix1_inv,1,n1,1,n1);
    free_dmatrix(extend_matrix2,1,n2,1,n2);
    free_dmatrix(extend_matrix2_inv,1,n2,1,n2);
    free_dmatrix(cross,1,n1,1,n2);
    free_dmatrix(doublekernal,1,n1,1,n2);
    free_dmatrix(quad1,1,n1,1,n1);
    free_dmatrix(quad2,1,n1,1,n2);
    free_dvector(tem_x,0,Npeo*Npeo-1);
    free_dvector(tem_m,0,Npeo*Npeo-1);
    free_dvector(new_nPlant,1,n1);
    free_dvector(new_nBac,1,n2);
    free_dmatrix(combY_initial,1,n1,1,n2);
    free_dmatrix(sigmacombY_initial,1,n1,1,n2);
    free_d3tensor(combX,1,n1,1,n2,1,Nbeta_start);
    free_dmatrix(temcombX,1,n1,1,Nbeta_start);
    free_dmatrix(temcombX2,1,n2,1,Nbeta_start);
    free_d3tensor(temcombX_acc,1,n1,Nbeta_start,Nbeta_start,1,acc_snp_total);
    free_d3tensor(temcombX2_acc,1,n2,Nbeta_start,Nbeta_start,1,acc_snp_total);
    free_d3tensor(temcombX_bac,1,n1,Nbeta_start-2,Nbeta_start-1,1,bac_snp_total);
    free_d3tensor(temcombX2_bac,1,n2,Nbeta_start-2,Nbeta_start-1,1,bac_snp_total);
    
    free_dmatrix(GG_tem,0,ngroup_start,1,Nbeta_start);
    free_dmatrix(GG,1,ngroup_start,1,ngroup_start);

}

    int inverse(double *params){
        
        double scaler;
        int i,j,k,t,q;
        int error=0;
        
        
        scaler=params[1]*params[2];
        
        
        //calcuate LU^T
        for (i=1;i<=n1;i++){
            for (j=1;j<=n2;j++){
                cross[i][j]=0;
                for (k=1;k<=n1;k++){
                    for (t=1;t<=n2;t++){
                        cross[i][j]+=new_nComb_nPair[k][t]*u1[i][k]*u2[j][t];
                    }
                    
                }
                
                
            }
            
        }
        
        
        for (i=1;i<=n1;i++){
            for (j=i;j<=n1;j++){
                extend_Matrix1[i][j]=0;
                for (k=1;k<=n1;k++){
                    extend_Matrix1[i][j]+=new_nPlant[k]*u1[i][k]*u1[j][k];
                }
                
                
            }
            
        }
        
        
        for (i=1;i<=n2;i++){
            for (j=i;j<=n2;j++){
                extend_Matrix2[i][j]=0;
                
                
                for (k=1;k<=n2;k++){
                    extend_Matrix2[i][j]+=new_nBac[k]*u2[i][k]*u2[j][k];
                    
                }
                
            }
            
        }
        
        //matrix N_a
        for (i=1;i<=n1;i++){
            for (j=i;j<=n1;j++){
                extend_matrix1[i][j]=extend_matrix1[j][i]=params[1]*extend_Matrix1[i][j]+(i==j);
                //  printf("%d\t%d\t%lf\n",i,j,extend_matrix1[i][j]);
            }
            
        }
        //printf("%lf\n",extend_matrix1[i][j]);
        
        
        //error+=solve(extend_matrix1_inv,extend_matrix1,n1,tem_m,tem_x,&info,&info2);
        
        matinvchol(extend_matrix1, n1, extend_matrix1_inv, &det_Na);
        
        //matrix N_b
        for (i=1;i<=n2;i++){
            for (j=i;j<=n2;j++){
                extend_matrix2[i][j]=extend_matrix2[j][i]=params[2]*extend_Matrix2[i][j]+(i==j);
                
                
            }
            
        }
        
      //  error+=solve(extend_matrix2_inv,extend_matrix2,n2,tem_m,tem_x,&info,&info2);
         matinvchol(extend_matrix2, n2, extend_matrix2_inv, &det_Na);
        
        
        
        // K_a
        for (i=1;i<=n1;i++){
            for (j=1;j<=n2;j++){
                quad2[i][j]=0;
                for (k=1;k<=n2;k++){
                    quad2[i][j]+=cross[i][k]*extend_matrix2_inv[k][j];
                }
                
            }
            
        }
        
        for (i=1;i<=n1;i++){
            for (j=i;j<=n1;j++){
                for (t=1;t<=n2;t++){
                    extend_matrix1[i][j]-=scaler*quad2[i][t]*cross[j][t];
                }
                extend_matrix1[j][i]=extend_matrix1[i][j];
                
            }
            
        }
        
        
        //K_b
        for (i=1;i<=n2;i++){
            for (j=i;j<=n2;j++){
                for (k=1;k<=n1;k++){
                    for (t=1;t<=n1;t++){
                        extend_matrix2[i][j]-=scaler*cross[k][i]*extend_matrix1_inv[k][t]*cross[t][j];
                    }
                }
                
                extend_matrix2[j][i]=extend_matrix2[i][j];
                
                
            }
            
            
        }
        
        // printf("this is extend_matrix1  %lf\n",extend_matrix1[1][1]);
        
        
        //error+=solve(extend_matrix1_inv,extend_matrix1,n1,tem_m,tem_x,&info,&info2);
        matinvchol(extend_matrix1, n1, extend_matrix1_inv, &det_Kb);
        

        
        // error+=solve(extend_matrix2_inv,extend_matrix2,n2,tem_m,tem_x,&info,&info2);
        matinvchol(extend_matrix2, n2, extend_matrix2_inv, &det_Kb);
        
        
        
        
        //printf("this is extend_matrix1 invv %lf\n",extend_matrix1_inv[1][1]);
        
        // M_a
        for (i=1;i<=n1;i++){
            for (j=1;j<=n1;j++){
                quad1[i][j]=0;
                for (k=1;k<=n1;k++){
                    quad1[i][j]+=u1[k][i]*extend_matrix1_inv[k][j];
                }
                
            }
            
            
        }
        
        
        
        
        for (i=1;i<=n1;i++){
            for (j=i;j<=n1;j++){
                extend_matrix1[i][j]=0;
                for (t=1;t<=n1;t++){
                    extend_matrix1[i][j]+=quad1[i][t]*u1[t][j];
                    
                }
                extend_matrix1[j][i]=extend_matrix1[i][j]=params[1]*extend_matrix1[i][j];
                
            }
            
        }
        
        //M_b
        for (i=1;i<=n2;i++){
            for (j=i;j<=n2;j++){
                extend_matrix2[i][j]=0;
                
                for (k=1;k<=n2;k++){
                    for (t=1;t<=n2;t++){
                        extend_matrix2[i][j]+=u2[k][i]*extend_matrix2_inv[k][t]*u2[t][j];
                    }
                }
                extend_matrix2[j][i]=extend_matrix2[i][j]=params[2]*extend_matrix2[i][j];
                
            }
            
        }
        
        //printf("iver? %lf\t%lf\n",u2[2][2],u2[1][1]);
        
        
        //M_c
        for (i=1;i<=n1;i++){
            for (j=1;j<=n2;j++){
                doublekernal[i][j]=0;
                
                for (k=1;k<=n1;k++){
                    for (t=1;t<=n2;t++){
                        doublekernal[i][j]+=quad1[i][k]*quad2[k][t]*u2[t][j];
                    }
                }
                doublekernal[i][j]=scaler*doublekernal[i][j];
                
            }
            
        }
        
        
        
        //  printf("inverse is %lf\n",extend_matrix1[1][1]+extend_matrix2[1][1]-2*doublekernal[1][1]);
        
        return(error);
        
    }

    
    
    
 void initial(void){
        
        int i,j,t,k;
        int tem;
        
        
        jacobi(matrix2, n2, w2, u2_tem,&tem);
        
        
        
        for (i=1;i<=n2;i++){
            for (j=1;j<=n2;j++){
                if(w2[i]>=0)
                    u2[i][j]=sqrt(w2[i])*u2_tem[j][i];
                else
                    u2[i][j]=0;
                
            }
        }
        
        
        jacobi(matrix1, n1, w1,u1_tem,&tem);
        
        for (i=1;i<=n1;i++){
            for (j=1;j<=n1;j++){
                if(w1[i]>=0)
                    u1[i][j]=sqrt(w1[i])*u1_tem[j][i];
                else
                    u1[i][j]=0;
            }
        }
        
        
        
        
        
        
        
        
    }
 
    
    
void readin(){
    
    
        int i,j,id,id2;
        int ind1, ind2,rep, leaf;
        double value, value2, value3, value4;
    
   // printf("%d\n",n1);
        for (i=1;i<=n1;i++){
            for (j=1;j<=n1;j++){
                fscanf(matrixf1,"%lf", &value);
                matrix1[i][j]=value;
                phi1[i][j]=value;
                //printf("%lf\t",value);
            }
          // printf("\n");
        }
    
   // printf("okay%d\t?\n",n2);
        
        for (i=1;i<=n2;i++){
            for (j=1;j<=n2;j++){
                //printf("%d\t%d\n",i,j);
                fscanf(matrixf2,"%lf", &value);
                matrix2[i][j]=value;
                phi2[i][j]=value;
              // printf("%lf\t",matrix2[i][j]);
            }
           //printf("\n");
        }
        
        
         //printf("okay%d\t?\n",n2);
        // printf("okay?\n");
        
        double valuenew=0;
        double total;
        // printf("sucess so far?\n");
        while(fscanf(covf,"%d %d %d %d %lf %lf %lf %lf %lf", &id, &id2, &rep,&leaf,&value,&value3,&value4,&value2,&valuenew)==9){
            
            
            // Plant[ind1]=id;
            nPlant[id]++;
            nBac[id2]++;
            nComb[id][id2]++;
            nTriple[id][id2][rep]++;
            // printf("%d\t%d\t%d\t%lf\n",id,id2,rep,Cov[id][id2][rep][1]);
            Cov_read[id][id2][rep][1]=value;
            Cov_read[id][id2][rep][2]=value3;
            Cov_read[id][id2][rep][3]=value4;
            value2=value2/4;
            Y_read[id][id2][rep]+=value2;
            Y_sq_read[id][id2][rep]+=value2*value2;
            /*if(value21==1){*/
            //printf("%d %d %d %d \n",  id, id2, rep,leaf);
            // }*/
        }
        
}
    
    void jacobi(double **a, int n, double d[], double **v, int *nrot)
    {
        // jacobi(matrix2, n2, w2, u2_tem,&tem);
        
        /* int i,j,k;
         
         double *S=dvector(0,n-1);
         double *U=dvector(0,n*n-1);
         double *VT=dvector(0,n*n-1);
         
         integer nn=n;
         integer nn2=5*n;
         
         double *work=dvector(1,nn2);
         
         int con;
         
         for (i=1;i<=n;i++){
         for (j=1;j<=n;j++)
         tem_m[(i-1)*n+j-1]=a[i][j];
         }
         integer info;
         
         dgesvd_("A", "A", &nn, &nn, tem_m, &nn, S, U, &nn, VT, &nn, work,&nn2,&info);
         
         con=S[0];
         for(i=1;i<=n;i++){
         // printf("%lf\t",S[i-1]);
         if (S[i-1]<=tol*con)
         d[i]=0;
         else
         d[i]=1/S[i-1];
         //printf("%lf\t",S[i-1]);
         }
         
         for (i=0;i<=(n-1);i++){
         for (j=0;j<=(n-1);j++){
         v[i+1][j+1]=U[i*n+j];
         //v[i+1][j+1]=VT[i*n+j];
         // fprintf(outf3,"%lf\t",phi1_G_inv_tem1[i+1][j+1]);
         }
         // fprintf(outf3,"\n");
         }
         }*/
        
        
        
        
        int j,iq,ip,i;
        double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
        
        
        b=dvector(1,n);
        z=dvector(1,n);
        for (ip=1;ip<=n;ip++) {
            for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
            v[ip][ip]=1.0;
        }
        for (ip=1;ip<=n;ip++) {
            b[ip]=d[ip]=a[ip][ip];
            z[ip]=0.0;
        }
        *nrot=0;
        for (i=1;i<=50;i++) {
            sm=0.0;
            for (ip=1;ip<=n-1;ip++) {
                for (iq=ip+1;iq<=n;iq++)
                    sm += fabs(a[ip][iq]);
            }
            if (sm == 0.0) {
                free_dvector(z,1,n);
                free_dvector(b,1,n);
                return;
            }
            if (i < 4)
                tresh=0.2*sm/(n*n);
            else
                tresh=0.0;
            for (ip=1;ip<=n-1;ip++) {
                for (iq=ip+1;iq<=n;iq++) {
                    g=100.0*fabs(a[ip][iq]);
                    if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
                        && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
                        a[ip][iq]=0.0;
                    else if (fabs(a[ip][iq]) > tresh) {
                        h=d[iq]-d[ip];
                        if ((double)(fabs(h)+g) == (double)fabs(h))
                            t=(a[ip][iq])/h;
                        else {
                            theta=0.5*h/(a[ip][iq]);
                            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                            if (theta < 0.0) t = -t;
                        }
                        c=1.0/sqrt(1+t*t);
                        s=t*c;
                        tau=s/(1.0+c);
                        h=t*a[ip][iq];
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        a[ip][iq]=0.0;
                        for (j=1;j<=ip-1;j++) {
                            ROTATE(a,j,ip,j,iq)
                        }
                        for (j=ip+1;j<=iq-1;j++) {
                            ROTATE(a,ip,j,j,iq)
                        }
                        for (j=iq+1;j<=n;j++) {
                            ROTATE(a,ip,j,iq,j)
                        }
                        for (j=1;j<=n;j++) {
                            ROTATE(v,j,ip,j,iq)
                        }
                        ++(*nrot);
                    }
                }
            }
            for (ip=1;ip<=n;ip++) {
                b[ip] += z[ip];
                d[ip]=b[ip];
                z[ip]=0.0;
            }
        }
        nrerror("Too many iterations in routine jacobi");
    }

    
    
    int matinvchol(double **a, int n, double **inv, double *det)
    /*
     Get the inverse of a matrix A by taking the Cholesky decomposition
     then solving for C^-T. Finally taking the producy C^-T * C^-1
     to get A^-1.
     
     Arguments:
     
     a    (input/output) double n x n array. Holds the matrix to invert on
     entry. On exit holds C^-T.// change inside the function so that a keeps unchanged.
     n    (input) integer. Number of rows and columns of a and inv.
     inv  (output) double n x n array. On output holds the inverted matrix.
     det  (output) pointer to double. On output holds the log of determinant of a.
     return 0 if non positive definite matrix; return 1 if okay
     
     */
    {
        double d, ldet, *col, *diag, sum;
        double **a_temp;
        a_temp=dmatrix(1,n,1,n);
        int i, j, k;
        int choldc(double **a, int n, double *diag);
        double error;
        for(i=1;i<=n;i++){
            for(j=1;j<=n;j++){
                a_temp[i][j]=a[i][j];
            }
        }
        
        diag = dvector(1,n);
        
        /* Do cholesky decomposition */
        error=choldc(a_temp,n,diag);
        
        
        /* Get the inverse of the cholesky decomp */
        for (i=1; i<=n; i++){
            a_temp[i][i] = 1.0/diag[i];
            for (j=i+1; j<=n; j++) {
                sum = 0.0;
                for (k=i; k<j; k++) sum -= a_temp[j][k] * a_temp[k][i];
                a_temp[j][i] = sum/diag[j];
            }
        }
        
        /* Get inverse of a */
        for (i=1; i<=n; i++) {
            for (j=1; j<=i; j++) {
                inv[i][j] = 0.0;
                for (k=i; k<=n; k++) inv[i][j] += a_temp[k][i] * a_temp[k][j];
                inv[j][i] = inv[i][j];
            }
        }
        
        /* Calculate log determinant */
        /*for (j=1; j<=n; j++) d *= a[j][j];
         *det = d;*/
        for (ldet=0.0, j=1; j<=n; j++) {
            if (diag[j] < 0) {
                d *= -1;
                diag[j] *= -1;
            }
            ldet += log(diag[j]);
        }
        *det = 2 * ldet;
        // printf("\nlog determinant= %lf\n",*det);
        fflush(stdout);
        
        free_dvector(diag,1,n);
        free_dmatrix(a_temp,1,n,1,n);
        return(error);
    }
    
    
    int choldc(double **a, int n, double p[])
    /* From Numerical Recipes */
    {
        int i,j,k, retval=1;
        double sum;
        
        for (i=1;i<=n;i++) {
            for (j=i;j<=n;j++) {
                for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
                if (i == j) {
                    if (sum <= 0.0){
                        fprintf(stderr,"\n WARNING: Non positive definite matrix.\n");
                        retval = 0;
                    }
                    p[i]=sqrt(sum);
                } else a[j][i]=sum/p[i];
            }
        }
        return(retval);
    }
    
    
/*    int solve(double **y, double **a,int n, double *m,double *x,int *info,int *info2){
     
     
     int i,j;
     
     
     
     for (i=1;i<=n;i++){
     for (j=1;j<=n;j++)
     m[(i-1)*n+j-1]=a[i][j];
     }
     
     for (i=1;i<=n;i++){
     for (j=1;j<=n;j++)
     x[(i-1)*n+j-1]=(i==j);
     }
     
     
     
     
     // integer nn=n;
     
     
     int *ipiv;
     ipiv=malloc(sizeof(int)*n);
     
     //dposv_("U",&nn,&nn,m,&nn,x,&nn,info);
     
     *info=clapack_dposv(CblasRowMajor, CblasUpper, n, n, m, n, x, n);
     
     
     *info2=0;
     
     if(*info!=0){
     
     for (i=1;i<=n;i++){
     for (j=1;j<=n;j++)
     m[(i-1)*n+j-1]=a[i][j];
     }
     
     for (i=1;i<=n;i++){
     for (j=1;j<=n;j++)
     x[(i-1)*n+j-1]=(i==j);
     }
 
     //dgesv_(&nn,&nn,m,&nn,ipiv,x,&nn,info2);
     *info2=clapack_dgesv(CblasRowMajor, n, n, m, n, ipiv,x, n);
     }
     
     // k=clapack_dposv(CblasRowMajor, CblasUpper, n, n, m, n, x, tmp);
     
     
     
     for (i=0;i<=(n-1);i++){
     for (j=0;j<=(n-1);j++)
     y[i+1][j+1]=x[i*n+j];
     }
     
     
     return ((*info!=0)&&(*info2!=0));
     
     }*/
    
 /*   int solve(double **y, double **a,int n, double *m,double *x,integer *info,integer *info2){
        int i,j;
        
        
        
        for (i=1;i<=n;i++){
            for (j=1;j<=n;j++)
                m[(i-1)*n+j-1]=a[i][j];
        }
        
        for (i=1;i<=n;i++){
            for (j=1;j<=n;j++)
                x[(i-1)*n+j-1]=(i==j);
        }
        
        
        
        
        integer nn=n;
        
        
        integer ipiv[n];
        
        dposv_("U",&nn,&nn,m,&nn,x,&nn,info);
        
        *info2=0;
        
        if(*info>0){
            
            
            for (i=1;i<=n;i++){
                for (j=1;j<=n;j++)
                    m[(i-1)*n+j-1]=a[i][j];
            }
            
            for (i=1;i<=n;i++){
                for (j=1;j<=n;j++)
                    x[(i-1)*n+j-1]=(i==j);
            }
            
            dgesv_(&nn,&nn,m,&nn,ipiv,x,&nn,info2);
        }
        
        // k=clapack_dposv(CblasRowMajor, CblasUpper, n, n, m, n, x, tmp);
        
        
        
        for (i=0;i<=(n-1);i++){
            for (j=0;j<=(n-1);j++)
                y[i+1][j+1]=x[i*n+j];
        }
        
        
        return ((*info!=0)&&(*info2!=0));
  
    }*/
    
    
    
    double pythag(double a, double b)
    {
        double absa,absb;
        absa=fabs(a);
        absb=fabs(b);
        if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
        else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
    }
    


