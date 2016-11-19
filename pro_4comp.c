#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string.h>
#include "nrutil.h"
#include "nrutil.c"
//#include <clapack.h>
#include "f2c.h"
#include "cblas.h"
#include "clapack.h"
#include <time.h>
//#include <omp.h>

#define MIN 10e-5
#define MAX 10e5
#define MAXLEN 2048
#define NTIMES 1 // change later
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);

double *result_opt;
double freq1,freq2;
double **G_bac;

double **Y_no;
double **Y_newTriple;
double **new_nComb_nPair;
double **Y_pair;
double **GG;
double *YG;
double **YG2;
double **GG_inv;
double mean;
double VV_general;
double tol;
/*const n1=20;
const n2=24;
const Npeo=33163;*/
int snp;

int nrep=3;
int n1;
//const int n1=4;
 //const n1=114;
//const n2=24;
//const n2=22;
//const int n222;
int n2;
 //const Npeo=30462;
//const Npeo=2hn627;
//const Npeo=28678;
//const Npeo=15204;
//const int Npeo=30469;
//const int Npeo=15204;//france
//const int Npeo=15265;//usa
//const int Npeo=30469;
int Npeo;//const int Npeo=1027;

//const Npeo=20484;

//const Npeo=15269;

const int nmember=1;
const int bio=4;
/*const n1=2;
const n2=3;
const Npeo=7;*/
double **phi1_G_inv_member;
double **member_phi1_G_inv_member;
double **member_phi1_G_inv_member_inv;
double **phi1_inv_general;
double **phi2_inv_general;
double **phi1_G_inv;
double **totalFisher;
double **totalFisher_inv;
double **sandwich;
double **sandwich_inv;

const int Nparams=4;
//const int Nbeta=6;
//const int Nbeta=0;
const int Nbeta=3;
int *membership;
double *freqmembership;
int *countmembership;
//const int Nparams=4;
//const Nbeta=2620;

//const ngenes=124963;
//const ngenes=1419007;
//const ngenes=4;
//const ngenes=127314;
//const ngenes=1225496;
//const ngenes=1;
char matrixfile[MAXLEN]="matrix",covfile[MAXLEN]="covaraite",matrixfile2[MAXLEN]="matrix2",genofile[MAXLEN]="genofile",outfile1[MAXLEN]="output",outfile2[MAXLEN]="output2",outfile3[MAXLEN]="chisq_bac",outfile4[MAXLEN]="chisq",countfile[MAXLEN]="countfile",memberfile[MAXLEN]="memberfile",genofile_bac[MAXLEN]="genotype_bac.txt";
double *tem_x;
double *tem_m;
double new_N;
double new_N_add;
int *nonmissing;
double freq;
double freqblue;

double *params;
//int *Plant;
//int *Bac;
clock_t begin, end;
double **matrix1;
double **matrix2;
double **phi1;
double **phi1test;
double **phi2;
FILE *matrixf1,*matrixf2,*covf,*outf1,*outf2,*outf3,*outf4,*genof,*countf,*memberf,*genof_bac;

int *nPlant;
int **nComb;
int *nBac;
double ***new_nTriple;
double **new_nComb;
double **new_nComb_nComb;
double **new_nComb_nComb1;
double **new_nComb_nComb2;
double **new_nComb_nComb3;
double **new_nComb_nComb4;
double *new_nBac;
double *new_nBac_1;
double *new_nPlant;




double time_spent;
double V1;
double **extend_Matrix1;
double **extend_Matrix2;
double *w1;
double **u1;
double **u1_tem;
double **extend_matrix1;
double **extend_matrix1_inv;
double *w2;
double **u2;
double **u2_tem;
double **extend_matrix2;
double **extend_matrix2_inv;
double **cross;
double **doublekernal;
double **quad1;
double **quad2;
double *tem_x1;
double *tem_m1;



double *G1;
double **GG_general;
double ***Y;
double ***Y_read;
double *G;
double Y_sq;
double ***Y_sq_read;
double ***Cov_read[MAXLEN];
double ***Cov_pair;
int ***nTriple;
double ***Cov[MAXLEN];
double **Sigmamat_inv;
double **Fisher;
double **Fisher_inv;
double **tem;
double sigma_t;
double a,b;
double *beta;
double *beta_new;
double ***mu;
double ***mu_read;
double ***M;
double **M_trans;
double *tvector;
double **tmatrix;
double *tvector2;
double **tmatrix2;
int *groupmember;


double **MsigmacombM_bac;
double **XsigmacombX_bac;
double *MsigmacombY_bac;

double ****MsigmacombM;
double **combM;
double ***MsigmacombX;
double **MsigmacombY;
double **XsigmacombY;
double ****XsigmacombX;

double Ftol;    
double *Point;
double *Delta;
double Output;

double *beta_add;
double ***combX;
double **combY;
double **combY_sq;
double **sigmacombY;
double **temcombX;
double **temcombX2;
integer info;
integer info2;
/*int info;
int info2;*/
double **matrix_phi1;
double **matrix_phi2;
double **matrix_phi3;
double **matrix_phi4;
double **Fisher_cov;
double **Fisher_cov_inv;
double quad;
double trace_inv;
double trace_invphi1;
double trace_invphi2;
double *Vvector;
double **phi1_inv;
double *temtest;
double **tem_R;

//double *test;


//change
int ngroup=22;
//int ngroup=1;
double log_det;
double det_Ka;
double det_Kb;
double det_Na;
double det_Nb;
void matinvchol(double **a, int n, double **inv, double *det);

double chisqtest_bac_retro(void);
double chisqtest_bac(void);
double chisqtest_new(void);
double chisqtest_strains(void);
int solve_general(double **y, double **a,int n, double *m,double *x,integer *info,integer *info2);
double chisqtest_retro(void);
//int solve(double **y, double **a,int n, double *m,double *x,int *info, int *info2);
int solve(double **y, double **a,int n, double *m,double *x,integer *info, integer *info2);
void readin();
int inverse(double *params);
void initial(void);
void jacobi(double **a, int n, double d[], double **v, int *nrot);

int coefficient(double *params);

double function(double *params); 
void minimize(double *start, double (*func)(double *));
double chisqtest(void);

int main(int argc, char *argv[])
{ 
  
   // printf("okay?\n");
     //tol=0.0001;
    int i,j,k,t,n,m,p;
    int mfile=0,pfile=0,ofile=0,cfile=0;
    if (argc>1){
        for (i=1;i<argc && argv[i][0]=='-';i++){
            switch(argv[i][1])
            {
                case 'n':
                    strncpy(countfile,argv[++i],MAXLEN);
                    printf("user specified count: %s\n",countfile);
                    mfile=1;
                    break;
                case 'm':
                    strncpy(matrixfile,argv[++i],MAXLEN);
                    printf("user specified matrix: %s\n",matrixfile);
                    mfile=1;
                    break;
                case 'l':
                    strncpy(memberfile,argv[++i],MAXLEN);
                    printf("user specified matrix test: %s\n",memberfile);
                    mfile=1;
                    break;
                case 'p':
                    strncpy(matrixfile2,argv[++i],MAXLEN);
                    printf("user specified matrix: %s\n",matrixfile2);
                    pfile=1;
                    break;
                case 'c':
                    strncpy(covfile,argv[++i],MAXLEN);
                    printf("user specified covaraite: %s\n",covfile);
                    cfile=1;
                    break;
                case 'g':
                    strncpy(genofile,argv[++i],MAXLEN);
                    printf("user specified genofile: %s\n",genofile);
                    ofile=1;
                    break;
                case 's':
                    strncpy(genofile_bac,argv[++i],MAXLEN);
                    printf("user specified bacteria genofile: %s\n",genofile_bac);
                    ofile=1;
                    break;

                case 'o':
                    strncpy(outfile1,argv[++i],MAXLEN);
                    printf("user specified output: %s\n",outfile1);
                    ofile=1;
                    break;
                case 'b':
                    strncpy(outfile2,argv[++i],MAXLEN);
                    printf("user specified output: %s\n",outfile2);
                    ofile=1;
                    break;
                case 'r':
                    strncpy(outfile3,argv[++i],MAXLEN);
                    printf("user specified output: %s\n",outfile3);
                    ofile=1;
                    break;
                case 't':
                    strncpy(outfile4,argv[++i],MAXLEN);
                    printf("user specified output: %s\n",outfile4);
                    ofile=1;
                    break;
                default:
                    printf ("Unknown option");
                    exit(1);
            }    
        }
    }

    
    
    if((countf=fopen(countfile,"r"))==NULL){
        printf("Cannot open count file.\n");
        exit(1);
    }
    
    fscanf(countf,"%d %d %d",&Npeo,&n1,&n2);
    //printf("%d\t%d\t%d\n",Npeo,n1,n2);
    
    int nmax=fmax(n1,n2);
    int error;
    info=info2=0;
   
    G1=dvector(1,2);
    GG_general=dmatrix(1,2,1,2);
    result_opt=dvector(1,4);
    Y_no=dmatrix(1,n1,1,n2);
    Y_newTriple=dmatrix(1,n1,1,n2);
    M=d3tensor(1,n1,1,n2,1,nrep);
    GG=dmatrix(1,ngroup,1,ngroup);
    YG=dvector(1,ngroup);
    YG2=dmatrix(1,n1,1,ngroup);
    GG_inv=dmatrix(1,ngroup,1,ngroup);
    groupmember=ivector(1,n2);
    totalFisher=dmatrix(1,ngroup+Nbeta,1,ngroup+Nbeta);
    totalFisher_inv=dmatrix(1,ngroup+Nbeta,1,ngroup+Nbeta);
    
    
    sandwich=dmatrix(1,ngroup+Nbeta,1,ngroup+Nbeta);
    sandwich_inv=dmatrix(1,ngroup+Nbeta,1,ngroup+Nbeta);
    
    //change
 for(i=1;i<=n2;i++){
  //groupmember[i]=1+(i==5)+(i==6)+(i==7)+(i==8);
     //groupmember[i]=i;
   groupmember[i]=i;
    }
    
    /* groupmember[1]=2;
     groupmember[2]=1;
     groupmember[3]=1;
     groupmember[4]=1;
     groupmember[5]=2;
     groupmember[6]=1;
     groupmember[7]=1;
     groupmember[8]=1;
     groupmember[9]=2;
     groupmember[10]=2;*/
    
  /*  groupmember[1]=3;
    groupmember[2]=2;
    groupmember[3]=2;
    groupmember[4]=2;
    groupmember[5]=3;
    groupmember[6]=1;
    groupmember[7]=2;
    groupmember[8]=2;
    groupmember[9]=1;
    groupmember[10]=1;
    groupmember[11]=2;
    groupmember[12]=1;
    groupmember[13]=3;
    groupmember[14]=3;*/
    
    
    //interaction clade 34 (acc * bac top1 SNP)
   /* groupmember[1]=1;
    groupmember[2]=1;
    groupmember[3]=2;
    groupmember[4]=1;
    groupmember[5]=2;
    groupmember[6]=1;
    groupmember[7]=2;
    groupmember[8]=1;
    groupmember[9]=1;*/
    
    
     //interaction clade 12 (acc * bac top1 SNP)
    /* groupmember[1]=1;
     groupmember[2]=2;
     groupmember[3]=2;
     groupmember[4]=1;
     groupmember[5]=1;
     groupmember[6]=1;
     groupmember[7]=2;
     groupmember[8]=1;
     groupmember[9]=2;
     groupmember[10]=1;
     groupmember[11]=1;
     groupmember[12]=1;
     groupmember[13]=1;*/
      
    

    
    
   /* groupmember[1]=4;
    groupmember[2]=1;
    groupmember[3]=1;
    groupmember[4]=1;
    groupmember[5]=1;
    groupmember[6]=3;
    groupmember[7]=3;
    groupmember[8]=3;
    groupmember[9]=4;
    groupmember[10]=2;
    groupmember[11]=3;
    groupmember[12]=3;
    groupmember[13]=2;
    groupmember[14]=2;
    groupmember[15]=3;
    groupmember[16]=2;
    groupmember[17]=4;
    groupmember[18]=1;
    groupmember[19]=1;
    groupmember[20]=1;
    groupmember[21]=1;
    groupmember[22]=1;*/
    
   /* groupmember[1]=2;
    groupmember[2]=1;
    groupmember[3]=1;
    groupmember[4]=1;
    groupmember[5]=1;
    groupmember[6]=2;
    groupmember[7]=2;
    groupmember[8]=2;
    groupmember[9]=2;
    groupmember[10]=1;
    groupmember[11]=2;
    groupmember[12]=2;
    groupmember[13]=1;
    groupmember[14]=1;
    groupmember[15]=2;
    groupmember[16]=1;
    groupmember[17]=2;
    groupmember[18]=1;
    groupmember[19]=2;
    groupmember[20]=1;
    groupmember[21]=1;
    groupmember[22]=1;*/
    
    
    membership=ivector(1,n1);
    countmembership=ivector(1,nmember);
    freqmembership=dvector(1,nmember);

    
   // test=dvector(1,n1);
    nonmissing=ivector(1,n1);
    new_nTriple=d3tensor(1,n1,1,n2,1,3);
     new_nComb_nPair=dmatrix(1,n1,1,n2);
    new_nComb=dmatrix(1,n1,1,n2);
    new_nComb_nComb=dmatrix(1,n1,1,n2);
    new_nComb_nComb1=dmatrix(1,n1,1,n2);
    new_nComb_nComb2=dmatrix(1,n1,1,n2);
    new_nComb_nComb3=dmatrix(1,n1,1,n2);
    new_nComb_nComb4=dmatrix(1,n1,1,n2);
    new_nPlant=dvector(1,n1);


    new_nBac=dvector(1,n2);
    new_nBac_1=dvector(1,n2);
    beta_add=dvector(1,Nbeta);
    beta_new=dvector(1,Nbeta);
    combX=d3tensor(1,n1,1,n2,1,Nbeta);
    combY=dmatrix(1,n1,1,n2);
    combY_sq=dmatrix(1,n1,1,n2);
    sigmacombY=dmatrix(1,n1,1,n2);
    double **matrix_phi1=dmatrix(1,n1,1,n1);
    double **matrix_phi2=dmatrix(1,n1,1,n2);
    
    double **matrix_phi3=dmatrix(1,n2,1,n2);
    double **matrix_phi4=dmatrix(1,n2,1,n1);
    double **Fisher_cov=dmatrix(1,Nparams+1,1,Nparams+1);
    double **Fisher_cov_inv=dmatrix(1,Nparams+1,1,Nparams+1);
    
    Vvector=dvector(1,n1);
    phi1_inv=dmatrix(1,n1,1,n1);
    temtest=dvector(1,n1);
    tem_R=dmatrix(1,n1,1,n2);
    M_trans=dmatrix(1,n1,1,n2);
    
    Cov_pair=d3tensor(1,n1,1,n2,1,Nbeta);
    
    for (i=1;i<=n1;i++){
        Cov[i]=d3tensor(1,n2,1,3,1,Nbeta);
        Cov_read[i]=d3tensor(1,n2,1,3,1,Nbeta);
    }
    
   phi1_inv_general=dmatrix(1,n1,1,n1);
    phi2_inv_general=dmatrix(1,n2,1,n2);
phi1_G_inv=dmatrix(1,n1,1,n1);
   phi1_G_inv_member=dmatrix(1,n1,1,nmember);
    member_phi1_G_inv_member=dmatrix(1,nmember,1,nmember);
 member_phi1_G_inv_member_inv=dmatrix(1,nmember,1,nmember);
    
   // Ftol=1e-5;//change!!
    Ftol=1e-2;
    Point=dvector(1,Nparams);
    Delta=dvector(1,Nparams);
    extend_Matrix1=dmatrix(1,n1,1,n1);
    extend_Matrix2=dmatrix(1,n2,1,n2);
    w1=dvector(1,n1);
    u1=dmatrix(1,n1,1,n1);
    u1_tem=dmatrix(1,n1,1,n1);
    extend_matrix1=dmatrix(1,n1,1,n1);
    extend_matrix1_inv=dmatrix(1,n1,1,n1);
    w2=dvector(1,n2);
    u2=dmatrix(1,n2,1,n2);
    u2_tem=dmatrix(1,n2,1,n2);
    extend_matrix2=dmatrix(1,n2,1,n2);
    extend_matrix2_inv=dmatrix(1,n2,1,n2);
    cross=dmatrix(1,n1,1,n2);
    doublekernal=dmatrix(1,n1,1,n2);
    quad1=dmatrix(1,n1,1,n1);
    quad2=dmatrix(1,n1,1,n2);
    matrix1=dmatrix(1,n1,1,n1);
    matrix2=dmatrix(1,n2,1,n2);
    phi1=dmatrix(1,n1,1,n1);
    phi1test=dmatrix(1,n1,1,n1);
    phi2=dmatrix(1,n2,1,n2);
    tem_x=dvector(0,Npeo*Npeo-1);
    tem_m=dvector(0,Npeo*Npeo-1);
    Y=d3tensor(1,n1,1,n2,1,3);
    Y_pair=dmatrix(1,n1,1,n2);
    Y_read=d3tensor(1,n1,1,n2,1,3);
    Y_sq_read=d3tensor(1,n1,1,n2,1,3);
    tem_x1=dvector(0,Nbeta*Nbeta-1);
    tem_m1=dvector(0,Nbeta*Nbeta-1);
    Fisher=dmatrix(1,Nbeta,1,Nbeta);
    Fisher_inv=dmatrix(1,Nbeta,1,Nbeta);
    tem=dmatrix(1,Npeo,1,Nbeta);
    beta=dvector(1,Nbeta);
    mu=d3tensor(1,n1,1,n2,1,3);
    mu_read=d3tensor(1,n1,1,n2,1,3);
    M=d3tensor(1,n1,1,n2,1,3);
    tvector=dvector(1,n1);
    tmatrix=dmatrix(1,n1,1,n1);
    tvector2=dvector(1,n2);
    tmatrix2=dmatrix(1,n2,1,n2);
    nPlant=ivector(1,n1);
    nBac=ivector(1,n2);
    nComb=imatrix(1,n1,1,n2);
    nTriple=i3tensor(1,n1,1,n2,1,3);
    params=dvector(1,Nparams);
    temcombX=dmatrix(1,n1,1,Nbeta);
    temcombX2=dmatrix(1,n2,1,Nbeta);
    
    G=dvector(1,n1);
    
    G_bac=dmatrix(1,2,1,n2);
    
   
    
    combM=dmatrix(1,n1,1,n2);

    MsigmacombM=d4tensor(1,n1,1,n2,1,n1,1,n2);
    MsigmacombX=d3tensor(1,n1,1,n2,1,Nbeta);
    MsigmacombY=dmatrix(1,n1,1,n2);
    XsigmacombX=d4tensor(1,n1,1,n2,1,n1,1,n2);
    XsigmacombY=dmatrix(1,n1,1,Nbeta);

    
    
    
    MsigmacombM_bac=dmatrix(1,n2,1,n2);
    XsigmacombX_bac=dmatrix(1,n2,1,n2);
    MsigmacombY_bac=dvector(1,n2);
    
   
    
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
    
    

    
    if((matrixf1=fopen(matrixfile,"r"))==NULL){
        printf("Cannot open matrix file.\n");
        exit(1);
    }
    
    if((memberf=fopen(memberfile,"r"))==NULL){
        printf("Cannot open matrix test file.\n");
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
    
    
    if((outf1=fopen(outfile1,"w"))==NULL){
        printf("Cannot open output file.\n");
        exit(1);
    }
   if((outf2=fopen(outfile2,"w"))==NULL){
        printf("Cannot open output file.\n");
        exit(1);
    }
    
    if((outf3=fopen(outfile3,"w"))==NULL){
        printf("Cannot open output file.\n");
        exit(1);
    }
    
    if((outf4=fopen(outfile4,"w"))==NULL){
        printf("Cannot open output file.\n");
        exit(1);
    }

    
  // printf("okay?\n");
    
    
    

    
readin();    
    
    
    
    
 // printf("okay?\n");
  
    Point[1]=0.0731253432;//acc
    Point[2]=1.5707643397;//bac
    Point[3]=0.2516760282;//plant
    Point[4]=0.0629946988;//interaction
    
    
    Point[1]=0.0862060002;//acc
    Point[2]=1.2212138002;//bac
    Point[3]=0.2481019143;//plant
    Point[4]=0.0621779429;//interaction
    
    
    
    Delta[1]=0.001;
    Delta[2]=0.1;
    Delta[3]=0.01;
    Delta[4]=0.001;
    initial();
    
   // printf("this is extend matrix 2 again %lf\n",extend_matrix2[1][1]);
     begin=clock();
    
//minimize(Point,function);
    
    //clade1234 logistic no any covaraites (3 regular covaraites, w/o green house fixed effect)
    
    
   /* params[1]=0.0859660544;
     params[2]=1.6701898035;
     params[3]=0.2480102491;
     params[4]=0.0619736787;*/
    
    
    //clade1234 linear no any covaraites (3 regular covaraites, w/o green house fixed effect)
     params[1]=0.0864937724;
     params[2]=1.8199200899;
     params[3]=0.2409041527;
     params[4]=0.0645489199;
    
    //new kinship for bac (more genotypes)
    params[1]=0.0865004127;
    params[2]=1.2614736635;
    params[3]=0.2409026570;
    params[4]=0.0645161743;
    

    
    //clade3_4 no any covaraites (3 regular covaraites, w/o green house fixed effect)
   /* params[1]=0.0739264784;
    params[2]=0.0572234669;
    params[3]=0.2124000060;
    params[4]=0.0457784886;*/
    
    //clade3_4 no any covaraites (4 regular covaraites, no green house fixed effect, b/c confounded with bac. effect, top 1 bac SNP)
   /* params[1]=0.0740864627;
    params[2]=0.0072782063;
    params[3]=0.2128151825;
    params[4]=0.0459052656;*/
    
    //clade3_4 no any covaraites (4 regular covaraites, no green house fixed effect, b/c confounded with bac. effect, top 2 bac SNP)
   /*params[1]=0.0739969514;
     params[2]=0.0186223155;
     params[3]=0.2129354605;
     params[4]=0.0457075262;*/
    
    
    
    //clade3_4 no any covaraites 3 regular covaraites, w/o green house fixed effect
   /*  params[1]=0.0740352190;
    params[2]=0.0656683453;
     params[3]=0.2128472810;
    params[4]=0.0457016078;*/
    
    
     //clade1_2 no any covaraites (3 regular covaraites, no green house fixed effect, b/c confounded with bac. effect); minimize score equation.
    /* params[1]=0.0919693192;
     params[2]=0.3918694517;
     params[3]=0.2711417998;
     params[4]=0.0799565132;*/
    
    
    //clade1_2 plus clade covaraites (4 regular covaraites, no green house fixed effect, b/c confounded with bac. effect); minimize score equation.
    /*params[1]=0.0904600076;//acc
    params[2]=2.0538597004;//bac
    params[3]=0.2711211701;//plant
    params[4]=0.0798219434;//interaction*/
   
    

    
    //clade1 no any covaraites (3 regular covaraites, w.o green house fixed effect
   /* params[1]=0.0904702474;
    params[2]=2.3674847394;
    params[3]=0.2711116177;
    params[4]=0.0798836104;*/

    
    
    //clade1_2 no any covaraites (4 regular covaraites, w./ green house fixed effect
   // params[1]=0.0904600076;
    //params[2]=2.0538597004;
    //params[3]=0.2711211701;
    //params[4]=0.0798219434;
    
    
    //clade1234 + clade 1, 2 indicator+ 3 covaraites.
  /*  params[1]=0.0739201564;
    params[2]=0.4749336845;
    params[3]=0.2441764386;
    params[4]=0.0637307989;*/
   
   
    //- log likelihood = -10804.559859;
    
     /* params[1]=0.055;*/
    //function(params);
minimize(Point,function);
/*for(i=1;i<=500;i++){
        params[2]=1+i*0.0001;
        fprintf(outf4,"%lf\t%lf\n",params[1],function(params));
        fflush(outf4);
    }
  
    fclose(outf4);*/
    
    //rm 6 top country interactions effects
    /*params[1]=0.0805884062;
     params[2]=0.0187878643;
     params[3]=0.2293074392;
     params[4]=0.0593120181;*/
    
    
    
    //rm 1 main effect
  /* params[1]=0.0862517183;
  params[2]=0.0206891709;
    params[3]=0.2293908493;
    params[4]=0.0606406387;*/
    
   // rm 1 main effect + 2 interactions
    /* params[1]=0.0603347901;
     params[2]=0.0203508276;
     params[3]=0.2294869194;
     params[4]=0.0583639910;*/
     //params[4]=0;

    
    // rm only 2 interactions
  /*  params[1]=0.0794730359;
    params[2]=0.0204274073;
    params[3]=0.2294961598;
    params[4]=0.0582670605;*/
    
    
    /*params[3]=0.232676;
    */
   //gwas matrix
    /*params[1]=0.070742;
    params[2]=0.019931;
    params[3]=0.289686;*/
 
 /* params[1]=0;
    params[2]=0;
    params[3]=0;*/

    
   //params[1]=0.0987500000;
/*params[1]=0.0910545426;
     params[2]=0;
     params[3]=	0.1266947147;*/
    
  /* params[1]=0;
    params[2]=0;
    //params[3]=0.3348871;
   // params[3]=0;
    //params[3]=0.384203;
    params[3]=0;*/

   /* params[1]=0;
    params[2]=0.1;
    params[3]=0;
    params[4]=0.1;*/
    

    
function(params);
    
    
    
        
    
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("time spent in total %lf\n", time_spent);
    
      
    for(i=1;i<=n1;i++){
        for(j=1;j<=n1;j++){
            matrix_phi1[i][j]=phi1[i][j];
            for(k=1;k<=n1;k++){
                matrix_phi1[i][j]-=nPlant[k]*phi1[i][k]*extend_matrix1[k][j];
            }   
        }   
    }
    
    
    
    for(i=1;i<=n2;i++){
        for(j=1;j<=n2;j++){
            matrix_phi3[i][j]=phi2[i][j];
            for(k=1;k<=n2;k++){
                matrix_phi3[i][j]-=nBac[k]*phi2[i][k]*extend_matrix2[k][j];
            }
        }     
    }
    
    
    for(i=1;i<=n1;i++){
        for(j=1;j<=n2;j++){
            matrix_phi2[i][j]=0;
            for(k=1;k<=n1;k++){
                for(t=1;t<=n2;t++){
                    matrix_phi2[i][j]-=nComb[k][t]*phi1[i][k]*extend_matrix2[t][j];
                }
            }
        }
    }
    
    
    
    
    for(i=1;i<=n2;i++){
        for(j=1;j<=n1;j++){
            matrix_phi4[i][j]=0;
            for(k=1;k<=n2;k++){
                for(t=1;t<=n1;t++){
                    matrix_phi4[i][j]-=nComb[t][k]*phi2[i][k]*extend_matrix1[t][j];
                }
            }
        }
    }
    
    
    for(i=1;i<=n1;i++){
        for(j=1;j<=n2;j++){
            for(k=1;k<=n1;k++){
                matrix_phi2[i][j]+=nPlant[k]*phi1[i][k]*doublekernal[k][j];
            }      
        }
    }
    
   
    
    
    for(i=1;i<=n2;i++){
        for(j=1;j<=n1;j++){
            for(k=1;k<=n2;k++){
                matrix_phi4[i][j]+=nBac[k]*phi2[i][k]*doublekernal[j][k];
            }     
        }
    }
    
    
    
    
    for(i=1;i<=n1;i++){
        for(j=1;j<=n1;j++){
            for(k=1;k<=n1;k++){
                for(t=1;t<=n2;t++){
                    matrix_phi1[i][j]+=phi1[i][k]*doublekernal[j][t]*nComb[k][t];
                }
            }
        }
    }
    
    
    for(i=1;i<=n2;i++){
        for(j=1;j<=n2;j++){
            
            for(k=1;k<=n2;k++){
                for(t=1;t<=n1;t++){
                    matrix_phi3[i][j]+=phi2[i][k]*doublekernal[t][j]*nComb[t][k];
                }
            } 
        }
    }
    
    
    
    Fisher_cov[1][1]=0;
    

    for(n=1;n<=n1;n++){
        for(i=1;i<=n1;i++){
            Fisher_cov[1][1]+=nPlant[i]*nPlant[n]*matrix_phi1[i][n]*matrix_phi1[n][i];
            for(m=1;m<=n2;m++){
                Fisher_cov[1][1]+= 2*nPlant[n]*nComb[i][m]*matrix_phi1[i][n]*matrix_phi2[n][m];
                for(j=1;j<=n2;j++){
                    Fisher_cov[1][1]+= nComb[i][j]*nComb[n][m]*matrix_phi2[i][m]*matrix_phi2[n][j];
                    
                }
            }
        }
    }     
 
    
    Fisher_cov[2][2]=0;
    for(i=1;i<=n2;i++){
        for(n=1;n<=n2;n++){
            Fisher_cov[2][2]+=nBac[i]*nBac[n]*matrix_phi3[i][n]*matrix_phi3[n][i];
            for(m=1;m<=n1;m++){
                Fisher_cov[2][2]+=2*nBac[n]*nComb[m][i]*matrix_phi3[n][i]*matrix_phi4[n][m];
                for(j=1;j<=n1;j++){
                    Fisher_cov[2][2]+= nComb[j][n]*nComb[m][i]*matrix_phi4[n][m]*matrix_phi4[i][j]; 
                }
            }
        }
    }  
    
    Fisher_cov[3][3]=0;
    for(i=1;i<=n1;i++){
        for(n=1;n<=n1;n++){
            for(j=1;j<=n2;j++){
                for(m=1;m<=n2;m++){
                    Fisher_cov[3][3]+=nComb[i][j]*nComb[n][m]*(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[n][j]-doublekernal[i][m])*(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[n][j]-doublekernal[i][m]);
                }
            }
        }
    }  
    
    Fisher_cov[3][3]+=-2*trace_inv-Npeo;

    
    Fisher_cov[1][2]=0;
    for(i=1;i<=n1;i++){
        for(n=1;n<=n1;n++){
            for(m=1;m<=n2;m++){
                for(j=1;j<=n2;j++){
                    Fisher_cov[1][2]+=nComb[i][j]*nComb[n][m]*(matrix_phi1[i][n]+matrix_phi2[i][m])*(matrix_phi3[m][j]+matrix_phi4[m][i]);
                }
            }
            
        }
    } 
    Fisher_cov[2][1]=Fisher_cov[1][2];
    
   
    
    Fisher_cov[1][3]=trace_invphi1;
    for(i=1;i<=n1;i++){
        for(n=1;n<=n1;n++){
            for(m=1;m<=n2;m++){
                for(j=1;j<=n2;j++){
                    Fisher_cov[1][3]-=nComb[i][j]*nComb[n][m]*(matrix_phi1[i][n]+matrix_phi2[i][m])*(extend_matrix1[n][i]+extend_matrix2[m][j]-doublekernal[n][j]-doublekernal[i][m]);
                }
            }
            
        }
    } 
    Fisher_cov[3][1]= Fisher_cov[1][3];
    
    
    Fisher_cov[2][3]=trace_invphi2;
    for(i=1;i<=n2;i++){
        for(n=1;n<=n2;n++){
            for(m=1;m<=n1;m++){
                for(j=1;j<=n1;j++){
                    Fisher_cov[2][3]-=nComb[j][i]*nComb[m][n]*(matrix_phi3[i][n]+matrix_phi4[i][m])*(extend_matrix1[m][j]+extend_matrix2[n][i]-doublekernal[m][i]-doublekernal[j][n]);
                    
                }
            }
            
            
        }
    }
    
    Fisher_cov[3][2]=Fisher_cov[2][3];
   
    
    
    
    
  /*  error=solve(Fisher_cov_inv,Fisher_cov,3,tem_m,tem_x,&info,&info2);
    if(info!=0){
        printf("fisher non positive!\n");
    }
    if(error!=0){
        printf("fisher non invertable!\n");
    }

    
    //not sure!!
    for(i=1;i<=nrep;i++){
        for(j=1;j<=nrep;j++){
              printf("%.10lf\t",2*sigma_t*Fisher_cov_inv[i][j]);
        
        }
        printf("\n");
    }*/

    printf("here?");
    fprintf(outf1,"%.10lf\t%.10lf\t%.10lf\t",params[1]/(1+params[1]+params[2]+params[3]+params[4]),params[2]/(1+params[1]+params[2]+params[3]+params[4]),params[3]/(1+params[1]+params[2]+params[3]+params[4]));
    for(i=1;i<=Nbeta;i++){
        fprintf(outf1,"%lf\t%lf\t",beta[i],sqrt(Fisher_inv[i][i]*sigma_t/(1+params[1]+params[2]+params[3]+params[4])));
    }
    fprintf(outf1,"%lf\t%lf\n",sigma_t,Output);
  /*fprintf(outf1,"%lf\t%lf\t%lf\t%lf\t%lf\n",bio*sigma_t,2*16*sigma_t*sigma_t/((1+params[1]+params[2])*(1+params[1]+params[2]))*Fisher_cov_inv[1][1],2*16*sigma_t*sigma_t/((1+params[1]+params[2])*(1+params[1]+params[2]))*Fisher_cov_inv[2][2],2*16*sigma_t*sigma_t/((1+params[1]+params[2])*(1+params[1]+params[2]))*Fisher_cov_inv[3][3],Output);*/
    
    //start
   fflush(outf1);
    
    
    

printf("chisq test\n");
        
   // printf("%lf\t%lf\t%lf\n",mu[1][1][1],mu[1][1][2],mu[1][1][3]);
    
    //start
    
    
    
    for(i=1;i<=n1;i++){
        
        for(p=1;p<=Nbeta;p++){
          XsigmacombY[i][p]=0;
           //printf("%d\t%d\t%lf\n",i,p,XsigmacombY[i][p]);
        }
        for(j=1;j<=n2;j++){
            MsigmacombY[i][j]=0;
        }
    }
    
    
       /* for(i=1;i<=n1;i++){
            for(j=1;j<=n2;j++){
                printf("%d\t%d\t%lf\n",i,j,MsigmacombY[i][j]);
            }
            }
      printf("okay again ?\n");*/
    
 for(i=1;i<=n1;i++){
        for(j=1;j<=n2;j++){
          MsigmacombY[i][groupmember[j]]+=sigmacombY[i][j]*M_trans[i][j]/(1+params[4]*new_nComb[i][j])-params[4]*M_trans[i][j]*combY[i][j];
                        
            for(k=1;k<=nrep;k++){
              MsigmacombY[i][groupmember[j]]+=(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k]);
                
                
                for(p=1;p<=Nbeta;p++){
                    
                    
              //   printf("kay%d\t %d\t  %d\t %d\t %d\t %lf\t %lf\t %lf\t %lf\t %lf\n",n1,n2,i,j,p,Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k],Cov_read[i][j][k][p],sigmacombY[i][j],M[i][j][k],XsigmacombY[i][p]);
                    
              XsigmacombY[i][p]+=(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])*Cov_read[i][j][k][p]+sigmacombY[i][j]*M[i][j][k]*Cov[i][j][k][p];
                    
                  //  printf("okay%d\t%d\t%d\t%d\t%d\n",n1,n2,i,j,p);
                    

                }
              
                
                
            }
            

        }
     
    }
            
  //  double *test2=dvector(1,n2);
        
 
    
    
    
   // printf("test is test\n");
    // printf("%lf\t%lf\t%lf\n",test[1],test[2],test[3]);
    /*for(j=1;j<=n2;j++){
        test[j]=0;
        test2[j]=0;
        for(i=1;i<=n1;i++){
            test[j]+=combY[i][j];
            test2[j]+=MsigmacombY[i][j];
            
        }
    }*/

    

    

double **tvectorcross=dmatrix(1,n1,1,ngroup);
    //prepare for retrospective test
  for(i=1;i<=n1;i++){
        tvector[i]=0;
         for(j=1;j<=ngroup;j++){
             tvectorcross[i][j]=MsigmacombY[i][j];
            tvector[i]+=MsigmacombY[i][j];
         }
    }
   
    for(i=1;i<=n1;i++){
        Vvector[i]=0;
      for(j=1;j<=n1;j++){
        Vvector[i]+=phi1[i][j]*tvector[j]; //enrichment vector
    }
    }
    
    //tol=0.0001;
    solve(phi1_inv,phi1, n1,tem_m,tem_x,&info,&info2);
    int q;
    
    
        mean=0;
    VV_general=0;
    
    
    //orthogonal to 1
    for(i=1;i<=n1;i++){
        mean+=tvector[i];
    }
    
    for(i=1;i<=n1;i++){
        tvector[i]-=mean/n1;
    }
    
   
    for(i=1;i<=n1;i++){
        for(j=1;j<=n1;j++){
         VV_general+=tvector[i]*phi1[i][j]*tvector[j];
            phi1_inv_general[i][j]=0;
        }
    }
    
    for(i=1;i<=n2;i++){
        for(j=1;j<=n2;j++){
            phi2_inv_general[i][j]=0;
        }
    }
    
    tol=0.000000001;
    solve_general(phi1_inv_general,phi1, n1,tem_m,tem_x,&info,&info2);
 //////finish
    
    tol=0.000001;
    solve_general(phi2_inv_general,phi2, n2,tem_m,tem_x,&info,&info2);
        
    
       //prepare for prospective test
    
  
         for(n=1;n<=n1;n++){
                for(i=1;i<=n1;i++){
                      for(j=1;j<=ngroup;j++){
                          for(k=1;k<=ngroup;k++){
                              MsigmacombM[n][j][i][k]=0;
                          }
    
                    for(p=1;p<=Nbeta;p++){
                    MsigmacombX[n][j][p]=0;
                    }
    
                      }
                }
         }
        

    for(n=1;n<=n1;n++){
        for(m=1;m<=n2;m++){
            
            combM[n][m]=M_trans[n][m]/(1+params[4]*new_nComb[n][m]);
            
           MsigmacombM[n][groupmember[m]][n][groupmember[m]]-=params[4]/(1+params[4]*new_nComb[n][m])*M_trans[n][m]*M_trans[n][m];
            
            for(p=1;p<=Nbeta;p++){
                MsigmacombX[n][groupmember[m]][p]-=params[4]*M_trans[n][m]*combX[n][m][p];
            }
            
            for(k=1;k<=nrep;k++){
                
                MsigmacombM[n][groupmember[m]][n][groupmember[m]]+=new_nTriple[n][m][k]*M[n][m][k]*M[n][m][k];
                
                
                
                for(p=1;p<=Nbeta;p++){
                    MsigmacombX[n][groupmember[m]][p]+=new_nTriple[n][m][k]*M[n][m][k]*Cov_read[n][m][k][p]*M[n][m][k];
                    
                    
                }
            }
            
             
            
        }
    }
    
    
    double tem_term;
        
    
    for(n=1;n<=n1;n++){
        for(i=1;i<=n1;i++){
            for(m=1;m<=n2;m++){
                for(j=1;j<=n2;j++){
                    
                    tem_term=combM[n][m]*(extend_matrix1[n][i]+extend_matrix2[m][j]-doublekernal[n][j]-doublekernal[i][m]);
                    
                  MsigmacombM[n][groupmember[m]][i][groupmember[j]]-=tem_term*combM[i][j];
                   
                    
                    for(p=1;p<=Nbeta;p++){
                        MsigmacombX[n][groupmember[m]][p]-=tem_term*combX[i][j][p];
                    }
                }
            }
        }
    }
   
    

        

    
    for(i=1;i<=n1;i++){
        for(n=1;n<=n1;n++){
            for(j=1;j<=ngroup;j++){
                for(m=1;m<=ngroup;m++){
    
            XsigmacombX[i][j][n][m]=0;
                    for(k=1;k<=Nbeta;k++){
                        for(p=1;p<=Nbeta;p++){
                            XsigmacombX[i][j][n][m]+=MsigmacombX[n][j][k]*Fisher_inv[k][p]*MsigmacombX[i][m][p];
                        }
                    }
                }
            }
        }
    }
    //#######finish 
      
    
    
    
    
   int gene;
    int chr;
    int pos;
    
 printf("bac chisq test\n");
    
    char letter[MAXLEN];
    char firstallele[MAXLEN];
    
    //// prepare for prospective test on bacteria snp
for(j=1;j<=n2;j++){
        MsigmacombY_bac[j]=0;
        
        for(i=1;i<=n1;i++){
            MsigmacombY_bac[j]+=MsigmacombY[i][j];
        }
        
        for(m=1;m<=n2;m++){
            MsigmacombM_bac[j][m]=0;
            XsigmacombX_bac[j][m]=0;
            for(i=1;i<=n1;i++){
                for(n=1;n<=n1;n++){
                    MsigmacombM_bac[j][m]+=MsigmacombM[i][j][n][m];
                    XsigmacombX_bac[j][m]+=XsigmacombX[i][j][n][m];
                }
            }
        }
      printf("%d\t%lf\n",j,MsigmacombY_bac[j]);

  }
    
    
    double value=0;
    for(j=1;j<=n2;j++){
         for(i=1;i<=n1;i++){
       value+= MsigmacombY[i][j];
    }
    }
    printf("this is the sum %lf\n",value);
    

    
    
  for(i=1;i<=n2+1;i++){
        fscanf(genof_bac,"%s ", letter);
        
    }
    
   int nonmiss_n2;
    
    //fprintf(outf3,"count \t freq1\t freq2\t stat_pro_num\t stat_pro_dnum\t stat_pro\t stat_retro \t df\t pattern\n");
    fprintf(outf3,"count \t freq1\t freq2\t stat_pro_num\t stat_pro_dnum\t stat_pro\t stat_retro_num\t stat_retro_dnum\t stat_retro \n");
    
    /*for(j=1;j<=2939;j++){
        fscanf(genof_bac,"%d ",&pos);
        for(i=1;i<=n2;i++){
            //below is for letter coding of genotype (original pattern)
            fscanf(genof_bac,"%s ", letter);
        }
    }*/
    
    
   
while(fscanf(genof_bac,"%d ",&pos)==1){
//for(j=1;j<=1;j++){
//fscanf(genof_bac,"%d ",&pos);
        freq1=freq2=0;
    nonmiss_n2=n2;
    
        for(i=1;i<=n2;i++){
            //below is for letter coding of genotype (original pattern)
       /* fscanf(genof_bac,"%s ", letter);
           // printf("%s\n",letter);
            
           if((strcmp(letter,"-")!=0)&(freq2==0)&(strcmp(letter,"N")!=0)){
                strcpy(firstallele,letter);
                gene=2;
            }
            else if (strcmp(letter,"-")==0){
                gene=1;//deletin
            }
            else if (strcmp(letter,firstallele)==0){
                gene=2;
            }
            else if (strcmp(letter,"N")==0){
                 nonmiss_n2--;
                gene=-9;
            }
            else {
                gene=0;
                //nonmiss_n2--;
            }*/
 
            //printf("%d\t",gene);
            //below is for interger coding of gentype (nonduplicate/distinct pattern)
           fscanf(genof_bac,"%d ",&gene);
            gene=gene;
            G_bac[1][i]=(gene==2);//deletion
             G_bac[2][i]=(gene==1);//mutation*/
            
         /*G_bac[1][i]=(gene==1);//deletion
            G_bac[2][i]=(gene==2);//mutation*/
            
                  
          freq1+=G_bac[1][i];//deletion
            freq2+=G_bac[2][i];//mutation
            
        }
        
    
        freq1=freq1/nonmiss_n2;
        freq2=freq2/nonmiss_n2;
        //printf("%d\t%lf\t%lf\n",pos,freq1,freq2);
    
  // if((fmin(freq2,1-freq2)>=0.1)*((fmin(freq1,1-freq1)>=0.1)+(freq1==0))){
        fprintf(outf3,"%d\t%lf\t%lf\t",pos,freq1,freq2);
  //  printf("%lf\n",nonmissing_n2);
        
 chisqtest_bac();
    
chisqtest_bac_retro();
 //   }
    
      //  fflush(outf3);
    
}
    
    fclose(outf3);
    
    /*##################### association test of acc snps ###########*/
    
    
    
    

    
    
    printf("acc chisq test\n");

   
    double cor=0;

    //#pragma omp parallel for
  
    /*for(i=1;i<=n1+1;i++){
        
        fscanf(genof,"%s ", letter);
    }
    
    printf("okay\n");*/
    fprintf(outf4,"chr\t pos\t chisq_pro\t chisq_retro\n");
    
  while(fscanf(genof,"%d %d ",&chr,&pos)==2){
//        fscanf(genof,"%d %d ",&chr,&pos);
//for(j=1;j<=1;j++){
//fscanf(genof,"%d %d ",&chr,&pos);
   // while(fscanf(genof,"%d ",&pos)==1){
     freq=0;
    // printf("okay\n");
       
     for(i=1;i<=n1;i++){
         
         fscanf(genof,"%d ", &gene);
         G[i]=0.5*gene;
         freq+=G[i];
         
     }
    
    fscanf(genof,"\n ");
    freq=freq/n1;
   if(freq>0.5)freq=1-freq;

 /*   cor=G[127]/n1;
  
    cor-=freq*1/n1;
  
  
    cor=cor/(sqrt(freq-freq*freq)*sqrt((double)1/n1-1/n1*1/n1));*/
    
   
 
    
 fprintf(outf4,"%d\t%d\t",chr,pos);
     // fprintf(outf4,"%lf\n",cor);
  fprintf(outf4,"%lf\t",chisqtest_strains());
    fprintf(outf4,"%lf\n",chisqtest_retro());
        fflush(outf4);
    
    }

    
    
    
    fclose(outf4);
    
    free_ivector(membership,1,n1);
    free_dvector(beta_add,1,Nbeta);
    free_dvector(beta_new,1,Nbeta);
    free_d3tensor(combX,1,n1,1,n2,1,Nbeta);
    free_dmatrix(combY,1,n1,1,n2);
    free_dmatrix(combY_sq,1,n1,1,n2);
    free_dmatrix(sigmacombY,1,n1,1,n2);
    for (i=1;i<=n1;i++){
        free_d3tensor(Cov[i],1,n2,1,3,1,Nbeta);
         free_d3tensor(Cov_read[i],1,n2,1,3,1,Nbeta);
    }
    
    free_dvector(Point,1,Nparams);
    free_dvector(Delta,1,Nparams);
    free_dmatrix(extend_Matrix1,1,n1,1,n1);
    free_dmatrix(extend_Matrix2,1,n2,1,n2);
    free_dvector(w1,1,n1);
    free_dmatrix(u1,1,n1,1,n1);
    free_dmatrix(u1_tem,1,n1,1,n1);
    free_dmatrix(extend_matrix1,1,n1,1,n1);
    free_dmatrix(extend_matrix1_inv,1,n1,1,n1);
    free_dvector(w2,1,n2);
    free_dmatrix(u2,1,n2,1,n2);
    free_dmatrix(u2_tem,1,n2,1,n2);
    free_dmatrix(extend_matrix2,1,n2,1,n2);
    free_dmatrix(extend_matrix2_inv,1,n2,1,n2);
    free_dmatrix(cross,1,n1,1,n2);
    free_dmatrix(doublekernal,1,n1,1,n2);
    free_dmatrix(quad1,1,n1,1,n1);
    free_dmatrix(quad2,1,n1,1,n2);
    free_dmatrix(matrix1,1,n1,1,n1);
    free_dmatrix(matrix2,1,n2,1,n2);
    free_dmatrix(phi1,1,n1,1,n1);
    free_dmatrix(phi2,1,n2,1,n2);
    free_dvector(tem_x,0,Npeo*Npeo-1);
    free_dvector(tem_m,0,Npeo*Npeo-1);
    free_d3tensor(Y,1,n1,1,n2,1,3);
    free_d3tensor(Y_read,1,n1,1,n2,1,3);
    free_d3tensor(Y_sq_read,1,n1,1,n2,1,3);
    free_dvector(tem_x1,0,Nbeta*Nbeta-1);
    free_dvector(tem_m1,0,Nbeta*Nbeta-1);
    free_dmatrix(Fisher,1,Nbeta,1,Nbeta);
    free_dmatrix(Fisher_inv,1,Nbeta,1,Nbeta);
    free_dmatrix(tem,1,Npeo,1,Nbeta);
    free_dvector(beta,1,Nbeta);
    free_d3tensor(mu,1,n1,1,n2,1,3);
    free_d3tensor(M,1,n1,1,n2,1,3);
    free_dvector(tvector,1,n1);
    free_dmatrix(tmatrix,1,n1,1,n1);
    free_dvector(tvector2,1,n2);
    free_dmatrix(tmatrix2,1,n2,1,n2);
    free_ivector(nPlant,1,n1);
    free_ivector(nBac,1,n2);
    free_imatrix(nComb,1,n1,1,n2);
    free_i3tensor(nTriple,1,n1,1,n2,1,3);
    free_dvector(params,1,Nparams);
    free_dmatrix(temcombX,1,n1,1,Nbeta);
    free_dmatrix(temcombX2,1,n2,1,Nbeta);
    free_ivector(groupmember,1,ngroup);

    fclose(outf1);
    fclose(outf2);
    fclose(outf3);
            
    fclose(matrixf1);
    fclose(matrixf2);
    fclose(covf);
    
}




double chisqtest_new(void){
    
 int i,j,n,k,p;
    double statistics=0;
    double error=0;
    
    for(p=1;p<=ngroup;p++){
        YG[p]=0;
        for(i=1;i<=n1;i++){
            YG[p]+=MsigmacombY[i][p]*G[i];
            YG2[i][p]=MsigmacombY[i][p]*G[i];
        }
        //printf("%d\t%lf\n",p,YG[p]);
    }

   // printf("new\n");
for(p=1;p<=ngroup;p++){
    
     for(k=1;k<=ngroup;k++){
         totalFisher[p][k]=0;
         sandwich[p][k]=0;
         for(i=1;i<=n1;i++){
        sandwich[p][k]+=YG2[i][p]*YG2[i][k];
         }
         //sandwich[p][k]=Y_sq_read[i][j][k]-2*mu[i][j][k]*Y_read[i][j][k]+nTriple[i][j][k]*mu[i][j][k];
        for(n=1;n<=n1;n++){
            for(i=1;i<=n1;i++){
                totalFisher[p][k]+=G[n]*MsigmacombM[n][p][i][k]*G[i];
            }
        }
     }
    
        for(k=1;k<=Nbeta;k++){
            totalFisher[p][ngroup+k]=0;
            sandwich[p][ngroup+k]=0;
            for(i=1;i<=n1;i++){
            sandwich[p][ngroup+k]+=YG2[i][p]*XsigmacombY[i][k];
            }
            sandwich[ngroup+k][p]=sandwich[p][ngroup+k];
            for(n=1;n<=n1;n++){
                totalFisher[p][ngroup+k]+=G[n]*MsigmacombX[n][p][k];
            }
            totalFisher[ngroup+k][p]=totalFisher[p][ngroup+k];
            
        }
    

}

    
    
    error=solve(totalFisher_inv,totalFisher,ngroup+Nbeta,tem_m,tem_x,&info,&info2);
    
    
   // printf("sandwich\n");
  

    
    
    for(p=1;p<=ngroup+Nbeta;p++){
        for(k=1;k<=ngroup+Nbeta;k++){
            sandwich_inv[p][k]=0;
            for(n=1;n<=ngroup+Nbeta;n++){
                for(j=1;j<=ngroup+Nbeta;j++){
            sandwich_inv[p][k]+=totalFisher_inv[p][n]*sandwich[n][j]*totalFisher_inv[j][k];
                }
            }
        }
    }
    


    for(p=1;p<=ngroup;p++){
        for(k=1;k<=ngroup;k++){
            statistics+=YG[p]*sandwich_inv[p][k]*YG[k];
           // printf("%lf\t",totalFisher_inv[p][k]);
        }
   // printf("\n");
    }

   return(statistics);

}
    
    



double chisqtest_bac(void){
    

    
    
    double statistics=0;
    int i,n,k,p;
    double error;
    
    
    for(k=1;k<=2;k++){
        for(p=1;p<=2;p++){
            GG_inv[k][p]=0;
            GG[k][p]=0;
            for(n=1;n<=n2;n++){
                for(i=1;i<=n2;i++){
                    GG[k][p]+=G_bac[k][n]*(MsigmacombM_bac[n][i]-XsigmacombX_bac[n][i])*G_bac[p][i];                }
            }
        }
    }
    
        
    
    for(k=1;k<=2;k++){
        YG[k]=0;
        for(i=1;i<=n2;i++){
            YG[k]+=MsigmacombY_bac[i]*G_bac[k][i];
        }
        
    }
    
    
    if(fmax(fmin(freq1,1-freq1),fmin(freq2,1-freq2))>MIN){//some variation
        
        if(fmin(freq1+freq2,1-freq1-freq2)<MIN){// 2 categories (NA+allele)
            statistics=YG[2]*YG[2]/GG[2][2]*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
           // fprintf(outf3,"%.10lf\t%d\tD_a\t",statistics,1);
        }
        else if (fmin(freq1,1-freq1)<MIN){// 2 categories (two alleles)
            statistics=YG[2]*YG[2]/GG[2][2]*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
         // fprintf(outf3,"%.10lf\t%d\ta_a\t",statistics,1);
        }
        else if (fmin(freq2,1-freq2)<MIN){// 2 categories (NA+allele)
            statistics=YG[1]*YG[1]/GG[1][1]*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
           // fprintf(outf3,"%.10lf\t%d\tD_a\t",statistics,1);
            
        }
        else{//three categories
            
            error=solve(GG_inv,GG,2,tem_m,tem_x,&info,&info2);
            
            if(error!=0){
                fprintf(outf3,"%d\t%d\tN\t",0,0);
                return(0);
            }
            else
                
                statistics=0;
            for(p=1;p<=2;p++){
                for(k=1;k<=2;k++){
                    statistics+=YG[p]*GG_inv[p][k]*YG[k]*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
                    
                }
            };
           // fprintf(outf3,"%.10lf\t%d\tD_a_a\t",statistics,2);
        }
    }
    else
        
    {
        fprintf(outf3,"%d\t%d\tN\t",0,0);
        statistics=0;
    }
    
    
//printf("this is stat %lf\t%lf\t%lf\t%.15lf\n",YG[2],GG[2][2],sigma_t/(1+params[1]+params[2]+params[3]+params[4]),statistics);
   fprintf(outf3,"%lf\t%lf\t%lf\t",YG[2],GG[2][2],statistics);
    return(statistics);
    
}


double chisqtest_bac_retro(void){
    
    double statistics=0;
    int i,n,k,p;
    double error;
        
    for(k=1;k<=2;k++){
        YG[k]=0;
        
        G1[k]=0;
        for(p=1;p<=2;p++){
            GG_general[k][p]=0;
            GG_inv[k][p]=0;
        }
        for(i=1;i<=n2;i++){
            YG[k]+=MsigmacombY_bac[i]*G_bac[k][i];
        }
        
    }
    
    
    double retro=0;
    for(i=1;i<=n2;i++){
        for(k=1;k<=n2;k++){
            retro+=MsigmacombY_bac[i]*phi2[i][k]*MsigmacombY_bac[k];
        }
    }
    
    
    double norm=0;
    int j,t;
    
    
  
  for(i=1;i<=n2;i++){
        for(j=1;j<=n2;j++){
            norm+=phi2_inv_general[i][j];
            
            for(p=1;p<=2;p++){
            G1[p]+=G_bac[p][i]*phi2_inv_general[i][j];
                           for(t=1;t<=2;t++){
            GG_general[p][t]+=G_bac[p][i]*phi2_inv_general[i][j]*G_bac[t][j];
            
        }
    }
          
      }
    }

       for(p=1;p<=2;p++){
                for(t=1;t<=2;t++){
                    GG_general[p][t]-=G1[p]*G1[t]/norm;
                }
       }
    
    //retro 2
    
    
   /*for(p=1;p<=2;p++){
        for(t=1;t<=2;t++){
            GG_general[p][t]= G1[t]/norm*(p==t)-G1[p]/norm*G1[t]/norm;
        }
    }*/
    
    /*for(i=1;i<=n2;i++){
        printf("%d\t%lf\t%lf\n",i,G_bac[1][i],G_bac[2][i]);
    }*/
    
    //printf("%lf\t%lf\t%lf\t%lf\n",G1[1]/norm,G1[2]/norm,freq1,freq2);

   // printf("%lf\t%lf\n",freq1,freq2);
    
    if(fmax(fmin(freq1,1-freq1),fmin(freq2,1-freq2))>MIN){//some variation
        
        if(fmin(freq1+freq2,1-freq1-freq2)<MIN){// 2 categories (NA+another allele)
            
            //retro
        statistics=YG[2]*YG[2]/(retro*GG_general[2][2])*(n2-1);
            // fprintf(outf3,"%.10lf\t%d\tD_a\n",statistics,1);
           
           
             //retro 2
            /* for(p=1;p<=2;p++){
                    for(t=1;t<=2;t++){
                        GG_general[p][t]= G1[t]/norm*(p==t)-G1[p]/norm*G1[t]/norm;
                    }
                }
        
          statistics=YG[2]*YG[2]/(retro* G1[2]/norm*(1-G1[2]/norm));
            fprintf(outf3,"%.10lf\t%d\tD_a\n",statistics,1);*/
            

     
        }
        else if (fmin(freq1,1-freq1)<MIN){// 2 categories (two alleles)
            
               
                 //retro
               statistics=YG[2]*YG[2]/(retro*GG_general[2][2])*(n2-1);
           // fprintf(outf3,"%.10lf\t%d\ta_a\n",statistics,1);
            

            
            //retro 2
            /* for(p=1;p<=2;p++){
                    for(t=1;t<=2;t++){
                        GG_general[p][t]= G1[t]/norm*(p==t)-G1[p]/norm*G1[t]/norm;
                    }
                }
            
            
            statistics=YG[2]*YG[2]/(retro*G1[2]/norm*(1-G1[2]/norm));
             fprintf(outf3,"%.10lf\t%d\ta_a\n",statistics,1);*/
            
            
        }
        else if (fmin(freq2,1-freq2)<MIN){// 2 categories (NA+another allele)
            
             
           
                //retro
            statistics=YG[1]*YG[1]/(retro*GG_general[1][1])*(n2-1);
            // fprintf(outf3,"%.10lf\t%d\tD_a\n",statistics,1);
           
                //retro 2
              /*for(p=1;p<=2;p++){
                    for(t=1;t<=2;t++){
                        GG_general[p][t]= G1[t]/norm*(p==t)-G1[p]/norm*G1[t]/norm;
                    }
                }
            
          statistics=YG[2]*YG[2]/(retro*G1[1]/norm*(1-G1[1]/norm));
            fprintf(outf3,"%.10lf\t%d\tD_a\n",statistics,1);*/
            
            
        }
        else{//three categories
            
           error=solve(GG_inv,GG_general,2,tem_m,tem_x,&info,&info2);
            
            
            if(error!=0){
                fprintf(outf3,"%d\t%d\tN\n",0,0);
                return(0);
            }
            else statistics=0;
            
            
            for(p=1;p<=2;p++){
                for(k=1;k<=2;k++){
                    //retro
                   statistics+=YG[p]*GG_inv[p][k]*YG[k]/retro*(n2-1);
                }
            }
       //  fprintf(outf3,"%.10lf\t%d\tD_a_a\n",statistics,2);
            
            //retro 2
           /*  statistics=0;
                
                for(p=1;p<=2;p++){
                    for(t=1;t<=2;t++){
                        GG_general[p][t]= G1[t]/norm*(p==t)-G1[p]/norm*G1[t]/norm;
                    }
                }
                
                error=solve(GG_inv,GG_general,2,tem_m,tem_x,&info,&info2);
                for(p=1;p<=2;p++){
                    for(k=1;k<=2;k++){

                    //retro 2
                     statistics+=YG[p]*GG_inv[p][k]*YG[k]/retro;
              
            }
            
                }
          
          //fprintf(outf3,"%.10lf\t%d\tD_a_a\n",statistics,2);*/
            
        }
    }
    else
        
    {
        fprintf(outf3,"%d\t%d\tN\n",0,0);
        statistics=0;
    }
    
    
    // statistics=YG[2]*YG[2]/(retro*G1[2]/norm*(1-G1[2]/norm));
//printf("this is retro stat %lf\t%lf\t%lf\t%.15lf\n",YG[2],GG_general[2][2],retro,statistics);
  //  printf("this is stat %.15lf\n",statistics);
   fprintf(outf3,"%lf\t%lf\t%lf\n",YG[2],GG_general[2][2]/(n2-1),statistics);
    return(statistics);
    
}



double chisqtest_strains(void){
  //  printf("hello\n");
    double statistics=0;
   // double GG=0,YG=0;
    int i,n,k,p;
    double error;
    
   
    // printf("\n");
    for(p=1;p<=ngroup;p++){
        YG[p]=0;
         for(k=1;k<=ngroup;k++){
             GG[p][k]=0;
             GG_inv[p][k]=0;
             for(n=1;n<=n1;n++){
                 for(i=1;i<=n1;i++){
              GG[p][k]+=G[n]*(MsigmacombM[n][p][i][k]-XsigmacombX[n][p][i][k])*G[i];
                 //  GG[p][k]+=G[n]*(XsigmacombX[n][p][i][k])*G[i];
                    // printf("%d\t%d\t%d\t%d\t%lf\n",n,p,i,k,MsigmacombM[n][p][i][k]);
                 }
             }
     //  printf("%d\t%d\t%d\t%d\t%lf\n",1,p,1,k,MsigmacombM[1][p][1][k]);
        }
   
      for(i=1;i<=n1;i++){
            YG[p]+=MsigmacombY[i][p]*G[i];
        }
     
    }
    
  
    
    
    
    
   
    error=solve(GG_inv,GG,ngroup,tem_m,tem_x,&info,&info2);
    
 
     for(p=1;p<=ngroup;p++){
                 for(k=1;k<=ngroup;k++){
         statistics+=YG[p]*GG_inv[p][k]*YG[k]*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
                     }
        
     }
    
  //printf("%lf\t%lf\t%lf\n",YG[1],GG_inv[1][1],statistics);
      return(statistics);
    
}



/*double chisqtest(void){
    
    
    double GG=0,YG=0;
    int i,n,k,p;
    
    //double freq=0;
    
    for(n=1;n<=n1;n++){
         //printf("%lf\t",G[n]);
        for(i=1;i<=n1;i++){
           GG+=G[n]*(MsigmacombM[n][i]-XsigmacombX[n][i])*G[i];

            // freq+=G[t][n];
        }
    }
    
   // printf("this is GG %lf\n",GG);
    //freq
    
    
    if(GG<1e-5){
        printf("error!not enough information!  %lf\n",GG);
        
        for(k=1;k<=Nbeta;k++){
            for(i=1;i<=Nbeta;i++){
                printf("%lf\t",Fisher_inv[k][i]);
            }
            printf("\n");
        }
        return(-1);
    }
    
    
    for(i=1;i<=n1;i++){
        YG+=MsigmacombY[i]*G[i];
        
    }
    
     // printf("YG is %lf\n",YG);
    
    
     //printf("YG is %lf\t",YG);
    
    //printf("final are %lf\n",YG);
    
   // printf("%lf\t%lf\t%lf\t%lf\n",YG,GG,sigma_t/(1+params[1]+params[2]+params[3]),YG*YG*(1+params[1]+params[2]+params[3])/(GG*sigma_t));
    
   // fflush(outf3);
    
    
    YG=YG*YG*(1+params[1]+params[2]+params[3])/(GG*sigma_t);
    
    return(YG);

}*/
    

  /*   freq=0;
    for(i=1;i<=nmember;i++){
       freqmembership[i]=0;
       countmembership[i]=0;
    }
    
    
    for(i=1;i<=n1;i++){
        if(G[i]!=-4.5){
            ncount++;
            nonmissing[ncount]=i;
            freq+=G[i];
            freqmembership[membership[i]]+=G[i];
            countmembership[membership[i]]++;
        }
    }
    
    freq=freq/ncount;
    for(i=1;i<=nmember;i++){
        freqmembership[i]=(freqmembership[i]+0.5)/(countmembership[i]+1);
        //printf("%d\t%lf\n", i,freqmembership[i]);
    }
    
   for(i=1;i<=ncount;i++){
        for(j=1;j<=ncount;j++){
            phi1_G_inv[i][j]=phi1_inv[nonmissing[i]][nonmissing[j]];
        }
    }
    
    
    //printf("this is GG %lf\n",phi1_G_inv[1][1]);
   for(i=1;i<=ncount;i++){
        for(j=1;j<=ncount;j++){
           G1+=(G[nonmissing[i]]-freqmembership[membership[nonmissing[i]]])*phi1_G_inv[i][j];
            GG+=(G[nonmissing[i]]-freqmembership[membership[nonmissing[i]]])*phi1_G_inv[i][j]*(G[nonmissing[j]]-freqmembership[membership[nonmissing[j]]]);
            
            norm+=phi1_G_inv[i][j];
        }
    }
    
   
     GG=GG-G1*G1/norm;
    GG=GG/(ncount-nmember);
      
  
    for(i=1;i<=ncount;i++){
     test[nonmissing[i]]=tvector[nonmissing[i]];
    }

  
    
        for(i=1;i<=ncount;i++){
       VG+=(G[nonmissing[i]]-freqmembership[membership[nonmissing[i]]])*test[nonmissing[i]];
            for(j=1;j<=ncount;j++){

                VV+=test[nonmissing[i]]*phi1[nonmissing[i]][nonmissing[j]]*test[nonmissing[j]];
            }
    }
    

    
//printf("%lf\t %lf\t %lf\t%lf\n",VG,VV,GG,VG*VG/(GG*VV));
    
 VG=VG*VG/(GG*VV);
    
    return(VG);   */

   
  
//}


double chisqtest_retro(void){
    
    
    int i,j,n,m,k,p,q;
    
    double VV=0,GG=0,VG=0,G1=0,V1=0,norm=0,VG_general=0,GG_general=0;
    int ncount=0;
    
    
    freq=0;    
    
    
    for(i=1;i<=n1;i++){
        freq+=G[i];
        VG_general+=G[i]*tvector[i];
        for(j=1;j<=n1;j++){
             
           norm+=phi1_inv_general[i][j];
             G1+=G[i]*phi1_inv_general[i][j];
            GG_general+=G[i]*phi1_inv_general[i][j]*G[j];         
           
            
            //test
       /*   norm+=(i==j);
             G1+=G[i]*(i==j);
             GG_general+=G[i]*(i==j)*G[j];*/
           
            
        }
    }
    
   /* for(i=1;i<=5;i++){
        for(j=1;j<=5;j++){
            printf("%lf\t",phi1_inv_general[i][j]);
        }
        printf("\n");
    }*/
    
   // printf("negative %lf\t%lf\n?",GG_general,G1*G1/norm);
    
    GG_general-=G1*G1/norm;
// printf("VG general is %lf\t%lf\t%lf\t%lf\t%lf\n",VG_general,mean,VG_general-mean*G1/norm,GG_general/(n1-nmember),VG_general*VG_general*(n1-nmember)/(GG_general*VV_general));
    freq=freq/n1;
    
    
    if(freq>0.5) freq=1-freq;
    //VG=VG_general*VG_general*(n1-nmember)/(GG_general*VV_general);
    //if(freq*n1>=5){
    VG=VG_general*VG_general*(n1-nmember)/(GG_general*VV_general);
    //}
    //else{
      //  VG=VG_general*VG_general/(G1/norm*(1-G1/norm)*VV_general);
   // }

   // fprintf(outf4,"%lf\t%lf\t%lf\t",freq,GG_general/(n1-nmember),G1/norm*(1-G1/norm));
    
    //printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",freq,GG_general/(n1-nmember),G1,norm,G1/norm,G1/norm*(1-G1/norm));
    return(VG);
    

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
        
        error+=solve(extend_matrix2_inv,extend_matrix2,n2,tem_m,tem_x,&info,&info2);
        
        
        
        
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
        
        
        error+=solve(extend_matrix1_inv,extend_matrix1,n1,tem_m,tem_x,&info,&info2);
             
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


void minimize(double *start, double (*func)(double *)){
    register int i, j, k;
    int neval=0;	/*Number of function evaluations */
    double **point;	/* Each row is a point in the parameter space. The 
                     number of points must be one more than the 
                     dimension of the space */
    double *functionvalue;	/* Objective function for each point in parameter 
                         space */
    void amoeba(double **p, double y[], int ndim, double ftol,
                double (*funk)(double []), int *nfunk);
  
    point=dmatrix(1, Nparams+1, 1, Nparams);
    functionvalue= dvector(1, Nparams+1);
    
    /* Do amoeba NTIMES */
    for (k=1; k<=NTIMES; k++) {
        
        
        
        /* Find starting points for the amoeba */
        for (j=1; j<=Nparams; j++) {
            point[1][j] = start[j];
            //printf("%lf\n",point[1][j]);
            fflush(stdout);
        }
        
        for (i=2; i<=Nparams+1; i++){
            for (j=1; j<=Nparams; j++){
                point[i][j] = point[1][j];
                if (j == i-1) {
                    //printf("this is j==i-1 and Delta j %d, %lf\n, point %lf",j, Delta[j],point[i][j]); 
                    point[i][j] += Delta[j];
                    //printf("after point %lf",point[i][j]);
                }
                //printf("first %lf\t", point[i][j]);
            }
            //printf("\n");
            fflush(stdout);
        }

        
        for (i=1; i<=Nparams+1; i++)
            functionvalue[i] = func(point[i]);
        
        /* Run the amoeba */
        /*printf("amoeba\n");
         fflush(stdout);*/
        amoeba(point, functionvalue, Nparams, Ftol, func, &neval);
        

        /*Reinitialize starting point to current minimization point */
        for (j=1; j<=Nparams; j++) {
            start[j] = point[1][j];
        }
        
       fprintf(outf2,"neval %d start over\n\n", neval);
        
    }
    
    
    
    Output = functionvalue[1];
    params[1]=point[1][1];
    params[2]=point[1][2];
    params[3]=point[1][3];
    params[4]=point[1][4];
    free_dvector(functionvalue,1,Nparams+1);
    free_dmatrix(point,1, Nparams+1, 1, Nparams);
}



#define NMAX 5000
#define GET_PSUM \
for (j=1;j<=ndim;j++) {\
for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void amoeba(double **p, double y[], int ndim, double ftol,
            double (*funk)(double []), int *nfunk)
{
	double amotry(double **p, double y[], double psum[], int ndim,
                  double (*funk)(double []), int ihi, double fac);
	int i,ihi,ilo,inhi,j,mpts=ndim+1;
	double rtol,sum,swap,ysave,ytry,*psum;
    
	psum=dvector(1,ndim);
	*nfunk=0;
	GET_PSUM
	for (;;) {
        //printf("%lf\t%lf\t%lf\n",y[1],y[2],y[3]);
		ilo=1;
		ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
		for (i=1;i<=mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
        
       // if ((fabs(y[ihi])+fabs(y[ilo]))==0)
         //   rtol=0;
        //else{
		//rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
            rtol=fabs(y[ihi]-y[ilo]);
       // }
       
		fprintf(stderr,"\nTest tolerance.\n");
       // if (y[ihi]>1e98)
		fprintf(stderr,"y_hi=%E\ny_lo=%.12f\nrtol=%E\n\n",y[ihi], y[ilo], rtol);

       //else
            fprintf(stderr,"y_hi=%.12f\ny_lo=%.12f\nrtol=%.12f\n",y[ihi], y[ilo], rtol);
    
	fflush(stderr);
        
		if (rtol < ftol) {
			SWAP(y[1],y[ilo])
			for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
                break;
		}
		if (*nfunk >= NMAX) nrerror("NMAX exceeded");
		*nfunk += 2;
		ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
		if (ytry <= y[ilo])
			ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
			if (ytry >= ysave) {
				for (i=1;i<=mpts;i++) {
					if (i != ilo) {
						for (j=1;j<=ndim;j++)
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						y[i]=(*funk)(psum);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		} else --(*nfunk);
	}
	free_dvector(psum,1,ndim);
}
#undef SWAP
#undef GET_PSUM 
#undef NMAX



double amotry(double **p, double y[], double psum[], int ndim,
              double (*funk)(double []), int ihi, double fac)
{
	int j;
	double fac1,fac2,ytry,*ptry;
    
	ptry=dvector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_dvector(ptry,1,ndim);
	return ytry;
}
#undef NRANSI



double function(double *params){
        
       
        int error,i,j,n,m,k,p;
        double result,delta,delta2,delta_new,delta2_new;
    double original_N=0;
    
    
    //restriction
// params[1]=0;
 //params[2]=0;
 // params[3]=0;
    
   // params[4]=0;
   
  if (params[1]<0) {
      //if (params[1]<0) {
        
      
      printf("hit the boundary! parameter is %lf\t%lf\t%lf\t%lf\n",params[1],params[2],params[3],params[4]);
      return(MAX);
      //params[1]=-params[1];
    }
    
   
   if (params[2] < 0) {
       // if (params[2]<0) {

       
       printf("hit the boundary! parameter is %lf\t%lf\t%lf\t%lf\n",params[1],params[2],params[3],params[4]);
       return(MAX);
        params[2]=-params[2];
      
    }
    
    if (params[3] < 0) {
        // if (params[2]<0) {
        printf("hit the boundary! parameter is %lf\t%lf\t%lf\t%lf\n",params[1],params[2],params[3],params[4]);
         return(MAX);
        params[3]=-params[3];
        
    }
    
    if (params[4] < 0) {
        // if (params[2]<0) {
     printf("hit the boundary! parameter is %lf\t%lf\t%lf\t%lf\n",params[1],params[2],params[3],params[4]);
         return(MAX);
        params[4]=-params[4];
        
    }

    log_det=0;
    new_N=0;
    new_N_add=0;
    for (i=1;i<=n1;i++){
        new_nPlant[i]=0;

        for (j=1;j<=n2;j++){
            new_nComb[i][j]=0;
            new_nComb_nComb1[i][j]=0;
            new_nComb_nComb2[i][j]=0;
            new_nComb_nComb3[i][j]=0;
          
            for (k=1;k<=nrep;k++){
               
                log_det+=log(1+params[3]*nTriple[i][j][k]);
                
                new_nTriple[i][j][k]=nTriple[i][j][k]/(1+params[3]*nTriple[i][j][k]);
                
                new_nComb[i][j]+= new_nTriple[i][j][k];
                
                new_nComb_nComb1[i][j]+=new_nTriple[i][j][k]/(1+params[3]*nTriple[i][j][k]);
                
                new_nComb_nComb2[i][j]+=nTriple[i][j][k]*new_nTriple[i][j][k];
                
                new_nComb_nComb3[i][j]+=new_nTriple[i][j][k]*new_nTriple[i][j][k];
                
                new_nComb_nComb4[i][j]+=new_nTriple[i][j][k]*new_nTriple[i][j][k]*new_nTriple[i][j][k];
           
               
            }
            
            log_det+=log(1+params[4]*new_nComb[i][j]);
            
            original_N+=new_nComb[i][j];
            
            
            
            new_nComb_nPair[i][j]=new_nComb[i][j]/(1+params[4]*new_nComb[i][j]);// 1^T J^{-1} 1
            
            
            new_nPlant[i]+=new_nComb_nPair[i][j];


        }
        new_N+=new_nPlant[i];
      
    }
    // printf("\n");
    
    for (j=1;j<=n2;j++){
        new_nBac[j]=0;
        new_nBac_1[j]=0;
        for (i=1;i<=n1;i++){
            new_nBac[j]+=new_nComb_nPair[i][j];
            new_nBac_1[j]+=new_nComb_nComb2[i][j];
            
        }
    }

   
    error=inverse(params);
     
    log_det+=det_Kb+det_Na;
        
    if (error!=0){
        result=1e99;
        printf("%lf\t%lf\t%lf\t%lf\n",params[1],params[2],params[3],params[4]);
        printf("Sigmamat not inveratable %d\n", error);
        return result;
    }
    
        
   fprintf(outf2,"%lf\t%lf\t%lf\t%lf\n",params[1],params[2],params[3],params[4]);
    
        
    
                                                                
                                                                
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            for (k=1;k<=nrep;k++){
                
                Y[i][j][k]=1/(1+params[3]*nTriple[i][j][k])*Y_read[i][j][k];
                
            
                
                for(p=1;p<=Nbeta;p++){
                    Cov[i][j][k][p]=1/(1+params[3]*nTriple[i][j][k])*Cov_read[i][j][k][p];      
                }
            }
        }
    }
    
                                                                
error=coefficient(params);
    
    
    printf("this is log_det %lf\n",log_det);
    

    log_det+=Npeo*log(sigma_t/(1+params[1]+params[2]+params[3]+params[4]));
    printf("this is log_likelihood %lf\t%lf\n",log_det,sigma_t/(1+params[1]+params[2]+params[3]+params[4]));
    
    // not derivative //logistic
//return (log_det);

    
       if(error!=0){
        result=1e99;
        printf("%lf\t%lf\t%lf\t%lf\n",params[1],params[2],params[3],params[4]);
        printf("Fisher matrix on X not positive definite %d\n", error);
    }
    
  
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            Y_newTriple[i][j]=Y_no[i][j]=0;
            
            for (k=1;k<=nrep;k++){
                Y_newTriple[i][j]+=new_nTriple[i][j][k]*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/M[i][j][k];
                Y_no[i][j]+=(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/M[i][j][k];
            }
        }
    }
    
      Y_sq=0;
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
          Y_sq+=Y_no[i][j]*Y_no[i][j]*params[4]/(1+params[4]*new_nComb[i][j])*(params[4]*new_nComb_nComb1[i][j]/(1+params[4]*new_nComb[i][j])-2)+2*params[3]*params[4]/(1+params[4]*new_nComb[i][j])*Y_newTriple[i][j]*Y_no[i][j];
            for (k=1;k<=nrep;k++){
            Y_sq+=-params[3]*(2+params[3]*nTriple[i][j][k])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/(M[i][j][k]*M[i][j][k]);
                
            }
        }
    }
    

    
    
    delta=0;
    delta_new=0;
    
 
        
    quad=Y_sq;
    trace_inv=0;
    
   
   
  
    for (i=1;i<=n1;i++){
        tvector[i]=0;
       
        for (j=1;j<=n2;j++){
            tvector[i]+=combY[i][j]+new_nComb[i][j]/(1+params[4]*new_nComb[i][j])*sigmacombY[i][j];
            
            for(k=1;k<=nrep;k++){
                
                quad+=(Y_sq_read[i][j][k]+mu_read[i][j][k]*(nTriple[i][j][k]*mu_read[i][j][k]-2*Y_read[i][j][k]))/(M[i][j][k]*M[i][j][k]);
                
            quad+=2*(1/(1+params[3]*nTriple[i][j][k])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/M[i][j][k]-new_nTriple[i][j][k]/(1+params[3]*nTriple[i][j][k])*params[4]*combY[i][j])*sigmacombY[i][j]/(1+params[4]*new_nComb[i][j]);
                
             quad+=new_nTriple[i][j][k]/((1+params[4]*new_nComb[i][j])*(1+params[4]*new_nComb[i][j])*(1+params[3]*nTriple[i][j][k]))*sigmacombY[i][j]*sigmacombY[i][j];
                              
        delta_new+=((Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/M[i][j][k]-params[4]*new_nTriple[i][j][k]*combY[i][j]+new_nTriple[i][j][k]/(1+params[4]*new_nComb[i][j])*sigmacombY[i][j])*((Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/M[i][j][k]-params[4]*new_nTriple[i][j][k]*combY[i][j]+new_nTriple[i][j][k]/(1+params[4]*new_nComb[i][j])*sigmacombY[i][j]);
                

            }
        
        
        }
        
        
        
    }
    
   
    //printf("this is delta_new is %lf\t quad is %lf \n",delta_new,quad);
    delta-=quad;
    delta2=0;
    
        
    
    for(i=1;i<=n1;i++){
          delta2+=new_nPlant[i]*phi1[i][i];
               
       for(j=1;j<=n2;j++){
           

           trace_inv-=(new_nComb_nComb1[i][j]+params[4]/(1+params[4]*new_nComb[i][j])*(params[4]*new_nComb_nComb1[i][j]/(1+params[4]*new_nComb[i][j])-2)*new_nComb[i][j]*new_nComb[i][j]+2*params[4]*params[3]/(1+params[4]*new_nComb[i][j])*new_nComb_nComb3[i][j]*new_nComb[i][j])*(extend_matrix1[i][i]+extend_matrix2[j][j]-2*doublekernal[i][j]);
           
           trace_inv-=params[3]*new_nComb[i][j]+params[4]/(1+params[4]*new_nComb[i][j])*new_nComb_nComb1[i][j];
            }
        
     
        for(n=i;n<=n1;n++){
             delta+=(2-(i==n))*phi1[i][n]*tvector[i]*tvector[n];
             tmatrix[i][n]=0;
            
            for(j=1;j<=n2;j++){
                 for(m=1;m<=n2;m++){
                    tmatrix[i][n]+=(new_nComb_nPair[i][j]*new_nComb_nPair[n][m])*(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[i][m]-doublekernal[n][j]);
                 }
            }
                    
                   delta2-=(2-(i==n))*phi1[n][i]*tmatrix[i][n];
        }
        
        
    }

    //printf("delta 2 is %lf\n",delta2);
    
    printf("delta is %lf\t%lf\t%lf\n",delta,(1+params[1]+params[2]+params[3]+params[4]),sigma_t);
    

    
    
    delta=delta*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
    
    

    
    trace_inv+=Npeo;
    
   // printf("trace_inv is %lf\n",trace_inv);
    
    
    trace_invphi1=delta2;

    delta2-=trace_inv;
        
   result=fabs(delta2-delta)*(1+params[1]+params[2]+params[3]+params[4]);
    
    result_opt[1]=abs(-delta+delta2);
    printf("first result %lf\t%lf\t%lf\n",delta,delta2,-delta+delta2);
    
    //return(-delta+delta2);
    //restriction
    
   // return(log_det);
 //return result;
    
      
  delta=0;
    
delta2_new=original_N;
    //printf("oringal %lf\n",original_N);
  for (j=1;j<=n2;j++){
        tvector2[j]=0;
       
        for (i=1;i<=n1;i++){
            tvector2[j]+=combY[i][j]+new_nComb[i][j]/(1+params[4]*new_nComb[i][j])*sigmacombY[i][j];
           
            
            delta2_new+=-params[4]/(1+params[4]*new_nComb[i][j])*new_nComb_nComb3[i][j];
            
        delta2_new-=new_nComb_nComb3[i][j]/((1+params[4]*new_nComb[i][j])*(1+params[4]*new_nComb[i][j]))*(extend_matrix1[i][i]+extend_matrix2[j][j]-2*doublekernal[i][j]);
            
            
            
            
            

        }
  }
    
   // printf("delta2 new %lf\n",delta2_new);
    delta-=quad;

    
    
    
    
    delta2=0;
   

    for (j=1;j<=n2;j++){
        delta2+=new_nBac[j]*phi2[j][j];
 
        for (m=j;m<=n2;m++){
            
           delta+=(2-(j==m))*phi2[j][m]*tvector2[j]*tvector2[m];
            //delta+=phi2[j][m]*tvector2[j]*tvector2[m];
            tmatrix2[j][m]=0;
    
            for (i=1;i<=n1;i++){
                for(n=1;n<=n1;n++){
                    tmatrix2[j][m]+=(new_nComb_nPair[i][j]*new_nComb_nPair[n][m])*(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[i][m]-doublekernal[n][j]);
                }
    
            }
           delta2-=(2-(j==m))*phi2[m][j]*tmatrix2[j][m];
           // delta2-=phi2[m][j]*tmatrix2[j][m];
        }
    }
    
    
    
    
    delta=delta*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
    
    

     trace_invphi2=delta2;
     delta2-=trace_inv;
    
    //restriction
  result+=fabs(delta2-delta)*(1+params[1]+params[2]+params[3]+params[4]);
    
   
        
    printf("second result %lf\t%lf\t%lf\n",delta,delta2,-delta+delta2);
    
    result_opt[2]=abs(-delta+delta2);
    //restriction
// return (-delta+delta2);
     // printf("delta new before is %lf\n",delta_new);

     delta_new-=quad;
   
   // printf("delta new is %lf\t%lf\n",delta_new,quad);
    
     delta_new=delta_new*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
         
      delta2_new-=trace_inv;
    
    
    
    
     result+=fabs(delta2_new-delta_new)*(1+params[1]+params[2]+params[3]+params[4]);
     
   printf("third result %lf\t%lf\t%lf\n",delta_new,delta2_new,delta_new-delta2_new);
    
    result_opt[3]=abs(delta_new-delta2_new);
    //restriction
   // return(result);
    
    //printf("result %lf\n",result);
    //restriction
   // return(result);
   

  /* fprintf(outf2,"%lf\t",(delta2_new-delta_new)*(1+params[1]+params[2]+params[3]+params[4]));
    fflush(outf2);*/
    
    
    
    double delta_new_new=0;
    
    for (i=1;i<=n1;i++){
        
        for (j=1;j<=n2;j++){
            
            delta_new_new+=(combY[i][j]+new_nComb_nPair[i][j]*sigmacombY[i][j])*(combY[i][j]+new_nComb_nPair[i][j]*sigmacombY[i][j]);
        }
    }
    
    
    
    delta_new_new-=quad;
    
   // printf("delta_new_new is %lf\n",delta_new_new);
    
    delta_new_new=delta_new_new*(1+params[1]+params[2]+params[3]+params[4])/sigma_t;
    
    
 
    double delta2_new_new=new_N;
    
  

    for (j=1;j<=n2;j++){
               for (i=1;i<=n1;i++){            
                  delta2_new_new-=new_nComb_nPair[i][j]*new_nComb_nPair[i][j]*(extend_matrix1[i][i]+extend_matrix2[j][j]-2*doublekernal[i][j]);
                   
        }
    }
    
    delta2_new_new-=trace_inv;
    
    
    
 result+=fabs(delta2_new_new-delta_new_new)*(1+params[1]+params[2]+params[3]+params[4]);
    
 
    printf("fourth result %lf\t%lf\t%lf\n",delta_new_new,delta2_new_new,delta_new_new-delta2_new_new);
    
    result_opt[4]=abs(delta_new_new-delta2_new_new);
    
    printf("result %lf\n",result);
    

    
    

    fprintf(outf2,"%lf\n",result);
    fflush(outf2);
    
// return (log_det);
    // derivative //logistic
    return result;
    
}


int coefficient(double *params){
    
    double delta=0.5;
    double thresh=1.0e-5;
    // double thresh=1.0e-2;
    
    int i,j,k,p,t,m,n;
    
    
    
    for (p=1;p<=Nbeta;p++){
        beta[p]=0;
       
    }
    
           
        
    
    int times=0;
    int atleast=0;
    while((delta>thresh)||(atleast<2)){
       // delta=0;
        //printf("this is atleast %d\t",atleast);
        //if(atleast2)
        atleast++;
     //  for(times=1;times<=20;times++){
        if(times>=20) {
            printf("no convergence in %d iterations\n",times);
         //   fprintf(outf3,"no convergence in %d iterations\n",times);
            break;
        }
        times++;
    
        for (i=1;i<=n1;i++){
            for (j=1;j<=n2;j++){
                
                
                M_trans[i][j]=0;
                combY[i][j]=0;
                combY_sq[i][j]=0;
                for (k=1;k<=nrep;k++){
                    
                    
                    
                    mu[i][j][k]=0;
                    mu_read[i][j][k]=0;
                   
                    for (p=1;p<=Nbeta;p++){
                        mu_read[i][j][k]+=Cov_read[i][j][k][p]*beta[p];
                     // printf("this is beta %lf\n",beta[p]);

                    }
                    
                     M[i][j][k]=1;
                   
                   
                     /*if(mu[i][j][k]>10){
                     mu[i][j][k]=1;
                     M[i][j][k]=0.01;
                     }
                     if(mu[i][j][k]<-10){
                     mu[i][j][k]=0;
                     M[i][j][k]=0.01;
                     }
                     else{*/
                    //logistic
  // mu_read[i][j][k]=exp(mu_read[i][j][k])/(1+exp(mu_read[i][j][k]));
    //         M[i][j][k]=sqrt(mu_read[i][j][k]*(1-mu_read[i][j][k]));
                    /* }*/
                    
                    //below is newly added
                    
                   mu[i][j][k]=1/(1+params[3]*nTriple[i][j][k])*mu_read[i][j][k];
                       
                    
                    combY[i][j]+=1/(1+params[4]*new_nComb[i][j])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/M[i][j][k];
                    
                     M_trans[i][j]+=new_nTriple[i][j][k]*M[i][j][k];
                    
                }
                
                for(p=1;p<=Nbeta;p++){
                    combX[i][j][p]=0;
                    for (k=1;k<=nrep;k++){
                        combX[i][j][p]+=1/(1+params[4]*new_nComb[i][j])*nTriple[i][j][k]*M[i][j][k]*Cov[i][j][k][p];
                    }
                }
                
            }
        }
        
        for (i=1;i<=n1;i++){
            for (j=1;j<=n2;j++){
                sigmacombY[i][j]=0;
                
                for(n=1;n<=n1;n++){
                    for(m=1;m<=n2;m++){
                        
                        sigmacombY[i][j]-=(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[i][m]-doublekernal[n][j])*combY[n][m];
                        
                    }
                }
                // printf("%d\t%d\t%lf\n",i,j,sigmacombY[i][j]);
            }
            
        }
        
        
        
        for (p=1;p<=Nbeta;p++){
            beta_add[p]=0;
            
            for (i=1;i<=n1;i++){
                temcombX[i][p]=0;
                for (j=1;j<=n2;j++){
                     temcombX[i][p]+=combX[i][j][p];
                    for (k=1;k<=nrep;k++){
                        beta_add[p]+=Cov_read[i][j][k][p]*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k]);
                    }
                    //beta_add[p]+=Cov_read[i][j][k]
                    beta_add[p]-=params[4]*combY[i][j]*combX[i][j][p]*(1+params[4]*new_nComb[i][j]);
            
                   beta_add[p]+=combX[i][j][p]*sigmacombY[i][j];
                    
                }
                
            }
         //printf("beta add %d \t %lf\n",p,beta_add[p]);
            
            for(j=1;j<=n2;j++){
                temcombX2[j][p]=0;
                for(i=1;i<=n1;i++){
                    temcombX2[j][p]+=combX[i][j][p];
                }
                
            }

            
        }
    
  
    
     /*  for (i=1;i<=n1;i++){
            for (j=1;j<=n2;j++){
                for (k=1;k<=nrep;k++){
                    printf("%d\t%d\t%d\t",i,j,k);
                    for(p=1;p<=Nbeta;p++){
                        printf("%d\t%lf\t",p,Cov[i][j][k][p]);
        }
                    printf("\n");
            }
        }
        }*/
       /* double total=0;
        printf("this is okay\n");
           for (i=1;i<=n1;i++){
                for (j=1;j<=n2;j++){
         for (k=1;k<=nrep;k++){
             total+=Cov[i][j][k][5]*(Y[i][j][k]);
            printf("%lf\t%lf\n",Cov[i][j][k][5]*(Y[i][j][k]),Cov[i][j][k][5],Y[i][j][k],total);
       // printf("%lf\n",temcombX2[j][5]);
        }
                }
           }
        printf("this is okay\n");*/
          
        
       // printf("this is Nbeta %d\n",Nbeta);
       
       
     
       
        for (p=1;p<=Nbeta;p++){
            for (t=p;t<=Nbeta;t++){
                
                Fisher[p][t]=0;
                
                
                for (j=1;j<=n2;j++){
                    for(m=1;m<=n2;m++){
                        Fisher[p][t]-=extend_matrix2[j][m]*temcombX2[j][p]*temcombX2[m][t];
                        
                    }
                }
                
              // printf("this is Fisher again %lf\n",Fisher[1][21]);
                
                for (i=1;i<=n1;i++){
                  
                    for(n=1;n<=n1;n++){
                      Fisher[p][t]-=extend_matrix1[i][n]*temcombX[i][p]*temcombX[n][t];
                    }
                    
                    
                    for (j=1;j<=n2;j++){
                        for (k=1;k<=nrep;k++){
                            Fisher[p][t]+=new_nTriple[i][j][k]*M[i][j][k]*Cov_read[i][j][k][p]*M[i][j][k]*Cov_read[i][j][k][t];
                        }
                        
                        Fisher[p][t]-=params[4]*(1+params[4]*new_nComb[i][j])*combX[i][j][p]*combX[i][j][t];
                        
                        
                      Fisher[p][t]-=-doublekernal[i][j]*temcombX[i][p]*temcombX2[j][t]-doublekernal[i][j]*temcombX2[j][p]*temcombX[i][t];
                        
                    }
                    
                }
                Fisher[t][p]=Fisher[p][t];
            // printf("%lf\t",Fisher[p][t]);
            }
        // printf("\n");
        }
       
      
       /* for(i=1;i<=n1;i++){
            for(j=1;j<=n2;j++){
                for(k=1;k<=nrep;k++){
            //    printf("%d\t%d\t%d\t%lf\t%lf\n",i,j,k,M[i][j][k],M[i][j][k]*Cov[i][j][k][1]);
            }
            //printf("\n");
        }
        }*/
      /* printf("Fisher matrix\n");
        for (p=1;p<=Nbeta;p++){
            for (t=1;t<=Nbeta;t++){
                printf("%lf\t",Fisher[p][t]);
            }
            printf("\n");
        }*/

        
    
        i=solve(Fisher_inv,Fisher,Nbeta,tem_m1,tem_x1,&info,&info2);
        
    
   
        
        
        if(i!=0){
            printf("serious warning: Fisher information is not invertable\t %lf\t%lf\t%lf\t%lf\n",params[1],params[2],params[3],params[4]);
        }
        
        
        
        delta=0;
        
        
       // only for test!!
        for (i=1;i<=Nbeta;i++){
            beta_new[i]=beta[i];
            for (j=1;j<=Nbeta;j++){
                beta_new[i]+=Fisher_inv[i][j]*beta_add[j];
            }
            
            delta=fmax(delta,fabs(beta_new[i]-beta[i]));// infinity norm
            //printf("what is the difference %lf\n",fabs(beta_new[i]-beta[i]));
            
            beta[i]=beta_new[i];
            //printf("this is beta %lf\n",beta_new[i]);
       //printf("this is beta new %lf\n",beta[i]);
        }
        
    }
    
          //printf("mu is %lf\t%lf\t%lf\n",M[1][1][1],mu_read[1][1][2],mu_read[1][1][3]);
        sigma_t=0;
        for (i=1;i<=n1;i++){
            for (j=1;j<=n2;j++){
                for (k=1;k<=nrep;k++){
                    
                  //printf("%d\t%d\t%d\t%lf\n",i,j,k,mu[i][j][k]);
                    
                sigma_t+=(Y_sq_read[i][j][k]+mu_read[i][j][k]*(nTriple[i][j][k]*mu_read[i][j][k]-2*Y_read[i][j][k]))/(M[i][j][k]*M[i][j][k]);
                  
                    //add
                   sigma_t-=params[3]*(1+nTriple[i][j][k]*params[3])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/(M[i][j][k]*M[i][j][k]);
                   
                    
              /*  sigma_t+=(Y_sq[i][j][k]+mu[i][j][k]*(nTriple[i][j][k]*mu[i][j][k]-2*Y[i][j][k]))/(M[i][j][k]*M[i][j][k]);
                    
                    
                    sigma_t-=params[3]*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/(M[i][j][k]*M[i][j][k]);*/
                    
                //    printf("%lf\t%lf\t%lf\n",params[3]*(1+nTriple[i][j][k]*params[3])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/(M[i][j][k]*M[i][j][k]),params[3]*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/(M[i][j][k]*M[i][j][k]),params[3]*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/(M[i][j][k]*M[i][j][k])*(1+nTriple[i][j][k]*params[3]));
                           

                     
                     
                }
                //add
               sigma_t-=params[4]*combY[i][j]*combY[i][j]*(1+params[4]*new_nComb[i][j]);
                
                //add
               sigma_t+=combY[i][j]*sigmacombY[i][j];
                
            }
        }
        //printf("sigma %lf\n",sigma_t);
    

        sigma_t=sigma_t*(1+params[1]+params[2]+params[3]+params[4])/Npeo;
        printf("sigma %lf\n",sigma_t);
    
        
        

    
   /* for(i=1;i<=n1;i++){
        printf("%d\t%lf\n",i,extend_matrix1[i][i]);
    }*/
    
   /*for(p=1;p<=Nbeta;p++){
        for(t=1;t<=Nbeta;t++){
            fprintf(outf3,"%lf\t",Fisher[p][t]);
        }
        fprintf(outf3,"\n");
    }
    fflush(outf3);*/
    
    //update everything after coefficient!!
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            M_trans[i][j]=0;
            combY[i][j]=0;
            combY_sq[i][j]=0;
            for (k=1;k<=nrep;k++){
                mu[i][j][k]=0;
                mu_read[i][j][k]=0;
                
                for (p=1;p<=Nbeta;p++){
                    mu_read[i][j][k]+=Cov_read[i][j][k][p]*beta[p];
                }
                
                M[i][j][k]=1;
                
                //logistic
      //  mu_read[i][j][k]=exp(mu_read[i][j][k])/(1+exp(mu_read[i][j][k]));
  //  M[i][j][k]=sqrt(mu_read[i][j][k]*(1-mu_read[i][j][k]));
                               
                //below is newly added
                
                mu[i][j][k]=1/(1+params[3]*nTriple[i][j][k])*mu_read[i][j][k];
                
                
                combY[i][j]+=1/(1+params[4]*new_nComb[i][j])*(Y[i][j][k]-nTriple[i][j][k]*mu[i][j][k])/M[i][j][k];
                
                M_trans[i][j]+=new_nTriple[i][j][k]*M[i][j][k];
                
            }
            
            for(p=1;p<=Nbeta;p++){
                combX[i][j][p]=0;
                for (k=1;k<=nrep;k++){
                    combX[i][j][p]+=1/(1+params[4]*new_nComb[i][j])*nTriple[i][j][k]*M[i][j][k]*Cov[i][j][k][p];
                }
            }
            
        }
    }
    
    for (i=1;i<=n1;i++){
        for (j=1;j<=n2;j++){
            sigmacombY[i][j]=0;
            
            for(n=1;n<=n1;n++){
                for(m=1;m<=n2;m++){
                    
                    sigmacombY[i][j]-=(extend_matrix1[i][n]+extend_matrix2[j][m]-doublekernal[i][m]-doublekernal[n][j])*combY[n][m];
                    
                }
            }
        }
        
    }
    
    
    
    for (p=1;p<=Nbeta;p++){
        
        for (i=1;i<=n1;i++){
            temcombX[i][p]=0;
            for (j=1;j<=n2;j++){
                temcombX[i][p]+=combX[i][j][p];
            }
            
        }
        
        for(j=1;j<=n2;j++){
            temcombX2[j][p]=0;
            for(i=1;i<=n1;i++){
                temcombX2[j][p]+=combX[i][j][p];
            }
            
        }
        
        
    }
        
    
    
    for (p=1;p<=Nbeta;p++){
        for (t=p;t<=Nbeta;t++){
            
            Fisher[p][t]=0;
            
            
            for (j=1;j<=n2;j++){
                for(m=1;m<=n2;m++){
                    Fisher[p][t]-=extend_matrix2[j][m]*temcombX2[j][p]*temcombX2[m][t];
                    
                }
            }
            
            // printf("this is Fisher again %lf\n",Fisher[1][21]);
            
            for (i=1;i<=n1;i++){
                
                for(n=1;n<=n1;n++){
                    Fisher[p][t]-=extend_matrix1[i][n]*temcombX[i][p]*temcombX[n][t];
                }
                
                
                for (j=1;j<=n2;j++){
                    for (k=1;k<=nrep;k++){
                        Fisher[p][t]+=new_nTriple[i][j][k]*M[i][j][k]*Cov_read[i][j][k][p]*M[i][j][k]*Cov_read[i][j][k][t];
                    }
                    
                    Fisher[p][t]-=params[4]*(1+params[4]*new_nComb[i][j])*combX[i][j][p]*combX[i][j][t];
                    
                    
                    Fisher[p][t]-=-doublekernal[i][j]*temcombX[i][p]*temcombX2[j][t]-doublekernal[i][j]*temcombX2[j][p]*temcombX[i][t];
                    
                }
                
            }
            Fisher[t][p]=Fisher[p][t];
            // printf("%lf\t",Fisher[p][t]);
        }
        // printf("\n");
    }
    
    
    i=solve(Fisher_inv,Fisher,Nbeta,tem_m1,tem_x1,&info,&info2);
   
    

    
    
    
    
    
    
    
    
    
    
    printf("parameter is %.10lf\t%.10lf\t%.10lf\t%.10lf\n",params[1],params[2],params[3],params[4]);
    printf("sigma_t_sq is %lf\n",sigma_t);
    printf("coefficient: estimate   \t   sd\n");
    for (i=1;i<=Nbeta;i++){
        printf("\t \t %lf\t%lf\n",beta[i],sqrt(Fisher_inv[i][i]*sigma_t/(1+params[1]+params[2]+params[3]+params[4])));
        fprintf(outf2,"%lf\t%lf\t",beta[i],sqrt(Fisher_inv[i][i]*sigma_t/(1+params[1]+params[2]+params[3]+params[4])));
    }
    fprintf(outf2,"%lf\t",sigma_t);
    
    return (int)info;
}



void readin(){

    int i,j,id,id2;
    int ind1, ind2,rep, leaf;
    double value, value2, value3, value4,value5,value6, value7, value8, value9,value10,value11,value12,value13,value14,value15,value16,value17,value18,value19,value20,value21,value22,value23,value24,value25,value26,value27;
    
       //printf("okay%d\t?\n",n1);
         // fscanf(matrixf1,"%lf", &value);
    // printf("%lf\t",value);
   
    
    /*for(i=1;i<=n2;i++){
        printf("test1 %d\t%d\n",i,groupmember[i]);
    }
    printf("this is n2 %d\n",n2);*/

    //  printf("this is groupmember %d\n",groupmember2[10]);
    
for (i=1;i<=n1;i++){
            fscanf(memberf,"%d", &id);
            //membership[i]=id;
    membership[i]=id;
 }
      
 /* for(i=1;i<=n2;i++){
        printf("test2 %d\t%d\n",i,groupmember[i]);
    }*/
    
   

   for (i=1;i<=n1;i++){
        for (j=1;j<=n1;j++){
            fscanf(matrixf1,"%lf", &value);
            matrix1[i][j]=value;
            phi1[i][j]=value;
            // printf("%lf\t",value);
        }
       // printf("\n");
    } 

    
  /*  double **matrix_test=dmatrix(1,354,1,354);
    double **matrix_test_inv=dmatrix(1,354,1,354);
    
   for (i=1;i<=354;i++){
        for (j=1;j<=354;j++){
            fscanf(matrixf2,"%lf", &value);
            matrix_test[i][j]=value;
        }
    }

    matinvchol(matrix_test, 354, matrix_test_inv, &det_Na);
    //965.066849
     printf("log determinant%lf \n",det_Na);*/
    
  for (i=1;i<=n2;i++){
        for (j=1;j<=n2;j++){
            fscanf(matrixf2,"%lf", &value);
            matrix2[i][j]=value;
            phi2[i][j]=value;
           // printf("%lf\t",matrix2[i][j]);
        }
     // printf("\n");
    }
    
  
 // printf("okay?\n");
  
      double valuenew=0;
    double total;
  // printf("sucess so far?\n");
while(fscanf(covf,"%d %d %d %d %lf %lf %lf %lf %lf", &id, &id2, &rep,&leaf,&value,&value3,&value4,&value2,&valuenew)==9){
    

   // for(i=1;i<=10;i++){
//while(fscanf(covf,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &id, &id2, &rep,&leaf,&value,&value3,&value4,&value2,&valuenew,&value5,&value6,&value7)==12){
    
      //  fscanf(covf,"%d %d %d %d %d %lf %lf %lf %lf", &ind1, &id, &id2, &rep,&leaf,&value,&value3,&value4,&value2);
   // total+=value5;
    
    //printf("total is %lf\n",total);
//while(fscanf(covf,"%d %d %d %d %lf %lf %lf %lf %lf  ", &id, &id2, &rep,&leaf,&value,&value3,&value4,&value2,&valuenew)==9){
        
       // Plant[ind1]=id;
        nPlant[id]++;
        nBac[id2]++;
        nComb[id][id2]++;
        nTriple[id][id2][rep]++;
    // printf("%d\t%d\t%d\t%lf\n",id,id2,rep,Cov[id][id2][rep][1]);
      Cov_read[id][id2][rep][1]=value;
    Cov_read[id][id2][rep][2]=value3;
    Cov_read[id][id2][rep][3]=value4;
   //Cov_read[id][id2][rep][4]=valuenew;
  // Cov_read[id][id2][rep][4]=(id2==5)+(id2==6)+(id2==7)+(id2==8);//distingish clade 1 vs 2
  
    
// Cov_read[id][id2][rep][4]=(id2==10)+(id2==13)+(id2==14)+(id2==16);//distingish clade 1 vs 2
    
 //Cov_read[id][id2][rep][5]=(id2==6)+(id2==7)+(id2==8)+(id2==11)+(id2==12)+(id2==15)+(id2==1)+(id2==9)+(id2==17);//distingish clade 1 vs 34*/
    
    
    // top 1 bac SNP
    // Cov_read[id][id2][rep][4]=(id2==1)+(id2==7)+(id2==10)+(id2==11)+(id2==12)+(id2==13)+(id2==14)+(id2==15)+(id2==17);
    
    
    // Cov_read[id][id2][rep][5]=(id2==1)+(id2==5)+(id2==6)+(id2==10)+(id2==11)+(id2==13);
 
    
    //TOP1 bac SNP (Clade 34)
 // Cov_read[id][id2][rep][4]=(id2==3)+(id2==5)+(id2==7);
    
    //TOP2 bac SNP (Clade 34)
  // Cov_read[id][id2][rep][4]=(id2==5)+(id2==7);
    
    
    //TOP bac SNP (Clade 12)
    // Cov_read[id][id2][rep][5]=(id2==2);//mutation
    //Cov_read[id][id2][rep][6]=(id2==4)+(id2==12)+(id2==13);//deletion
    
    
    
   /* Cov_read[id][id2][rep][5]=(id==127);*/
    
 /* Cov_read[id][id2][rep][5]=(id==84)*(groupmember[id2]==2);
  Cov_read[id][id2][rep][6]=(id==88)*(groupmember[id2]==1);*/
 
    
    /*Cov_read[id][id2][rep][5]=(membership[id]==7)*(groupmember[id2]==4);//FRA:clade4
    
     Cov_read[id][id2][rep][6]=(membership[id]==7)*(groupmember[id2]==1);//FRA:clade1
    
        Cov_read[id][id2][rep][7]=(membership[id]==16)*(groupmember[id2]==1);//USA:clade1
    
    Cov_read[id][id2][rep][8]=(membership[id]==2)*(groupmember[id2]==1);//BEL:clade1

    Cov_read[id][id2][rep][9]=(membership[id]==11)*(groupmember[id2]==1);//POL:clade1
    
    Cov_read[id][id2][rep][10]=(membership[id]==8)*(groupmember[id2]==1);//GER:clade1*/
    

 /*   Cov_read[id][id2][rep][5]=value5*(groupmember[id2]==1);
    
    
    Cov_read[id][id2][rep][6]=value5*(groupmember[id2]==2);
    
    Cov_read[id][id2][rep][7]=value5;*/
    
    /*Cov_read[id][id2][rep][8]=*/
    
  /* Cov_read[id][id2][rep][5]=value5;
    Cov_read[id][id2][rep][6]=value6;
     Cov_read[id][id2][rep][7]=value7;
  Cov_read[id][id2][rep][8]=value8;
     Cov_read[id][id2][rep][9]=value9;
     Cov_read[id][id2][rep][10]=value10;
     Cov_read[id][id2][rep][11]=value11;
     Cov_read[id][id2][rep][12]=value12;
     Cov_read[id][id2][rep][13]=value13;
     Cov_read[id][id2][rep][14]=value14;
     Cov_read[id][id2][rep][15]=value15;
     Cov_read[id][id2][rep][16]=value16;
     Cov_read[id][id2][rep][17]=value17;
     Cov_read[id][id2][rep][18]=value18;
    Cov_read[id][id2][rep][19]=value19;
   Cov_read[id][id2][rep][20]=value20;
     Cov_read[id][id2][rep][21]=value21;
     Cov_read[id][id2][rep][22]=value22;*/
    
    /*Cov_read[id][id2][rep][5]=value8;
    Cov_read[id][id2][rep][6]=value9;
    Cov_read[id][id2][rep][7]=value10;
    Cov_read[id][id2][rep][8]=value11;
    Cov_read[id][id2][rep][9]=value12;
    Cov_read[id][id2][rep][10]=value13;
    Cov_read[id][id2][rep][11]=value14;
    Cov_read[id][id2][rep][12]=value15;
    Cov_read[id][id2][rep][13]=value16;
    Cov_read[id][id2][rep][14]=value17;
    Cov_read[id][id2][rep][15]=value18;
    Cov_read[id][id2][rep][16]=value19;
    Cov_read[id][id2][rep][17]=value20;
    Cov_read[id][id2][rep][18]=value21;
    Cov_read[id][id2][rep][19]=value22;*/
    /*Cov_read[id][id2][rep][23]=value23;
     Cov_read[id][id2][rep][24]=value24;
    Cov_read[id][id2][rep][25]=value25;*/
     /*Cov[id][id2][rep][26]=value26;
     Cov[id][id2][rep][27]=value27;*/
    
    value2=value2/4;
        Y_read[id][id2][rep]+=value2;
        Y_sq_read[id][id2][rep]+=value2*value2;
    /*if(value21==1){*/
//printf("%d %d %d %d \n",  id, id2, rep,leaf);
   // }*/
    }
    
        
   //printf("okay?\n");
   // printf("okay2?\n");
   /* double *freq=dvector(1,ngenes);
    double *sd=dvector(1,ngenes);
    
    for(i=1;i<=ngenes;i++){
        freq[i]=0;
        sd[i]=0;
        for(j=1;j<=n1;j++){
            freq[i]+=G[i][j];
        }
        freq[i]=freq[i]/n1;
        for(j=1;j<=n1;j++){
            sd[i]+=(G[i][j]-freq[i])*(G[i][j]-freq[i]);
        }
        
        sd[i]=sqrt(sd[i]/(n1-1));
        
        for(j=1;j<=n1;j++){
            G[i][j]=(G[i][j]-freq[i])/sd[i];
        }
    }*/
    
    
       }
/*int solve(double **y, double **a,int n, double *m,double *x,int *info,int *info2){
    
 
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

int solve(double **y, double **a,int n, double *m,double *x,integer *info,integer *info2){
    
   /* int i,j,k;
    
    double *S=dvector(0,n-1);
    double *U=dvector(0,n*n-1);
    double *VT=dvector(0,n*n-1);
    integer nn=n;
    integer nn2=5*n;
    
    double *work=dvector(1,nn2);
    
    int con;
    double **phi1_G_inv_tem1=dmatrix(1,n,1,n);
    double **phi1_G_inv_tem2=dmatrix(1,n,1,n);
    
    for (i=1;i<=n;i++){
        for (j=1;j<=n;j++)
            m[(i-1)*n+j-1]=a[i][j];
    }
    
    dgesvd_("A", "A", &nn, &nn, m, &nn, S, U, &nn, VT, &nn, work,&nn2,info);
    
    
    
    con=S[0];
    for(i=1;i<=n;i++){
        // printf("%lf\t",S[i-1]);
        if (S[i-1]<=tol*con)
            S[i-1]=0;
        else
            S[i-1]=1/S[i-1];
        //printf("%lf\t",S[i-1]);
    }
    
    
    
    
    for (i=0;i<=(n-1);i++){
        for (j=0;j<=(n-1);j++){
            phi1_G_inv_tem1[i+1][j+1]=U[i*n+j];
            phi1_G_inv_tem2[i+1][j+1]=VT[i*n+j];
            // fprintf(outf3,"%lf\t",phi1_G_inv_tem1[i+1][j+1]);
        }
        // fprintf(outf3,"\n");
    }
    
    for (i=1;i<=n;i++){
        for (j=1;j<=n;j++){
            y[i][j]=0;
            for(k=1;k<=n;k++){
                y[i][j]+= phi1_G_inv_tem1[k][i]*S[k-1]*phi1_G_inv_tem2[j][k];
            }
            // y[j][i]=y[i][j];
            // fprintf(outf3,"%lf\t",phi1_G_inv[i][j]);
        }
        //fprintf(outf3,"\n");
    }
    
    free_dmatrix(phi1_G_inv_tem1,1,n,1,n);
    free_dmatrix(phi1_G_inv_tem2,1,n,1,n);
    free_dvector(S,0,n-1);
    free_dvector(U,0,n*n-1);
    free_dvector(VT,0,n*n-1);
    free_dvector(work,1,nn2);
    
    return (*info);
}*/

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
 
}



double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
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



int solve_general(double **y, double **a,int n, double *m,double *x,integer *info,integer *info2){
    
    int i,j,k;
    
    double *S=dvector(0,n-1);
    double *U=dvector(0,n*n-1);
    double *VT=dvector(0,n*n-1);
    integer nn=n;
    integer nn2=5*n;
    
    double *work=dvector(1,nn2);
    
    int con;
    double **phi1_G_inv_tem1=dmatrix(1,n,1,n);
    double **phi1_G_inv_tem2=dmatrix(1,n,1,n);
    
    for (i=1;i<=n;i++){
        for (j=1;j<=n;j++)
            m[(i-1)*n+j-1]=a[i][j];
    }
    
    dgesvd_("A", "A", &nn, &nn, m, &nn, S, U, &nn, VT, &nn, work,&nn2,info);
    
    
    
    con=S[0];
    for(i=1;i<=n;i++){
        // printf("%lf\t",S[i-1]);
        if (S[i-1]<=tol*con)
            S[i-1]=0;
        else
            S[i-1]=1/S[i-1];
        //printf("%lf\t",S[i-1]);
    }
    
    
    
    
    for (i=0;i<=(n-1);i++){
        for (j=0;j<=(n-1);j++){
            phi1_G_inv_tem1[i+1][j+1]=U[i*n+j];
            phi1_G_inv_tem2[i+1][j+1]=VT[i*n+j];
            // fprintf(outf3,"%lf\t",phi1_G_inv_tem1[i+1][j+1]);
        }
        // fprintf(outf3,"\n");
    }
    
    for (i=1;i<=n;i++){
        for (j=1;j<=n;j++){
            y[i][j]=0;
            for(k=1;k<=n;k++){
                y[i][j]+= phi1_G_inv_tem1[k][i]*S[k-1]*phi1_G_inv_tem2[j][k];
            }
            // y[j][i]=y[i][j];
            // fprintf(outf3,"%lf\t",phi1_G_inv[i][j]);
        }
        //fprintf(outf3,"\n");
    }
    
    free_dmatrix(phi1_G_inv_tem1,1,n,1,n);
    free_dmatrix(phi1_G_inv_tem2,1,n,1,n);
    free_dvector(S,0,n-1);
    free_dvector(U,0,n*n-1);
    free_dvector(VT,0,n*n-1);
    free_dvector(work,1,nn2);
    
    return (*info);
    
}

void matinvchol(double **a, int n, double **inv, double *det)
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
 
 */
{
    double d, ldet, *col, *diag, sum;
    double **a_temp;
    a_temp=dmatrix(1,n,1,n);
    int i, j, k;
    int choldc(double **a, int n, double *diag);
    
    for(i=1;i<=n;i++){
          for(j=1;j<=n;j++){
              a_temp[i][j]=a[i][j];
          }
    }
    
    diag = dvector(1,n);
    
    /* Do cholesky decomposition */
    choldc(a_temp,n,diag);
    
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
                    fprintf(stderr,"\nWARNING: Non positive definite matrix.\n");
                    retval = 0;
                }
                p[i]=sqrt(sum);
            } else a[j][i]=sum/p[i];
        }
    }
    return(retval);
}

