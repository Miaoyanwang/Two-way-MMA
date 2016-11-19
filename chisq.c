#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "nrutil.c"

#define MAXLEN 2048
char coeffile[MAXLEN]="coeffile.txt",kinshipfile[MAXLEN]="kinship.txt",genofile[MAXLEN]="genotype.txt",invkinshipfile[MAXLEN]="inv_kinship.txt",outfile[MAXLEN]="out.txt";
FILE *coef,*kinshipf,*genof,*invkinshipf,*outf;


int n1=129;
int n2=22;

int main(int argc, char *argv[])
{
    int i,j,n;
    if (argc>1){
        for (i=1;i<argc && argv[i][0]=='-';i++){
            switch(argv[i][1])
            {
                case 'c':
                    strncpy(coeffile,argv[++i],MAXLEN);
                    printf("user specified coef: %s\n",coeffile);
                break;
                case 'm':
                    strncpy(kinshipfile,argv[++i],MAXLEN);
                    printf("user specified matrix: %s\n",kinshipfile);
                    break;
                case 'g':
                    strncpy(genofile,argv[++i],MAXLEN);
                    printf("user specified genofile: %s\n",genofile);
                    break;
                case 'i':
                    strncpy(invkinshipfile,argv[++i],MAXLEN);
                    printf("user specified invkinship: %s\n",invkinshipfile);
                    break;
                case 'o':
                    strncpy(outfile,argv[++i],MAXLEN);
                    printf("user specified out: %s\n",outfile);
                    break;
                    
            }
        }
    }
    if((coef=fopen(coeffile,"r"))==NULL){
        printf("Cannot open coef file.\n");
        exit(1);
    }
    
    if((kinshipf=fopen(kinshipfile,"r"))==NULL){
        printf("Cannot open kinship file.\n");
        exit(1);
    }
    
    if((invkinshipf=fopen(invkinshipfile,"r"))==NULL){
        printf("Cannot open inv kinship file.\n");
        exit(1);
    }
    
    if((genof=fopen(genofile,"r"))==NULL){
        printf("Cannot open genotype file.\n");
        exit(1);
    }
    
    if((outf=fopen(outfile,"w"))==NULL){
        printf("Cannot open out file.\n");
        exit(1);
    }

    double **coefvector=dmatrix(1,n1,1,n2);
    double value;
    double **kinship=dmatrix(1,n1,1,n1);
    double *G=dvector(1,n1);
    int id,id2;
    
  
    while(fscanf(coef,"%d\t %d\t %lf\n",&id,&id2,&value)==3){
        coefvector[id][id2]=value;
    }
  

    for (i=1;i<=n1;i++){
        for (j=1;j<=n1;j++){
            fscanf(kinshipf,"%lf", &value);
            kinship[i][j]=value;
           // printf("%lf\t",kinship[i][j]);
        }
       // printf("\n");
    }

     
    int chr, pos;
    
    double GG_general=0,G1=0,norm=0;
    double *VG=dvector(1,n2),*VV_general=dvector(1,n2);
    
    double **invkinship;
    invkinship=dmatrix(1,n1,1,n1);
    
    
    for (i=1;i<=n1;i++){
    for (j=1;j<=n1;j++){
        fscanf(invkinshipf,"%lf", &value);
        invkinship[i][j]=value;
        // printf("%lf\t",kinship[i][j]);
    }
    }
    
    
    
    while(fscanf(genof,"%d %d",&chr,&pos)==2){
        
        for(i=1;i<=n1;i++){
        fscanf(genof,"%lf",&value);
            G[i]=value/2;
            //printf("%lf\t",G[i]);
        }
        
        
        GG_general=0,G1=0,norm=0;
        for(n=1;n<=n2;n++){
            VG[n]=0;
            VV_general[n]=0;
    }
        fprintf(outf,"%d\t%d\t",chr,pos);
        
         for(i=1;i<=n1;i++){
        for(j=1;j<=n1;j++){
            GG_general+=G[i]*invkinship[i][j]*G[j];
            G1+=G[i]*invkinship[i][j];
            norm+=invkinship[i][j];
        }
         }
        
        GG_general-=G1*G1/norm;
  
        for(n=1;n<=n2;n++){
    for(i=1;i<=n1;i++){
         VG[n]+=G[i]*coefvector[i][n];
     for(j=1;j<=n1;j++){
            VV_general[n]+=coefvector[i][n]*kinship[i][j]*coefvector[j][n];
        }
    }
        

    
     // printf("%lf\t%lf\t%lf\t",VG,GG_general/(n1-1),VV_general);
    
    VG[n]=VG[n]*VG[n]*(n1-1)/(GG_general*VV_general[n]);
   fprintf(outf,"%lf\t",VG[n]);

    }
        fprintf(outf,"\n");
        fflush(outf);
    }

}






