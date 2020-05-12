#include <stdio.h>       
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SEG 15
#define NLAYER1 330
#define NLAYER2 150
#define NLAYER3 2
#define MLAYER1 60
#define MLAYER2 45
#define MLAYER3 2
#define ALPHA 0.001
#define LAMDA 0.9
#define REPEAT 40
#define MAX_LEN 5000
#define LenScale 1000

void Test(char sprffile[20]);
int ReadSprf(char file[20]);
int Readmat(char file[20]);
void initNet();
void FeedNet(int I);
void Netout();
void trainNet();
void changeWgt();
void result(char *file);
void initNet1();
void FeedNet1(int I);
void Netout1();
void trainNet1();
void changeWgt1();
void result1(char *file);
void change(float *w,float *d,float *D);
int Max(float *p,int N);

char StrA1[MAX_LEN],StrA2[MAX_LEN],prdA2[MAX_LEN];
float Score[MAX_LEN][21];
float nout[MAX_LEN][NLAYER3],pout[MAX_LEN][NLAYER3];
int Len;
float w1[NLAYER1][NLAYER2],Dw1[NLAYER1][NLAYER2];
float w2[NLAYER2][NLAYER3],Dw2[NLAYER2][NLAYER3];
float in1[NLAYER1],in2[NLAYER2],out2[NLAYER2],in3[NLAYER3],out[NLAYER3];
float tgtout[NLAYER3];
float v1[MLAYER1][MLAYER2],Dv1[MLAYER1][MLAYER2];
float v2[MLAYER2][MLAYER3],Dv2[MLAYER2][MLAYER3];
float im1[MLAYER1],im2[MLAYER2], mout2[MLAYER2],im3[MLAYER3],mout[MLAYER3];

float dw1[NLAYER1][NLAYER2],dw2[NLAYER2][NLAYER3];
float dv1[MLAYER1][MLAYER2],dv2[MLAYER2][MLAYER3]; // Derivatives
 
int n[NLAYER3][NLAYER3], m[MLAYER3][MLAYER3];
int Ntotal, Ngood, Mgood;
int Mmax=0;
 
int main(int argu, char *argv[]){
  int II,i,j,k,count;

  char wghtfile[30];
  FILE *fpout;
   

  //needed files: wgtfile, pdb.mat3
  if(argu<2) {
    fprintf(stderr,"USAGE: solve wgtfile pdb\n"); 
    return(0);
  }
  
  strcpy(wghtfile,argv[1]);
  if((fpout=fopen(wghtfile,"r"))==NULL) {
    printf("cannot open wghtfile %s\n", wghtfile);
    exit(0);
  }
  
  for(i=0;i<NLAYER1;i++){
    for(j=0;j<NLAYER2;j++){
      fscanf(fpout,"%f\n",&w1[i][j]); //read trained parameters
      Dw1[i][j]=0.0;
    }
  }
  for(i=0;i<NLAYER2;i++){
    for(j=0;j<NLAYER3;j++){
      fscanf(fpout,"%f\n",&w2[i][j]); //read trained parameters
      Dw2[i][j]=0.0;
    }
  }
  for(i=0;i<MLAYER1;i++){
    for(j=0;j<MLAYER2;j++){
      fscanf(fpout,"%f\n",&v1[i][j]); //read trained parameters
      Dv1[i][j]=0.0;
    }
  }
  for(i=0;i<MLAYER2;i++){
    for(j=0;j<MLAYER3;j++){
      fscanf(fpout,"%f\n",&v2[i][j]); //read trained parameters
      Dv2[i][j]=0.0;
    }
  }
  fclose(fpout);
  
  Test(argv[2]);
}

/******************************************************/
void Test(char sprffile[20]){
  int i,j;
  
  Readmat(sprffile); //read Score from protein.mat
  
  for(i=1;i<=Len;i++){
    FeedNet(i); //Score->in1
    Netout(); //in1->out, use w1[][], w2[][]
    for(j=0;j<NLAYER3;j++)
      nout[i][j]=out[j];
  }
  
  for(i=1;i<=Len;i++)  {
    FeedNet1(i); //nout->im1
    Netout1(); //im1->mout, use v1[][], v2[][]
    for(j=0;j<MLAYER3;j++)
      pout[i][j]=mout[j];
  }
  result1(sprffile); //pout->prediction
  
  Ntotal+=Len;
} 
  
/*************************************************************/
int Readmat(char file[20]){
     int i,j,k;
     FILE *fp;
     int tmp[20];
     char ss[1000], filenm[45]="./";

     strcat(filenm, file);
     strcat(filenm,".mat3");
     if((fp=fopen(filenm,"r"))==NULL){
           fprintf(stderr,"cannot open %s\n",file);
	   exit(-1);
      }

     while(fgetc(fp)!='\n'){}
     while(fgetc(fp)!='\n'){}
     while(fgetc(fp)!='\n'){}
     Len=0;
     while(1){
	   fgets(ss,1000,fp);
	   if(strlen(ss)<50) break;
	    Len++;
	   sscanf(ss+6,"%c",&StrA1[Len]);
	   /*      printf("%c",StrA1[Len]); */
	         for(i=0;i<20;i++){
	        sscanf(ss+10+3*i,"%3f",&Score[Len][i]);
	      Score[Len][i]=1.0/(1.0+exp(-1.0*Score[Len][i]));
	     /*        printf("%f ",Score[Len][i]); */
	       }
	  /*      printf("\n"); */
    }

        fclose(fp);

	  return(0);
 }

/****************************************************************/
int ReadSprf(char file[20]){
  int i,j,k;
  float SecScore[7];
  FILE *fp;
  int tmp[20];
  char filenm[45]="PRF/";
  
  strcat(filenm,file);
  strcat(filenm,".prf");
  if((fp=fopen(filenm,"r"))==NULL)
    {  printf("cannot open file %s\n",filenm);
       return(-1);
    }

  while(fgetc(fp)!='\n'){}
  while(fgetc(fp)!='\n'){}
   Len=0;
    for(i=1;i<=10000;i++) {
      if(fscanf(fp,"  %d %c%c",&k,&StrA1[i],&StrA2[i])!=3)
        break;
      Len++;
    fgetc(fp);

    for(j=0;j<20;j++){
      fscanf(fp," %2d",tmp+j);
    }
    Score[i][0]=(float)tmp[0]; Score[i][1]=(float)tmp[17];
    Score[i][2]=(float)tmp[13]; Score[i][3]=(float)tmp[3];
    Score[i][4]=(float)tmp[2]; Score[i][5]=(float)tmp[16];
    Score[i][6]=(float)tmp[4]; Score[i][7]=(float)tmp[6];
    Score[i][8]=(float)tmp[7]; Score[i][9]=(float)tmp[8];
    Score[i][10]=(float)tmp[11]; Score[i][11]=(float)tmp[10];
    Score[i][12]=(float)tmp[12]; Score[i][13]=(float)tmp[5];
    Score[i][14]=(float)tmp[15]; Score[i][15]=(float)tmp[18];
    Score[i][16]=(float)tmp[19]; Score[i][17]=(float)tmp[9];
    Score[i][18]=(float)tmp[14]; Score[i][19]=(float)tmp[1];
       for(j=0;j<20;j++)
	  Score[i][j]=1.0/(1.0+exp(-1.0*Score[i][j])); 
    for(j=0;j<6;j++){
      fscanf(fp,"%f",&SecScore[j]);
    }
  }

  fclose(fp);

  return(0);
}
  
/****************************************************************/
void FeedNet(int I){
  int i,j,k;
  
  k=0;
  for(i=(I-SEG/2);i<(I-SEG/2)+SEG;i++){
    if(i<1 | i>Len){
      for(j=0;j<20;j++){
	in1[k]=0.0;
	k++;
      }
      in1[k]=1.0;
      k++;
      in1[k]=(float)Len/LenScale;
      k++;
    }
    else{
      for(j=0;j<20;j++){
	in1[k]=Score[i][j];
	k++;
      }
      in1[k]=0.0;
      k++;
      in1[k]=(float)Len/LenScale;
      k++;
    }   
  }
  for(i=0;i<NLAYER3;i++)
    tgtout[i]=0.0;
  tgtout[(StrA2[I]-48)%2]=1.0; 
  
}

/**************************************************************************/
void initNet(){
  int i,j;
  
  for(i=0;i<NLAYER1;i++){
    for(j=0;j<NLAYER2;j++){
      w1[i][j]=(0.5-drand48())*0.2;
      dw1[i][j]=0.0;
      Dw1[i][j]=0.0;
    }
  }
  
  for(i=0;i<NLAYER2;i++){
    for(j=0;j<NLAYER3;j++){
      w2[i][j]=(0.5-drand48())*0.2;
      dw2[i][j]=0.0;
      Dw2[i][j]=0.0;
    }
  }
}

/***************************************************************/
void Netout(){
  int i,j;

  for(i=0;i<NLAYER2;i++){
    in2[i]=0.0;
    for(j=0;j<NLAYER1;j++){
      in2[i]+=w1[j][i]*in1[j];
    }
    out2[i]=1.0/(1.0+exp(-1.0*in2[i]));
  }

  for(i=0;i<NLAYER3;i++){
    in3[i]=0.0;
    for(j=0;j<NLAYER2;j++){
      in3[i]+=w2[j][i]*out2[j];
    }
    out[i]=1.0/(1.0+exp(-1.0*in3[i]));
  }
}

/***************************************************************/
void trainNet(){
  int i,j;
  float dwO[NLAYER3];
  float dwH[NLAYER2];

  for(i=0;i<NLAYER2;i++) dwH[i]=0.0;
  
  for(j=0;j<NLAYER3;j++){
    dwO[j]=(tgtout[j]-out[j])*out[j]*(1.0-out[j]);   
    for(i=0;i<NLAYER2;i++){
      dw2[i][j]+=dwO[j]*out2[i];
      dwH[i]+=dwO[j]*w2[i][j];
    }
  }

  for(j=0;j<NLAYER2;j++){
    dwH[j]=dwH[j]*out2[j]*(1.0-out2[j]);
    for(i=0;i<NLAYER1;i++){
      dw1[i][j]+=dwH[j]*in1[i];
      }     
    }
}


/*****************************************************************/

void changeWgt()
{
  int i,j;

  for(j=0;j<NLAYER2;j++)
    for(i=0;i<NLAYER1;i++)
       change(&w1[i][j],&dw1[i][j],&Dw1[i][j]);           
     
  for(j=0;j<NLAYER3;j++)
    for(i=0;i<NLAYER2;i++)
       change(&w2[i][j],&dw2[i][j],&Dw2[i][j]);

 for(j=0;j<NLAYER2;j++)
    for(i=0;i<NLAYER1;i++)
       dw1[i][j]=0.0;

  for(j=0;j<NLAYER3;j++)
    for(i=0;i<NLAYER2;i++)
       dw2[i][j]=0.0;
  
}
/*****************************************************************/
void result(char *file){
  int i,j;

  j=0; 
  for(i=1;i<=Len;i++){
    
    prdA2[i]=Max(nout[i],NLAYER3)+48;
    if(prdA2[i]==(StrA2[i]-48)%2+48) 
      j++;
    
    n[(StrA2[i]-48)%2][prdA2[i]-48]++;
  }
  
  /*  printf("%s  Result1: %5.3f %4d\n", file,
      (float)j/(float)Len,Len);*/
  
  Ngood+=j;
}

/*****************************************************************/
void FeedNet1(int I){
  int i,j,k;
  
  k=0;
  for(i=(I-SEG/2);i<(I-SEG/2)+SEG;i++){
    if(i<1 | i>Len){
      for(j=0;j<NLAYER3;j++){
	im1[k]=0.0;
	k++;
      }
      im1[k]=1.0;
      k++;
      im1[k]=(float)Len/LenScale;
      k++;

    }
    else{
      for(j=0;j<NLAYER3;j++){
	im1[k]=nout[i][j];
	k++;
      }
      im1[k]=0.0;
      k++;
      im1[k]=(float)Len/LenScale;
      k++;
    }
  } 
  for(i=0;i<MLAYER3;i++)
    tgtout[i]=0.0;
  
  tgtout[(StrA2[I]-48)%2]=1.0; 
}


/**************************************************************************/
void initNet1(){
  int i,j;
 
  for(i=0;i<MLAYER1;i++){
    for(j=0;j<MLAYER2;j++){
      v1[i][j]=(0.5-drand48())*0.2;
      Dv1[i][j]=0.0;
    }
  }
 
  for(i=0;i<MLAYER2;i++){
    for(j=0;j<MLAYER3;j++){
      v2[i][j]=(0.5-drand48())*0.2;
      Dv2[i][j]=0.0;
    }
  }
}

/***************************************************************/
void Netout1(){
  int i,j;

  for(i=0;i<MLAYER2;i++){
    im2[i]=0.0;
    for(j=0;j<MLAYER1;j++){
      im2[i]+=v1[j][i]*im1[j];
    }
    mout2[i]=1.0/(1.0+exp(-1.0*im2[i]));
  }
  
  for(i=0;i<MLAYER3;i++){
    im3[i]=0.0;
    for(j=0;j<MLAYER2;j++){
      im3[i]+=v2[j][i]*mout2[j];
    }
    mout[i]=1.0/(1.0+exp(-1.0*im3[i]));
  }
}

/***************************************************************/
void trainNet1(){
  int i,j;
  float dwO[MLAYER3];
  float dwH[MLAYER2];

  for(i=0;i<MLAYER2;i++) dwH[i]=0.0;

  for(j=0;j<MLAYER3;j++){
    dwO[j]=(tgtout[j]-mout[j])*mout[j]*(1.0-mout[j]);
    for(i=0;i<MLAYER2;i++){
      dv2[i][j]+=dwO[j]*mout2[i];
      dwH[i]+=dwO[j]*v2[i][j];
    }
  }

  for(j=0;j<MLAYER2;j++){
    dwH[j]=dwH[j]*mout2[j]*(1.0-mout2[j]);
    for(i=0;i<MLAYER1;i++){
      dv1[i][j]+=dwH[j]*im1[i];
      }
    }
}
/*****************************************************************/

void changeWgt1()
{
  int i,j;

  for(j=0;j<MLAYER2;j++)
    for(i=0;i<MLAYER1;i++)
       change(&v1[i][j],&dv1[i][j],&Dv1[i][j]);

  for(j=0;j<MLAYER3;j++)
    for(i=0;i<MLAYER2;i++)
       change(&v2[i][j],&dv2[i][j],&Dv2[i][j]);
  
  for(j=0;j<MLAYER2;j++)
    for(i=0;i<=MLAYER1;i++)
       dv1[i][j]=0.0;

  for(j=0;j<MLAYER3;j++)
    for(i=0;i<=MLAYER2;i++)
       dv2[i][j]=0.0;
  
}


/***************************************************************/

void result1(char *file){
  int i,j;
  FILE *fp;
  char filenm[25]="./", suffix[6]=".neu";
  char true_exp;
  
  //strncat(filenm,file,5);
  strcat(filenm,file);
  strcat(filenm,suffix);
  if((fp=fopen(filenm,"w"))==NULL) {
    printf("cannot open file %s\n",filenm);
    return;
  }
  fprintf(fp,"%5d\n",Len);
  
  j=0; 
  for(i=1;i<=Len;i++){
    prdA2[i]=Max(pout[i],MLAYER3)+48; //MLAYER3=2
    fprintf(fp,"%5d %c  %5.2f  %5.2f  %c\n",i,StrA1[i],pout[i][0],
	    pout[i][1],prdA2[i]);
  }
  
  fclose(fp);
}

/*************************************************************/
void change(float *w,float *d,float *D)
{
   float ww,dd,DD;

   ww=*w; dd=*d; DD=*D;

   DD=LAMDA*DD + ALPHA*dd;
   ww+=DD;

  *w=ww;  *D=DD;
}

/***************************************************************/

int Max(float *p, int N){
  int i,n;
  float tmp=-9999.9;

  for(i=0;i<N;i++)
    if(*(p+i)>tmp){
      tmp=*(p+i);
      n=i;
    }
  return(n);
}





















