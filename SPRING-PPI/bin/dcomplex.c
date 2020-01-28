


/**********************************
 *  
 *  force centroid: all heavy atom
 *  area: interface
 *  bin: hongyi
 *  ref stat: hongyi
 *  
 *
 *  by Z.C. @ 3/21/2003
 *
 **********************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



#define PDB_LIST_LENG	2850
#define RES_TYPE	20
#define MAX_ATOMS_IN_ONE_RESIDUE	1
#define MAXA		10000  //total atom number in one protein.
#define BINSIZE		20



struct residueinf{
       char location;
       int residue;
       int atomnum;
       int atoms[20];
       int oninterface;
};
typedef struct residueinf RESIDUEINF;

struct atominf{
	double X;
	double Y;
	double Z;
        //char pdbline[90];
};
typedef struct atominf ATOMINF;



//*********************** function list ***********

//read charge file
void read_charge_file(char * filename);
//hash function
int hashfunc(char * resnm);
//read a wholeline from a file
int readline(FILE *fp, char *line);
// exit when meet an error
void stop(char *s );
void readPDBfile(FILE * fp, int *ires1, int *ires2 );

//*********************** END **********************


//********************** GLobal data  **************
int  residhash[101];  //residue hashtable
char resnm[RES_TYPE][4];  //residue name
int  resatmnum[RES_TYPE];
char cind[RES_TYPE][MAX_ATOMS_IN_ONE_RESIDUE][5];
float chg[RES_TYPE][MAX_ATOMS_IN_ONE_RESIDUE];
char residname3[20][4]={"CYS", "MET", "PHE", "ILE", "LEU",
	                "VAL", "TRP", "TYR", "ALA", "GLY",
			"THR", "SER", "GLN", "ASN", "GLU",
			"ASP", "HIS", "ARG", "LYS", "PRO" };
char residname1[21]="CMFILVWYAGTSQNEDHRKP";

RESIDUEINF * iResidue1;
RESIDUEINF * iResidue2;
ATOMINF * iAtom;

int RecptChainNum; 
char RecptChain[10];

int LegandChainNum; 
char LegandChain[10];

//********************** END  **********************





int main(int argc, char *argv[]){
   

    
    if(argc!=4){
       printf("usage: %s filename chain_IDs_1 chain_IDs_2\n", argv[0]);
       exit(2);
    }

    FILE *fp;
    char filename[80], dum[100];
    char afil[PDB_LIST_LENG][40];
     int    ret, 
	nfil, // total number of PDB files
	resid_id, mm ;  
        
    // read the charge information from file
    // buid residue hash table and
    // atom name dictionary for future using.
    //scanf("%s", filename);
    read_charge_file("/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/bin/charge_inp.dat");



    //read chain and decide the group of chain
    //scanf("%s %s", RecptChain, LegandChain);
    //RecptChainNum=1;
    //LegandChainNum=1;
    //printf("%c %c\n", RecptChain[0], LegandChain[0]);
    //
    //read chain and decide the group of chain
    strcpy(RecptChain, argv[2]);
    RecptChainNum=strlen(argv[2]);
    for(int i=0; i<RecptChainNum; i++){
      if(RecptChain[i]=='_') RecptChain[i]=' ';
        //printf("%c\n", RecptChain[i]);
    }
    
    strcpy(LegandChain, argv[3]);
    LegandChainNum=strlen(argv[3]);
    for(int i=0; i<LegandChainNum; i++){
      if(LegandChain[i]=='_') LegandChain[i]=' ';
       // printf("%c\n", LegandChain[i]);
    }
    

    // read potential file
    double * pot; // contact number, i.e. N(i,j,r).
    pot=  (double  *) malloc(27*20*15*20*15*sizeof(double));
    memset(pot, 0, 27*20*15*20*15*sizeof(double));
    //fp=(FILE *)fopen("fort.21_alla_hy","r");
    //fp=(FILE *)fopen("fort.21_alla_hy_1.61","r"); //dimeric 1.61
    //fp=(FILE *)fopen("fort.21","r");
    fp=(FILE *)fopen("/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/bin/fort.21_alla","r"); //monomic 1.61
    if(fp==NULL) stop("Could not find potential file: fort.21_alla");
    char ctmp1[5], ctmp2[5], ctmp3[5], ctmp4[5];
    double U;
    int ll=0;
    do{
	int  m, i,j, k,l;
        ret=fscanf(fp, "%s %s %s %s %lf\t%d\t%d %d %d %d\n",
		ctmp1, ctmp2, ctmp3, ctmp4, &U, &m, &i, &j, &k, &l);
	// printf("%s %s %s %s %f %d %d %d %d %d\n",
	// ctmp1, ctmp2, ctmp3, ctmp4, U, m, i, j, k, l);
	if(ret==EOF) break;
	int D1=20*15*20*15;
	int D2=20*15*15;
	int D3=15*15;
	int D4=15;
	pot[m*D1+i*D2+k*D3+j*D4+l]=U;
	pot[m*D1+k*D2+i*D3+l*D4+j]=U;
	//printf("%s %s %s %s %f %d %d %d %d %d\n",
	// ctmp1, ctmp2, ctmp3, ctmp4, pot[m*D1+i*D2+k*D3+j*D4+l], m, i, j, k, l);
	ll++; 
    }while(ret!=EOF );
    // printf("Potential Lines=%d\n", ll);
    
                                    
    // make a map array for future using.
    // map arry decides bin location.
    // refs. Samudrala and Moult (1998)
    //       Lu and Skolnick (2001)
    int *map;
    map=(int *) malloc(700*sizeof(int));
    for(int i=1; i<500; i++){
        if(i<4) map[i]=1;
        if(i>=4 && i<16) map[i]=i-2;
        if(i>=16 && i< 50) map[i]=(int)((i)*0.5)+6;
        if(i>=50 ) map[i]=-1;
    }
    

    /* read pdb file 
     * if a protein has more than one chain, we read 
     * one chain once time. For the next chain, reopen the
     * pdb  file and read again.
     */
     iResidue1=(RESIDUEINF *) malloc( 2000*sizeof(RESIDUEINF));
     if(iResidue1==NULL) stop("iResidue");
     memset(iResidue1, 0, 2000*sizeof(RESIDUEINF));

     iResidue2=(RESIDUEINF *) malloc( 2000*sizeof(RESIDUEINF));
     if(iResidue2==NULL) stop("iResidue");
     memset(iResidue2, 0, 2000*sizeof(RESIDUEINF));


     iAtom=(ATOMINF *) malloc(MAXA*sizeof(ATOMINF));
     if(iAtom==NULL) stop("iAtom");
     memset(iAtom, 0, MAXA*sizeof(ATOMINF));

     double fTotalEnergy=0;
     int    iNativeIndex;

     fp=(FILE *)fopen(argv[1], "r");
     if(fp==NULL){
         perror("PBD Files");
         exit(1);
     }


     int ires1, ires2;
     readPDBfile(fp, &ires1, &ires2 );

     //printf("ires1=%d , ires2=%d\n", ires1, ires2);
     fclose(fp);
	
	
	int ri, rj, i0,i1, j0,j1;
        int jj; //bin index
	double rd;

	//find the residues on interface
	int find=0;
	for(int i=0; i<ires1; i++){
	for(int atm1=0; atm1<iResidue1[i].atomnum ; atm1++){
	    for(int j=0; j<  ires2; j++){
	    for(int atm2=0; atm2<iResidue2[j].atomnum ; atm2++){
	      //printf("%d,%d\n",iResidue1[i].residue,iResidue2[j].residue);
	        int id1=iResidue1[i].atoms[atm1];
                int id2=iResidue2[j].atoms[atm2];		 

		//printf("%d,%d\n",id1,id2);
		//calculate the distance based on the coordinates.
		rd=sqrt( (iAtom[id1].X-iAtom[id2].X)
		         *(iAtom[id1].X-iAtom[id2].X)
		         +(iAtom[id1].Y-iAtom[id2].Y)
		         *(iAtom[id1].Y-iAtom[id2].Y)
		         +(iAtom[id1].Z-iAtom[id2].Z)
		         *(iAtom[id1].Z-iAtom[id2].Z));
		if(rd < 10.0){
		//if(rd < 10000000){
		   find=1;
		   iResidue1[i].oninterface=8;
		   iResidue2[j].oninterface=8;
		   //printf("%d %d %s %s\n", i,j, resnm[iResidue1[i].residue],
		   //resnm[iResidue2[j].residue]);
		}

	    }
	    }
	}
	}
	
	


	//count the contacts
	for(int i=0; i< ires1; i++){
	if(iResidue1[i].oninterface!=8) continue;
	for(int ia=0; ia <iResidue1[i].atomnum; ia++){
	    int id1=iResidue1[i].atoms[ia];
            for(int j=0; j<  ires2; j++){
            if(iResidue2[j].oninterface!=8) continue;
            for(int ib=0; ib <iResidue2[j].atomnum; ib++){
                //id of alpha carbone
	        int id2=iResidue2[j].atoms[ib];

		//printf("X: %.3f %.3f\nY: %.3f %.3f\nZ: %.3f %.3f\n",iAtom[id1].X,iAtom[id2].X,iAtom[id1].Y,iAtom[id2].Y,iAtom[id1].Z,iAtom[id2].Z);
	        //calculate the distance based on the coordinates.
                rd=sqrt(
	         (iAtom[id1].X-iAtom[id2].X)
	        *(iAtom[id1].X-iAtom[id2].X)
		+(iAtom[id1].Y-iAtom[id2].Y)
		*(iAtom[id1].Y-iAtom[id2].Y)
		+(iAtom[id1].Z-iAtom[id2].Z)
	        *(iAtom[id1].Z-iAtom[id2].Z));
		//printf("%.3f\n",rd);
		
	        if(rd > 30.0) continue;
	        jj=map[(int)(rd*2.0)];
	        if(jj>0 && jj <= 20){
	           i0=iResidue1[i].residue;
	           j0=iResidue2[j].residue;
	           int D0=20*20*15*15;
	           int D1=20*15*15;
		   int D2=15*15;
		   int D3=15;
	           double U=  pot[jj*D0+i0*D1+j0*D2+ia*D3+ib];
		   //printf("%d %d %d %d %d %.3f %.5f\n",jj,i0,j0,ia,ib,rd,U);
		   //if(U > 8) U=2.0; //NOTE here!!!!
		   fTotalEnergy+=U;
                   //          if (jj<16)  printf("m=%d, %d, %d, %d, %d, POT=%f\n",
		//		   jj,i0, ia, j0, ib, U);
	        } 
	    } 
	    } 
        }
	} 


     printf("%f\n", fTotalEnergy*0.0157-4.7);


    free(iResidue1);
    free(iResidue2);
    free(iAtom);

    return(0); 
}


void readPDBfile(FILE * fp, int *ires1, int *ires2 ){
	char wholeline[90];
	int  lines=0, imdl=0, id=0, ime=-1;
        char dupl, ares1[6], ares0[6];
	char temp[10]; // for tempery use
	int ret;
	int index1, index2;
		
	char CurntChain='1';
	int CurtChainStat=-1;
	int avalible=0;

        memset(ares1, 0, sizeof(ares1));
        memset(ares0, 0, sizeof(ares0));
        *ires1=-1;
        *ires2=-1;
        id=0;
	int atomnum=0;

	//read lines
	//we need take care the PDB file format.
	do{
	    lines++;
	    memset(wholeline, 0, sizeof(wholeline));
	    ret=readline(fp, wholeline);

	    // check if there are more than two MODEL 
	    if(strncmp(wholeline, "MODEL", 5)==0) imdl++;
	    if(imdl >= 2) break;

	    //we only need atom line.
	    if(strncmp(wholeline, "ATOM", 4)!=0) continue;

	    //Check the current line is for heavy atom
	    dupl= wholeline[13];
	    if(dupl!='C' && dupl!='O' && dupl!='N' && dupl!='S') continue;


	    //check if the current atom belong to the chain we need.
	    //if is not ATOM, continue.
	    if(wholeline[21]==CurntChain){
	        if( avalible!=1){
	              continue;
	        }
	    } else{  //not equal current chain type
                memset(ares0, 'X', sizeof(ares0));

	        int find=0;
	        for(int i=0; i< RecptChainNum; i++ ){
	          if(RecptChain[i]==wholeline[21]){
                     if(CurtChainStat==1) {
                        *ires2=*ires2+1;
                        iResidue1[*ires1-1].atomnum=atomnum;
                     }
	             find=1;
	             CurntChain=RecptChain[i];
	             CurtChainStat=0;
	             avalible=1;
	             //printf("find chain %c, ires=%d\n", CurntChain,*ires1); 
		  } 
		}
	        if(find==0){
	           for(int i=0; i< LegandChainNum; i++){
	               if(LegandChain[i]==wholeline[21] ){
                          if(CurtChainStat==0) {
                             *ires1=*ires1+1;
                             iResidue1[*ires1-1].atomnum=atomnum;
                          }
	                  find=1;
	                  CurntChain=LegandChain[i];
	                  CurtChainStat=1;
	                  avalible=1;
	                  //printf("find chain %c, ires=%d\n", CurntChain,*ires2); 
		       }
		   } 
		} 
		if(find==0){
	           CurntChain==wholeline[21];
                   CurtChainStat=-1;
	           avalible=0;
	           continue;
	        } 
	    }
	    


	    //now we take care the alternative location of an atom.
	    double prob=1.0;
            if(wholeline[16]!=' ') {
	       memset(temp, 0, sizeof(temp));
	       strncpy(temp, &wholeline[57], 4);
               prob=atof(temp);
	    }

	    
	    if(wholeline[16]!=' ')
		  if( wholeline[16]!='A' &&
		      wholeline[16]!='1' ) continue;
	    
	    //if(wholeline[16]!=' ' && wholeline[16]!='A' ) continue;
	    // if(wholeline[16]=='B' && fabs(prob-0.5) < 0.0001) continue;

	    char residname[4];
	    memset(residname, 0, sizeof(residname));
	    strncpy(residname, &wholeline[17], 3);
	    int position=hashfunc(residname);
	    index1=residhash[position];
	    if(index1<0){
	         // printf(" this residue is not in hash table\n");
		 continue;  // this residue is not in hash table.
	    }
	    if(strncmp(residname,"ASX",3)==0) continue;
         

	    // check if the current atom belong to a
	    // new residue
	    // ares1: the residue name in which the current atom is
	    // ares0: the residue name in which last atom is
	    strncpy(ares1, &wholeline[22], 5);
	    if(strncmp(ares1, ares0, 4)==0 &&
	       ares1[4]!=ares0[4] ) continue;
	    if(strncmp(ares1, ares0, 4)!=0){
		switch(CurtChainStat){
		    case 0: if(*ires1>=0){
		               *ires1=*ires1+1;
	                       iResidue1[*ires1].residue=index1;
		               iResidue1[*ires1-1].atomnum=atomnum;
		            }else{
		               *ires1=*ires1+1;
			       iResidue1[*ires1].residue=index1;
		            }
		            atomnum=0;
			    break;
		    case 1: if(*ires2>=0){
			       *ires2=*ires2+1;
			       iResidue2[*ires2].residue=index1;
			       iResidue2[*ires2-1].atomnum=atomnum;
			    }else{ 
			       *ires2=*ires2+1; 
			       iResidue2[*ires2].residue=index1;
			    }
			    atomnum=0;
			    break;
		}
	        strncpy(ares0, ares1, 5);
	    }
             
            //read residue name and atom name
            //and then retrieve their indexes.
	    char atomname[5];
            memset(atomname, 0, sizeof(atomname));
            strncpy(atomname, &wholeline[13], 3);
	    if(strncmp(atomname, "OT", 2)==0 ||
	       strncmp(atomname, "OXT", 3)==0) strncpy(atomname, "O   ",4);
            if(strncmp(residname, "ILE", 3)==0 &&
               strncmp(atomname, "CD1", 3)==0){
		    strncpy(atomname, "CD  ",4);
	    }

            int gotit=0;
	    for(int i=0; i<MAX_ATOMS_IN_ONE_RESIDUE; i++){
	       if(strncmp(atomname, cind[index1][i],3)==0){
		       index2=i;
		       gotit=1;
	       }
	    }
	    if(!gotit) {
	      printf("Can not find atom %s in residue %s\n", atomname, residname);
	      printf("%s\n", wholeline);
	      exit(3);
	    }

	    //if(index2==1 ) printf("%d %d %d %s\n",*ires, id, index2, atomname);
	    //printf("%d %d %d %d\n",CurtChainStat,*ires1,index2,id);
	    switch(CurtChainStat){
	           case 0: iResidue1[*ires1].atoms[index2]=id; break;
	           case 1: iResidue2[*ires2].atoms[index2]=id; break;
	    }
	    // read coordinates of the current atomes
	    memset(temp, 0, sizeof(temp));
	    strncpy(temp, &wholeline[30],8);
	    iAtom[id].X=atof(temp);
	    memset(temp, 0, sizeof(temp));
	    strncpy(temp, &wholeline[38],8);
	    iAtom[id].Y=atof(temp);
	    memset(temp, 0, sizeof(temp));
	    strncpy(temp, &wholeline[46],8);
	    iAtom[id].Z=atof(temp);
	    //strncpy(iAtom[id].pdbline, wholeline, 90);
	    //printf("%.3f,%.3f,%.3f\n",iAtom[id].X,iAtom[id].Y,iAtom[id].Z);
	    id++;
	    atomnum++;
	}while(ret!=EOF);

	 switch(CurtChainStat){

	    case 0: *ires1=*ires1+1;
		    iResidue1[*ires1-1].atomnum=atomnum;
		    atomnum=0;
		    break;
	    case 1: *ires2=*ires2+1;
		    iResidue2[*ires2-1].atomnum=atomnum;
		    atomnum=0;
		    break;
	 }
}



int readline(FILE *fp, char *line){
	int i, ret;

	i=-1;
	do{
	    i++;
	    ret=fscanf(fp, "%c", &line[i]);
	}while(ret!=EOF && line[i]!='\n');
	line[i]='\0';
	return(ret);
}


int hashfunc(char * resnm){
    int position;
        position=30*resnm[0]
		 +70*resnm[1]
		+ 7*resnm[2];
	position=position%101;
    return position;
}

void stop(char * s){
   perror(s);
   exit(3);
}

// read the charge information from file 
// buid residue hash table and
// atom name dictionary for future using.
void read_charge_file(char * filename){

    FILE * fp=(FILE *)fopen(filename, "r");
    if(fp==NULL)  stop("Charge Information");
    int position;
    memset(residhash, -1, sizeof(residhash)); //initial hash table.
    int resid_id=0;
    do{ 
      
       int ret, mm;
       ret=fscanf(fp, "%s %d\n", resnm[resid_id], &mm);
       if(ret==EOF) break;
       //printf("%s %d %s %c\n", resnm[resid_id], 
       //resid_id, residname3[resid_id], residname1[resid_id] );
       resatmnum[resid_id]=mm;
       position=hashfunc(resnm[resid_id]);
       residhash[position]=resid_id;
       for(int i=0; i< mm; i++){
          fscanf(fp, "%c", &cind[resid_id][i][0]);
          fscanf(fp, "%c", &cind[resid_id][i][1]);
          fscanf(fp, "%c", &cind[resid_id][i][2]);
          fscanf(fp, "%c", &cind[resid_id][i][3]);
          fscanf(fp, "%c", &cind[resid_id][i][4]);
	  cind[resid_id][i][4]='\0';
          fscanf(fp, "%f\n", &chg[resid_id][i]);
          //printf("%d %s %f\n",i, cind[resid_id][i], chg[resid_id][i]);
       } 
       resid_id++;
    }while(1);
    fclose(fp);
    return;
}
