/**************************************************
MIT License

Copyright (c) 2019 Julie C. Mitchell and Oak Ridge National Laboratory

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

**************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define MAXATOM 35000
#define MAXRES 15000
#define MAXCON 2250000
#define DEBUG 0

float  xyz[2][MAXATOM][3];
char   res[2][MAXATOM][20];
char   atm[2][MAXATOM][20];
char   ren[2][MAXATOM][20];
char   atn[2][MAXATOM][20];
char   chn[2][MAXATOM][20];
int    num[2];

int    minatom[MAXATOM];

int    connum[MAXCON][2];
float  convec[MAXCON][4];
float  dmin, dmax;

float  mindist[MAXRES][MAXRES],smindist[MAXRES][MAXRES];

 void ExtractString(int len, char *src, char *dst )
{
    register char *ptr;
    register char ch;
    register int i;

    ptr = dst;
    for( i=0; i<len; i++ )
    {   if( *src )
	{   ch = *src++;
            *dst++ = ch;
            if( ch != ' ' ) 
		ptr = dst;
	} else break;
    }
    *ptr = 0;
    
}

/* read a pdb file into global vars */

int pdbread(char *pdbfile, int which)
{
  FILE *fp;
  int i, n;
  char sNewString[132];
  char sCoord[5][20];

  /* Open input file. */
  
  fp=fopen(pdbfile,"r");
  
  if (fp == NULL){
  	fprintf(stderr,"Problem opening file %s.\n",pdbfile);
  	exit(1);
  }
  
  /* Loop over input lines. */

  n = 0;
  
   while (fgets(sNewString,132,fp) != NULL){
 	 if ( (strncmp(sNewString,"ATOM",4) == 0) || (strncmp(sNewString,"HETATM",6) == 0) ) {

		ExtractString(8,&sNewString[30],sCoord[0]);
	 	ExtractString(8,&sNewString[38],sCoord[1]);
		ExtractString(8,&sNewString[46],sCoord[2]);
 
 		ExtractString(5,&sNewString[6],atn[which][n]);
		ExtractString(2,&sNewString[13],atm[which][n]);
		ExtractString(3,&sNewString[17],res[which][n]);
		ExtractString(1,&sNewString[21],chn[which][n]);
 		ExtractString(4,&sNewString[22],ren[which][n]);
                
               if (DEBUG)
                    fprintf(stdout,"%5s %2s %3s %1s\n",ren[which][n],atm[which][n],res[which][n],chn[which][n]); 

                for (i=0;i<3;i++) {
			xyz[which][n][i] = atof(sCoord[i]);
                        if (DEBUG) fprintf(stdout,"%f ",xyz[which][n][i]);
                }
                if (DEBUG) fprintf(stdout,"\n");
                
		n++;
                
		if (n >= MAXATOM) {
			fprintf(stderr,"You have exceeded the maximum allowed number of atoms.\n");
			exit(1);
		}
	}	 
  }
  
  fclose(fp);   
  
return n;}

/* get automatic naming convention based on molecule names */

void getOutName(char *outName, char *sName, char *mName)
{
    char sNameT[300], mNameT[300], sTemp[300];
    int period, slash;
    
/*
    strcpy(sNameT,sName);
    period = strcspn(sNameT,".");
    if (period < strlen(sNameT)) strcpy(&sNameT[period],"");
    while (strcspn(sNameT,"/") < strlen(sNameT)) {
        slash = strcspn(sNameT,"/");
        strcpy(sTemp,&sNameT[slash+1]);
        strcpy(sNameT,sTemp);
    }
    
    strcpy(outName,sNameT);
    
    strcpy(mNameT,mName);
    period = strcspn(mNameT,".");
    if (period < strlen(mNameT)) strcpy(&mNameT[period],"");
    while (strcspn(mNameT,"/") < strlen(mNameT)) {
        slash = strcspn(mNameT,"/");
        strcpy(sTemp,&mNameT[slash+1]);
        strcpy(mNameT,sTemp);
    }
    
    strcat(outName,"_");
    strcat(outName,mNameT);
*/ 

	strcpy(outName,mName);

return;}

/* MAIN PROGRAM */

int main(int argc, char **argv){

    int i,j,k,ii,jj;
    int nco;
    int smash=0;
    float dx,dy,dz,d;
    char outName[300], tempName[300];
    FILE *mFile, *muFile, *fFile, *fuFile, *bFile, *buFile, *smFile;
    
    /* process command line call */
    
    if (argc < 5){
        fprintf(stderr,"Usage contacts m1.pdb m2.pdb rmin rmax \n");
        exit(1);
    }
    
    num[0] = pdbread(argv[1],0);
    num[1] = pdbread(argv[2],1);
    dmin = atof(argv[3]);
    dmax = atof(argv[4]);
    
    if (argc == 6) smash = 1;
    
    /* open output files */
    
    getOutName(outName,argv[1],argv[2]);
    
//    strcpy(tempName,outName);    strcat(tempName,".mc");     mFile = fopen(tempName,"w");
//    strcpy(tempName,outName);    strcat(tempName,".muc");    muFile = fopen(tempName,"w");
    strcpy(tempName,outName);    strcat(tempName,".fc");     fFile = fopen(tempName,"w");    
    strcpy(tempName,outName);    strcat(tempName,".fuc");    fuFile = fopen(tempName,"w");
//    strcpy(tempName,outName);    strcat(tempName,".bc");     bFile = fopen(tempName,"w");    
//    strcpy(tempName,outName);    strcat(tempName,".buc");    buFile = fopen(tempName,"w");
	if (smash){
			strcpy(tempName,outName);    strcat(tempName,".sm");     smFile = fopen(tempName,"w");
	}
        
    /* initialize some things */
    
    nco = 0;    
    for (i=0;i<num[0];i++){
        for (j=0;j<num[1];j++){
            mindist[i][j] = dmax + 1.0;
        }
    }

    for (k=0;k<MAXATOM;k++){
        minatom[k] = 0;
    }
    
    /* compute distances, store information and print full contact lists */
    
    for (i=0;i<num[0];i++){
        for (j=0;j<num[1];j++){
        
            dx = (xyz[0][i][0] - xyz[1][j][0]);
            dy = (xyz[0][i][1] - xyz[1][j][1]);
            dz = (xyz[0][i][2] - xyz[1][j][2]);
           
            d = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (d <= dmax){
                connum[nco][0] = i;
                connum[nco][1] = j;
                convec[nco][0] = dx;
                convec[nco][1] = dy;
                convec[nco][2] = dz;
                convec[nco][3] = d;
                nco++;

                ii = atoi(ren[0][i]);
                jj = atoi(ren[1][j]);
                
                if (d < mindist[ii][jj]){
                    mindist[ii][jj] = d;
                }
                
                fprintf(fFile,"%5s %5s %5s %5s %5s   %5s %5s %5s %5s %5s   %8.3f\n", \
                    res[0][i],ren[0][i],chn[0][i],atm[0][i],atn[0][i],res[1][j],\
                    ren[1][j],chn[1][j],atm[1][j],atn[1][j],d);
/*
				if (d <= dmin){
                    fprintf(bFile,"%5s %5s %5s %5s %5s   %5s %5s %5s %5s %5s   %8.3f\n", \
                        res[0][i],ren[0][i],chn[0][i],atm[0][i],atn[0][i],res[1][j],\
                        ren[1][j],chn[1][j],atm[1][j],atn[1][j],d);
                } else {
                    fprintf(mFile,"%5s %5s %5s %5s %5s   %5s %5s %5s %5s %5s   %8.3f\n", \
                        res[0][i],ren[0][i],chn[0][i],atm[0][i],atn[0][i],res[1][j],\
                        ren[1][j],chn[1][j],atm[1][j],atn[1][j],d);
                }
 */

            }
        }    
    }
    
    /* determine unique contacts and print those with smallest distance */
    
    for (k=0;k<nco;k++){
        i = connum[k][0];
        j = connum[k][1];
 	d = convec[k][3];
        ii = atoi(ren[0][i]);
        jj = atoi(ren[1][j]);

        if (d == mindist[ii][jj])  {

            if (d <= dmin) {
                minatom[i] = 1;
                minatom[j] = 1;
            }
            
            fprintf(fuFile,"%5s %5s %5s %5s %5s   %5s %5s %5s %5s %5s   %8.3f\n", \
                res[0][i],ren[0][i],chn[0][i],atm[0][i],atn[0][i],res[1][j],\
                ren[1][j],chn[1][j],atm[1][j],atn[1][j],d);

/*
			if (d <= dmin) {
                fprintf(buFile,"%5s %5s %5s %5s %5s   %5s %5s %5s %5s %5s   %8.3f\n", \
                    res[0][i],ren[0][i],chn[0][i],atm[0][i],atn[0][i],res[1][j],\
                    ren[1][j],chn[1][j],atm[1][j],atn[1][j],d);
            } else {
                fprintf(muFile,"%5s %5s %5s %5s %5s   %5s %5s %5s %5s %5s   %8.3f\n", \
                    res[0][i],ren[0][i],chn[0][i],atm[0][i],atn[0][i],res[1][j],\
                    ren[1][j],chn[1][j],atm[1][j],atn[1][j],d);
            }
 */
        }
    }

	
    /* set min distance to be the minimum of those with minatom[i]=minatom[j]=0 */
    
    for (i=0;i<num[0];i++){
        for (j=0;j<num[1];j++){
            smindist[i][j] = dmax + 1.0;
        }
    }

    for (k=0;k<nco;k++){
    
        i = connum[k][0];
        j = connum[k][1];
 	d = convec[k][3];
        ii = atoi(ren[0][i]);
        jj = atoi(ren[1][j]);
        
        if ( (dmin < mindist[ii][jj]) && (d < smindist[ii][jj]) && (!minatom[i]) && (!minatom[j]) ){
            smindist[ii][jj] = d;
        }
        
    }
    
   /* find ways to adjust atom positions */

/*
    for (k=0;k<nco;k++){
        i = connum[k][0];
        j = connum[k][1];
 	d = convec[k][3];
        ii = atoi(ren[0][i]);
        jj = atoi(ren[1][j]);

        if ( (d > dmin) && (d <= dmax) && (!minatom[i]) && (!minatom[j]) )  {
        
            if (d == mindist[ii][jj]) {
                fprintf(smFile,"%5s %5s %5s %5s %5s   %5s %5s %5s %5s %5s   %8.3f %8.3f %8.3f %8.3f\n", \
                    res[0][i],ren[0][i],chn[0][i],atm[0][i],atn[0][i],res[1][j],\
                    ren[1][j],chn[1][j],atm[1][j],atn[1][j],d,convec[k][0],convec[k][1],convec[k][2]);
            }
            
        }
    }
*/    
    /* close files ... we are done */
    
//    fclose(mFile);
//    fclose(muFile);
    fclose(fFile);
    fclose(fuFile);
//    fclose(bFile);
//    fclose(buFile);
    
    if (smash) fclose(smFile);
    
return 0;}








