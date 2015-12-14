/*********************************************************************
*                              Q V O L                               *
**********************************************************************
* Program is part of the Situs package (c) Willy Wriggers, 1998-2012 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    * 
* Vector quantization of volumetric maps.                            *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

       
#include "situs.h"
#include "lib_std.h"
#include "lib_pio.h"
#include "lib_vio.h"
#include "lib_vwk.h"
#include "lib_rnd.h"

#define NCONSTR 1000                 /* max. number of constraints  */
#define SEED 7566                    /* for random number generator */
#define SHKPRE 10                    /* number of shake pre steps   */
#define SHKMAX 20000                 /* max. number of SHAKE iter.  */
#define FLENGTH 1000  /* file name length */

/* function declarations */

static int compar (int *, int *);

/* global variables */
static RseqNN dist;
static Rseq3 vglob;
static Rmat3NN wglob;

/* the main-program */

int main(int argc, char *argv[])

{ int cycle, skip, psfmode, done, conncount, shakemaxflag = 0; 
  char psf_file [FLENGTH];
  int ch1, ch2;
  RseqNN effrad;
  FILE *fout;
  unsigned char cconn[NNMAX][NNMAX];
  int numshake, nn;
  int pair1shake[NCONSTR];
  int pair2shake[NCONSTR];
  double distshake[NCONSTR];
  double currm, currnm;
  int opt, nextopt;
  unsigned numP;
  PDB *pdbP;
  int shakecount;
  double *phi;
  int *pid;
  int indx, indy, indz;
  unsigned long indv, count, nvox, extxy;  
  double dist0;
  double dtemp;
  double epsilon, lambda;
  IseqNN tk, order;
  double ei, ef, li, lf;
  double maxstep;
  int maxdata=0;
  unsigned extx, exty, extz;
  double origx, origy, origz, width;  
  double cutoff, randoff;
  int menumode, distmode, methmode, connmode;
  double mindensity, maxdensity;
  double wdiff;
  int i,j,k;
  int found, itercount;
  double mindist, currdist;
  int currindex=0;
  double xcum, ycum, zcum, rmscum, currmass, eeps;
  double xdev, ydev, zdev;
  Rmat3DT vecdata;  
  ImatNNDT list;
  IseqNN listcount;
  double shakescale, condiff;
  double w10,w11,w12,w20,w21,w22,w30,w31,w32,w40,w41,w42,w50,w51,w52,w60,w61,w62; 
  Rmat3NN average, wlast;
  RseqNN mcl, rmsd;
  Rseq3 com;
  double rgyr;
  int sprec;
  double distpre[NCONSTR];
  FILE *fin;  
  char con_file [FLENGTH];
  char line[101];

  if (argc < 3 || argc > 4) {
    fprintf(stderr, "...qvol> Usage: qvol input-map [optional: input-vectors] output-vectors \n");
    exit (1);
  }
  
  read_vol(argv[1], &width, &origx, &origy, &origz, &extx, &exty, &extz, &phi);
  nvox = extx*exty*extz;
  extxy = extx*exty;
  pid = (int *) malloc(nvox * sizeof(int));
  if (pid == NULL) {
    fprintf(stderr, "...qvol> Error: Unable to satisfy memory allocation request [e.c. 30210]\n"); 
    exit(30210);
  }

  printf ("...qvol> Density values below a user-defined cutoff value will not be considered\n");
  printf ("...qvol> Do you want to inspect the input density values before entering the cutoff value?\n");
  printf ("...qvol> Choose one of the following three options -\n");
  printf ("...qvol>      1: No (continue)\n");
  printf ("...qvol>      2: Show me the minimum and maximum density values only \n");
  printf ("...qvol>      3: Show me the voxel histogram \n");
  printf ("...qvol> ");
  menumode = readln_int();
 
  maxdensity = calc_max(phi, nvox);
  mindensity = calc_min(phi, nvox);

  switch (menumode) {

  case 1:
    
    break;
    
  case 2:

    printf ("...qvol> Min. / max. density values: %f / %f\n", mindensity, maxdensity);
    break;

  case 3:

    if (maxdensity <= mindensity) {
      printf ("...qvol> Min. / max. density values: %f / %f\n", mindensity, maxdensity);
      break;
    } else print_histogram(&extx, &exty, &extz, &phi, 40);
    break;
    
  default:
    fprintf(stderr, "...qvol> Error: Unable to identify option [e.c. 30230]\n"); 
    exit(30230);   
  }

  printf ("...qvol> Now enter the cutoff density value: "); 
  cutoff = readln_double();
  
  indv = nvox;
  for (count=0;count<nvox;count++) if (*(phi+count) < cutoff) {
    *(phi+count)=0;
    --indv;
  }
  printf("...qvol> Cutting off density values < %f, remaining occupied volume: %ld voxels (%e Angstrom^3)\n", cutoff, indv, (indv*width*width*width));  

  /* normalizing (dividing) by maxdensity */
  for (count=0;count<nvox;count++) *(phi+count) /= maxdensity;
   
  if (argc == 4) { 
    read_pdb(argv[2], &numP, &pdbP);
    nn = numP;
  } else {
    printf ("...qvol> Enter desired number of codebook vectors: "); 
    nn = readln_int();
    if (nn>NNMAX) {
      nn=NNMAX; 
      printf ("...qvol> Note: maximum number of codebook vectors (defined in situs.h): %d\n",nn); 
    }
    if (nn<NNMIN) {
      nn=NNMIN; 
      printf ("...qvol> Note: minimum number of codebook vectors (defined in situs.h): %d\n",nn); 
    }
  }


  /* optimization method to be used (methmode): */
  /* 1: Use TRN (global search) method (no distance constraints) */
  /* 2: Use LBG (gradient descent) method (distance constraints optional) */
  /* 3: Proceed directly to analysis/output */

  if (argc == 4 ) {
    printf ("...qvol> Do you want to optimize the start vectors or skip and proceed to the connectivity analysis?\n");
    printf ("...qvol> Choose one of the following two options -\n");
    printf ("...qvol>      1: Optimize start vectors with LBG\n");
    printf ("...qvol>      2: Skip and proceed directly to connectivity analysis \n");
    printf ("...qvol> ");
    connmode = readln_int();
    switch (connmode) {
    case 1:
      methmode = 2;
      break;
    case 2:  
      methmode = 3;
      break;
    default:
      fprintf(stderr, "...qvol> Error: Unable to identify option [e.c. 30231]\n"); 
      exit(30231);   
    }
  }  
  else methmode = 1;
 
  /* initializing start vectors */
  printf("...qvol> \n");
  if (argc == 4 ) {
    printf("...qvol> Using start vectors from file %s.\n",argv[2]);
    for (i=0; i<nn; ++i) {
      wglob[i][0] = (pdbP[i].x-origx)/width;
      wglob[i][1] = (pdbP[i].y-origy)/width;
      wglob[i][2] = (pdbP[i].z-origz)/width;
    }
    for (i=0; i<nn; ++i) {
      if (wglob[i][0]<0 || wglob[i][0] > extx || wglob[i][1]<0 || wglob[i][1] > exty || wglob[i][2]<0 || wglob[i][2] > extz) {
	fprintf(stderr, "...qvol> Error: Start vectors from file %s are not compatible with map from file %s. [e.c. 30260]\n",argv[2],argv[1]); 
	exit(30260);
      }
    }
  } else {
    printf("...qvol> Using random start vectors.\n");
    for (i=0; i<nn; ++i) {
      found = 0;
      do {
	indv = genrand()*nvox;
	if (*(phi+indv) > genrand()) found = 1; }
      while (found != 1); 
      randoff = 0.1-0.2*genrand();
      indz = indv / extxy;
      indv -= indz * extxy;
      indy = indv / extx;
      indx = indv - indy * extx;
      wglob[i][0] = indx+randoff;
      wglob[i][1] = indy+randoff;
      wglob[i][2] = indz+randoff;
    }
  }

  /* enter constraints, if any, and set numshake */
  if (methmode == 2) {
    printf ("...qvol> \n");
    printf ("...qvol> Vector distance constraints restrict undesired degrees of freedom.\n");
    printf ("...qvol> Do you want to add distance constraints?\n");
    printf ("...qvol> Choose one of the following three options -\n");
    printf ("...qvol>      1: No\n");
    printf ("...qvol>      2: Yes. I want to enter them manually\n");
    printf ("...qvol>      3: Yes. I want to read connectivities from a PSF file and use start vector distances\n");   
    printf ("...qvol>      4: Yes. I want to read them from a Situs constraints file\n");
    printf ("...qvol> ");
    distmode = readln_int();
  } else distmode = 1;
  
  switch (distmode) {

  case 1:
    
    numshake=0;
    break;

  case 2:

    printf ("...qvol> Enter a maximum of %d constraints.\n",NCONSTR);
    
    for (numshake=0;numshake<NCONSTR;++numshake) {
      printf ("...qvol> Distance constraint %d. Enter index (1-%d) of first vector (0 to finish, -n to backtrack n constraints): ",numshake+1,nn);
      pair1shake[numshake] = readln_int();
      if (pair1shake[numshake] == 0) break;
      else if (pair1shake[numshake] < 0) {
	numshake += pair1shake[numshake]-1;
	if (numshake < -1) numshake = -1;
	continue;
      }
      else if (pair1shake[numshake] > nn) {
	printf ("...qvol> Index must be within (1-%d). Please repeat.\n",nn);
	--numshake;
	continue;
      }
      printf ("...qvol> Distance constraint %d. Enter index (1-%d) of second vector: ",numshake+1,nn);
      pair2shake[numshake] = readln_int();      
      if (pair2shake[numshake] <= 0 || pair2shake[numshake] > nn) {
	printf ("...qvol> Index must be within (1-%d). Please repeat.\n",nn);
	--numshake;
	continue;
      }
      printf ("...qvol> Distance constraint %d. Enter distance between vectors %d and %d in Angstrom: ",numshake+1,pair1shake[numshake],pair2shake[numshake]);
      distshake[numshake] = readln_double();      
      if (distshake[numshake] <= 0) {
	printf ("...qvol> Distance must be > 0. Please repeat.\n");
	--numshake;
	continue;
      }
    }      
    for (i=0; i<numshake; ++i) {
      distshake[i] *= distshake[i]/(width*width);
      pair1shake[i] -= 1;
      pair2shake[i] -= 1;
    }
    break;

  case 3:

    /* read filename and open */ 

    printf("...qvol> Enter filename: ");
    if(fgets(con_file,FLENGTH,stdin)==NULL) {
      fprintf(stderr, "...qvol> Error: Can't read filename [e.c. 30400]\n"); 
      exit(30400);   
    }
    removespaces(con_file,FLENGTH);
    fin = fopen(con_file, "r"); 
    if(fin == NULL) {
      fprintf(stderr, "...qvol> Error: Can't open file %s [e.c. 30410]\n",con_file); 
      exit(30410);   
    }

    for(;;) {
      fgets(line, 100, fin);
      if (feof(fin)) break;
      if(strstr(line,"NBOND")!=NULL) { 
      sscanf(line,"%d",&numshake);
      break;
      }
    }
 
    for(j=0;j<numshake;++j) {
      if(j == NCONSTR || fscanf (fin, "%d", (pair1shake+j)) != 1) {
	fprintf(stderr, "...qvol> Error: Can't initiate reading %d. connectivity entry in file %s [e.c. 30419]\n",j+1, con_file);
	if (j==NCONSTR) fprintf(stderr, "...qvol> Increase NCONSTR in qvol.c and recompile qvol.\n");
	exit(30419);
      }
      if (fscanf (fin, "%d", (pair2shake+j)) != 1) {
	fprintf(stderr, "...qvol> Error: Can't complete reading %d. connectivity entry in file %s [e.c. 30420]\n",numshake+1, con_file); 
	exit(30420);
      }
       --pair1shake[j]; --pair2shake[j]; 
      distshake[j] = ((pdbP[pair1shake[j]].x-pdbP[pair2shake[j]].x)*
		      (pdbP[pair1shake[j]].x-pdbP[pair2shake[j]].x)+
		      (pdbP[pair1shake[j]].y-pdbP[pair2shake[j]].y)*
		      (pdbP[pair1shake[j]].y-pdbP[pair2shake[j]].y)+
		      (pdbP[pair1shake[j]].z-pdbP[pair2shake[j]].z)*
		      (pdbP[pair1shake[j]].z-pdbP[pair2shake[j]].z)) / (width*width);
    }

    printf("...qvol> %d connectivities read from file %s\n", numshake, con_file);
    printf("...qvol> The corresponding distances were assigned from file %s\n", argv[2]);
    fclose(fin);
    break;

  case 4:

    /* read filename and open */ 

    printf("...qvol> Enter filename: ");
    if(fgets(con_file,FLENGTH,stdin)==NULL) {
      fprintf(stderr, "...qvol> Error: Can't read filename [e.c. 30400]\n"); 
      exit(30400);   
    }
    removespaces(con_file,FLENGTH);
    fin = fopen(con_file, "r"); 
    if(fin == NULL) {
      fprintf(stderr, "...qvol> Error: Can't open file %s [e.c. 30410]\n",con_file); 
      exit(30410);   
    }

    numshake = 0;
    for(;;) {
      if(numshake == NCONSTR || fscanf (fin, "%d", (pair1shake+numshake)) != 1) break;
      if (fscanf (fin, "%d %le", (pair2shake+numshake), &dtemp) != 2) {
	fprintf(stderr, "...qvol> Error: Can't complete reading %d. constraint entry in file %s [e.c. 30420]\n",numshake+1, con_file); 
	exit(30420);
      }
      distshake[numshake] = dtemp * dtemp / (width*width);
      --pair1shake[numshake]; --pair2shake[numshake]; 
      numshake++;
    }
    printf("...qvol> %d constraints read from file %s\n", numshake, con_file);
    rewind(fin); for (i=0;i<numshake;i++) fscanf(fin, "%d %d %le", &j, &k, &dtemp);
    if (fscanf (fin, "%d %d %le", &j, &k, &dtemp) == 3 && numshake == NCONSTR) {
      fprintf(stderr, "...qvol> Error: There are additional (unread) constraints in file %s [e.c. 30430]\n", con_file);
      fprintf(stderr, "...qvol> Edit file %s to remove any redundant constraints and restart program.\n", con_file);
      fprintf(stderr, "...qvol> Constraints are redundant e.g. if the connectivity network is over-determined.\n");
      fprintf(stderr, "...qvol> If you rule out this possibility, try to increase the NCONSTR constant in qvol.c.\n");
      exit(30430);
    }
    fclose(fin);
    break;

  default:
    fprintf(stderr, "...qvol> Error: Unable to identify option [e.c. 30500]\n"); 
    exit(30500);   
  }
 

  if (numshake > 0) { /* SHAKE preconditioning */

    for (sprec=1; sprec<=SHKPRE; ++sprec) {
      for (k=0; k<numshake; ++k) {
	currdist = sqrt((wglob[pair1shake[k]][0]-wglob[pair2shake[k]][0])*(wglob[pair1shake[k]][0]-wglob[pair2shake[k]][0]) + (wglob[pair1shake[k]][1]-wglob[pair2shake[k]][1])*(wglob[pair1shake[k]][1]-wglob[pair2shake[k]][1]) + (wglob[pair1shake[k]][2]-wglob[pair2shake[k]][2])*(wglob[pair1shake[k]][2]-wglob[pair2shake[k]][2])); 
	distpre[k] = ((SHKPRE-sprec) / (double)SHKPRE) * currdist + (sprec/(double)SHKPRE) * distshake[k];
      }

      shakecount = 0;
      do {   /* loop until constraints are converged or SHKMAX reached */ 
	condiff = 0;
	for (k=0; k<numshake; ++k) {
	  i = pair1shake[k];
	  j = pair2shake[k];
	  w10 = wglob[i][0];
	  w11 = wglob[i][1];
	  w12 = wglob[i][2];
	  w20 = wglob[j][0];
	  w21 = wglob[j][1];
	  w22 = wglob[j][2];
	  w50 = w20-w10; w51 = w21-w11; w52 = w22-w12;
	  shakescale = (distpre[k] - w50*w50-w51*w51-w52*w52) /
	    (4.0 * (w50*w50+w51*w51+w52*w52));
	  wglob[i][0] -= shakescale * w50;
	  wglob[i][1] -= shakescale * w51;
	  wglob[i][2] -= shakescale * w52;
	  wglob[j][0] += shakescale * w50;
	  wglob[j][1] += shakescale * w51;
	  wglob[j][2] += shakescale * w52;
	  shakescale = fabs(distpre[k] - w50*w50-w51*w51-w52*w52) / distpre[k];
	  if (shakescale > condiff) condiff = shakescale;
	}
	++ shakecount;
      } while (condiff > 1e-3 && shakecount < SHKMAX);
      if (shakecount == SHKMAX && sprec==SHKPRE) shakemaxflag = 1;
      printf("...qvol> Distance preconditioning step %d -- %d SHAKE distance iterations \n",sprec,shakecount);
    }

    if (shakemaxflag){
      printf("...qvol> \n");
      printf("...qvol> Warning: SHKMAX parameter exceeded. Distances did not converge!\n"); 
      printf("...qvol> Try to remove any redundant constraints and restart program.\n");
      printf("...qvol> Constraints are redundant e.g. if the connectivity network is\n");
      printf("...qvol> over-determined. This case results in competition among the constraints.\n");
      printf("...qvol> If you rule out this possibility, try to increase the SHKMAX constant in qvol.c.\n");
      printf("...qvol> If this doesn't work either, start with a smaller number of well-chosen\n");
      printf("...qvol> constraints and add constraints incrementally until you've identified the cause\n"); 
      printf("...qvol> of the convergence problem.\n");
      printf("...qvol> \n");
      shakemaxflag = 0;
    }
  }
  

  switch (methmode) { /* choose optimization algorithm */
    
  case 3:  /* none */
    
    maxdata = nn;
    for (i=0; i<nn; ++i) for (j=0; j<3; ++j) vecdata[i][j] = wglob[i][j];
    break;
    
  case 1:  /* TRN */
    
    maxstep = SMAXS;
    maxdata = nn*MAXCYCLE;
    ei = 0.1;      
    ef = 0.001;   
    li = 0.2 * nn;
    lf = 0.02;    
    sgenrand(SEED);     

    printf("...qvol> Computing %d datasets, %d iterations each... \n", MAXCYCLE, (int)maxstep); 
    
    /* tessellate and store found codebook vectors in vecdata array */
    
    for (cycle=0; cycle<MAXCYCLE; ++cycle) {     
      printf("...qvol> Now producing dataset %d\n",cycle+1);
      
      for (i=0; i<nn; ++i) {
	found = 0;
	do {
	  indv = genrand()*nvox; 
	  if (*(phi+indv) > genrand()) found = 1; }
	while (found != 1); 
	randoff = 0.1-0.2*genrand();
	indz = indv / extxy;
	indv -= indz * extxy;
	indy = indv / extx;
	indx = indv - indy * extx;
	wglob[i][0] = indx+randoff;
	wglob[i][1] = indy+randoff;
	wglob[i][2] = indz+randoff;
      }
      
      for (count=0; count<maxstep; ++count) {
	
	found = 0;
	do {
	  indv = genrand()*nvox; 
	  if (*(phi+indv) > genrand()) found = 1; }
	while (found != 1); 
	
	indz = indv / extxy;
	indv -= indz * extxy;
	indy = indv / extx;
	indx = indv - indy * extx;
	vglob[0] = indx;
	vglob[1] = indy;
	vglob[2] = indz;
	
	for (i=0; i<nn; ++i) order[i]=i;
	qsort(order,nn,sizeof(int),(int (*)(const void *, const void *))compar); 
	for (i=0; i<nn; ++i) tk[order[i]] = i;	
	epsilon = ei * exp((count/maxstep)*log(ef/ei));
	lambda =  li * exp((count/maxstep)*log(lf/li));
	for (i=0; i<nn; ++i) for (j=0; j<3; ++j)
	  wglob[i][j] = wglob[i][j] + epsilon * exp (-tk[i]/lambda) * (vglob[j] - wglob[i][j]); 
      }   

      for (i=0; i<nn; ++i) for (j=0; j<3; ++j) vecdata[nn*cycle+i][j] = wglob[i][j];
    }
    break;

  case 2: /* LBG with optional constraints */

    maxdata = nn;

    printf("...qvol> Starting standard LBG vector quantization.\n");
    eeps = 0.1; /* plasticity, soft LBG, empirical */
    itercount = 0;
    
    do {   /* loop until LBG converged */ 
      
      ++itercount;
      
      /* determine clusters */ 
      for (count=0; count<nvox; ++count) if ((*(phi+count))>0) {
	indv = count;
	indz = indv / extxy;
	indv -= indz * extxy;
	indy = indv / extx;
	indx = indv - indy * extx;
	mindist = 1e20;
	for (i=0; i<nn; ++i) {
	  currdist = ((wglob[i][0]-indx)*(wglob[i][0]-indx)+(wglob[i][1]-indy)*(wglob[i][1]-indy)+(wglob[i][2]-indz)*(wglob[i][2]-indz));
	  if (currdist < mindist) {
	    mindist = currdist;
	    currindex = i;
	  }
	}
	*(pid+count)= currindex;
      }
      
      /* compute cluster centers */
      for (i=0; i<nn; ++i) for (j=0; j<3; ++j) average[i][j] = 0;
      for (i=0; i<nn; ++i) mcl[i] = 0;
      for (count=0; count<nvox; ++count) if ((*(phi+count))>0) {
	indv = count;
	indz = indv / extxy;
	indv -= indz * extxy;
	indy = indv / extx;
	indx = indv - indy * extx;
	currindex = (*(pid+count)); 
	currmass = (*(phi+count));
	average[currindex][0] += currmass * indx;
	average[currindex][1] += currmass * indy;
	average[currindex][2] += currmass * indz;
	mcl[currindex] += currmass;
      }
      
      for (i=0; i<nn; ++i) for (j=0; j<3; ++j) {
	if (mcl[i]==0) {
	  fprintf(stderr, "...qvol> Error: Voronoi cell %d contains no density [e.c. 30520]\n",i+1);
	  if (numshake>0) fprintf(stderr, "...qvol> Check start vectors and constraint distances.\n");
	  else if (argc==4) fprintf(stderr, "...qvol> Check start vectors.\n");
	  exit(30520);   
	} else average[i][j] /= mcl[i];
      }   
      
      /* update wglob[i][j] without constraints */
      for (i=0; i<nn; ++i) for (j=0; j<3; ++j) {
	wlast[i][j] = wglob[i][j];
	wglob[i][j] += eeps * (average[i][j] - wglob[i][j]);
      } 
      
      if (numshake > 0) {
	
	shakecount = 0;
	do {   /* loop until constraints are converged or SHKMAX reached */ 
	  condiff = 0;
	  for (k=0; k<numshake; ++k) {
	    i = pair1shake[k];
	    j = pair2shake[k]; 
	    w10 = wglob[i][0];
	    w11 = wglob[i][1];
	    w12 = wglob[i][2];
	    w20 = wglob[j][0];
	    w21 = wglob[j][1];
	    w22 = wglob[j][2];
	    w30 = wlast[i][0];
	    w31 = wlast[i][1];
	    w32 = wlast[i][2];
	    w40 = wlast[j][0];
	    w41 = wlast[j][1];
	    w42 = wlast[j][2];
	    w50 = w20-w10; w51 = w21-w11; w52 = w22-w12;
	    w60 = w40-w30; w61 = w41-w31; w62 = w42-w32;
	    shakescale = (distshake[k] - w50*w50-w51*w51-w52*w52) /
	      (4.0 * (w50*w60+w51*w61+w52*w62));
	    wglob[i][0] -= shakescale * w60;
	    wglob[i][1] -= shakescale * w61;
	    wglob[i][2] -= shakescale * w62;
	    wglob[j][0] += shakescale * w60;
	    wglob[j][1] += shakescale * w61;
	    wglob[j][2] += shakescale * w62;
	    shakescale = fabs(distshake[k] - w50*w50-w51*w51-w52*w52) / distshake[k];
	    if (shakescale > condiff) condiff = shakescale;
	  }
	  ++ shakecount;
	} while (condiff > 1e-4 && shakecount < SHKMAX);
	if (shakecount == SHKMAX) shakemaxflag = 1;
	printf("...qvol> It. %d -- %d SHAKE distance iterations \n",itercount,shakecount);    
      } 
  
      /* compute stopping criterion */
      wdiff = 0;
      for (i=0; i<nn; ++i) for (j=0; j<3; ++j) {
	wdiff += ((wglob[i][j]-wlast[i][j])*(wglob[i][j]-wlast[i][j]));
      } 
      wdiff = width*sqrt(wdiff/(double)nn);
      fprintf(stderr, "...qvol> It. %d -- Average vector update: %e Angstrom \n",itercount,wdiff); 
    } while (wdiff > 1e-2);   

    if (shakemaxflag){
      printf("...qvol> \n");
      printf("...qvol> Warning: SHKMAX parameter exceeded. Distances did not converge!\n"); 
      printf("...qvol> Try to remove any redundant constraints and restart program.\n");
      printf("...qvol> Constraints are redundant e.g. if the connectivity network is\n");
      printf("...qvol> over-determined. This case results in competition among the constraints.\n");
      printf("...qvol> If you rule out this possibility, try to increase the SHKMAX constant in qvol.c.\n");
      printf("...qvol> If this doesn't work either, start with a smaller number of well-chosen\n");
      printf("...qvol> constraints and add constraints incrementally until you've identified the cause\n"); 
      printf("...qvol> of the convergence problem.\n");
      printf("...qvol> \n");
      shakemaxflag = 0;
    }

    fprintf(stderr, "...qvol> \n");
    for (i=0; i<nn; ++i) for (j=0; j<3; ++j) vecdata[i][j] = wglob[i][j];
    break;
    
  } /* switch (methmode) */

  for (j=0; j<3; ++j) for (i=0; i<nn; ++i)  wglob[i][j] = vecdata[i][j];

  do { /* loop until clustering converged */ 
    
    for (i=0; i<nn; ++i) listcount[i] = 0;
    
    /* determine clusters */
    for (j=0; j<maxdata; ++j) { 
      mindist = 1e20;
      for (i=0; i<nn; ++i) {
	currdist = (wglob[i][0]-vecdata[j][0])*(wglob[i][0]-vecdata[j][0]) + 
	  (wglob[i][1]-vecdata[j][1])*(wglob[i][1]-vecdata[j][1]) + 
	  (wglob[i][2]-vecdata[j][2])*(wglob[i][2]-vecdata[j][2]);
	if (currdist < mindist) {
	  mindist = currdist;
	  currindex = i;
	}
      }
      list[currindex][listcount[currindex]] = j;
      ++listcount[currindex];
    }

    /* compute cluster center */
    for (currindex=0; currindex<nn; ++currindex) {
      xcum = 0; ycum = 0; zcum = 0;
      for (i=0; i<listcount[currindex]; ++i) {
	xcum += vecdata[list[currindex][i]][0];
	ycum += vecdata[list[currindex][i]][1];
	zcum += vecdata[list[currindex][i]][2];
      }
      if (listcount[currindex]>0) { 
	average[currindex][0] = xcum / (double)listcount[currindex];
	average[currindex][1] = ycum / (double)listcount[currindex];    
	average[currindex][2] = zcum / (double)listcount[currindex];
      }
    }

    /* update */
    for (i=0; i<nn; ++i) {
      for (j=0; j<3; ++j) {
	wlast[i][j] = wglob[i][j];
	wglob[i][j] = average[i][j];
      }
    }
 
    if (numshake > 0) {
      
      shakecount = 0;
      do {   /* loop until constraints are converged or SHKMAX reached */ 
	condiff = 0;
	for (k=0; k<numshake; ++k) {
	  i = pair1shake[k];
	  j = pair2shake[k];
	  w10 = wglob[i][0];
	  w11 = wglob[i][1];
	  w12 = wglob[i][2];
	  w20 = wglob[j][0];
	  w21 = wglob[j][1];
	  w22 = wglob[j][2];
	  w30 = wlast[i][0];
	  w31 = wlast[i][1];
	  w32 = wlast[i][2];
	  w40 = wlast[j][0];
	  w41 = wlast[j][1];
	  w42 = wlast[j][2];
	  w50 = w20-w10; w51 = w21-w11; w52 = w22-w12;
	  w60 = w40-w30; w61 = w41-w31; w62 = w42-w32;
	  shakescale = (distshake[k] - w50*w50-w51*w51-w52*w52) /
	    (4.0 * (w50*w60+w51*w61+w52*w62));	
	  wglob[i][0] -= shakescale * w60;
	  wglob[i][1] -= shakescale * w61;
	  wglob[i][2] -= shakescale * w62;
	  wglob[j][0] += shakescale * w60;
	  wglob[j][1] += shakescale * w61;
	  wglob[j][2] += shakescale * w62;
	  shakescale = fabs(distshake[k] - w50*w50-w51*w51-w52*w52) / distshake[k];
	  if (shakescale > condiff) condiff = shakescale;
	}
	++shakecount;
      } while (condiff > 1e-5 && shakecount < SHKMAX);
      if (shakecount == SHKMAX) shakemaxflag = 1;
      printf("...qvol> Final clustering -- %d SHAKE distance iterations \n",shakecount); 
    }

    /* compute stopping criterion */
    wdiff = 0;
    for (i=0; i<nn; ++i) for (j=0; j<3; ++j) {
      wdiff += ((wglob[i][j]-wlast[i][j])*(wglob[i][j]-wlast[i][j]));
    } 
    wdiff = width*sqrt(wdiff/(double)nn);
    fprintf(stderr, "...qvol> Final clustering -- Average vector update: %e Angstrom \n",wdiff);
  } while (wdiff > 1e-5);   
  
  if (shakemaxflag){
    printf("...qvol> \n");
    printf("...qvol> Warning: SHKMAX parameter exceeded. Distances did not converge!\n"); 
    printf("...qvol> Try to remove any redundant constraints and restart program.\n");
    printf("...qvol> Constraints are redundant e.g. if the connectivity network is\n");
    printf("...qvol> over-determined. This case results in competition among the constraints.\n");
    printf("...qvol> If you rule out this possibility, try to increase the SHKMAX constant in qvol.c.\n");
    printf("...qvol> If this doesn't work either, start with a smaller number of well-chosen\n");
    printf("...qvol> constraints and add constraints incrementally until you've identified the cause\n"); 
    printf("...qvol> of the convergence problem.\n");
    printf("...qvol> \n");
    shakemaxflag = 0;
  }
  
  /* variability within each cluster */
  
  for (currindex=0; currindex<nn; ++currindex) {
    rmscum = 0;
    for (i=0; i<listcount[currindex]; ++i) {
      xdev = vecdata[list[currindex][i]][0] - wglob[currindex][0]; 
      ydev = vecdata[list[currindex][i]][1] - wglob[currindex][1]; 
      zdev = vecdata[list[currindex][i]][2] - wglob[currindex][2]; 
      rmscum += xdev*xdev+ydev*ydev+zdev*zdev;
    }
    if (listcount[currindex] > 0) {
      rmscum /=  (double)listcount[currindex];
      rmsd[currindex] = sqrt(rmscum);
    } else {  /* rec. field contains zero */
      rmsd[currindex] = -1; 
      printf("...qvol> Warning: Voronoi cell %d contains no data points.\n",currindex+1); 
      printf("...qvol> Vector variability in output file %s will be set to 99.99.\n",argv[argc-1]); 
    }
  }
  
  /* write output to pdb file */

  fout = fopen(argv[argc-1], "w");
  if( fout == NULL ) {
    fprintf(stderr, "...qvol> Error: Can't open file! %s  [e.c. 30700]\n", argv[argc-1]); 
    exit(30700);   
  }
  
  /* compute equivalent spherical radius */
      
  /* determine clusters */ 
  for (count=0; count<nvox; ++count) if ((*(phi+count))>0) {
    indv = count;
    indz = indv / extxy;
    indv -= indz * extxy;
    indy = indv / extx;
    indx = indv - indy * extx;
    mindist = 1e20;
    for (i=0; i<nn; ++i) {
      currdist = ((wglob[i][0]-indx)*(wglob[i][0]-indx)+(wglob[i][1]-indy)*(wglob[i][1]-indy)+(wglob[i][2]-indz)*(wglob[i][2]-indz));
      if (currdist < mindist) {
	mindist = currdist;
	currindex = i;
      }
    }
    (*(pid+count))= currindex;
  }
  
  /* compute effective radius of each cluster */
  for (i=0; i<nn; ++i) effrad[i] = 0;
  for (i=0; i<nn; ++i) mcl[i] = 0;
  for (count=0; count<nvox; ++count) if ((*(phi+count))>0) {
    indv = count;
    indz = indv / extxy;
    indv -= indz * extxy;
    indy = indv / extx;
    indx = indv - indy * extx;
    currindex = (*(pid+count)); 
    currmass = (*(phi+count));
    effrad[currindex] += currmass * ( (indx-wglob[currindex][0])*(indx-wglob[currindex][0]) + 
				      (indy-wglob[currindex][1])*(indy-wglob[currindex][1]) + 
				      (indz-wglob[currindex][2])*(indz-wglob[currindex][2]) );
    mcl[currindex] += currmass;     
  }

  for (i=0; i<nn; ++i) {
    if (mcl[i]==0) {
      fprintf(stderr, "...qvol> Error: Voronoi cell %d contains no density [e.c. 30800]\n",i+1); 
      if (numshake>0) fprintf(stderr, "...qvol> Check start vectors and constraint distances.\n");
      else if (argc==4) fprintf(stderr, "...qvol> Check start vectors.\n");
      exit(30800);   
    }	else effrad[i] = sqrt(5.0*effrad[i]/(3.0*mcl[i])); 
  }   
  
  skip = 0;
  for(i=0;i<nn;++i) {    
    if (rmsd[i] < 0) {
      skip = 1;
      fprintf(fout, "ATOM  %5d QVOL QVOL%5d %11.*f%8.*f%8.*f %5.2f %5.2f      QVOL\n",i+1,i+1,coord_precision(wglob[i][0]*width+origx),wglob[i][0]*width+origx,coord_precision(wglob[i][1]*width+origy),wglob[i][1]*width+origy,coord_precision(wglob[i][2]*width+origz),wglob[i][2]*width+origz,99.99,effrad[i]*width);
    } else {
      fprintf(fout, "ATOM  %5d QVOL QVOL%5d %11.*f%8.*f%8.*f %5.2f %5.2f      QVOL\n",i+1,i+1,coord_precision(wglob[i][0]*width+origx),wglob[i][0]*width+origx,coord_precision(wglob[i][1]*width+origy),wglob[i][1]*width+origy,coord_precision(wglob[i][2]*width+origz),wglob[i][2]*width+origz,rmsd[i]*width,effrad[i]*width);
    }
  }
  fclose (fout);
  
  printf("...qvol>\n");
  printf("...qvol> Codebook vectors have been written to file %s\n", argv[argc-1]);  
  printf("...qvol> The PDB B-factor field contains the equivalent spherical radii\n");
  printf("...qvol> of the corresponding Voronoi cells (in Angstrom).\n");
  if (maxdata>nn) {
    printf("...qvol> Cluster analysis of the %d independent runs:\n",MAXCYCLE);
    printf("...qvol> The PDB occupancy field in %s contains the rms variabilities of the vectors.\n",argv[argc-1]);
    if (skip == 0) {
      rmscum = 0;
      for(i=0;i<nn;++i) rmscum += rmsd[i];
      rmscum /= (double)nn;
      printf("...qvol> Average rms fluctuation of the %d codebook vectors: %6.3f Angstrom\n", nn, rmscum*width);
    } else printf("...qvol> Unable to compute average rms fluctuation, cluster analysis did not converge\n");
  }

  /* write statistics of cluster size deviations*/
  k = 0;
  for(i=0;i<nn;++i) k += abs(listcount[i]-maxdata/nn);  
  if (k > 0) {
    printf("...qvol> \n");
    printf("...qvol> Warning: There are cluster size deviations:\n");
    printf("...qvol> Expected cluster size: %d\n",maxdata/nn);
    printf("...qvol> Actual cluster sizes: ");
    for(i=0;i<(nn-1);++i) printf("%d,",listcount[i]);
    printf("%d.\n",listcount[nn-1]);
  }
 
  /* compute com and rgyr of the codebook vectors */
  for (j=0;j<3;++j) com[j]=0; 
  for (i=0;i<nn;++i) for (j=0;j<3;++j) com[j]+=wglob[i][j];
  for (j=0;j<3;++j) com[j] /= (double)nn;
  rgyr=0;
  for (i=0;i<nn;++i) for (j=0;j<3;++j) rgyr += (wglob[i][j]-com[j])*(wglob[i][j]-com[j]);
  rgyr /= (double)nn;
  rgyr = sqrt(rgyr);
  printf("...qvol> Radius of gyration of the %d codebook vectors: %6.3f Angstrom\n", nn, rgyr*width);

  /* compute or write connectivities if desired */
  
  printf ("...qvol> \n");
  printf ("...qvol> Do you want to update or save the input connectivities?\n");
  printf ("...qvol> Choose one of the following options -\n");
  printf ("...qvol>      1: No. I'm done \n");
  printf ("...qvol>      2: Update and save to a PSF file\n");
  printf ("...qvol>      3: Update and save to a constraints file\n");
  printf ("...qvol>      4: Update and save to both PSF and constraints files\n");
  if (numshake > 0) printf ("...qvol>      5: Just save (don't update) to a PSF file \n");
  
  printf ("...qvol> ");
  psfmode = readln_int();
  
  /* compute connectivities, if necessary */
  switch (psfmode) {
  
  case 2: case 3: case 4:

    for (j=0; j<nn; ++j) for (i=0; i<nn; ++i) cconn[i][j] = 0;      
    for (count=0; count<nvox; ++count) {
      if (*(phi+count) > 0) {
	indv = count;
	indz = indv / extxy;
	indv -= indz * extxy;
	indy = indv / extx;
	indx = indv - indy * extx;
	for (i=0; i<nn; ++i) dist[i] = (indx-wglob[i][0])*(indx-wglob[i][0]) +
			       (indy-wglob[i][1])*(indy-wglob[i][1]) +
			       (indz-wglob[i][2])*(indz-wglob[i][2]);
	currm = 1e20;
	currnm = 1e20;
	nextopt = 0;
	opt = 0;
	for (i=0;i<nn;i++) {
	  if (dist[i]<=currm) {
	  currnm = currm;
	  nextopt = opt;
	  currm = dist[i];
	  opt = i;
	  } else if (dist[i]<=currnm) {
	    currnm = dist[i];
	    nextopt = i;
	  }   
	}
	cconn[opt][nextopt] = 1; 
	cconn[nextopt][opt] = 1;
      }
    }
    break;
  }


  /* write files, if necessary */
  switch (psfmode) {
    
  case 1: default:
    
    printf ("...qvol> Bye bye!\n");
    break;
    
  case 2: case 4:

    /* read filename and open */
    for (done=0;done==0;) {      
      printf("...qvol> Enter PSF filename: ");
      if(fgets(psf_file,FLENGTH,stdin)==NULL) {
	fprintf(stderr, "...qvol> Error: Can't read filename [e.c. 30900]\n"); 
	exit(30900);   
      }
      removespaces(psf_file,FLENGTH);
      fout = fopen(psf_file, "r"); 
      if(fout != NULL) {
	printf("...qvol> Warning: File exists: %s \n", psf_file); 
	printf("...qvol> Do you want to overwrite? (yes/no) "); 
	ch1=getchar();
	for (;;)  {
	  ch2 = getchar(); 
	  if (ch2 == EOF ) {
	    fprintf(stderr, "...qvol> Error: EOF while reading input [e.c. 30910]\n");
	    exit (30910);
	  }
	  if (ch2 == '\n') break;
	} 
	if (ch1=='y' || ch1=='Y') done = 1;
      } else {
	done = 1;
      }
    }

    fout = fopen(psf_file, "w");
    fprintf(fout, "PSF \n");
    fprintf(fout, " \n");
    fprintf(fout, "       1 !NTITLE\n");
    fprintf(fout, " REMARK Vector connectivity psf file, computed with Competitive Hebbian Rule \n");
    fprintf(fout, " \n");
    fprintf(fout, "%8d !NATOM\n",nn);
    for(i=0;i<nn;++i) fprintf(fout, "%8d QVOL %-4d QVOL QVOL QVOL   0.000000       0.00000           0\n", i+1, i+1);
    fprintf(fout, " \n");
    conncount=0;
    for (j=0; j<nn; ++j) for (i=0; i<j; ++i) if (cconn[i][j] == 1) ++conncount; 
    fprintf(fout, "%8d !NBOND: bonds\n", conncount);
    k = 0;
    for (j=0; j<nn; ++j) for (i=0; i<j; ++i) if (cconn[i][j] == 1) {
      if (k == 4) {
        k = 0;
        fprintf(fout, " \n");
      }
      ++k;
      fprintf(fout, "%8d%8d",i+1,j+1);
    }
    if (k > 0) fprintf(fout, "\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NTHETA: angles\n");
    fprintf(fout, "\n");   
    fprintf(fout, "       0 !NPHI: dihedrals\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NIMPHI: impropers\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NDON: donors\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NACC: acceptors\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NNB\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0       0 !NGRP\n");
    fprintf(fout, "\n");
    fclose (fout);
    printf("...qvol> Connectivity data written to PSF file %s.\n",psf_file);
    if (psfmode == 2) break;
    
  case 3:
    
    /* read filename and open */
    for (done=0;done==0;) {      
      printf("...qvol> Enter connectivity filename: ");
      if(fgets(psf_file,FLENGTH,stdin)==NULL) {
	fprintf(stderr, "...qvol> Error: Can't read filename [e.c. 30930]\n"); 
	exit(30930);   
      }
      removespaces(psf_file,FLENGTH);
      fout = fopen(psf_file, "r"); 
      if(fout != NULL) {
	printf("...qvol> Warning: File exists: %s \n", psf_file); 
	printf("...qvol> Do you want to overwrite? (yes/no) "); 
	ch1=getchar();
	for (;;)  {
	  ch2 = getchar(); 
	  if (ch2 == EOF ) {
	    fprintf(stderr, "...qvol> Error: EOF while reading input [e.c. 30940]\n");
	    exit (30940);
	  }
	  if (ch2 == '\n') break;
	} 
	if (ch1=='y' || ch1=='Y') done = 1;
      } else {
	done = 1;
      }
    }

    fout = fopen(psf_file, "w");

    for (j=0; j<nn; ++j) for (i=0; i<j; ++i) if (cconn[i][j] == 1) {
      dist0 = (wglob[i][0]-wglob[j][0])*(wglob[i][0]-wglob[j][0]) +
	(wglob[i][1]-wglob[j][1])*(wglob[i][1]-wglob[j][1]) +
	(wglob[i][2]-wglob[j][2])*(wglob[i][2]-wglob[j][2]);
      dist0 = sqrt(dist0)*width;
      fprintf(fout, "%d %d %f\n",i+1,j+1,dist0);
    }
    fclose (fout);
    printf("...qvol> Connectivity data / distances for all vectors written to file %s.\n",psf_file);
    break;    
  
  case 5:
    
    if (numshake == 0) { printf ("...qvol> No constraints loaded. Bye bye!\n"); exit(1); }
    /* read filename and open */
    for (done=0;done==0;) {      
      printf("...qvol> Enter PSF filename: ");
      if(fgets(psf_file,FLENGTH,stdin)==NULL) {
	fprintf(stderr, "...qvol> Error: Can't read filename [e.c. 30950]\n"); 
	exit(30950);   
      }
      removespaces(psf_file,FLENGTH);
      fout = fopen(psf_file, "r"); 
      if(fout != NULL) {
	printf("...qvol> Warning: File exists: %s \n", psf_file); 
	printf("...qvol> Do you want to overwrite? (yes/no) "); 
	ch1=getchar();
	for (;;)  {
	  ch2 = getchar(); 
	  if (ch2 == EOF ) {
	    fprintf(stderr, "...qvol> Error: EOF while reading input [e.c. 30960]\n");
	    exit (30960);
	  }
	  if (ch2 == '\n') break;
	} 
	if (ch1=='y' || ch1=='Y') done = 1;
      } else {
	done = 1;
      }
    }

    fout = fopen(psf_file, "w");
    fprintf(fout, "PSF \n");
    fprintf(fout, " \n");
    fprintf(fout, "       1 !NTITLE\n");
    fprintf(fout, " REMARK Vector connectivity psf file, converted from constraints \n");
    fprintf(fout, " \n");
    fprintf(fout, "%8d !NATOM\n",nn);
    for(i=0;i<nn;++i) fprintf(fout, "%8d QVOL %-4d QVOL QVOL QVOL   0.000000       0.00000           0\n", i+1, i+1);
    fprintf(fout, " \n");
    fprintf(fout, "%8d !NBOND: bonds\n", numshake);
    k = 0;
    for (i=0; i<numshake; ++i) {
      if (k == 4) {
        k = 0;
        fprintf(fout, " \n");
      }
      ++k;
      fprintf(fout, "%8d%8d",pair1shake[i]+1,pair2shake[i]+1);
    }
    if (k > 0) fprintf(fout, "\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NTHETA: angles\n");
    fprintf(fout, "\n");   
    fprintf(fout, "       0 !NPHI: dihedrals\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NIMPHI: impropers\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NDON: donors\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NACC: acceptors\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0 !NNB\n");
    fprintf(fout, "\n");
    fprintf(fout, "       0       0 !NGRP\n");
    fprintf(fout, "\n");
    fclose (fout);
    printf("...qvol> Connectivity data written to PSF file %s.\n",psf_file);
    break;
  }
  return (0);
}

static int compar (int *e1, int *e2)
/* distance comparison function for qsort */

{ double ax, ay, az, d1, d2;

  ax = vglob[0]-wglob[*e1][0];
  ay = vglob[1]-wglob[*e1][1];
  az = vglob[2]-wglob[*e1][2];
  d1 = ax*ax + ay*ay + az*az;
  ax = vglob[0]-wglob[*e2][0];
  ay = vglob[1]-wglob[*e2][1];
  az = vglob[2]-wglob[*e2][2];
  d2 = ax*ax + ay*ay + az*az;
  if (d1<d2) return -1;
  else return 1;
}

