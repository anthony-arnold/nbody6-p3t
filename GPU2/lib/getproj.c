#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include "fieldinfo.h"

#define PI 3.1415927

#define NMAX  1400000

#define GND  (90-27.128251)/180.0*PI
#define GNR  (90+192.859481)/180.0*PI
#define L0   32.931918/180.0*PI

#define GCDIST 8178.0    // Gravity collaboration

#define LSR_U   11.10    // Values from Schoenreich et al. (2010)
#define LSR_V   12.24
#define LSR_W    7.25

//#define LSR_U    9.44    // Values from Piskunov et al. 2006
//#define LSR_V   11.90
//#define LSR_W    7.20

#define VSUN  242.0        // km/sec, Irrgang

//NGC 3201
//#define RACL3201 154.40343
//#define DECL3201 -46.41248

// NGC 4590
//#define RACLCL 189.86658
//#define DECLCL -26.74406

//NGC 1261
#define RACLCL 48.06754
#define DECLCL -55.21622

void discard_len(FILE* fd) {
   static int LEN = 4;
   int i, c;
   for (i=0;i<LEN;i++) {
      if ((c=fgetc(fd))==EOF) {
         exit(-1);
      }
   }
}

int cmpmy(double *x1, double *x2) {
 if(*x1<*x2) return(-1);
 return(1);
}

double **t,**t2;

double *mass,x[NMAX],y[NMAX],z[NMAX],u[NMAX],v[NMAX],w[NMAX];
int nstar,*kstar;

int main(int argc, char *argv[]) {
 int i,k;
 char sgn[2];
 void show_data(),read_data(),get_tform_mat();
 void invert_tform_mat();
 double parallax,vr,vx,vy,vz,vn,ve;
 double pmra,pmdec;
 double b,l,l2,dis,ra,dec,ras,decs;
 int rah,ram,decg,decm;
 char fname[200];
 FILE *dat;
 int isin[NFIELD],findin();
 int nfieldtrue,assignf[20];
 double rashift, decshift ;
 void get_cen();
 int shift=0;

 if (argc<2) {
   printf("Usage: getproj filename time <shift?>\n");
   exit(-1);
 }

 if (argc==2) show_data(argv[1]);

 if (argc==4) shift=1;

 t = malloc(3*sizeof(double *));
 for (i=0;i<3;i++)
    t[i] = malloc(3*sizeof(double));

 get_tform_mat();

 invert_tform_mat(t);

 t2 = malloc(3*sizeof(double *));
 for (i=0;i<3;i++)
    t2[i] = malloc(3*sizeof(double));

 read_data(argv[1],argv[2]);

 get_cen(nstar, x, y, z, &rashift, &decshift, shift);

 strcpy(fname,argv[1]);
 if (!strcmp((fname+strlen(fname)-3),"POS")) fname[strlen(fname)-4]=0;

// Discriminate between different observations and split fields
  nfieldtrue=1;
  assignf[0]=0;

  for (i=1;i<NFIELD;i++)
   if (strcmp(namef[i],namef[i-1])) {
       nfieldtrue++;
       assignf[i] = assignf[i-1]+1;
   } else
       assignf[i] = assignf[i-1];

 strcat(fname,".out");

 dat = fopen(fname,"w");

 for (i=0;i<nstar;i++) {
    v[i] -= VSUN;
    x[i] = -(x[i]-GCDIST);

    x[i] /= 1000.0; y[i] /= 1000.0; z[i]/= 1000.0;

// Turn u around so it points towards GC
     u[i] = -u[i];

// Subtract solar motion relative to LSR
     u[i] -= LSR_U;
     v[i] -= LSR_V;
     w[i] -= LSR_W;

     vx = u[i]*t[0][0] + v[i]*t[0][1]+w[i]*t[0][2];
     vy = u[i]*t[1][0] + v[i]*t[1][1]+w[i]*t[1][2];
     vz = u[i]*t[2][0] + v[i]*t[2][1]+w[i]*t[2][2];

     dis  = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
     b = asin(z[i]/dis);
     l = acos(x[i]/dis/cos(b));
     l2= asin(y[i]/dis/cos(b));

     if (l2<0.0) l=2.0*PI-l;

     dec = asin(cos(b)*sin(l-L0)*sin(GND)+sin(b)*cos(GND));
     double x1 = cos(b)*cos(l-L0);
     double x2 = sin(b)*sin(GND)-cos(b)*cos(GND)*sin(l-L0);

     if (x2>=0.0)
        ra = (GNR-PI/2.0)+atan(x1/x2);
     else if (x1>=0.0)
        ra = (GNR-PI/2.0)+atan(x1/x2)+PI;
     else
        ra = (GNR-PI/2.0)+atan(x1/x2)-PI;

     if (ra>2.0*PI) ra=ra-2.0*PI;

     t2[0][0] =-cos(ra)*sin(dec);  t2[1][0] =-sin(ra);  t2[2][0] = cos(ra)*cos(dec);
     t2[0][1] =-sin(ra)*sin(dec);  t2[1][1] = cos(ra);  t2[2][1] = sin(ra)*cos(dec);
     t2[0][2] = cos(dec);          t2[1][2] = 0.0;      t2[2][2] = sin(dec);

     invert_tform_mat(t2);

     ra  = ra*180.0/PI;
     if (dec<0.0) strcpy(sgn,"-"); else strcpy(sgn," ");
     dec *= 180.0/PI;
     double decabs = fabs(dec);

     rah=floor(ra/15.0);
     ram=floor((ra/15.0-rah)*60.0);
     ras=((ra/15.0-rah)*60.0-ram)*60.0;

     decg=floor(decabs);
     decm=floor((decabs-decg)*60.0);
     decs=((decabs-decg)*60.0-decm)*60.0;

     vn = t2[0][0]*vx+t2[1][0]*vy+t2[2][0]*vz;
     ve = t2[0][1]*vx+t2[1][1]*vy+t2[2][1]*vz;
     vr = t2[0][2]*vx+t2[1][2]*vy+t2[2][2]*vz;

     parallax = 1.0/(1000.0*dis);

     pmra  = 1000.0*ve/4.7406*parallax;
     pmdec = 1000.0*vn/4.7406*parallax;

// Apply shift in RA/DEC or no shift
     ra += rashift;
     dec += decshift;

     for (k=0;k<nfieldtrue;k++) isin[k]=0;

     for (k=0;k<NFIELD;k++)
       isin[assignf[k]]=isin[assignf[k]]+findin(ra,dec,nlin[k],&ragrid[k][0],&decgrid[k][0]);

//   fprintf(dat,"%12i   %02i %02i %6.3f %1s%02i %02i %6.3f %7.1lf %7.2lf  0.00  %6.2lf  0.00  %6.2lf  0.00\n",i,rah,ram,ras,sgn,decg,decm,decs,dis*1000.0,vr,pmra,pmdec);
     fprintf(dat,"%12i %12.7f %12.7f %12.3lf %9.3lf %9.3lf %9.3lf %9.5lf %3i",i,ra,dec,dis*1000.0,pmra,pmdec,vr,mass[i],kstar[i]);
     for (k=0;k<nfieldtrue;k++) fprintf(dat,"%2i",isin[k]);
     fprintf(dat,"\n");
 }
 fclose(dat);
 exit(0);
}

void get_cen(int nstar, double *x, double *y, double *z, double *rashift, double *decshift, int shift) {
  double *rast,*decst,*disij,dens;
  double xst,yst,zst,dis,b,l,l2,x1,x2;
  double ramean,decmean,weight;
  int i,j,k,merk;
  double mind;

  rast  = malloc(nstar*sizeof(double));
  decst = malloc(nstar*sizeof(double));

  disij = malloc(nstar*sizeof(double));

  for (i=0;i<nstar;i++) {
    xst = -(x[i]-GCDIST)/1000.0;
    yst = y[i]/1000.0;
    zst = z[i]/1000.0;

     dis  = sqrt(xst*xst+yst*yst+zst*zst);
     b = asin(zst/dis);
     l = acos(xst/dis/cos(b));
     l2= asin(yst/dis/cos(b));

     if (l2<0.0) l=2.0*PI-l;

     decst[i] = asin(cos(b)*sin(l-L0)*sin(GND)+sin(b)*cos(GND));
     x1 = cos(b)*cos(l-L0);
     x2 = sin(b)*sin(GND)-cos(b)*cos(GND)*sin(l-L0);

     if (x2>=0.0)
        rast[i] = (GNR-PI/2.0)+atan(x1/x2);
     else if (x1>=0.0)
        rast[i] = (GNR-PI/2.0)+atan(x1/x2)+PI;
     else
        rast[i] = (GNR-PI/2.0)+atan(x1/x2)-PI;

     if (rast[i]>2.0*PI) rast[i]=rast[i]-2.0*PI;

     rast[i]  *= 180.0/PI;
     decst[i] *= 180.0/PI;
  }

  ramean  = 0.0;
  decmean = 0.0;
  weight  = 0.0;

  for (i=0;i<2000;i++) {
    for (j=0;j<nstar;j++)
       disij[j]=sqrt((rast[i]-rast[j])*(rast[i]-rast[j])+(decst[i]-decst[j])*(decst[i]-decst[j]));
    disij[i] = 1.E10;

    for (k=0;k<10;k++) {
       mind = disij[0];
       merk = 0;
       for (j=0;j<nstar;j++) {
         if (disij[j]<mind) {
           mind=disij[j];
           merk = j;
         }
       }
       disij[merk]=1E10;
     }

     dens = 1.0/(mind*mind*mind);

     ramean  += rast[i]*dens;
     decmean += decst[i]*dens;
     weight  += dens;
  }

  ramean /= weight;
  decmean /= weight;

  printf("Cluster center at ra=%lf dec=%lf\n",ramean,decmean);

  if (shift==1) {
    printf("Shifting to observed center of NGC 3201 !\n");
    *rashift  = RACLCL-ramean;
    *decshift = DECLCL-decmean;
  } else {
    printf("Keeping this center !\n");
    *rashift  = 0.0;
    *decshift = 0.0;
  }
}



void show_data(char *fnameu) {
   static const int KMAX = 20;
 FILE *dat;
 char fname[200];
 double *bla,tmyr,as[KMAX];
 int i,c,ntot,nk,buf[4];

 bla = malloc(NMAX*sizeof(double));

 strcpy(fname,fnameu);
 if (strcmp((fnameu+strlen(fnameu)-3),"POS") && strncmp((fnameu+strlen(fnameu)-4),"POS",3)) strcat(fname,".POS");
 dat = fopen(fname,"r");

 do {
    discard_len(dat); // Read 1st record length of 1st record
  fread(buf,4,4,dat);
  ntot=buf[0];
  if (ntot<1000 || ntot>1E6) {
     fprintf(stderr, "Read failed %i\n",ntot);
     exit(-1);
  }
  nk = buf[3];
  if (nk > KMAX) {
     fprintf(stderr, "Read failed nk=%i\n",nk);
     exit(-1);
  }

// Read 2nd record length of 1st record and 1st of 2nd
  discard_len(dat);
  discard_len(dat);

// Read data
  fread(as,8,nk,dat);         // Header data
  printf("%lf %i\n",as[0],ntot);

  fread(bla,8,ntot,dat);     // Masses
  fread(bla,8,3*ntot,dat);   // Positions
  fread(bla,8,3*ntot,dat);   // Velocities
  fread(bla,4,ntot,dat);     // Names

  discard_len(dat); // Read 2nd record length of 2nd record
 } while (1);
}


void read_data(char *fnameu, char *tsnaps) {
 FILE *dat;
 double tsnap;
 char fname[200];
 int ntot,ntot0=0,buf[3],c,i;
 float *lsev,*rsev;
 int *name;
 double *phi;
 double tmyr,as[30];
 double rbar,zmbar,vstar,dr[3],dv[3];
 double xbuf[3],rgal[3],vgal[3];
 double xx;
 double *xc,*yc,*zc,*uc,*vc,*wc;

 xc = malloc(NMAX*sizeof(double));
 yc = malloc(NMAX*sizeof(double));
 zc = malloc(NMAX*sizeof(double));

 uc = malloc(NMAX*sizeof(double));
 vc = malloc(NMAX*sizeof(double));
 wc = malloc(NMAX*sizeof(double));

 name=malloc(NMAX*sizeof(int));
 kstar=malloc(NMAX*sizeof(int));

 mass=malloc(NMAX*sizeof(double));
 phi=malloc(NMAX*sizeof(double));

 lsev=malloc(NMAX*sizeof(float));
 rsev=malloc(NMAX*sizeof(float));

 strcpy(fname,fnameu);
 if (strcmp((fnameu+strlen(fnameu)-3),"POS") && strncmp((fnameu+strlen(fnameu)-4),"POS",3))  strcat(fname,".POS");

 if ((dat=fopen(fname,"r"))==NULL) {
   printf("Opening POS-File Failed !\n");
   exit(-1);
 }

 tsnap=strtod(tsnaps,NULL);

 do {
  for (i=0;i<4;i++) if ((c=fgetc(dat))==EOF) {
     printf("No snapshot found for t=%lf !\n",tsnap);
     exit(-1);    // Read 1st record length of 1st record
  }
  fread(buf,4,3,dat);
  ntot=buf[0];
  if (ntot<1000 || ntot>1E6) {
     printf("Read failed %i\n",ntot);
     exit(-1);
  }

  if (ntot0==0) ntot0=ntot;

 for (i=0;i<8;i++) c=fgetc(dat);   // Read 2nd record length of 1st record and 1st of 2nd

// Read data
    fread(as,8,30,dat);         // Header data
    fread(mass,8,ntot,dat);     // Masses

    for (i=0;i<ntot;i++) {
       fread(xbuf,8,3,dat);     // Positions
       xc[i]=xbuf[0];
       yc[i]=xbuf[1];
       zc[i]=xbuf[2];
    }
    for (i=0;i<ntot;i++) {
       fread(xbuf,8,3,dat);     // Velocities
       uc[i]=xbuf[0];
       vc[i]=xbuf[1];
       wc[i]=xbuf[2];
    }

    if (as[29]<=0.5) {
      fread(name,4,ntot,dat);     // Names
    } else {
      fread(phi,8,ntot,dat);      // Potentials
      fread(name,4,ntot,dat);     // Names
      fread(kstar,4,ntot,dat);    // Stellar types
      fread(lsev,4,ntot,dat);     // Luminosities
      fread(rsev,4,ntot,dat);     // Radii
    }

    for (i=0;i<4;i++) c=fgetc(dat);   // Read 2nd record length of 2nd record
    tmyr = as[9];
  } while (fabs(tmyr-tsnap)>0.01);

  nstar = ntot-as[1];
  rbar  = as[2];
  zmbar = as[3];
  vstar = as[11];

  fprintf(stderr,"rbar  = %lf\n",rbar);
  fprintf(stderr,"zmbar = %lf\n",zmbar);

  dr[0] = as[6];
  dr[1] = as[7];
  dr[2] = as[8];
  dv[0] = as[26];
  dv[1] = as[27];
  dv[2] = as[28];

  rgal[0] = as[20]*rbar;
  rgal[1] = as[21]*rbar;
  rgal[2] = as[22]*rbar;

  vgal[0] = as[23]*vstar;
  vgal[1] = as[24]*vstar;
  vgal[2] = as[25]*vstar;

  for (i=0;i<nstar;i++) {
/*
      x[i] = (x[i]-dr[0])*rbar;
      y[i] = (y[i]-dr[1])*rbar;
      z[i] = (z[i]-dr[2])*rbar;
      u[i] = (u[i]-dv[0])*vstar;
      v[i] = (v[i]-dv[1])*vstar;
      w[i] = (w[i]-dv[2])*vstar;
*/
      xc[i] = xc[i]*rbar;
      yc[i] = yc[i]*rbar;
      zc[i] = zc[i]*rbar;
      uc[i] = uc[i]*vstar;
      vc[i] = vc[i]*vstar;
      wc[i] = wc[i]*vstar;
  }

  for (i=0;i<nstar;i++) {
      x[i] = rgal[0]-xc[i];
      y[i] = rgal[1]-yc[i];
      z[i] = rgal[2]-zc[i];
      u[i] = vgal[0]-uc[i];
      v[i] = vgal[1]-vc[i];
      w[i] = vgal[2]-wc[i];
  }

  for (i=0;i<nstar;i++) {
      xc[i] = xc[i]-dr[0]*rbar;
      yc[i] = yc[i]-dr[1]*rbar;
      zc[i] = zc[i]-dr[2]*rbar;
      uc[i] = uc[i]-dv[0]*vstar;
      vc[i] = vc[i]-dv[1]*vstar;
      wc[i] = wc[i]-dv[2]*vstar;
  }

  if (as[29]==0.0) {
    for (i=0;i<nstar;i++) {
       xx = 1E13*(mass[i])-floor(1E13*mass[i]+1.E-5);
       kstar[i]=floor(100*(xx+0.001));
    }
  }

  for (i=0;i<nstar;i++)        // Convert masses to Msun
    mass[i] *= zmbar;

 printf("Read %i stars at snapshot %lf\n",nstar,tmyr);

 fclose(dat);

 strcpy(fname,fnameu);
 if (!strcmp((fname+strlen(fname)-3),"POS")) fname[strlen(fname)-4]=0;

 strcat(fname,".nbform");

 if ((dat=fopen(fname,"w"))==NULL) {
   printf("Opening Failed !\n");
   exit(-1);
 }

 for (i=0;i<nstar;i++)
//   fprintf(dat,"%12i %12.3lf %12.3lf %12.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %3i\n",i,x[i],y[i],z[i],u[i],0.0,v[i],0.0,w[i],0.0,kstar[i]);
 fprintf(dat,"%8.5lf %12.3lf %12.3lf %12.3lf %10.3lf %10.3lf %10.3lf %12.4lf %12.4lf %12.4lf %10.4lf %10.4lf %10.4lf %3i\n",mass[i],x[i],y[i],z[i],u[i],v[i],w[i],xc[i],yc[i],zc[i],uc[i],vc[i],wc[i],kstar[i]);

 fclose(dat);
}


void invert_tform_mat(double **a) {
 double Determinant();
 void CoFactor(),Transpose();
 double det;
 double **t2;
 int i,j;

 t2 = malloc(3*sizeof(double *));
 for (i=0;i<3;i++)
    t2[i] = malloc(3*sizeof(double));

 det = Determinant(a,3);
 CoFactor(a,3,t2);
 Transpose(t2,3);

 for(i=0;i<3;i++)
    for(j=0;j<3;j++)
     a[i][j]=t2[i][j]/det;

 for(i=0;i<3;i++)
    free(t2[i]);
 free(t2);
}


/*
   Recursive definition of determinate using expansion by minors.
*/
double Determinant(double **a, int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}

/*
   Find the cofactor matrix of a square matrix
*/
void CoFactor(double **a,int n,double **b)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = malloc((n-1)*sizeof(double *));
   for (i=0;i<n-1;i++)
     c[i] = malloc((n-1)*sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}

/*
   Transpose of a square matrix, do it in place
*/
void Transpose(double **a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}


// Calculate transformation matrix from equatorial to galactic coordinates according to Johnson & Soderblom
void get_tform_mat() {

 double x11,x12,x13,x21,x22,x23,x31,x32,x33;

 x11 = -cos(PI/2.0+L0)*sin(PI/2.0-GND);
 x12 = -sin(PI/2.0+L0);
 x13 = cos(PI/2.0+L0)*cos(PI/2.0-GND);

 x21 = -sin(PI/2.0+L0)*sin(PI/2.0-GND);
 x22 = cos(PI/2.0+L0);
 x23 = sin(PI/2.0+L0)*cos(PI/2.0-GND);

 x31 = cos(PI/2.0-GND);
 x32 = 0.0;
 x33 = sin(PI/2.0-GND);

 t[0][0] = x11*cos(GNR-PI/2.0)+x12*sin(GNR-PI/2.0);
 t[0][1] = x11*sin(GNR-PI/2.0)-x12*cos(GNR-PI/2.0);
 t[0][2] = x13;

 t[1][0] = x21*cos(GNR-PI/2.0)+x22*sin(GNR-PI/2.0);
 t[1][1] = x21*sin(GNR-PI/2.0)-x22*cos(GNR-PI/2.0);
 t[1][2] = x23;

 t[2][0] = x31*cos(GNR-PI/2.0)+x32*sin(GNR-PI/2.0);
 t[2][1] = x31*sin(GNR-PI/2.0)-x32*cos(GNR-PI/2.0);
 t[2][2] = x33;
}



int findin(double ra, double dec, int nlinf, double *ragridf, double *decgridf) {
  double dirmat[16][2]={{0, 1}, {0.5, 1}, {1, 1}, {1, 0.5}, {1, 0}, {1, -0.5}, {1, -1}, {0.5, -1.0}, {0, -1}, {-0.5, -1.0}, {-1, -1}, {-1.0, -0.5}, {-1, 0}, {-1.0, 0.5}, {-1, 1}, {-0.5, 1.0}};
  int i,j,isect,intersect();

  for (j=0;j<16;j++) {
    isect=0;
    for (i=1;i<nlinf;i++) {
      if (intersect(ragridf[i-1],ragridf[i],decgridf[i-1],decgridf[i],ra,dec,ra+dirmat[j][0],dec+dirmat[j][1])==1) isect=1;
    }
    if (isect==0) return 0;
  }

  return 1;
}


int intersect(double x1, double x2, double y1, double y2, double x3, double y3, double x4, double y4) {
  double t,k,dx1,dy1,dx2,dy2;

  dx1 = x2-x1;
  dy1 = y2-y1;
  dx2 = x4-x3;
  dy2 = y4-y3;

  if (dx1*dy2-dy1*dx2==0.0) return 0;
  t = (x3*dy2-y3*dx2-x1*dy2+y1*dx2)/(dx1*dy2-dy1*dx2);
  if (dx2 != 0.0) k = (x1-x3+t*dx1)/dx2;
  else k = (y1-y3+t*dy1)/dy2;
  if (t<0 || t>1 || k<0) return 0;
  return 1;
}
