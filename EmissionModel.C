/* my first program in C++ */


#include <stdlib.h>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <fitsio.h>
#include <cmath>
#define pi 3.14159265358979323846L
using namespace std;
//Declaring some global variables
string outname, infile, outname2,rotcur,line,keystring;
char *commentstring;
double kpcscale, Distance;
double *rad,*rotvelor,*vv;
float ***Res;
float * Res_cpy;
float ** Diskdusttwodim;
char tmp2[16];
char * tmp3;
void Gauleg (double X1, double X2, double * X,double * W,int N);
void dimake (int NN,double td,double csttau,double Rdi,double y, double z,double sini,double cosi,double * xint,double * yint,double trRdi);
double Func3 (double t, double y,double z,double R0,double z0,double sini,double cosi,double td,double tau,double csttau,
               double * xint,double * yint,int NA,double * rad,double * rotvel,int nr,double dveldisp,double * vadd,
               double vlow,double vscale,int nv,double * boun,double bdrat,double b0,double bratio,double buldisp,
	      int p, int k,int l,double trRdi,double *bounem);
double dustint (double t,int NN,double * xint,double * yint,double * boun,double y,double z,double cosi,double sini);
void tplint(double * XA,double * YA,int N, double X,double& Y);
int main ()
{
  //Declaring Local variables for main
  //this sets the values for these two variables
  double inc,incout,R0,z0,Rdi,vlow,y,z, result;
  double lag,dveldisp,dveldispor,td,vscale,xs,zs,csttau;
  double dum,Xlo,Xhi,Zlo,Zhi,frac;
  double tau,rdust,cenvel,cosi,sini,xxx,Val,Val2;
  double dveldispchange,dveldispouter,trRdi,bounem[2];
  double boun[4],bdrat,b0,bratio,buldisp,f2,Sum3;
  double halfax1,halfax2,halfax3,cdeltscale1,cdeltscale2,xcrval,ycrval;
  double *xint,*yint,*vadd,*vprof;  
  double *XX,*W,*rotvel, *vprn, *vpro;
  int nx,nz,nv,N,nzn,nxn,p,f,k;
  int counter,intnum,nelements;
  int i,NN,NA,ii,j;
  int bitpix,naxis,nr,nvp2;
  int status[1],exists[1];
  long naxes[3],naxes2[2];
  long firstpix[3] = {1,1,1};
  bool simple,extend;
  bool Log,Fold,Unproj;
  fitsfile *fptr;
  
  cout << "What's is the input file?\n ";
  getline (cin, infile);
  // infile="Exampleparameter5.txt";
  ifstream inputfile;
  inputfile.open (infile.c_str(),ios::in);
  if (inputfile.is_open ()) {
    getline (inputfile,line);
    getline (inputfile,outname);
    getline (inputfile,line);
    getline (inputfile,outname2);
    getline (inputfile,line); 	    
    inputfile >> inc; 
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> R0;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> z0;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> Rdi;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> td;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> tau;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> rdust;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> bdrat;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> b0;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> bratio;
    getline (inputfile,line); getline (inputfile,line); 	
    getline (inputfile,line); 
    if (line == "N" || line == "n") { Log = false; } else { Log = true;}
    getline (inputfile,line); 	
    getline (inputfile,line); 
    if (line == "N" || line == "n") { Fold = false; } else { Fold = true;}
    getline (inputfile,line); 	
    inputfile >>  rotcur;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> lag;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> dveldisp;
    dveldispor =dveldisp;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> dveldispouter;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> trRdi;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> buldisp;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> vlow;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> vscale;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> cenvel;
    getline (inputfile,line); getline (inputfile,line); 	
    getline (inputfile,line);
    if (line == "N" || line == "n") { Unproj = false; } else { Unproj = true;}
    getline (inputfile,line); 	
    inputfile >> Xlo;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> Xhi;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> Zlo;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> Zhi;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> nx;	
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> nz;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> nv;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> N;
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> bounem[0];
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> bounem[1];
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> boun[0];
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> boun[1];
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> boun[2];
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> boun[3];
    getline (inputfile,line); getline (inputfile,line); 	
    inputfile >> Distance;
    inputfile.close();
  } else cout << "Unable to open file";

  Distance *=1E3;
  kpcscale = Distance*((1./(360.*60.))*2*pi);
  cout << "You have given a Distance of " << Distance <<  " kpc \n";
  if (Log) {
    if (Xlo <= 0.) { Xlo = 1.;}
    if (Zlo <= 0.) { Zlo = 1.;} 
    Xlo = log10(Xlo);
    Xhi = log10(Xhi);
    Zlo = log10(Zlo);
    Zhi = log10(Zhi);
    }
  nr = 0;
  ifstream rotationfile;
  rotationfile.open (rotcur.c_str(),ios::in);
  while ( !line.empty()){getline(rotationfile,line);nr++;}
  rotationfile.close();
  nr=nr-1;
  rad = new double[nr];
  rotvelor = new double[nr];
  counter=0;
  rotationfile.open (rotcur.c_str(),ios::in);
  if (rotationfile.is_open ()) {
    while ( rotationfile.good() &&  counter < nr){
	rotationfile >> rad[counter]; 
	rotationfile >> rotvelor[counter];
	counter++;
	}
  } else {cout << "The rotation file is missing \n"; return 0;}
  nvp2=nv+1;

  cout << "You have provided the following input\n";
  cout << "The dimensions of x=" << nx << "z=" << nz << " v=" << nv <<endl;
  cout << "Physically ranging in X from " <<  Xlo << " to " << Xhi <<  "kpc in X\n";
  cout << "Physically ranging in Z from " <<  Zlo << " to " << Zhi <<  "kpc in Z\n";
  cout << "Physically ranging in vel from " <<  vlow << " to " << vscale*nv+vlow  <<  " km/s in velocity\n";
  cout << "A Hole from " << bounem[0] << " to " << bounem[1] << " kpc\n";
  cout << "The emission disk has a scale length " << R0  << " and height " << z0  << endl;
  cout << "The emission disk has a inclination of " << inc << endl;
  cout << "The disk is truncated at " << trRdi  << " kpc\n";
  cout << "You have given a Distance of " <<  Distance <<  " kpc\n";
  cout << "The dust disk has a scale length " << Rdi  << " and height " <<  td << endl; 
  cout << "With an optical depth of " <<  tau << ";and a spatial extend of " <<  rdust << endl;
  cout << "The bulge to disk ratio is " << bdrat << " with an ellipticity of " <<  bratio << endl;
  cout << "With a scale length of " <<  b0 << endl;
  cout << "You have logarithmic binning? " <<  Log << endl;
  cout << "You fold the disk? "  << Fold << endl;
  cout << "You have a rotation curve \n";
  for (int i=0; i < nr; i++) {
    cout << rad[i] << " " << rotvelor[i] << endl;
  }
  cout << "The disk has a central velocity of " << cenvel << endl;
  cout << "With central dispersion " << dveldisp <<  " km/s declining to " << dveldispouter << endl;
  cout << "And a bulge dispersion of " <<  buldisp << endl;
  cout << "And is lagging with " <<  lag << " km/s/kpc " <<endl;
  cout << "The Boundaries are " << boun[0] << " "<< boun[1] << " " << boun[2] << " " << boun[3] << endl;
  cout << "The Number of steps in Gauss-legendre integration "  <<  N <<endl;
  cout << "The Cube will be written to the fits file" << outname << endl;
  cout << "The Disk will be written to the fits file " << outname2 << endl;
  incout = inc;
  inc = abs(inc);
  if (inc >90.){ inc=90.;}
  inc=inc*pi/180L;
  cosi=cos(inc);
  //  cout << cosi << endl;
  sini=sin(inc);
  // cout << sini << endl;
  if (Fold) { nxn=(nx)/2 ; nzn=(nz)/2;} else {
    nxn=nx+1;
    nzn=nz+1;
  }
  xs=(Xhi-Xlo)/(nxn-1L);
  zs=(Zhi-Zlo)/(nzn-1L);
  NN=128;
  csttau=1L;
  y=0L;
  z=0L;
  vv = new double[nv];
  for (int i=0 ; i < nv; i++) {
    vv[i]=vlow+(i-1)*vscale;
  }
  xint = new double[NN]; yint = new double[NN];
  dimake(NN,td,csttau,Rdi,y,z,sini,cosi,xint,yint,trRdi);
  NA=NN;
  csttau=tau/yint[1];
  delete [] yint;
  delete [] xint;
  if (Fold) {
    dveldispchange=(dveldisp-dveldispouter)/nxn;
    dveldisp=dveldisp+dveldispchange;
  } else{
    dveldispchange=((dveldisp-dveldispouter)/((nxn-1)/2L));
    dveldisp=dveldispouter-dveldispchange;
  }
  //float Res[nv][nz][nx];
  float ***Res = new float **[nv];
  for ( i = 0; i < nv; i++){
    Res[i] = new float *[nz];
    for (j = 0; j < nz; j++){
      Res[i][j] = new float[nx];
    }
  }
 
  //Res_cpy = new float[(nx * nv * nz)];

  float **Diskdusttwodim = new float *[nz];
  for ( i = 0; i < nz; i++){
    Diskdusttwodim[i] = new float[nx];
  }
  Diskdusttwodim[4][4]=0.;
  
  for (i=0; i< nxn-1; i++){
    if (i%5 == 0) {
      cout << "We are at x pixel " << i << endl;
    }
    if (Fold) {
      dveldisp=dveldisp-dveldispchange;
    }else{
      if (i < nxn/2.+1) {	     
	dveldisp=dveldisp+dveldispchange;
      }	else {
	dveldisp=dveldisp-dveldispchange;
      }
    }
    for (p=0; p < nzn-1; p++){ 
      y=Xlo+(i)*xs;
      z=Zlo+(p)*zs;
      rotvel = new double[nr];
      for (f=0 ; f < nr; f++){
	if (rotvelor[f] < 0) {
	  rotvel[f]=rotvelor[f]+(lag*(abs(z)));
	} else {
	  rotvel[f]=rotvelor[f]-(lag*(abs(z)));
	}
	if (rotvelor[f] == 0) {
	  rotvel[f]=0.;
	}
      }
      if (N == 0) {
	frac=(Zhi-z)/(Zhi-Zlo);
	NN=64*frac;
	if (NN < 16) {NN=16;}
      }else {NN=N;}
      
      xint = new double[NN]; yint = new double[NN];
      dimake(NN,td,csttau,Rdi,y,z,sini,cosi,xint,yint,trRdi);
      NA=NN;
      dum=pi/2L;
    
      XX=new double[NN]; W=new double[NN];
      Gauleg(-dum,dum,XX,W,NN);
     
      vprof = new double[nv]; 
      for (ii=0; ii< nv; ii++){ vprof[ii]=0;}
      Sum3=0L;
      for (k=0; k < NN; k++){
	f2=0.;
	vadd =new double[nv];
	f2=Func3 (XX[k],y,z,R0,z0,sini,cosi,td,tau,csttau,xint, yint,NN,rad,rotvel,nr,dveldisp,vadd, vlow,  vscale, nv,
		  boun,  bdrat,  b0,  bratio,  buldisp, p, k,i,  trRdi,  bounem);
	Sum3=Sum3+W[k]*f2;
        for (ii=0; ii< nv; ii++){
	  vprof[ii]=vprof[ii]+W[k]*vadd[ii];
	}
       	delete [] vadd;

      }
      if (Sum3 < 0.) { Sum3=0. ; }
      delete [] XX;
      delete [] W;
      delete [] xint;
      delete [] yint;
      delete [] rotvel;
      Val=0;
      Val2=Sum3;
      
      if (Fold){ 
	Diskdusttwodim[nz-nzn+p][nx-nxn+i]= Sum3;
	Diskdusttwodim[nzn-p+2][nx-nxn+i]= Sum3;
	Diskdusttwodim[nz-nzn+p][nxn-i+2]= Sum3;
	Diskdusttwodim[nzn-p+2][nxn-i+2]= Sum3;
      }else {
	if (Sum3 <= 0.) {Diskdusttwodim[p][i]=0.;} else {Diskdusttwodim[p][i]=(float)Sum3;}}
      vpro = new double[nv]; 
      for (ii=0; ii < nv; ii++){
	if (vprof[ii] < 0){vprof[ii] = 0.;} 
	vpro[ii]=vprof[ii];
      }
      if (Fold) {
	for (j=0;j<nv;j++){
	  // Res_cpy[((j*nz*nx) + ((nz-nzn+p)*nx) + (nx-nxn+i))];
	  // Res_cpy[((j*nz*nx) + ((nzn-p+2)*nx) + (nx-nxn+i))];
	  Res[j][nz-nzn+p][nx-nxn+i]=vprof[j];
	  Res[j][nzn-p+2][nx-nxn+i]=vprof[j];
        }
	vprn = new double[nv];
	for (ii=0; ii< nv; ii++){
	  xxx=0-vv[ii];
	  tplint(vv,vpro,nv,xxx,result);
	  vprn[ii]=result;
	}
	
	Val2=0;
	for (ii=0; ii<nv; ii++){ Val2=Val2+vprn[ii];}
	for (j=0; j<nv; j++){ 
	  //Res_cpy[((j*nz*nx) + ((nz-nzn+p)*nx) + (nx-i+2))];
	  //Res_cpy[((j*nz*nx) + ((nzn-p+2)*nx) + (nxn-i+2))]
	  Res[j][nz-nzn+p][nxn-i+2]=vprn[j];
	  Res[j][nzn-p+2][nxn-i+2]=vprn[j];
	}
	delete [] vprn;
      } else {	
	for (j=0; j< nv; j++){
	  if (vprof[j] <= 0.) {
	    Res[j][p][i]=0.;
	  } else {
	    Res[j][p][i]=(float)vprof[j]; 
	  }
	}
      }
      delete [] vpro;
      delete [] vprof;
    }
  } 

  delete [] vv;
  delete [] rad;	
  delete [] rotvelor;
  fits_file_exists (outname.c_str(), exists, status);
  cout << "Does " << outname  << " exists y=1 n=0 |" << exists[0] << endl;
  if (exists[0] == 1) {	  
    ffopen(&fptr,outname.c_str(),1,status);
    ffdelt(fptr, status);
    fits_file_exists (outname.c_str(), exists, status);
    cout <<  "Do we remove it? y=0 n=1 |"  << exists[0] << endl;
  }
  //	create a new empty fitsfile	
  ffinit(&fptr,outname.c_str(),status);
  cout << "Writing the file" << outname << " " << status[0] << endl;
  //   initialize parameters about the FITS image 
  simple=1;
  naxes[0]=nx;
  naxes[1]=nz;
  naxes[2]=nv;
  halfax1=(nx/2.)+1;
  halfax2=(nz/2.)+1;
  halfax3=(abs(vlow)/vscale)+1;
  cdeltscale1=(xs/kpcscale)*60;
  xcrval=0;
  cdeltscale2=(zs/kpcscale)*60;
  ycrval=0;
  intnum=abs(N);
  naxis=3;
  bitpix=-32;
  extend=1;
  keystring="CRPIX1";
  commentstring= new char[1];
  ffphpr(fptr,simple,bitpix,naxis,naxes,0,1,extend,status);
  ffdkey(fptr,"COMMENT",status);
  ffdkey(fptr,"COMMENT",status);
  ffuky(fptr,TDOUBLE,keystring.c_str(),&halfax1,commentstring,status);
  keystring="CRVAL1";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&xcrval,commentstring,status);
  keystring="CDELT1";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&cdeltscale1,commentstring,status);
  keystring="CTYPE1";
  strcpy(tmp2,"ARCSEC1");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);
  keystring="CUNIT1";
  strcpy(tmp2,"ARCSEC");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);
  keystring="CRPIX2";
 ffuky(fptr,TDOUBLE,keystring.c_str(),&halfax2,commentstring,status);
  keystring="CRVAL2";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&ycrval,commentstring,status);
  keystring="CDELT2";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&cdeltscale2,commentstring,status);
  keystring="CTYPE2";
  strcpy(tmp2,"ARCSEC2");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);
  keystring="CUNIT1";
  strcpy(tmp2,"ARCSEC");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);

  keystring="CRPIX3";
 ffuky(fptr,TDOUBLE,keystring.c_str(),&halfax3,commentstring,status);
  keystring="CRVAL3";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&cenvel,commentstring,status);
  keystring="CDELT3";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&vscale,commentstring,status);
  keystring="CTYPE3";
  strcpy(tmp2,"VELOCITY");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);
  keystring="CUNIT3";
  strcpy(tmp2,"KM/S");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);
  cout << "where" <<endl;
  tmp3 = new char[120];
  sprintf( tmp3 , "The Optical Depth = %5.5f",tau );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Inclination = %5.5f degrees",incout );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Number of Integrations =  %d",N );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Disk Scale Length = %5.5f kpc",R0 );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Disk Scale Height = %5.5f kpc",z0 );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Dust Scale Length = %5.5f kpc",Rdi );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Dust Scale Height = %5.5f kpc",td );
  ffpcom(fptr,tmp3,status);
   sprintf( tmp3 , "The Spatial Extent of the Dust = %5.5f kpc",rdust);
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Bulge to Disk Ratio = %5.5f",bdrat );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "B0 = %5.5f",b0 );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Ellipticity of the Bulge = %5.5f",bratio );
  ffpcom(fptr,tmp3,status);
  if (Log) {sprintf( tmp3 , "You Requested Logarithmic Binning" );} else
    {sprintf( tmp3 , "You Did Not Request Logarithmic Binning" );}
  ffpcom(fptr,tmp3,status);
  if (Fold) {sprintf( tmp3 , "You Requested Folding for Symmetry" );} else
    {sprintf( tmp3 , "You Did Not Request Folding for Symmetry" );}
  ffpcom(fptr,tmp3,status);
  if (Unproj) {sprintf( tmp3 , "You Requested a Test Parameter" );} else
    {sprintf( tmp3 , "You Did Not Request a Test Parameter" );}
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Rotation Curve Was Read From %s",rotcur.c_str());
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Lag Is = %5.5f km/s/kpc",lag );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Isotropic Central Velocity Dispersion Is = %5.5f km/s",dveldispor );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Isotropic Outer Velocity Dispersion Is = %5.5f km/s",dveldispouter );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Truncation Radius = %5.5f kpc",trRdi );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Velocity Dispersion of the Bulge = %5.5f km/s",buldisp );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Central Velocity of the Galaxy = %5.5f km/s",cenvel );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Upper Limit on the Xaxis = %5.5f kpc",Xhi );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Lower Limit on the Xaxis = %5.5f kpc",Xlo );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Upper Limit on the Yaxis = %5.5f kpc",Zhi );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Lower Limit on the Yaxis = %5.5f kpc",Zlo );
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "A Hole Starts At %5.5f And Ends At  %5.5f kpc",bounem[0],bounem[1]);
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Dust Runs From %5.5f To  %5.5f kpc",boun[0],boun[1]);
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Dust Runs From %5.5f To  %5.5f kpc",boun[2],boun[3]);
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Distance = %5.5f kpc",Distance);
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "The Disk File Belonging To This Cube = %s",outname2.c_str());
  ffpcom(fptr,tmp3,status);
  nelements=naxes[0]*naxes[1]*naxes[2];
  Res_cpy = new float[(nx * nv * nz)];
  for(k = 0; k < nv; k++){
    for(j = 0; j < nz; j++){
      for(i = 0; i < nx; i++){
  	Res_cpy[((k*nz*nx) + (j*nx) + i)] = Res[k][j][i];
      }
    }
  }
  
  fits_write_pix(fptr, TFLOAT, firstpix, nelements, Res_cpy, status);
  delete [] Res_cpy;
  for ( i = 0; i < nv; i++){
    for (int j = 0; j < nz; j++){
      delete [] Res[i][j];
    }
    delete [] Res[i];
  }
  delete [] Res;
  fits_close_file(fptr, status);
  status[0]=0;
  fits_file_exists (outname2.c_str(), exists, status);
  cout << "Does " << outname2  << " exists y=1 n=0 " << exists[0] << endl;
  if (exists[0] == 1) {	  
    ffopen(&fptr,outname2.c_str(),1,status);
    ffdelt(fptr, status);
    fits_file_exists (outname2.c_str(), exists, status);
    cout <<  "Do we remove it? y=0 n=1// "  << exists[0] << endl;
  }
  //	create a new empty fitsfile	
  ffinit(&fptr,outname2.c_str(),status);
  cout << "Writing the file" << outname2 << " " << status[0] << endl;
  //     initialize parameters about the FITS image 
  naxis=2;
  naxes2[0]=nx;
  naxes2[1]=nz;
  ffphpr(fptr,simple,bitpix,naxis,naxes2,0,1,extend,status);
  ffdkey(fptr,"COMMENT",status);
  ffdkey(fptr,"COMMENT",status);
  keystring="CRPIX1";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&halfax1,commentstring,status);
  keystring="CRVAL1";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&xcrval,commentstring,status);
  keystring="CDELT1";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&cdeltscale1,commentstring,status);
  keystring="CTYPE1";
  strcpy(tmp2,"ARCSEC1");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);
  keystring="CUNIT1";
  strcpy(tmp2,"ARCSEC");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);
  keystring="CRPIX2";
 ffuky(fptr,TDOUBLE,keystring.c_str(),&halfax2,commentstring,status);
  keystring="CRVAL2";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&ycrval,commentstring,status);
  keystring="CDELT2";
  ffuky(fptr,TDOUBLE,keystring.c_str(),&cdeltscale2,commentstring,status);
  keystring="CTYPE2";
  strcpy(tmp2,"ARCSEC2");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);
  keystring="CUNIT1";
  strcpy(tmp2,"ARCSEC");
  ffuky(fptr,TSTRING,keystring.c_str(),&tmp2,commentstring,status);
  delete [] commentstring;
  sprintf( tmp3 , "The Cube File Belonging To This Image = %s",outname.c_str());
  ffpcom(fptr,tmp3,status);
  sprintf( tmp3 , "Look there for all the input parameters ");
  ffpcom(fptr,tmp3,status);
  delete [] tmp3;
  nelements=naxes2[0]*naxes2[1];
  fits_write_img(fptr, TFLOAT, 1, nelements, Diskdusttwodim[0], status);
  fits_close_file(fptr, status);
  
  for ( i = 0; i < nz; i++){
    delete [] Diskdusttwodim[i];
  }
  delete []  Diskdusttwodim;

  return 0;
}


double Func3 (double t, double y,double z,double R0,double z0,double sini,double cosi,double td,double tau,double csttau,
               double * xint,double * yint,int NA,double * rad,double * rotvel,int nr,double dveldisp,double * vadd,
               double vlow,double vscale,int nv,double * boun,double bdrat,double b0,double bratio,double buldisp,
	      int p, int k,int l,double trRdi,double *bounem){
  double f2b,f2d,tt,fd,fb,vel;
  double bz0,radialpos,z1,y1,noemer,result,test;
  double *hulp;
  double d1,dd2,dd3,db2,db3,arg1,arg2,arg3,res1;
  int i;
  bz0=b0*bratio;
  radialpos=sqrt(y*y+(z/cosi+tan(t)*sini)*(z/cosi+tan(t)*sini));	
  if (t <= -pi/2L+1e-4L ||  t >= pi/2L-1e-4L || radialpos > trRdi || (radialpos >= bounem[1] && radialpos <= bounem[2])) {
    
    f2b=0L;
    f2d=0L;
  } else {
    tt=tan(t);
    if (abs(tau) < 1.e-4L){
      z1=z;
      fd=exp(-(sqrt(y*y+(z/cosi+tt*sini) *(z/cosi+tt*sini)))/R0-abs(tt*cosi/z0));
      fb=exp(-(sqrt(y*y+(z/cosi+tt*sini) *(z/cosi+tt*sini)))/b0-abs(tt*cosi/bz0));
      f2d=fd*(1L+tt*tt);
      f2b=bdrat*fb*(1L+tt*tt);
    } else {
      d1=dustint(t,NA,xint,yint,boun,y,z,cosi,sini);
      z1=z;
      dd2=abs(tt*cosi/z0);
      dd3=sqrt(y*y+(z/cosi+tt*sini)*(z/cosi+tt*sini))/R0;
      db2=abs(tt*cosi/bz0);
      db3=sqrt(y*y+(z/cosi+tt*sini)*(z/cosi+tt*sini))/b0;
      arg1 =  d1;
      arg2 = d1 + dd2 + dd3;
      arg3 = db2 + db3;
      f2d=exp(-arg2)*(1+tt*tt);
      f2b=exp(-arg1)*(bdrat*exp(-arg3)) * (1+tt*tt);
    }
  }
  y1=y;
  tplint(rad,rotvel,nr,y1,res1);
  noemer=y*y+pow((tt*sini+z/cosi),2);
  if (noemer == 0.) { result=0.;} else { result=res1*sini*abs(y)/sqrt(noemer);}
  test=0.;
  if (f2d < 1e-50L) { f2d=0.;}
  for (i=0; i<nv; i++){
    vel=vlow+(i)*vscale;
    vadd[i]=exp(-0.5L*(pow(((vel-result)/dveldisp),2.)));
  }
  test=0;
  for (i=0; i<nv; i++){test=test+vadd[i];}
  for (i=0; i<nv; i++){ vadd[i]=vadd[i]*f2d/test; } 	
  hulp = new double[nv];
  for (i=0; i<nv; i++){vel=vlow+(i)*vscale;  
    hulp[i]=exp(-0.5*(pow(((vel-result)/buldisp),2.))); }
  test=0;
  for (i=0L; i<nv; i++){test=test+hulp[i]; }
  for (i=0L; i<nv; i++){ vadd[i]=vadd[i]+((hulp[i]*f2b)/test);}
  delete [] hulp;
  return f2d;
}

void dimake (int NN,double td,double csttau,double Rdi,double y, double z,double sini,double cosi,double * xint,double * yint,double trRdi){
  int i,k;
  double tt,Sum2,xst;
  double *LL,*MM;
  double f2,xend;
 
 
  for (i=0; i < NN; i++){
    Sum2 = 0.;
    xst = -pi/2.+(i)*pi/NN;
    xend = xst+pi/NN;
    LL = new double[16];
    MM = new double[16];
    Gauleg(xst,xend,LL,MM,16);

 
    for ( k=0; k< 16L; k++){
      if ( LL[k] <= -pi/2L+1e-4L ||
	   LL[k] >= pi/2L-1e-4L  ||
	   sqrt(y*y+(z/cosi+tan(LL[k])*sini)*(z/cosi+tan(LL[k])*sini)) > trRdi) {Sum2=Sum2;
      } else {
	tt=tan(LL[k]);
	f2=exp(-(sqrt(y*y+(z/cosi+tt*sini) *
		      (z/cosi+tt*sini)))/Rdi -
	       abs(tt*cosi/td)) * 
	  (1L+tt*tt);
	Sum2=Sum2+MM[k]*f2;
      }
    }
    delete [] LL;
    delete [] MM;
    xint[i]=xst;
    yint[i]=Sum2*csttau;
  }
  for (i=NN-2; i >= 0;i--){
    yint[i]=yint[i]+yint[i+1];
    xint[i]=xint[i];
  }
  for (i=0; i < NN-1; i++){
    xint[i]=(xint[i]+xint[i+1])/2L;
  }

}
void Gauleg (double X1, double X2, double * X,double * W,int N) {
  int I,J,M;
  double EPS = 3E-14L,P3,P2,P1,PP;
  double XM,XL,Z,Z1;
  M=(N+1)/2;
  XM=0.5L*(X2+X1);
  XL=0.5L*(X2-X1);
  
  for ( I=0 ; I<M ; I++){
    Z=cos(pi*(I+1L-.25L)/(N+.5L));
    while (abs(Z-Z1) > EPS){
      P1=1L;
      P2=0L;
      for (J=1; J<=N; J++){
	P3=P2;
	P2=P1;
	P1=((2L*J-1L)*Z*P2-(J-1L)*P3)/J;
      }
      PP=N*(Z*P1-P2)/(Z*Z-1L);
     
      Z1=Z;
      Z=Z1-P1/PP;
    }
    X[I]=XM-XL*Z;
    X[N-I-1]=XM+XL*Z;
    W[I]=2L*XL/((1L-Z*Z)*PP*PP);
    W[N-I-1]=W[I];
  }

}



double dustint(double t,int NN,double * xint,double * yint,double * boun,double y,double z,double cosi,double sini){
  double dustint,rad;
  rad=sqrt(y*y+(z/cosi+tan(t)*sini)*(z/cosi+tan(t)*sini));
  if ((rad >= boun[0] && rad <= boun[1]) || (rad >= boun[2] && rad <= boun[3])) {
    tplint(xint,yint,NN,t,dustint);
    if (dustint < 0.) {dustint=0.;}
  } else { dustint=0L;}
  return dustint;
}


void tplint(double * XA,double * YA,int N, double X,double& Y){
  int K,KLO,KHI;
  KLO=0L;
  KHI=N-1L;
  while (KHI-KLO > 1L) {
    K=(KHI+KLO)/2L;
    if (XA[K] > X){KHI=K;} else {KLO=K;}
  }
  Y=YA[KLO]+(X-XA[KLO])*(YA[KHI]-YA[KLO])/(XA[KHI]-XA[KLO]);
}


