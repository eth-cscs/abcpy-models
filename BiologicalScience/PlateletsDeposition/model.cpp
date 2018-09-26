#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sstream>

#include <boost/random.hpp>

using namespace std;

class Random
{
  public:
  boost::mt19937 eng; 
  boost::random::uniform_real_distribution<double> random;
  boost::random::uniform_int_distribution<int> randint;

  Random()
  {
    eng.seed(time(NULL)); 
    random=boost::random::uniform_real_distribution<double>(0,1);
  }

  Random(int a, int b)
  {
    eng.seed(time(NULL)); 
    randint=boost::random::uniform_int_distribution<int>(a,b) ;
  }

  void init()
  {
   eng.seed(time(NULL)); 
  }

  void init(int seed)
  {
   eng.seed(seed); 
  }

  double get()
  {
    return random(eng);
  }

  int Get()
  {
    return randint(eng);
  }

};

// Random rng;  BC_


double Prob(double p, double maxAlbumine, double albumine)
{
  if(albumine>=maxAlbumine)return 0.;
  return p*(1 - albumine/maxAlbumine);
}


double Prob2(double p, double attenuation, double maxAlbumine, double albumine)
{
  if(albumine>=maxAlbumine)return 0.;
  return p*exp(-attenuation*(albumine/maxAlbumine));
}


class d1q3
{
public:
  int sizeX_;
  double v,dt,dx;
  double w,w0,cs2,tau;
  double* f0; double* f1; double* f2;
  double* feq0; double* feq1; double* feq2;
  double* rho; double* u;
  
  d1q3(int nx, double density, double Diffusion,double Dt, double Dx, double W)  // ----constructor
  {
    sizeX_=nx;  
    dx=Dx;//dx=2.e-6;
    dt=Dt;//dt=(4./3)*1e-4;
    v=dx/dt;

    w=W; w0=1.-2*w;
    cs2=(1-w0)*v*v;
    tau=Diffusion/(cs2*dt) + 0.5;
    //    cout<<"Diffusion model:"<<endl;
    //    cout<<"dx="<<dx<<" dt="<<dt<<" tau="<<tau<<" v="<<v<<" cs2="<<cs2<<endl;
    // cin.get();

    f1=new double[sizeX_]; f2=new double[sizeX_];  f0=new double[sizeX_]; 
    feq1=new double[sizeX_]; feq2=new double[sizeX_];  feq0=new double[sizeX_]; 
    rho=new double[sizeX_]; u=new double[sizeX_]; 

    for(int i=0; i<sizeX_; i++){   // initialize the fi at equilibrium
	f1[i]=w*density;   f2[i]=w*density;  	f0[i]=w0*density;
	rho[i]=f0[i]+f1[i]+f2[i];
	u[i]=v*(f1[i]-f2[i]);
    }
  } //----------------- end constructor ---------------


  ~d1q3()   //---------------destructeur
  {
      delete [] f1;      delete [] f2;  delete [] f0;
      delete [] feq1;      delete [] feq2;  delete [] feq0;
      delete [] rho;      delete [] u;
  } //------------------------  end destructor

  void init(double density)
  { 
    for(int i=0; i<sizeX_; i++){   // initialize the fi at equilibrium
	f1[i]=w*density;   f2[i]=w*density;  	f0[i]=w0*density;
	rho[i]=f0[i]+f1[i]+f2[i];
	u[i]=v*(f1[i]-f2[i]);
    }
  }


  double getTau(){return tau;}

  double getRho(int i)
  { return f0[i]+f1[i]+f2[i];}

  double getJ(int i)
  { return v*(1-1/(2*tau))*(f1[i]-f2[i]);}


  void computeRhoU(int i)
  {
      rho[i]=f0[i]+f1[i]+f2[i];
      u[i]=v*(f1[i]-f2[i]);
  }

  void computeRhoU()
  {
    for(int i=0; i<sizeX_; i++){  
      rho[i]=f0[i]+f1[i]+f2[i];
      u[i]=v*(f1[i]-f2[i]);
    }
  }

  double computeNoPart()
  {
    double n=0;
    for(int i=0; i<sizeX_; i++){  
      n+=f0[i]+f1[i]+f2[i];
    }
    return  n*dx;
  }

  double computeNoPart(int i0, int i1)
  {
    double N=0;

    for(int i=i0; i<i1; i++){  
      N+=f0[i]+f1[i]+f2[i];
    }
    return N*dx; 
  }



  void printRhoU()
  {
    for(int i=0; i<sizeX_; i++){  
      rho[i]=f0[i]+f1[i]+f2[i];
      u[i]=v*(f1[i]-f2[i])/rho[i];
      //cout <<i<<"   "<<rho[i]<<"   "<<u[i]<<endl;
    }
  }

  void printQ(int j, int k)
  {
    for(int i=j; i<=k; i++){  
      //cout <<v*(f1[i]-f2[i])<<" ";
    }
    //cout <<endl;
  }

  void printQLink(int j, int k)
  {
    for(int i=j; i<=k; i++){  
      //cout <<v*(f1[i+1]-f2[i])<<" ";
    }
    //cout <<endl;
  }


  void printRho()
  {
    for(int i=0; i<sizeX_; i++){  
      rho[i]=f0[i]+f1[i]+f2[i];
      //cout <<i<<"   "<<rho[i]<<endl;
    }
  }


  void printRho(int start, int end)
  {
    if(start<0)start=0;    if(end>=sizeX_)end=sizeX_;
    for(int i=start; i<=end; i++){  
      rho[i]=f0[i]+f1[i]+f2[i];
      //cout <<rho[i]<<"  ";
    }
    //cout <<endl;
  }


  void initLattice()  //-------------- initialize the lattice
  {
    double delta;
    for(int i=0; i<sizeX_; i++){
      //      delta=.006*exp( -((i-sizeX_/2.)*(i-sizeX_/2.))/5.); 
      delta=0;
      f1[i]=.01+delta;
      f2[i]=.01+delta;
      f0[i]=.01+delta;
    }
  }  //-------------------------- end initLattice
  
  
  void collision()  //-------------- collision--------------
  {
    for(int i=0; i<sizeX_; i++){  
      rho[i]=f0[i]+f1[i]+f2[i];

      feq0[i]=w0*rho[i];
      feq1[i]=w*rho[i];
      feq2[i]=w*rho[i];

      f0[i]=f0[i]+ (1/tau)*(feq0[i]-f0[i]);
      f1[i]=f1[i]+ (1/tau)*(feq1[i]-f1[i]);
      f2[i]=f2[i]+ (1/tau)*(feq2[i]-f2[i]);
    }
}  //----------------- end collision-----------------------

  void propagation(double& left, double& right) //--------------------------
  {

    //---------   f1    ---------------------------------
    // right cshift 
      right=f1[sizeX_-1];
      for(int i=sizeX_-1; i>0; i--)f1[i]=f1[i-1];
      f1[0]=right;  // this creates periodic boundary condition
    //---------   f2    ---------------------------------
    // leftt cshift 
      left=f2[0];
      for(int i=0; i<sizeX_-1; i++)f2[i]=f2[i+1];
      f2[sizeX_-1]=left;// this creates periodic boundary condition
  }  // end propagation-----------------------------------


  void leftRightBoundary(double f1_left, double f2_right)
  {
    f1[0]=f1_left;  f2[sizeX_-1]=f2_right;
  }

  void boundary(double rho, double f2_right)
  {
    f1[0]=rho-f0[0]-f2[0];  f2[sizeX_-1]=f2_right;
  }


  void leftBoundaryFixedDensity(double rho)
  {
    f1[0]=rho-f0[0]-f2[0];
  }
  void rightBoundaryFixedDensity(double rho)
  {
    f2[sizeX_-1]=rho-f0[sizeX_-1]-f1[sizeX_-1];
  }

  void rightBoundaryFixedQ(double Q)
  {
    f2[sizeX_-1]=f1[sizeX_-1]+Q/v;
  }



}; //================  End class d1q3  =========


class Particles
{
  public:
  boost::mt19937 eng2;
  boost::random::uniform_real_distribution<double> rng;// BC

  int nx, ny;
  int* sizeOfCluster;
  int* volumeOfCluster;
  int*  maxHeightOfCluster;
  int** deposit;
  int** label;
  double** avail;  // array to hold if site is available for adhesion

  int noCluster;
  int noAggreg;
  int noColl;
  int noPart;   // total no of particle created in space 
  int aggregateIndex;     // to enumerate the clusters
  double averageSize;
  double averageVol;
  double averageHeight;


  Particles(int nnx, int nny,int seed)  // --------- constructor BC_
  {
    rng=boost::random::uniform_real_distribution<double>(0,1); // BC
    eng2.seed(seed);  //BC_

    nx=nnx; ny=nny;
    sizeOfCluster=new int[nx*ny];  //  size of each the clusters
    volumeOfCluster=new int[nx*ny];  // volume of each the clusters
    maxHeightOfCluster=new int[nx*ny];  // max height of each clusters
    // allocation memory for the arrays
    deposit=new int*[nx];
    for(int i=0;i<nx;i++)deposit[i]=new int[ny];

    avail=new double*[nx];
    for(int i=0;i<nx;i++)avail[i]=new double[ny];

    label=new int*[nx];
    for(int i=0;i<nx;i++)label[i]=new int[ny];


    // initialized variables to zeros: 
    for(int i=0;i<nx;i++)
      for(int j=0;j<ny;j++){
	deposit[i][j]=0;
	avail[i][j]=0;
	label[i][j]=0;
      }
    noPart=0;

    noCluster=0;
    noAggreg=0;

  }// end constructor

  ~Particles()
  {
    for(int i=0;i<nx;i++){
      delete []  deposit[i];
      delete []  avail[i];
      delete []  label[i];
    }

    delete [] deposit;
    delete [] avail;
    delete [] label;
    delete [] sizeOfCluster; 
    delete [] volumeOfCluster;
    delete [] maxHeightOfCluster;
  }

  void randomDepositInit(double rhoDeposit, boost::mt19937 eng)
  { 
    noCluster=0;
    noAggreg=0;
    noColl=0;
    noPart=0;

    int part;
    for(int i=0;i<nx;i++)
      for(int j=0;j<ny;j++){
        part=(int) (rng(eng2) + rhoDeposit);
        deposit[i][j]=part;
        noCluster+=part;
      }
  }


  int countDeposit(){
    int b=0;
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){b+=deposit[i][j];}
    }
    return b;
  }

  int printDeposit(){
    ofstream fout("deposit.dat");

    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){fout<<deposit[i][j]<<" ";}
      fout<<endl;
    }
    fout.close();
  }

  void histoAlbumine(int* histo, int bin, int maxAlbumine){

    int binSize=maxAlbumine/bin;

    for(int i=0;i<=bin;i++)histo[i]=0;

    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
	//	if(deposit[i][j]>0)continue;
	
	for(int k=0;k<=bin;k++){
	  if(avail[i][j]>=maxAlbumine-k*binSize){histo[k]+=1; break;}
	}

      }  // end j loop
    }  // end i loop
  }


  double countDepositedAlbumine(){
    double b=0;
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){b+=avail[i][j];}
    }
    return b;
  }

  double varianceDepositedAlbumine(double average){
    double b=0.;
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
	b+=(avail[i][j]-average)*(avail[i][j]-average);
      }
    }
    return b/(double)(nx*ny);
  }


  void setDeposit(){


    for(int i=0;i<nx;i++)  // clear up label[i][j]
      for(int j=0;j<ny;j++){
	deposit[i][j]=0;//(int) (rng(eng2)+.0005);
      }

    for(int i=10;i<14;i++) 
      for(int j=10;j<14;j++)deposit[i][j]=1;

    for(int i=30;i<32;i++) 
      for(int j=30;j<32;j++)deposit[i][j]=1;


    int test=0;
    for(int i=0;i<nx;i++) 
      for(int j=0;j<ny;j++)if(deposit[i][j]==1)test++;
    //cout <<"test="<<test<<endl;
  }

  //----  Counts the no of clusters, label each pixel of each cluster with a
  // a number going from 1 to no of cluster.
  int countAggregate(){

    int ii,jj;
    int* myStack;
    myStack=new int[2*nx*ny]; // pointer to hold pixels (i,j) that still need 
                          // to be explored for a continuation of the cluster 
    int top=0;   // position in stack

    for(int i=0;i<nx;i++)  // clear up label[i][j]
      for(int j=0;j<ny;j++){
	label[i][j]=0;
      }

    aggregateIndex=0;     // to enumerate the clusters

    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
	if(deposit[i][j]==0)continue;  // this is not a cluster
	if(label[i][j]!=0) continue;   // this an already visited cluster

	myStack[top]=i; top++; myStack[top]=j; top++;
	// Get the same label as your neighbors, if any.
	// Else get a new label given by newIndex

	aggregateIndex++;     // set a unique number to all pixels belonging
                              // to the same cluster
        int clusterSize=0;
        int clusterVolume=0;
	int maxHeight=0.;

	while(top>0){
	  top--; jj=myStack[top]; top--; ii=myStack[top];
	  label[ii][jj]=aggregateIndex;
	  clusterSize++;
	  clusterVolume+=deposit[ii][jj];
	  if(deposit[ii][jj]>maxHeight)maxHeight=deposit[ii][jj];

	  // visit all 8 neighbors of (ii,jj) to see if the cluster continues
	  if(deposit[ii][(jj+ny-1)%ny]!=0 && label[ii][(jj+ny-1)%ny]==0){
	    label[ii][(jj+ny-1)%ny]=-1; // to avoid to push it again on stack
                                         // when seen by another neighbor
	    myStack[top]=ii; top++; 
	    myStack[top]=(jj+ny-1)%ny; top++;}

	  if(deposit[(ii+nx-1)%nx][(jj+1)%ny]!=0 && 
             label[(ii+nx-1)%nx][(jj+1)%ny]==0){
	    label[(ii+nx-1)%nx][(jj+1)%ny]=-1;
	    myStack[top]=(ii+nx-1)%nx; top++; 
	    myStack[top]=(jj+1)%ny; top++;}

	  if(deposit[(ii+nx-1)%nx][jj]!=0 && label[(ii+nx-1)%nx][jj]==0){
	    label[(ii+nx-1)%nx][jj]=-1;
	    myStack[top]=(ii+nx-1)%nx; top++; 
	    myStack[top]=jj; top++;}

	  if(deposit[(ii+nx-1)%nx][(jj+ny-1)%ny]!=0 && 
             label[(ii+nx-1)%nx][(jj+ny-1)%ny]==0){
	    label[(ii+nx-1)%nx][(jj+ny-1)%ny]=-1;
	    myStack[top]=(ii+nx-1)%nx; top++; 
	    myStack[top]=(jj+ny-1)%ny; top++;}

	  if(deposit[(ii+nx+1)%nx][(jj+ny-1)%ny]!=0 && 
             label[(ii+nx+1)%nx][(jj+ny-1)%ny]==0){
	    label[(ii+nx+1)%nx][(jj+ny-1)%ny]=-1;
	    myStack[top]=(ii+nx+1)%nx; top++; 
	    myStack[top]=(jj+ny-1)%ny; top++;}

	  if(deposit[(ii+nx+1)%nx][jj]!=0 && 
             label[(ii+nx+1)%nx][jj]==0){
	    label[(ii+nx+1)%nx][jj]=-1;
	    myStack[top]=(ii+nx+1)%nx; top++; 
	    myStack[top]=jj; top++;}

	  if(deposit[(ii+nx+1)%nx][(jj+ny+1)%ny]!=0 && 
             label[(ii+nx+1)%nx][(jj+ny+1)%ny]==0){
	    label[(ii+nx+1)%nx][(jj+ny+1)%ny]=-1;
	    myStack[top]=(ii+nx+1)%nx; top++; 
	    myStack[top]=(jj+ny+1)%ny; top++;}

	  if(deposit[ii][(jj+ny+1)%ny]!=0 && 
             label[ii][(jj+ny+1)%ny]==0){
	    label[ii][(jj+ny+1)%ny]=-1;
	    myStack[top]=ii; top++; 
	    myStack[top]=(jj+ny+1)%ny; top++;}

	}// end while

	sizeOfCluster[aggregateIndex]=clusterSize;
	volumeOfCluster[aggregateIndex]=clusterVolume;
	maxHeightOfCluster[aggregateIndex]=maxHeight;

      }// end j loop
    } // end i loop

    averageSize=0.;
    averageVol=0.;
    averageHeight=0.;

    if(aggregateIndex>0){
      // ofstream sizeHeightVol("size-height-vol.dat");
      // sizeHeightVol<<"% size, max-height, volume"<<endl;
      // sizeHeightVol<<"startReadData"<<endl;
      for(int i=1;i<=aggregateIndex;i++){
	averageSize+=sizeOfCluster[i];
	averageVol+=volumeOfCluster[i];
	averageHeight+=maxHeightOfCluster[i];
	//	sizeHeightVol<<sizeOfCluster[i]<<"  "<<maxHeightOfCluster[i]<<"  "<<volumeOfCluster[i]<<"  addData"<<endl;
	}
      //      sizeHeightVol<<"endReadData"<<endl;
      //      sizeHeightVol.close();
      
    averageSize=averageSize/(double)aggregateIndex;
    averageVol=averageVol/(double)aggregateIndex;
    //    averageHeight=3.*averageVol/averageSize; // 3 for a cone
    averageHeight=averageHeight/(double)aggregateIndex;
    }

    delete [] myStack;
    return aggregateIndex; // this is the no of aggregate
  } //end countAggregate


// adhesion process: only possible on available sites, with prob pAdhesion
// Particles arrive with probability flux on the surface
// The function returns the fraction of rejected particles: 
// the incoming fraction (flux)  minus those who could adhere
// flux=-pAdhesion*N*dt is the prob of deposition
  double adhesion(double pAdhesion, double attenuation, double maxAlbumine, double flux) 
  { 
    int platelets=0;
    int newdeposit;
    double p;

    // cout<<"adhesion:";
    // for(int i=0;i<20;i++)cout<<rng(eng2)<<" ";
    // cout<<endl;
    // cin.get();

    for(int i=0;i<nx;i++){
      for (int j=0;j<ny;j++){
	if (avail[i][j]>=maxAlbumine)continue;  // site not available

	// deposition takes place with prob pAdhesion on existing platelets
	if(deposit[i][j]!=0){
	  newdeposit=int(rng(eng2)+pAdhesion*flux);
	  //if(newdeposit>1)cout<<newdeposit<<endl;  // check consistency
	  deposit[i][j]+=newdeposit; 
	  platelets+=newdeposit;    // a platelet has been used
	  continue;
	}

	// deposition takes place with prob pAdhesion on a new site
	p=Prob2(pAdhesion,attenuation,maxAlbumine,avail[i][j]);

	newdeposit=0;
	if(rng(eng2)<sqrt(p*flux) && rng(eng2)<sqrt(p*flux))newdeposit=1;

	deposit[i][j]=newdeposit; 
	platelets+=newdeposit;    // a platelet has been used
	noCluster+=newdeposit;
     } // end of j loop
    } // end of i loop

    return (double)platelets/(double)(nx*ny);  // 
  }  // end adhesion



  // Aggregation process: particles can aggregate on top or on the side
  // of an existing cluster. The new particles arrive with probability flux.
  // The function returns the fraction of rejected particles (see above)
  double aggregation(double pAggreg, double attenuation,double pTop, double maxAlbumine, double flux)  
  {                                
    int local_density;
    int platelets=0;
    int newdeposit;
    double p;

    // cout<<"aggregation:";
    // for(int i=0;i<20;i++)cout<<rng(eng2)<<" ";
    // cout<<endl;
    // cin.get();

    for(int i=0;i<nx;i++){
      for (int j=0;j<ny;j++){
	
	// deposit on the side
	local_density=deposit[(i+1)%nx][j]+
                      deposit[(i-1+nx)%nx][j]+
                      deposit[i][(j-1+ny)%ny]+               
                      deposit[i][(j+1)%ny];
	              // deposit[(i+1)%nx][(j+1)%ny]+
	              // deposit[(i+1)%nx][(j-1+ny)%ny]+
	              // deposit[(i-1+nx)%nx][(j+1)%ny]+               
	              // deposit[(i-1+nx)%nx][(j-1+ny)%ny];               

	// adsorb next to an existing deposit
	if(local_density>0 & deposit[i][j]==0){
	  p=Prob2(pAggreg,attenuation,maxAlbumine,avail[i][j]);

	  newdeposit=0;
	  if(rng(eng2)<sqrt(p*flux) && rng(eng2)<sqrt(p*flux))newdeposit=1;

	  deposit[i][j]=newdeposit;
	  platelets+=newdeposit; // the platelets has been consumed
          noAggreg+=newdeposit;
	  continue;
	}

	// deposit on top
	if (deposit[i][j]>0){
	  newdeposit=0;
	  if(rng(eng2)<sqrt(pTop*flux) && rng(eng2)<sqrt(pTop*flux))newdeposit=1;

	  deposit[i][j]+=newdeposit;
	  platelets+=newdeposit; // the platelets has been consumed
	  noAggreg+=newdeposit;}
      }
    }
    return (double) platelets/(double) (nx*ny);
  }     // end Aggregation


  double competingAdhesion(double pfilling, double maxAlbumine, double flux)
  {
    double albumine=0;
    double newPart;
    double p;

    for(int i=0;i<nx;i++){
      for (int j=0;j<ny;j++){
	if (avail[i][j]>=maxAlbumine || deposit[i][j]>0)continue;   // site already filled
	p=Prob(pfilling,maxAlbumine,avail[i][j]);
	newPart=flux*p;
	avail[i][j]+=newPart;
	albumine+=newPart;
      }
    }

    return (double)albumine/(double)(nx*ny);
  }


}; // -----------------------------end class particles

class Parameters
{
public:
  int nx, ny, nz;
  double Activated, NonActivated, Albumine;
  double D;
  double DZ;
  double tau;
  double W;
  double dz;
  double dt;
  double pAdhesion, pAggreg, pTop, pfilling;
  double maxAlbumine;
  double attenuation;
  double Lz;

  void readParameterFile(string filename)
  {
    string comment;
    ifstream fin;

    fin.open(filename.c_str());
    fin>>nx>>ny;     getline(fin,comment);
    fin>>NonActivated;getline(fin,comment);
    fin>>Activated;   getline(fin,comment);
    fin>>Albumine;    getline(fin,comment);
    fin>>D;              getline(fin,comment);
    fin>>DZ;              getline(fin,comment);
    fin>>tau;              getline(fin,comment);
    fin>>W;              getline(fin,comment);
    fin>>pAdhesion;         getline(fin,comment);
    fin>>pAggreg;         getline(fin,comment);
    fin>>pTop;         getline(fin,comment);
    fin>>pfilling;         getline(fin,comment);
    fin>>maxAlbumine;         getline(fin,comment);
    fin>>attenuation;         getline(fin,comment);
    fin>>dt;         getline(fin,comment);
    fin>>Lz;         getline(fin,comment);
    fin.close();

    double w0=1.-2.*W;
    double cs2=(1-w0);
    dz=sqrt(D*dt/(tau-0.5)/cs2);
    nz=Lz/dz;
  }

  void printParameters()
  {
    //    cout <<"nx="<<nx<<"  ny="<<ny<<"  nz="<<nz<<endl;
    //cout <<"nx="<<nx<<"  ny="<<ny<<"  nz="<<nz<<endl; 
    //cout <<"non activated="<<NonActivated<<"  activated="<<Activated<<endl;
    //cout <<"albumine="<<Albumine<<endl;
    //cout <<"Diffusion Constant="<<D<<endl;
    //cout <<"Relaxation time="<<tau<<endl;
    //cout <<"Thrickness of boundary layer="<<DZ<<endl;
    //cout <<"prob Adhesion="<<pAdhesion<<endl;
    //cout <<"prob Aggregation="<<pAggreg<<endl;
    //cout <<"prob Aggregation on top="<<pTop<<endl;
    //cout <<"prob filling="<<pfilling<<endl;
    //cout <<"max number of albumin per cell="<<maxAlbumine<<endl;
    //cout <<"probability attenuation="<<attenuation<<endl;
    //cout <<"time step="<<dt<<endl;
    //cout <<"Blood thickness="<<Lz<<endl;
  }

  void writeParametersToFile(string filename)
  {
    ofstream fout;

    fout.open(filename.c_str());

    fout <<"/nx "<<nx<<" def /ny "<<ny<<" def /nz "<<nz<<" def"<<endl;
    fout <<"%non activated="<<NonActivated<<"  activated="<<Activated<<endl;
    fout <<"%albumine="<<Albumine<<endl;
    fout <<"/D "<<D<<" def"<<endl;
    fout <<"/DZ "<<DZ<<" def"<<endl;
    fout <<"/tau "<<tau<<" def"<<endl;
    fout <<"/pAdhesion "<<pAdhesion<<" def"<<endl;
    fout <<"/pAggreg "<<pAggreg<<"  def"<<endl;
    fout <<"/pTop "<<pTop<<" def"<<endl;
    fout <<"/pfilling "<<pfilling<<" def"<<endl;
    fout <<"/maxAlbumine "<<maxAlbumine<<" def"<<endl;
    fout <<"/attenuation "<<attenuation<<" def"<<endl;
    fout <<"/dt "<<dt<<" def"<<endl;
    fout <<"/dz "<<dz<<" def"<<endl;
    fout <<"/w "<<W<<" def"<<endl;

    fout.close();
  }
};

void Deposition(double* results, unsigned int nrows, unsigned int ncols, double pAd, double pAg, double pT, double pF, double aT, int seed) 
{

  // here we use the random engine eng from main() to 
  // seed the random engine of the deposition process,
  // as used by class Particles. Otherwise, if using
  // directly eng in Particle, the same sequence keep repeating
  boost::random::uniform_int_distribution<int> randint;
  randint=boost::random::uniform_int_distribution<int>(0,2147483647) ;

  double tMax=300; // the simulation will run for thant many seconds

  // units are sqrt(5) microns per lattice site 
  double ds=5.e-12;  // surface of a dx*dy cell in m^2

  //  rng.init();  BC_

  int nx, ny, nz;
  double rhoActivated, rhoNonActivated, rhoAlbumine;
  double D;
  double dt,dz;
  double f1_right, f2_left;
  double Lz;   // thickness of the bloodlayer
  double v;
  double w;    // weight for the d1q3
  double tau;  // relaxation time for the diffusion model
  double W;   // the weight for the distribution functions in the diffusion model

 // deposition substrate
  double DZ;   // thickness of the adsorption boundary layer

  double pAdhesion;
  double pAggreg;
  double pTop;
  double pfilling;  // prob to block a site for platelet adhesion
  double maxAlbumine;   // max number of albumin to saturate a cell
  double attenuation;   // attenuation factor for deposition Prob2

  int iter;
  int Clusters;

  double vol_microliter;
  double mm2;  

  double rejected;  // number of particle that could not deposit
  double accepted;  // number of particle that could deposit
  double J_al, J_ap, J_nap;
  double N_al, N_ap, N_nap;
  double A_al, A_ap, A_nap;

  // int bin=8;
  // int* histo=new int[bin+1];

  //----------- Reading data from file 
  // string filename;
  // filename="parameter.dat";
  // Parameters myParam;
  // myParam.readParameterFile(filename);
  // //  myParam.printParameters();

  // ----------- Initialization of model parameters
  //  nx=myParam.nx; ny=myParam.ny;  // size of the deposition substrate
  nx=447; ny=447;  // size of the deposition substrate

  D=1.3e-8;  // diffusion coefficient 1.e-8
  Lz=.82e-3;  // thickness of the blood layer
  dt=1.e-2;  // time step

  tau=0.8;  // relaxation time for the diffusion process
  W=1./3.;   // the weight for the distribution function
  double w0=1.-2.*W;
  double cs2=(1-w0);
  dz=sqrt(D*dt/(tau-0.5)/cs2); //spatial discretization 

  nz=Lz/dz;  // spatial discretization for the blood layer
  DZ=2.e-5;    // thickness of the boundary layer
  maxAlbumine=100000;

  // Experimental data are given per micro-liter = 10^9 (um)^3
  // One micro-liter = 10^{-9} m^3 = 1 (mm)^3 
  // the total volume is 130 ul
  // Density is computed as the number of particles per m^3, each cell is 
  // a cube of side ds*dz.
  // One then multiplies the number of part per ul by 10^9 to get it per m^3
  rhoNonActivated=172200*1.e9;
  rhoActivated=4808*1.e9;
  rhoAlbumine=2.69e13*1.e9;

  // values given to Deposition() as an argument
  pAdhesion=pAd; //myParam.pAdhesion;
  pAggreg=pAg;//myParam.pAggreg;
  pTop=pT;//myParam.pTop;
  pfilling=pF;//myParam.pfilling;
  attenuation=aT; // attenuation factor for deposition Prob2: attenuation=6


  d1q3 activatePlatelets(nz,rhoActivated,D,dt,dz,W);
  d1q3 nonActivatePlatelets(nz,rhoNonActivated,D,dt,dz,W);
  d1q3 albumine(nz,rhoAlbumine,D,dt,dz,W); 

  Particles deposition(nx,ny,seed); 


  //----------------------conversion to experimental measurement
  vol_microliter=(nx/1000.)*(ny/1000.)*ds*.82;
  mm2=nx*ny*ds*1.e6;  



  J_al=0. ;   J_ap=0. ;  J_nap=0. ;

  // initial values of the boundary layer
  N_ap=rhoActivated*ds*DZ;  
  N_al=rhoAlbumine*ds*DZ; 
  N_nap=rhoNonActivated*ds*DZ;

  A_al=0.;  A_ap=0.;  A_nap=0.;

  int nbIter=tMax/dt; // the number of iterations

  //------------------- time iterations
  // Open file to output data
  //  ofstream fout("output.dat");


  iter=0;
  Clusters=deposition.countAggregate();

  int outputIndex=0;
  results[nrows * outputIndex + 0]=iter*dt;
  results[nrows * outputIndex + 1]=Clusters/mm2;
  results[nrows * outputIndex + 2]=deposition.averageSize*ds*1e12;
  results[nrows * outputIndex + 3]=nonActivatePlatelets.computeNoPart()*1e-6/0.82;
  results[nrows * outputIndex + 4]=activatePlatelets.computeNoPart()*1e-6/0.82;

  // fout <<iter*dt<<" "<<Clusters/mm2<<"  "<<deposition.averageSize*ds*1e12<<" "
  //      <<nonActivatePlatelets.computeNoPart()*1e-6/0.82<<" "
  //      <<activatePlatelets.computeNoPart()*1e-6/0.82<<endl;

  while (iter<nbIter) {
    iter++;

    activatePlatelets.propagation(f2_left,f1_right);
    activatePlatelets.boundary(N_ap/(ds*DZ), f1_right);
    J_ap=activatePlatelets.getJ(0);

    accepted=deposition.adhesion(pAdhesion,attenuation,maxAlbumine,N_ap*dt);

    N_ap-=accepted;    // remove from the interface what has been accepted 
    N_ap+=(-J_ap)*dt*ds;  // number of activated platelets hitting the wall
    A_ap+=accepted;

    albumine.propagation(f2_left,f1_right);
    albumine.boundary(N_al/(ds*DZ), f1_right);
    J_al=albumine.getJ(0);

    accepted=deposition.competingAdhesion(pfilling,maxAlbumine,N_al*dt);

    N_al-=accepted;
    N_al+=(-J_al)*dt*ds;
    A_al+=accepted;

    nonActivatePlatelets.propagation(f2_left,f1_right);
    nonActivatePlatelets.boundary(N_nap/(ds*DZ), f1_right);
    J_nap=nonActivatePlatelets.getJ(0);
    accepted=deposition.aggregation(pAggreg,attenuation,pTop,maxAlbumine,N_nap*dt);
    N_nap-=accepted;
    N_nap+=(-J_nap)*dt*ds;
    A_nap+=accepted;

    activatePlatelets.collision();
    nonActivatePlatelets.collision();
    albumine.collision();

    //------------- measurement branch 
    if(iter*dt==20 || iter*dt==60 || iter*dt==120 || iter*dt==300){

      activatePlatelets.computeRhoU();
      nonActivatePlatelets.computeRhoU();
      albumine.computeRhoU();

      Clusters=deposition.countAggregate(); // we have to run this first

      outputIndex++;
      results[nrows * outputIndex + 0]=iter*dt;
      results[nrows * outputIndex + 1]=Clusters/mm2;
      results[nrows * outputIndex + 2]=deposition.averageSize*ds*1e12;
      results[nrows * outputIndex + 3]=nonActivatePlatelets.computeNoPart()*1e-6/0.82;
      results[nrows * outputIndex + 4]=activatePlatelets.computeNoPart()*1e-6/0.82;

      // fout <<iter*dt<<" "<<Clusters/mm2<<"  "<<deposition.averageSize*ds*1e12<<" "
      // 	   <<nonActivatePlatelets.computeNoPart()*1.e-6/0.82<<" "
      // 	   <<activatePlatelets.computeNoPart()*1.e-6/0.82<<endl;

      // cout <<iter*dt<<" "<<Clusters<<"  "<<deposition.averageSize<<" "
      // 	   <<deposition.averageHeight<<"  "
      // 	   <<nx*ny*ds*nonActivatePlatelets.computeNoPart()<<" "
      // 	   <<nx*ny*ds*activatePlatelets.computeNoPart()<<endl;

    }
    //---------- end measurmement  branch


  }
//------------ end time iteration

//  fout.close();


  //  delete  [] histo;
}// end Deposition


void model(double* results, unsigned int rsize, unsigned int k, double pAd, double pAg, double pT, double pF, double aT, int seed) {
  boost::mt19937 rng(seed);
  for (int i=0; i<k; ++i) {
    int msize = rsize / k;
    int pos = msize * i;
    double* ires = results+pos;
    Deposition(ires, 5, 5, pAd, pAg, pT, pF, aT, rng());
  }
}
