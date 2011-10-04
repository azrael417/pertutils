#include "mathutils.hpp"
#include "pertutils.hpp"

void alpha::set_betavec(){
  betavec.clear();

  betavec.push_back(betazero());
  if(looporder>1){
    betavec.push_back(betaone());
  }
  if(looporder>2){
    betavec.push_back(betatwo());
  }
  if(looporder>3){
    betavec.push_back(betathree());
  }
  betfunc.set(betavec);
}

alpha::alpha(double MUSTART, double ALPHASTART, double MUMIN, double MUMAX, int MUCOUNT, int NF, int LO) : mustart(MUSTART), alphastart(ALPHASTART), mumin(MUMIN), mumax(MUMAX), mucount(MUCOUNT), numflavours(NF), looporder(LO) {
  xvec.clear();
  yvec.clear();
  double rho, dummy;

  double rhostart=log(pow(mustart,2));
  double rhomin=log(pow(mumin,2));
  double rhomax=log(pow(mumax,2));
  double rhostep=(rhomax-rhomin)/((double)mucount);

  set_betavec();

  double alphamin=afunc(alphastart/pimath,rhostart,rhomin,mucount)*pimath;
  xvec.push_back(rhomin);
  yvec.push_back(alphamin);
  for(int i=1; i<=mucount; i++){
    rho=rhomin+(double)i*rhostep;
    xvec.push_back(rho);
    dummy=afunc(yvec[i-1]/pimath,xvec[i-1],rho,100)*pimath;
    yvec.push_back(dummy);
  }
  interpolation=new spline_interp(xvec,yvec);
}

void alpha::change_LO(int loopord){
  looporder=loopord;
  delete [] interpolation;

  xvec.clear();
  yvec.clear();
  double rho, dummy;
  double avector[2];

  mustart=0.9;
  alphastart=alpha_step_from_mz(mustart,10000,avector);

  double rhostart=log(pow(mustart,2));
  double rhomin=log(pow(mumin,2));
  double rhomax=log(pow(mumax,2));
  double rhostep=(rhomax-rhomin)/((double)mucount);

  double alphamin=afunc(alphastart/pimath,rhostart,rhomin,mucount)*pimath;
  xvec.push_back(rhomin);
  yvec.push_back(alphamin);
  for(int i=1; i<=mucount; i++){
    rho=rhomin+(double)i*rhostep;
    xvec.push_back(rho);
    dummy=afunc(yvec[i-1]/pimath,xvec[i-1],rho,100)*pimath;
    yvec.push_back(dummy);
  }
  interpolation=new spline_interp(xvec,yvec);
}

void alpha::change_NF(int nf){
  numflavours=nf;
  delete [] interpolation;

  xvec.clear();
  yvec.clear();
  double rho, dummy;
  double avector[2];

  mustart=0.9;
  alphastart=alpha_step_from_mz(mustart,10000,avector);

  double rhostart=log(pow(mustart,2));
  double rhomin=log(pow(mumin,2));
  double rhomax=log(pow(mumax,2));
  double rhostep=(rhomax-rhomin)/((double)mucount);

  double alphamin=afunc(alphastart/pimath,rhostart,rhomin,mucount)*pimath;
  xvec.push_back(rhomin);
  yvec.push_back(alphamin);
  for(int i=1; i<=mucount; i++){
    rho=rhomin+(double)i*rhostep;
    xvec.push_back(rho);
    dummy=afunc(yvec[i-1]/pimath,xvec[i-1],rho,100)*pimath;
    yvec.push_back(dummy);
  }
  interpolation=new spline_interp(xvec,yvec);
}

alpha::alpha(double MUMIN, double MUMAX, int MUCOUNT, int NF, int LO) : mumin(MUMIN), mumax(MUMAX), mucount(MUCOUNT), numflavours(NF), looporder(LO){
  xvec.clear();
  yvec.clear();
  double rho, dummy;
  double avector[2];

  mustart=0.9;
  alphastart=alpha_step_from_mz(mustart,10000,avector);

  double rhostart=log(pow(mustart,2));
  double rhomin=log(pow(mumin,2));
  double rhomax=log(pow(mumax,2));
  double rhostep=(rhomax-rhomin)/((double)mucount);

  double alphamin=afunc(alphastart/pimath,rhostart,rhomin,mucount)*pimath;
  xvec.push_back(rhomin);
  yvec.push_back(alphamin);
  for(int i=1; i<=mucount; i++){
    rho=rhomin+(double)i*rhostep;
    xvec.push_back(rho);
    dummy=afunc(yvec[i-1]/pimath,xvec[i-1],rho,100)*pimath;
    yvec.push_back(dummy);
  }
  interpolation=new spline_interp(xvec,yvec);
}

alpha::alpha(std::string filename){
  std::ifstream input;
  double token;

  input.open(filename.c_str());
  if(!input.good()){
    std::cerr << "alpha: error, cannot open " << filename << "for reading!" << std::endl;
  }
  else{
    xvec.clear();
    yvec.clear();
    input >> mustart;
    input >> alphastart;
    input >> mumin;
    input >> mumax;
    input >> mucount;
    for(int i=0; i<=mucount; i++){
      input >> token;
      xvec.push_back(token);
      input >> token;
      yvec.push_back(token);
    }
    input.close();
    interpolation=new spline_interp(xvec,yvec);
  }
}

//Print stored xvec,yvec into file:
int alpha::writefile(std::string filename){
  std::ofstream output;

  output.open(filename.c_str());
  if(!output.good()){
    std::cerr << "alpha: error, cannot open " << filename << " for writing!" << std::endl;
    return EXIT_FAILURE;
  }
  else{
    output << mustart << std::endl;
    output << alphastart << std::endl;
    output << mumin << std::endl;
    output << mumax << std::endl;
    output << mucount << std::endl;
    for(unsigned int i=0; i<xvec.size(); i++){
      output << xvec[i] << " " << yvec[i] << std::endl;
    }
  }
  output.close();
  return EXIT_SUCCESS;
}

double alpha::get_mumin(){
  return mumin;
}

double alpha::get_mumax(){
  return mumax;
}

double alpha::get_mustart(){
  return mustart;
}

double alpha::get_alphastart(){
  return alphastart;
}

double alpha::operator()(const double mu){
  double rho=log(pow(mu,2));
  double result=interpolation->interp(rho);
  return result;
}

double alpha::get(const double mu){
  double rho=log(pow(mu,2));
  double result=interpolation->interp(rho);
  return result;
}

alpha::~alpha(){
  xvec.clear();
  yvec.clear();
}

double alpha::alpha_step_from_mz(double muziel, int numsteps, double* avec){
  double must=mz;
  double alphast=alphas_mz;
  double alphaerror=alphas_mz_error;
  int nflavour=numflavours;

  double rhostart=log(must*must);
  double rhoziel=log(muziel*muziel);
  double rhozielsave=rhoziel;
  double astart=alphast/pimath;
  double aerror=alphaerror/pimath;
  double* avalues = new double[4];
  double asave, aerrorsave, rhostep;

  double* xvalues=new double[numsteps];
  double* yvalues=new double[numsteps];
  double* yerrorplus=new double[numsteps];
  double* yerrorminus=new double[numsteps];

  if(rhozielsave>=log(mb_msbar*mb_msbar)){
    numflavours=5;
    set_betavec();

    rhostep=(rhoziel-rhostart)/numsteps;
    avalues[0]=astart;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yvalues,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yvalues,numsteps);                                                                                                                                                                             
    asave=yvalues[numsteps-1];

    avalues[0]=astart+aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);                                                                                                                                                                          

    avalues[0]=astart-aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);                                                                                                                                                                         
    aerrorsave=(yerrorplus[numsteps-1]-yerrorminus[numsteps-1])/2.;
  }
  else if(rhozielsave>=log(mc_msbar*mc_msbar)){
    //Rennen auf m_b:                                                                                                                                                                                                                     
    numflavours=5;
    set_betavec();

    rhoziel=log(mb_msbar*mb_msbar);
    rhostep=(rhoziel-rhostart)/numsteps;

    //Central Value:                                                                                                                                                                                                                      
    avalues[0]=astart;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yvalues,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yvalues,numsteps);                                                                                                                                                                             
    asave=yvalues[numsteps-1];

    //Obere und untere Grenze:                                                                                                                                                                                                            
    avalues[0]=astart+aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);                                                                                                                                                                          

    avalues[0]=astart-aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);                                                                                                                                                                         
    aerrorsave=(yerrorplus[numsteps-1]-yerrorminus[numsteps-1])/2.;
    astart=asave;
    aerror=aerrorsave;

    std::cout << "alpha_s(m_b^2= " << mb_msbar << "^2) = " << astart*pimath << " +/- " << aerror*pimath << " " << std::endl;
    
    //Rennen vor m_c:                                                                                                                                                                                                                     
    numflavours=4;
    set_betavec();

    rhostart=log(mb_msbar*mb_msbar);
    rhoziel=rhozielsave;
    rhostep=(rhoziel-rhostart)/numsteps;

    //Central Value:                                                                                                                                                                                                                      
    avalues[0]=astart;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yvalues,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yvalues,numsteps);                                                                                                                                                                             
    asave=yvalues[numsteps-1];

    //Obere und untere Grenze:                                                                                                                                                                                                            
    avalues[0]=astart+aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);                                                                                                                                                                          

    avalues[0]=astart-aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);                                                                                                                                                                         
    aerrorsave=(yerrorplus[numsteps-1]-yerrorminus[numsteps-1])/2.;
  }
  else{
    //Rennen auf m_b:                                                                                                                                                                                                                     
    numflavours=5;
    set_betavec();

    rhoziel=log(mb_msbar*mb_msbar);
    rhostep=(rhoziel-rhostart)/numsteps;

    //Central Value:                                                                                                                                                                                                                      
    avalues[0]=astart;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yvalues,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yvalues,numsteps);                                                                                                                                                                             
    asave=yvalues[numsteps-1];

    //Obere und untere Grenze:                                                                                                                                                                                                            
    avalues[0]=astart+aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);                                                                                                                                                                          

    avalues[0]=astart-aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);                                                                                                                                                                         
    aerrorsave=(yerrorplus[numsteps-1]-yerrorminus[numsteps-1])/2.;
    astart=asave;
    aerror=aerrorsave;

    std::cout << "alpha_s(m_b^2= " << mb_msbar << "^2) = " << astart*pimath << " +/- " << aerror*pimath << std::endl;
    
    //Rennen auf m_c:                                                                                                                                                                                                                     
    numflavours=4;
    set_betavec();

    rhostart=log(mb_msbar*mb_msbar);
    rhoziel=log(mc_msbar*mc_msbar);
    rhostep=(rhoziel-rhostart)/numsteps;

    rhostart=log(mb_msbar*mb_msbar);
    rhoziel=log(mc_msbar*mc_msbar);
    rhostep=(rhoziel-rhostart)/numsteps;

    //Central Value:                                                                                                                                                                                                                      
    avalues[0]=astart;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yvalues,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yvalues,numsteps);                                                                                                                                                                             
    asave=yvalues[numsteps-1];

    //Obere und untere Grenze:                                                                                                                                                                                                            
    avalues[0]=astart+aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);                                                                                                                                                                          

    avalues[0]=astart-aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);                                                                                                                                                                         
    aerrorsave=(yerrorplus[numsteps-1]-yerrorminus[numsteps-1])/2.;
    astart=asave;
    aerror=aerrorsave;

    std::cout << "alpha_s(m_c^2= " << mc_msbar << "^2) = ) = " << astart*pimath << " +/- " << aerror*pimath << std::endl;

    //Rennen unter mc_c:                                                                                                                                                                                                                  
    numflavours=3;
    set_betavec();

    rhostart=log(mc_msbar*mc_msbar);
    rhoziel=rhozielsave;
    rhostep=(rhoziel-rhostart)/numsteps;

    //Central Value:                                                                                                                                                                                                                      
    avalues[0]=astart;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yvalues,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yvalues,numsteps);                                                                                                                                                                             
    asave=yvalues[numsteps-1];

    //Obere und untere Grenze:                                                                                                                                                                                                            
    avalues[0]=astart+aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorplus,numsteps);                                                                                                                                                                          

    avalues[0]=astart-aerror;
    avalues[1]=rhostep*betafunc(avalues[0]);
    avalues[2]=(rhostep*rhostep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(rhostep*rhostep*rhostep/6.)*d2betafunc(avalues[0]);

    mvalPC(betfunc,avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);
    //rk4(avalues,rhostart,rhostep,xvalues,yerrorminus,numsteps);                                                                                                                                                                         
    aerrorsave=(yerrorplus[numsteps-1]-yerrorminus[numsteps-1])/2.;
  }
  std::cout << "alpha_s(" << muziel <<"^2 GeV^2) = " << asave*pimath << " +/- " << aerrorsave*pimath << std::endl;

  numflavours=nflavour;
  set_betavec();

  delete [] avalues;
  delete [] xvalues;
  delete [] yvalues;
  delete [] yerrorplus;
  delete [] yerrorminus;

  avec[0]=asave;
  avec[1]=aerrorsave;

  return asave*pimath;
}

//QCD Beta-function at 4-loop (hep-ph/9910332,p.20, 72):
double alpha::betazero(){
  return (11.-2.*(double)numflavours/3.)/4.;
}

double alpha::betaone(){
  return (102.-38.*(double)numflavours/3.)/16.;
}

double alpha::betatwo(){
  return (2857./2.-5033.*(double)numflavours/18.+325.*(double)numflavours*(double)numflavours/54.)/64.;
}

double alpha::betathree(){
  return (149753./6.+3564.*zeta3-(1078361./162.+6508.*zeta3/27.)*(double)numflavours+(50065./162.+6472.*zeta3/81.)*(double)numflavours*(double)numflavours+1093./729.*(double)numflavours*(double)numflavours*(double)numflavours)/256.;
}

double alpha::betafunc(double a){
  double result=0.0;
  for(unsigned int i=0; i<betavec.size(); i++){
    result-=betavec[i]*pow(a,(i+2));
  }
  return result;
}

double alpha::pdbetafunc(double a){
  double result=0.0;
  for(unsigned int i=0; i<betavec.size(); i++){
    result-=(double)(i+2)*betavec[i]*pow(a,(i+1));
  }
  return result;
}

double alpha::dbetafunc(double a){
  double result=pdbetafunc(a);
  return result*betafunc(a);
}

double alpha::pd2betafunc(double a){
  double result=0.0;
  for(unsigned int i=0; i<betavec.size(); i++){
    result-=(double)(i+2)*(i+1)*betavec[i]*pow(a,i);
  }
  return result;
}

double alpha::d2betafunc(double a){
  double result=pdbetafunc(a)*dbetafunc(a);
  result+=pow(betafunc(a),2)*pd2betafunc(a);
  return result;
}

double alpha::afunc(double astart, double sstart, double sziel, int numsteps){ //liefert a(log(mu^2)) bei gegebenem a(log(mu_0^2)) und log(mu_0^2).                                                                                                  
  double* avalues;
  double sstep=(sziel-sstart)/numsteps;
  double result;
  double *xvalues, *yvalues;

  if(fabs(sziel-sstart)<1e-7){
    result=astart;
  }
  else{
    avalues = new double[4];
    xvalues=new double[numsteps];
    yvalues=new double[numsteps];

    avalues[0]=astart;
    avalues[1]=sstep*betafunc(avalues[0]);
    avalues[2]=(sstep*sstep/2.)*dbetafunc(avalues[0]);
    avalues[3]=(sstep*sstep*sstep/6.)*d2betafunc(avalues[0]);

    betfunc.set(betavec);
    mvalPC(betfunc,avalues,sstart,sstep,xvalues,yvalues,numsteps);
    result=yvalues[numsteps-1];

    delete[] avalues;
    delete [] xvalues;
    delete [] yvalues;
  }
  return result;
}

