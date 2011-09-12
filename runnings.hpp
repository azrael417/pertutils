#ifndef _RUNNINGS
#define _RUNNINGS

//Konstanten:
const double zeta3=1.20205690;
const double zeta4=1.08232323;
const double zeta5=1.03692776;

const double LQCD=0.24243; //hep-lat/9810063, (6.2) with r0=0.49fm                                                                                                                                                                            
const double alphas_mz=0.1184;
const double alphas_mz_error=0.0007;
const double mz=91.1876;
const double mz_error=0.0021;
const double mb_msbar=4.19;
const double mb_uperror_msbar=0.18;
const double mb_downerror_msbar=0.06;
const double mc_msbar=1.27;
const double mc_uperror_msbar=0.07;
const double mc_downerror_msbar=0.09;

const double scale_conversion=0.197327;

struct betafunction : public TFunctor{
  std::vector<double> betavec;
  betafunction(std::vector<double> betas) : betavec(betas) {}
  betafunction(){}
  void set(std::vector<double> betas){ betavec=betas; }
  double operator()(const double a){
    double result=0.0;
    for(int i=0; i<betavec.size(); i++){
      result-=betavec[i]*pow(a,(i+2));
    }
    return result; 
  }
};

class alpha{

 private:
  betafunction betfunc;
  double mumin, mumax, mustart, alphastart;
  int mucount, numflavours, looporder;
  std::vector<double> xvec, yvec, betavec;
  general_spline_interp* interpolation;
  void set_betavec();
  double afunc(double astart, double sstart, double sziel, int numsteps);
  double alpha_step_from_mz(double muziel, int numsteps, double* avec);
  double betazero();
  double betaone();
  double betatwo();
  double betathree();
  double betafunc(double a);
  double pdbetafunc(double a);
  double dbetafunc(double a);
  double pd2betafunc(double a);
  double d2betafunc(double a);

 public:
  alpha(double MUSTART, double ALPHASTART, double MUMIN, double MUMAX, int MUCOUNT, int NF, int LO);
  alpha(double MUMIN, double MUMAX, int MUCOUNT, int NF, int LO);
  alpha(std::string filename);
  void change_LO(int loopord);
  void change_NF(int nf);
  int writefile(std::string filename);
  double operator()(const double mu);
  double get(const double mu);
  double get_mumin();
  double get_mumax();
  double get_mustart();
  double get_alphastart();
  ~alpha();

};
#endif
