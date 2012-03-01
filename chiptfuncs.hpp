#ifndef _CHIPTFUNCS
#define _CHIPTFUNCS

//******************************************************************
//******************************************************************
//START chiPT FV functions
//******************************************************************
//******************************************************************
class Xifunc{
private:
  double MLsqmin, MLsqmax, MLcount, MLstep, s;
  std::vector<double> xvec, yvec;
  spline_interp* interpolation;

public:
  Xifunc(double MLSQMIN, double MLSQMAX, int MLCOUNT, double ss);
  Xifunc(std::string filename);
  int writefile(std::string filename);
  double operator()(const double Msq, const double L);
  double get(const double Msq, const double L);
  double get_MLsqmin();
  double get_MLsqmax();
  double get_s();
  ~Xifunc();
};

class Xifuncintegrandtransformed : public TFunctor{
 private:
  double MLsq, s;
  int flag;
  
 public:
 Xifuncintegrandtransformed(const double MLML, const double ss, int fflag) : MLsq(MLML), s(ss), flag(fflag) {};
  double operator()(const double t);
};

double Xifuncintegrand(double tau, double MLsq,  double s);
double Xifuncnorm(double,double,int flag=0);

//Colangelo+Stephan
const int multiplicities[301]={1,6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24,48,24,0,24,30,72,32,0,72,48,0,12,48,48,48,30,24,72,0,24,96,48,24,24,72,48,0,8,54,84,48,24,72,96,0,48,48,24,72,0,72,96,0,6,96,96,24,48,96,48,0,36,48,120,56,24,96,48,0,24,102,48,72,48,48,120,0,24,144,120,48,0,48,96,0,24,48,108,72,30,168,48,0,72,96,72,72,32,72,144,0,0,96,96,48,72,120,72,0,48,78,120,48,0,144,144,0,12,144,48,120,48,48,168,0,48,96,96,72,48,96,48,0,30,96,192,56,24,168,120,0,72,144,96,96,0,72,96,0,24,192,108,24,96,96,120,0,48,78,144,120,24,168,144,0,24,48,96,120,72,120,144,0,48,192,144,48,0,192,48,0,8,48,240,96,54,120,120,0,84,144,72,96,48,96,240,0,24,240,96,72,72,96,72,0,96,96,120,96,0,192,144,0,48,150,96,120,48,120,240,0,24,144,216,48,72,144,96,0,0,144,132,104,72,168,144,0,96,144,144,168,0,48,192,0,6,192,96,96,96,216,72,0,96,96,240,48,24,264,192,0,48,96,144,120,96,72,168,0,48,240,96,72,0,192,144,0,36,102,240,96,48,216,168,0,120,192,72,192,56};
double Ikfunc2(int n, double mpiL);
double Ipifunc2(int n, double mpiL);
double mpiLovermpi(double mpi, double fpi, double L);
double mkLovermk(double mpi, double mk, double fpi, double L);
double fkLoverfk(double mpi, double fpi, double fk, double L);
double fpiLoverfpi(double mpi, double fpi, double L);
//******************************************************************
//******************************************************************
//END chiPT FV functions
//******************************************************************
//******************************************************************


//******************************************************************
//******************************************************************
//START Zeta-function on a torus
//******************************************************************
//******************************************************************
//implementation of expression (5) in arXiv:1107.5023
class Zetafunc{
private:
	threevec<double> boost;
	double gamma,lambda, qsq, ghatwnorm;
	int l,m;
	int MAXRUN;
	double integrand2(const double t);
	dcomplex term1(double q2);
	double term1zeroboost(double q2);
	double term2(double q2);
	dcomplex term3(double q2);
	double term3zeroboost(double q2);
	
public:
	Zetafunc(double gammaa, int ll=0, int mm=0, threevec<double> boostvec=threevec<double>(0.,0.,0.), double lambdaa=1., int maxrun=1000) : boost(boostvec), gamma(gammaa), lambda(lambdaa), l(ll), m(mm), MAXRUN(maxrun) {};
	dcomplex operator()(const double q2);
};
//******************************************************************
//******************************************************************
//END Zeta-function on a torus
//******************************************************************
//******************************************************************
#endif
