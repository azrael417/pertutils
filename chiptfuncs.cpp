#include "mathutils.hpp"
#include "pertutils.hpp"

//******************************************************************
//******************************************************************
//START chiPT FV functions
//******************************************************************
//******************************************************************
Xifunc::Xifunc(double MLSQMIN, double MLSQMAX, int MLCOUNT, double ss) : MLsqmin(MLSQMIN), MLsqmax(MLSQMAX), MLcount(MLCOUNT), s(ss) {
    MLstep=(MLsqmax-MLsqmin)/((double)MLcount);
    
    xvec.clear();
    yvec.clear();
    double MLsq, dummy;
    for(int i=0; i<=MLcount; i++){
        MLsq=MLsqmin+(double)i*MLstep;
        xvec.push_back(MLsq);
        dummy=Xifuncnorm(MLsq,s);
        yvec.push_back(dummy);
    }
    interpolation=new spline_interp(xvec,yvec);
}

Xifunc::Xifunc(std::string filename){
    std::ifstream input;
    double token;
    
    input.open(filename.c_str());
    if(!input.good()){
        std::cerr << "Xifunc: error, cannot open " << filename << "for reading!" << std::endl;
    }
    else{
        xvec.clear();
        yvec.clear();
        input >> MLsqmin;
        input >> MLsqmax;
        input >> MLcount;
        input >> s;
        for(int i=0; i<=MLcount; i++){
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
int Xifunc::writefile(std::string filename){
    std::ofstream output;
	output.precision(10);
    output.setf(std::ios::scientific,std::ios::floatfield);
    
    output.open(filename.c_str());
    if(!output.good()){
        std::cerr << "Xifunc: error, cannot open " << filename << " for writing!" << std::endl;
        return EXIT_FAILURE;
    }
    else{
        output << MLsqmin << std::endl;
        output << MLsqmax << std::endl;
        output << MLcount << std::endl;
        output << s << std::endl;;
        for(unsigned int i=0; i<xvec.size(); i++){
            output << xvec[i] << " " << yvec[i] << std::endl;
        }
    }
    output.close();
    return EXIT_SUCCESS;
}

double Xifunc::get_MLsqmin(){
    return MLsqmin;
}

double Xifunc::get_MLsqmax(){
    return MLsqmax;
}

double Xifunc::get_s(){
    return s;
}

double Xifunc::operator()(const double Msq, const double L){
    double result=pow(L,(2.*s-3.))*(interpolation->interp(Msq*L*L));
    return result;
}

double Xifunc::get(const double Msq, const double L){
    double result=pow(L,(2.*s-3.))*(interpolation->interp(Msq*L*L));
    return result;
}

Xifunc::~Xifunc(){
    xvec.clear();
    yvec.clear();
}

double Xifuncintegrandtransformed::operator()(const double t){
    double tau;
    double dtaudt;
    double result;
    
    //Print XI(ML,s)
    tau=exp(t-exp(-t));
    dtaudt=tau*(1.+exp(-t));
    result=Xifuncintegrand(tau,MLsq,s)*dtaudt;
    
    //Print dXI(ML,s)/dML
    if(flag !=0 ){
        result*=((-1.)*tau*MLsq/(4.*pimath*pimath));
    }
    return result;
}

double Xifuncintegrand(double tau, double MLsq,  double s){
    double result;
    
    result=pow(EllipticTheta(tau),3)-pow(pimath/tau,3./2.);
    result*=pow(tau,(s-1.))*exp(-tau*MLsq/pow(2.*pimath,2.));
    
    return result;
}

double Xifuncnorm(double MLsq, double s, int flag){
    double result=0.;
    Xifuncintegrandtransformed Xif(MLsq,s,flag);
    
    result=qtrap(Xif, 0.0, 6.);
    result/=(pow(2.*pimath,(2.*s))*gammfunc(s));
    
    return result;
}

//Stephan Dürrs functions:
double Ikfunc2(int n, double mpiL){
    Bessik bessel;
    return -1.5*bessel.k1(sqrt((double)n)*mpiL);
}

double Ipifunc2(int n, double mpiL){
    Bessik bessel;
    return -4.*bessel.k1(sqrt((double)n)*mpiL);
}

double mpiLovermpi(double mpi, double fpi, double L){
    double eps=1e-10;
    int nmax=250;
    int nup=5;
    double beta=1.;
    Levin acc(nmax,eps);
    
    double sn=1.;
    for(int n=1; n<=nup; n++){
        sn+=(-1.)*(double)multiplicities[n]*mpi/(32.*fpi*fpi*pimath*pimath*sqrt((double)n)*L)*Ipifunc2(n,mpi*L);
    }
    
    //Levin u transformation:
    double result, omegan, an;
    int nrun=nup+1;
    while( (!acc.cnvgd) && ((nrun-(nup+1))<nmax) ){
        if(multiplicities[nrun]==0){
            nrun++;
            continue;
        }
        else{
            an=(-1.)*(double)multiplicities[nrun]*mpi/(32.*fpi*fpi*pimath*pimath*sqrt((double)nrun)*L)*Ipifunc2(nrun,mpi*L);
            sn+=an;
            omegan=(beta+(double)nrun)*an;
            result=acc.next(sn,omegan);
            nrun++;
        }
    }
    
    //Switching to direct resummation if accelleration fails:
    if(isnan(result)){
        result=1.;
        for(int n=1; n<=nmax; n++){
            result+=(-1.)*(double)multiplicities[n]*mpi/(32.*fpi*fpi*pimath*pimath*sqrt((double)n)*L)*Ipifunc2(n,mpi*L);
        }
    }
    return result;
}

double mkLovermk(double mpi, double mk, double fpi, double L){
    double eps=1e-10;
    int nmax=250;
    int nup=5;
    double beta=1.;
    Levin acc(nmax,eps);
    
    double sn=1.;
    for(int n=1; n<=nup; n++){
        sn+=(-1.)*(double)multiplicities[n]*mpi*mpi/(32.*mk*fpi*fpi*pimath*pimath*sqrt((double)n)*L)*Ikfunc2(n,mpi*L);
    }
    
    //Levin u transformation:
    double result, omegan, an;
    int nrun=nup+1;
    while( (!acc.cnvgd) && ((nrun-(nup+1))<nmax) ){
        if(multiplicities[nrun]==0){
            nrun++;
            continue;
        }
        else{
            an=(-1.)*(double)multiplicities[nrun]*mpi*mpi/(32.*mk*fpi*fpi*pimath*pimath*sqrt((double)nrun)*L)*Ikfunc2(nrun,mpi*L);
            sn+=an;
            omegan=(beta+(double)nrun)*an;
            result=acc.next(sn,omegan);
            nrun++;
        }
    }
    
    //Switching to direct resummation if accelleration fails:
    if(isnan(result)){
        result=1.;
        for(int n=1; n<=nmax; n++){
            result+=(-1.)*(double)multiplicities[n]*mpi*mpi/(32.*mk*fpi*fpi*pimath*pimath*sqrt((double)n)*L)*Ikfunc2(n,mpi*L);
        }
    }
    return result;
}

double fkLoverfk(double mpi, double fpi, double fk, double L){
    double eps=1e-10;
    int nmax=250;
    int nup=5;
    double beta=1.;
    Levin acc(nmax,eps);
    
    double sn=1.;
    for(int n=1; n<=nup; n++){
        sn+=(double)multiplicities[n]*mpi/(16*pimath*pimath*fk*fpi*sqrt((double)n)*L)*Ikfunc2(n,mpi*L);
    }
    
    //Levin u transformation:                                                                                                                                                        
    double result, omegan, an;
    int nrun=nup+1;
    while( (!acc.cnvgd) && ((nrun-(nup+1))<nmax) ){
        if(multiplicities[nrun]==0){
            nrun++;
            continue;
        }
        else{
            an=(double)multiplicities[nrun]*mpi/(16*pimath*pimath*fk*fpi*sqrt((double)nrun)*L)*Ikfunc2(nrun,mpi*L);
            sn+=an;
            omegan=(beta+(double)nrun)*an;
            result=acc.next(sn,omegan);
            nrun++;
        }
    }
    
    //Switching to direct resummation if accelleration fails:                                                                                                                        
    if(isnan(result)){
        result=1.;
        for(int n=1; n<=nmax; n++){
            result+=(double)multiplicities[n]*mpi/(16*pimath*pimath*fk*fpi*sqrt((double)n)*L)*Ikfunc2(n,mpi*L);
        }
    }
    return result;
}

double fpiLoverfpi(double mpi, double fpi, double L){
    double eps=1e-10;
    int nmax=250;
    int nup=2;
    double beta=1.;
    Levin acc(nmax,eps);
    
    double sn=1.;
    for(int n=1; n<=nup; n++){
        sn+=(double)multiplicities[n]*mpi/(16*pimath*pimath*fpi*fpi*sqrt((double)n)*L)*Ipifunc2(n,mpi*L);
    }
    
    //Levin u transformation:
    double result=0., omegan, an;
    int nrun=nup+1;
    while( (!acc.cnvgd) && ((nrun-(nup+1))<nmax) ){
        if(multiplicities[nrun]==0){
            nrun++;
            continue;
        }
        else{
            an=(double)multiplicities[nrun]*mpi/(16*pimath*pimath*fpi*fpi*sqrt((double)nrun)*L)*Ipifunc2(nrun,mpi*L);
            sn+=an;
            omegan=(beta+(double)nrun)*an;
            result=acc.next(sn,omegan);
            nrun++;
        }
    }
    
    //Switching to direct resummation if accelleration fails:
    if(isnan(result)){
        result=1.;
        for(int n=1; n<=nmax; n++){
            result+=(double)multiplicities[n]*mpi/(16*pimath*pimath*fpi*fpi*sqrt((double)n)*L)*Ipifunc2(n,mpi*L);
        }
    }
    return result;
}
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
double Zetafunc::integrand1(const double t){
	return exp(t*qsq)/sqrt(t);
}

//removd pi^{3/2+l} from integrand since it will be computed often in the integral and sum
double Zetafunc::integrand2(const double t){
	return pow(t,-(1.5+(double)l))*exp(t*qsq)*exp(-pimath*pimath*ghatwnorm*ghatwnorm/t);
}

dcomplex Zetafunc::term2(double q2){
	double dresult,tmp;
	qsq=q2;
	
	if(l!=0) dresult=0.;
	else{
		tmp=2*qsq*qromb(integrand1,0.,lambda);
		tmp-=2./sqrt(lambda)*exp(lambda*qsq);
		dresult=gamma*pimath/2.*tmp;
	}
	dcomplex result(dresult,0.);
	return result;
}

dcomplex Zetafunc::term3(double q2){
	double tmp;
	dcomplex result;
	threevec<double> w,ghatw,wpar,wperp;
	qsq=q2;
	for(int z=-MAXRUN; z<=MAXRUN; z++){
		w.set_coord(2)=(double)z;
		for(int y=-MAXRUN; y<=MAXRUN; y++){
			w.set_coord(1)=(double)y;
			for(int x=-MAXRUN; x<=MAXRUN; x++){
				w.set_coord(0)=(double)x;
				orthogonal_projection(w,boostvec,wpar,wperp);
				ghatw=gamma*wpar+wperp;
				ghatwnorm=ghatw.norm();
				tmp=qromb(integrand2,0.,lambda);
				tmp*=pow(ghatwnorm,l);
			}
		}
	}
	result*=gamma*pow(pimath,(double)(1.5+l));
	
	return result;
}
//******************************************************************
//******************************************************************
//END Zeta-function on a torus
//******************************************************************
//******************************************************************