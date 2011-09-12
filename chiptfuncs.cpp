#include "mathutils.hpp"
#include "pertutils.hpp"

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
    interpolation=new general_spline_interp(xvec,yvec);
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
        interpolation=new general_spline_interp(xvec,yvec);
    }
}

//Print stored xvec,yvec into file:
int Xifunc::writefile(std::string filename){
    std::ofstream output;
    
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
    double result, omegan, an;
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
