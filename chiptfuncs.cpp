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
void Zetafunc::set_numterms_sum(const int maxrun){
	MAXRUN=maxrun;
}

void Zetafunc::set_boost(const threevec<double> boostvec, const double gammaa){
	boost=boostvec;
	gamma=gammaa;
	if(boost.norm()<1.e-5){
		is_zeroboost=true;
		gamma=1.;
	}
	else is_zeroboost=false;
}

dcomplex Zetafunc::operator()(const double q2){
	dcomplex result;
	//special treatment since sums are real:
	if(is_improved){
		double dresult=term1improved(q2)+term3improved(q2);
		result(dresult,0.);
	}
	else{
		result=term1(q2)+term3(q2);
	}
	result+=term2(q2);
	
	return result;
}

//this is term1, which is the simple sum:
dcomplex Zetafunc::term1(const double q2){
	dcomplex result(0.,0.),tmpcomp;
	double fact,rsq,theta,phi;
	threevec<double> nvec,npar,nperp,rvec;
	
	for(int z=-MAXRUN; z<=MAXRUN; z++){
		for(int y=-MAXRUN; y<=MAXRUN; y++){
			for(int x=-MAXRUN; x<=MAXRUN; x++){
				nvec((double)x,(double)y,(double)z);
				
				//compute r:
				if(!is_zeroboost){
					orthogonal_projection(nvec,boost,npar,nperp);
					rvec=npar/gamma-boost/(2.*gamma)+nperp;
				}
				else{
					rvec=nvec;
				}
				rvec.get_spherical_coordinates(rsq,theta,phi);
				if(l!=0)fact=pow(rsq,l); //compute |r|^l
				else fact=1.;
				rsq*=rsq; //compute r^2
				
				//compute prefact:
				fact*=exp(-lambda*(rsq-q2))/(rsq-q2); //compute the ratio
				tmpcomp=fact*spherical_harmonicy(l,m,theta,phi); //multiply with spherical harmonics
				result+=tmpcomp;
			}
		}
	}
	
	return result;
}

double Zetafunc::term1improved(const double q2){
	double result=0.,rsq;
	threevec<double> nvec;
	
	//zero term:
	result=-exp(lambda*q2)/q2;
	
	//x!=0, y=z=0, 3 times for different combinations (replace x e.g. with y) and factor 2 for +/-:
	for(int x=1; x<=MAXRUN; x++){
		rsq=x*x;
		result+=6.*exp(-lambda*(rsq-q2))/(rsq-q2)*pow(fabs(x),l);
	}
	
	//x,y!=0,z=0, 3 times for different combinations and factor 4 for +/-:
	for(int x=1; x<=MAXRUN; x++){
		for(int y=1; y<=MAXRUN; y++){
			rsq=x*x+y*y;
			result+=12.*exp(-lambda*(rsq-q2))/(rsq-q2)*pow(sqrt(rsq),l);
		}
	}
	
	//x,y,z!=0 and factor 8 for +/-:
	for(int x=1; x<=MAXRUN; x++){
		for(int y=1; y<=MAXRUN; y++){
			for(int z=1; z<=MAXRUN; z++){
				rsq=x*x+y*y+z*z;
				result+=8.*exp(-lambda*(rsq-q2))/(rsq-q2)*pow(sqrt(rsq),l); //compute the ratio
			}
		}
	}
	result/=sqrt(4.*pimath);

	return result;
}

//this term is proportional to Dawson's Integral: term2=2/q *exp(lambda*q^2)*F(sqrt(lambda)q)
//where F(x)=exp(-x^2) int_0^x ds exp(-s^2):
double Zetafunc::term2(const double q2){
	double result=0.,tmp;
	
	if(l==0){
		tmp=4.*sqrt(q2)*exp(lambda*q2)*dawson(sqrt(lambda*q2));
		tmp-=2./sqrt(lambda)*exp(lambda*q2);
		result=gamma*pimath/2.*tmp;
	}
	return result;
}

dcomplex Zetafunc::term3(const double q2){
	double tmp,r,theta,phi,wdprod,ghatwnorm;
	dcomplex result(0.,0.),tmpcomp1,tmpcomp2,tmpcomp3;
	threevec<double> wvec,ghatw,wpar,wperp;
	
	for(int z=-MAXRUN; z<=MAXRUN; z++){
		for(int y=-MAXRUN; y<=MAXRUN; y++){
			for(int x=-MAXRUN; x<=MAXRUN; x++){
				if( (x==0) && (y==0) && (z==0) ) continue; //exclude zero!
				wvec((double)x,(double)y,(double)z);
				
				//compute scalar product and orthogonal projection:
				if(!is_zeroboost){
					orthogonal_projection(wvec,boost,wpar,wperp);
					ghatw=gamma*wpar+wperp;
				}
				else{
					ghatw=wvec;
				}
				ghatwnorm=ghatw.norm();
				
				//solve the integral:
				integrand2.set(q2,ghatwnorm,l);
				Midpnt<TFunctor> int2(integrand2,0.,lambda);
				tmp=qromo(int2);
				if(l!=0) tmp*=pow(ghatwnorm,l);
				
				//compute the complex parts:
				ghatw.get_spherical_coordinates(r,theta,phi);
				tmpcomp1=spherical_harmonicy(l,m,theta,phi);
				//exip(-ipi w*d)-term
				if(!is_zeroboost){
					wdprod=wvec*boost;
					tmpcomp2(cos(pimath*wdprod),-sin(pimath*wdprod));
				}
				else tmpcomp2(1.,0.);
				tmpcomp3=tmpcomp1*tmpcomp2*tmp;
				result+=tmpcomp3;
			}
		}
	}
	result*=gamma*pow(pimath,(double)(1.5+l));
	
	return result;
}

double Zetafunc::term3improved(const double q2){
	double result=0.,integral;
	threevec<double> wvec;
	
	//double wrapper since lower integration boundary is a singularity:
	integrand2.set(q2,0.,0,is_improved);
	Midpnt<TFunctor> int2(integrand2,0.,lambda);
	integral=qromo(int2);
	
	//multiply integral with stuff:	
	result=pimath*integral/2.;
	
	return result;
}
//******************************************************************
//******************************************************************
//END Zeta-function on a torus
//******************************************************************
//******************************************************************