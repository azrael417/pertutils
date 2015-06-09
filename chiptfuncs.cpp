#include "mathutils.hpp"
#include "pertutils.hpp"

namespace pertutils{
    
    //******************************************************************
    //******************************************************************
    //START chiPT FV functions
    //******************************************************************
    //******************************************************************
    Xifunc::Xifunc(double MLSQMIN, double MLSQMAX, int MLCOUNT, double ss) : MLsqmin(MLSQMIN), MLsqmax(MLSQMAX), MLcount(MLCOUNT), s(ss) {
        MLstep=(MLsqmax-MLsqmin)/(static_cast<double>(MLcount));
        
        double MLsq, dummy;
        std::vector<double> xvectmp, yvectmp;
        for(int i=0; i<=MLcount; i++){
            MLsq=MLsqmin+static_cast<double>(i)*MLstep;
            xvectmp.push_back(MLsq);
            dummy=Xifuncnorm(MLsq,s);
            yvectmp.push_back(dummy);
        }
        xvec.assign(xvectmp);
        yvec.assign(yvectmp);
        
        interpolation=new spline_interp(xvec,yvec);
    }
    
    Xifunc::Xifunc(::std::string filename){
        ::std::ifstream input;
        double token;
        
        input.open(filename.c_str());
        if(!input.good()){
            ::std::cerr << "Xifunc: error, cannot open " << filename << "for reading!" << ::std::endl;
        }
        else{
            std::vector<double> xvectmp, yvectmp;
            yvec.clear();
            input >> MLsqmin;
            input >> MLsqmax;
            input >> MLcount;
            input >> s;
            for(int i=0; i<=MLcount; i++){
                input >> token;
                xvectmp.push_back(token);
                input >> token;
                yvectmp.push_back(token);
            }
            input.close();
            xvec.assign(xvectmp);
            yvec.assign(yvectmp);
            
            interpolation=new spline_interp(xvec,yvec);
        }
    }
    
    //Print stored xvec,yvec into file:
    int Xifunc::writefile(::std::string filename){
        ::std::ofstream output;
        output.precision(10);
        output.setf(::std::ios::scientific,::std::ios::floatfield);
        
        output.open(filename.c_str());
        if(!output.good()){
            ::std::cerr << "Xifunc: error, cannot open " << filename << " for writing!" << ::std::endl;
            return EXIT_FAILURE;
        }
        else{
            output << MLsqmin << ::std::endl;
            output << MLsqmax << ::std::endl;
            output << MLcount << ::std::endl;
            output << s << ::std::endl;;
            for(unsigned int i=0; i<xvec.size(); i++){
                output << xvec[i] << " " << yvec[i] << ::std::endl;
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
        double result=::std::pow(L,(2.*s-3.))*(interpolation->interp(Msq*L*L));
        return result;
    }
    
    double Xifunc::get(const double Msq, const double L){
        double result=::std::pow(L,(2.*s-3.))*(interpolation->interp(Msq*L*L));
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
        tau=::std::exp(t-::std::exp(-t));
        dtaudt=tau*(1.+::std::exp(-t));
        result=Xifuncintegrand(tau,MLsq,s)*dtaudt;
        
        //Print dXI(ML,s)/dML
        if(flag !=0 ){
            result*=((-1.)*tau*MLsq/(4.*pimath*pimath));
        }
        return result;
    }
    
    double Xifuncintegrand(double tau, double MLsq,  double s){
        double result;
        
        result=::std::pow(EllipticTheta(tau),3)-::std::pow(pimath/tau,3./2.);
        result*=::std::pow(tau,(s-1.))*::std::exp(-tau*MLsq/::std::pow(2.*pimath,2.));
        
        return result;
    }
    
    double Xifuncnorm(double MLsq, double s, int flag){
        double result=0.;
        Xifuncintegrandtransformed Xif(MLsq,s,flag);
        
        result=qtrap(Xif, 0.0, 6.);
        result/=(::std::pow(2.*pimath,(2.*s))*gammfunc(s));
        
        return result;
    }
    
    //Stephan Dürrs functions:
    double Ikfunc2(int n, double mpiL){
        Bessik bessel;
        return -1.5*bessel.k1(::std::sqrt(static_cast<double>(n))*mpiL);
    }
    
    double Ipifunc2(int n, double mpiL){
        Bessik bessel;
        return -4.*bessel.k1(::std::sqrt(static_cast<double>(n))*mpiL);
    }
    
    double mpiLovermpi(double mpi, double fpi, double L){
        double eps=1e-10;
        int nmax=250;
        int nup=5;
        double beta=1.;
        Levin acc(nmax,eps);
        
        double sn=1.;
        for(int n=1; n<=nup; n++){
            sn+=(-1.)*static_cast<double>(multiplicities[n])*mpi/(32.*fpi*fpi*pimath*pimath*::std::sqrt(static_cast<double>(n))*L)*Ipifunc2(n,mpi*L);
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
                an=(-1.)*static_cast<double>(multiplicities[nrun])*mpi/(32.*fpi*fpi*pimath*pimath*::std::sqrt(static_cast<double>(nrun))*L)*Ipifunc2(nrun,mpi*L);
                sn+=an;
                omegan=(beta+static_cast<double>(nrun))*an;
                result=acc.next(sn,omegan);
                nrun++;
            }
        }
        
        //Switching to direct resummation if accelleration fails:
        if(isnan(result)){
            result=1.;
            for(int n=1; n<=nmax; n++){
                result+=(-1.)*static_cast<double>(multiplicities[n])*mpi/(32.*fpi*fpi*pimath*pimath*::std::sqrt(static_cast<double>(n))*L)*Ipifunc2(n,mpi*L);
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
            sn+=(-1.)*static_cast<double>(multiplicities[n])*mpi*mpi/(32.*mk*fpi*fpi*pimath*pimath*::std::sqrt(static_cast<double>(n))*L)*Ikfunc2(n,mpi*L);
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
                an=(-1.)*static_cast<double>(multiplicities[nrun])*mpi*mpi/(32.*mk*fpi*fpi*pimath*pimath*::std::sqrt(static_cast<double>(nrun))*L)*Ikfunc2(nrun,mpi*L);
                sn+=an;
                omegan=(beta+static_cast<double>(nrun))*an;
                result=acc.next(sn,omegan);
                nrun++;
            }
        }
        
        //Switching to direct resummation if accelleration fails:
        if(isnan(result)){
            result=1.;
            for(int n=1; n<=nmax; n++){
                result+=(-1.)*static_cast<double>(multiplicities[n])*mpi*mpi/(32.*mk*fpi*fpi*pimath*pimath*::std::sqrt(static_cast<double>(n))*L)*Ikfunc2(n,mpi*L);
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
            sn+=static_cast<double>(multiplicities[n])*mpi/(16*pimath*pimath*fk*fpi*::std::sqrt(static_cast<double>(n))*L)*Ikfunc2(n,mpi*L);
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
                an=static_cast<double>(multiplicities[nrun])*mpi/(16*pimath*pimath*fk*fpi*::std::sqrt(static_cast<double>(nrun))*L)*Ikfunc2(nrun,mpi*L);
                sn+=an;
                omegan=(beta+static_cast<double>(nrun))*an;
                result=acc.next(sn,omegan);
                nrun++;
            }
        }
        
        //Switching to direct resummation if accelleration fails:
        if(isnan(result)){
            result=1.;
            for(int n=1; n<=nmax; n++){
                result+=static_cast<double>(multiplicities[n])*mpi/(16*pimath*pimath*fk*fpi*::std::sqrt(static_cast<double>(n))*L)*Ikfunc2(n,mpi*L);
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
            sn+=static_cast<double>(multiplicities[n])*mpi/(16*pimath*pimath*fpi*fpi*::std::sqrt(static_cast<double>(n))*L)*Ipifunc2(n,mpi*L);
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
                an=static_cast<double>(multiplicities[nrun])*mpi/(16*pimath*pimath*fpi*fpi*::std::sqrt(static_cast<double>(nrun))*L)*Ipifunc2(nrun,mpi*L);
                sn+=an;
                omegan=(beta+static_cast<double>(nrun))*an;
                result=acc.next(sn,omegan);
                nrun++;
            }
        }
        
        //Switching to direct resummation if accelleration fails:
        if(isnan(result)){
            result=1.;
            for(int n=1; n<=nmax; n++){
                result+=static_cast<double>(multiplicities[n])*mpi/(16*pimath*pimath*fpi*fpi*::std::sqrt(static_cast<double>(n))*L)*Ipifunc2(n,mpi*L);
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
        double tmp;
        //special treatment since sums are real:
        if(is_improved && is_zeroboost){
            tmp=term1spherical(q2);
            term1val(tmp,0.);
            tmp=term3spherical(q2);
            term3val(tmp,0.);
        }
        else if(is_zeroboost){
            term1val=term1noboost(q2);
            term3val=term3noboost(q2);
        }
        else{
            term1val=term1full(q2);
            term3val=term3full(q2);
        }
        tmp=term2full(q2);
        term2val(tmp,0.);
        
        //result:
        result=term1val+term2val+term3val;
        
        return result;
    }
    
    void Zetafunc::print_contributions(){
        ::std::cout << "Printing contributions from individual terms!" << ::std::endl;
        ::std::cout << "term1= " << term1val << ::std::endl;
        ::std::cout << "term2= " << term2val << ::std::endl;
        ::std::cout << "term3= " << term3val << ::std::endl;
    }
    
    //*************************************************TERM 1*************************************************
    //brute force sum
    dcomplex Zetafunc::term1full(const double q2){
        dcomplex result(0.,0.),tmpcomp;
        double fact,rsq,theta,phi;
        threevec<double> nvec,npar,nperp,rvec;
            
        double resultre=0.,resultim=0.;
#pragma omp parallel for reduction(+:resultre) reduction(+:resultim) private(rvec,rsq,theta,phi,nvec,fact,tmpcomp,npar,nperp)
        for(int z=-MAXRUN; z<=MAXRUN; z++){
            for(int y=-MAXRUN; y<=MAXRUN; y++){
                for(int x=-MAXRUN; x<=MAXRUN; x++){
                    nvec(static_cast<double>(x),static_cast<double>(y),static_cast<double>(z));
                    
                    //compute r:
                    if(!is_zeroboost){
                        orthogonal_projection(nvec,boost,npar,nperp);
                        rvec=npar/gamma-boost/(2.*gamma)+nperp;
                    }
                    else{
                        rvec=nvec;
                    }
                    rvec.get_spherical_coordinates(rsq,theta,phi);
                    if(l!=0)fact=::std::pow(rsq,l); //compute |r|^l
                    else fact=1.;
                    rsq*=rsq; //compute r^2
                    
                    //compute prefact:
                    fact*=::std::exp(-lambda*(rsq-q2))/(rsq-q2); //compute the ratio
                    tmpcomp=fact*spherical_harmonicy(l,m,theta,phi); //multiply with spherical harmonics
                    resultre+=tmpcomp.re();
                    resultim+=tmpcomp.im();
                    //result+=tmpcomp;
                }
            }
        }
        result=dcomplex(resultre,resultim);
        
        return result;
    }
    
    //can be used if boost is zero:
    dcomplex Zetafunc::term1noboost(const double q2){
        dcomplex result(0.,0.),sphfact,tmpcomp;
        double fact,rsq,r,theta,phi;
        threevec<double> nvec;
        
        //zero term is zero:
        result=dcomplex(0.,0.);
        
        //line terms
        //x>0, y=z=0 and y>0, x=z=0 and z>0, x=y=0:
        //note that z->-z is equivalent to theta->pi-theta, similar to this: x->-x yields phi->pi-phi and y->-y is phi->2pi-phi
        sphfact=    spherical_harmonicy(l,m,pimath/2.   ,   0.);                                //x>0
        sphfact+=   spherical_harmonicy(l,m,pimath/2.   ,   pimath);                            //x<0
        sphfact+=   spherical_harmonicy(l,m,pimath/2.   ,   pimath/2.);                         //y>0
        sphfact+=   spherical_harmonicy(l,m,pimath/2.   ,   3.*pimath/2.);                      //y<0
        sphfact+=   spherical_harmonicy(l,m,0.          ,   0.);                                //z>0
        sphfact+=   spherical_harmonicy(l,m,pimath      ,   0.);                                //z<0
        for(int x=1; x<=MAXRUN; x++){
            rsq=x*x;
            result+=dcomplex(::std::exp(-lambda*(rsq-q2))*::std::pow(static_cast<double>(x),l)/(rsq-q2),0.)*sphfact;
        }
        
        
        //plane terms
        //x,y!>0, z=0:
        //the four ylm account for the possibilities +/+, +/-, -/+, -/-:
        for(int y=1; y<=MAXRUN; y++){
            for(int x=1; x<=MAXRUN; x++){
                threevec<int>(x,y,0).get_spherical_coordinates(r,theta,phi);
                rsq=r*r;
                //z=0:
                sphfact=    spherical_harmonicy(l,m,theta   ,   phi);                           //x,y>0
                sphfact+=   spherical_harmonicy(l,m,theta   ,   pimath-phi);                    //x<0,y>0
                sphfact+=   spherical_harmonicy(l,m,theta   ,   pimath+phi);                    //x,y<0
                sphfact+=   spherical_harmonicy(l,m,theta   ,   2.*pimath-phi);                 //x>0,y<0
                
                //result
                result+=dcomplex(::std::exp(-lambda*(rsq-q2))*::std::pow(r,l)/(rsq-q2),0.)*sphfact;
            }
        }
        
        //x,z>0, y=0 and y,z>0, x=0:
        for(int z=1; z<=MAXRUN; z++){
            for(int x=1; x<=MAXRUN; x++){
                threevec<int>(x,0,z).get_spherical_coordinates(r,theta,phi);
                rsq=r*r;
                
                //y=0:
                sphfact=    spherical_harmonicy(l,m,theta          ,    0);                    //x,z>0
                sphfact+=   spherical_harmonicy(l,m,pimath-theta   ,    0);                    //x>0,z<0
                sphfact+=   spherical_harmonicy(l,m,theta          ,    pimath);               //x<0,z>0
                sphfact+=   spherical_harmonicy(l,m,pimath-theta   ,    pimath);               //x,z<0
                //x=0:
                sphfact+=   spherical_harmonicy(l,m,theta          ,    pimath/2.);            //y,z>0
                sphfact+=   spherical_harmonicy(l,m,pimath-theta   ,    pimath/2.);            //y>0,z<0
                sphfact+=   spherical_harmonicy(l,m,theta          ,    3.*pimath/2.);         //y<0,z>0
                sphfact+=   spherical_harmonicy(l,m,pimath-theta   ,    3.*pimath/2.);         //y,z<0
                
                //result:
                result+=dcomplex(::std::exp(-lambda*(rsq-q2))*::std::pow(r,l)/(rsq-q2),0.)*sphfact;
            }
        }
        
        
        //cubic terms
        //x,y,z>0
        double resultre=0.,resultim=0.;
#pragma omp parallel for reduction(+:resultre) reduction(+:resultim) private(nvec,rsq,r,theta,phi,fact,tmpcomp)
        for(int z=1; z<=MAXRUN; z++){
            for(int y=1; y<=MAXRUN; y++){
                for(int x=1; x<=MAXRUN; x++){
                    nvec(static_cast<double>(x),static_cast<double>(y),static_cast<double>(z));
                    nvec.get_spherical_coordinates(r,theta,phi);
                    fact=::std::pow(r,l); //compute |r|^l
                    rsq=r*r; //compute r^2
                    fact*=::std::exp(-lambda*(rsq-q2))/(rsq-q2); //compute the ratio
                    
                    //compute sphfact: account for all possible orientations
                    sphfact=dcomplex(0.,0.);
                    for(unsigned int thetaa=0; thetaa<2; thetaa++){
                        sphfact+=spherical_harmonicy(l,m,static_cast<double>(thetaa)*pimath-theta,phi);            //x,y>0
                        sphfact+=spherical_harmonicy(l,m,static_cast<double>(thetaa)*pimath-theta,pimath-phi);     //x<0,y>0
                        sphfact+=spherical_harmonicy(l,m,static_cast<double>(thetaa)*pimath-theta,pimath+phi);     //x,y<0
                        sphfact+=spherical_harmonicy(l,m,static_cast<double>(thetaa)*pimath-theta,2.*pimath-phi);  //x>0,y<0
                    }

                    //compute the result:
                    tmpcomp=fact*sphfact; //ratio * ylm factor
                    resultre+=tmpcomp.re();
                    resultim+=tmpcomp.im();
                }
            }
        }
        result+=dcomplex(resultre,resultim);
        
        return result;
    }

    double Zetafunc::term1spherical(const double q2){
        double result=0.,rsq;
        threevec<double> nvec;
        
        //zero term:
        result=-::std::exp(lambda*q2)/q2;
        
        //x!=0, y=z=0, 3 times for different combinations (replace x e.g. with y) and factor 2 for +/-: makes an overall factor of 6:
        for(int x=1; x<=MAXRUN; x++){
            rsq=x*x;
            result+=6.*::std::exp(-lambda*(rsq-q2))/(rsq-q2);
        }
        
        //x,y!=0,z=0, 3 times for different combinations and factor 4 for +/-:
        for(int x=1; x<=MAXRUN; x++){
            for(int y=1; y<=MAXRUN; y++){
                rsq=x*x+y*y;
                result+=12.*::std::exp(-lambda*(rsq-q2))/(rsq-q2);
            }
        }
        
        //x,y,z!=0 and factor 8 for +/-:
        for(int x=1; x<=MAXRUN; x++){
            for(int y=1; y<=MAXRUN; y++){
                for(int z=1; z<=MAXRUN; z++){
                    rsq=x*x+y*y+z*z;
                    result+=8.*::std::exp(-lambda*(rsq-q2))/(rsq-q2); //compute the ratio
                }
            }
        }
        result/=::std::sqrt(4.*pimath);
        
        return result;
    }
    
    //*************************************************TERM 2*************************************************
    //this term is proportional to Dawson's Integral: term2=2/q *exp(lambda*q^2)*F(sqrt(lambda)q)
    //where F(x)=exp(-x^2) int_0^x ds exp(s^2) for q^2>0. For q^2<0, one can express is as the error Function: sqrt(pi)/|q| Erf(sqrt(lambda*|q|^2)):
    double Zetafunc::term2full(const double q2){
        double result=0.,tmp;
        
        if(l==0 && q2>=0.){
            tmp=4.*::std::sqrt(q2)*::std::exp(lambda*q2)*dawson(::std::sqrt(lambda*q2));
            tmp-=2./::std::sqrt(lambda)*::std::exp(lambda*q2);
            result=gamma*pimath/2.*tmp;
        }
        else if(l==0){
            Erf erf;
            tmp=-2.*::std::sqrt(std::abs(q2))*::std::sqrt(pimath)*erf.erf(::std::sqrt(lambda*std::abs(q2)));
            tmp-=2.*::std::exp(lambda*q2)/::std::sqrt(lambda);
            result=gamma*pimath/2.*tmp;
        }
        
        return result;
    }
    
    //*************************************************TERM 3*************************************************
    dcomplex Zetafunc::term3full(const double q2){
        double resultre=0.,resultim=0.;
#pragma omp parallel
        {
            double tmp,r,theta,phi,wdprod,ghatwnorm;
            dcomplex result(0.,0.),tmpcomp1,tmpcomp2,tmpcomp3;
            threevec<double> wvec,ghatw,wpar,wperp;
            
            Zetafuncint integrand2;
#pragma omp parallel for reduction(+:resultre) reduction(+:resultim)
            for(int z=-MAXRUN; z<=MAXRUN; z++){
                for(int y=-MAXRUN; y<=MAXRUN; y++){
                    for(int x=-MAXRUN; x<=MAXRUN; x++){
                        if( (x==0) && (y==0) && (z==0) ) continue; //exclude zero!
                        wvec(static_cast<double>(x),static_cast<double>(y),static_cast<double>(z));
                        
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
                        if(l!=0) tmp*=::std::pow(ghatwnorm,static_cast<double>(l));
                        
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
                        //result+=tmpcomp3;
                        resultre+=tmpcomp3.re();
                        resultim+=tmpcomp3.im();
                    }
                }
            }
        }
        dcomplex result(resultre,resultim);
        result*=gamma*::std::pow(pimath,static_cast<double>(1.5+l));
        
        return result;
    }
    
    dcomplex Zetafunc::term3noboost(const double q2){
        dcomplex result(0.,0.),sphfact,tmpcomp;
        double tmp,r,theta,phi;
        
        //zero term is zero:
        result=dcomplex(0.,0.);
        
        //line terms
        //x>0, y=z=0 and y>0, x=z=0 and z>0, x=y=0:
        //note that z->-z is equivalent to theta->pi-theta, similar to this: x->-x yields phi->pi-phi and y->-y is phi->2pi-phi
        sphfact=    spherical_harmonicy(l,m,pimath/2.   ,   0.);                                //x>0
        sphfact+=   spherical_harmonicy(l,m,pimath/2.   ,   pimath);                            //x<0
        sphfact+=   spherical_harmonicy(l,m,pimath/2.   ,   pimath/2.);                         //y>0
        sphfact+=   spherical_harmonicy(l,m,pimath/2.   ,   3.*pimath/2.);                      //y<0
        sphfact+=   spherical_harmonicy(l,m,0.          ,   0.);                                //z>0
        sphfact+=   spherical_harmonicy(l,m,pimath      ,   0.);                                //z<0
        
        double resultre=0., resultim=0.;
#pragma omp parallel
        {
            Zetafuncint integrand2;
#pragma omp parallel for reduction(+:resultre) private(tmp)
            for(int x=1; x<=MAXRUN; x++){
                integrand2.set(q2,x,l);
                Midpnt<TFunctor> int2(integrand2,0.,lambda);
                tmp=qromo(int2);
                if(l!=0) tmp*=::std::pow(x,static_cast<double>(l));
                resultre+=tmp;
            }
        }
        result+=resultre*sphfact;
        
        
        //plane-terms
        resultre=0.;
        resultim=0.;
        double* integrals=new double[MAXRUN*MAXRUN];
#pragma omp parallel
        {
            Zetafuncint integrand2;
            //do integrals first:
#pragma omp parallel for shared(integrals) private(r,theta,phi,tmp)
            for(int y=1; y<=MAXRUN; y++){
                for(int x=y; x<=MAXRUN; x++){
                    threevec<int>(x,y,0).get_spherical_coordinates(r,theta,phi);
                    
                    //integral
                    integrand2.set(q2,r,l);
                    Midpnt<TFunctor> int2(integrand2,0.,lambda);
                    tmp=qromo(int2);
                    if(l!=0) tmp*=::std::pow(r,static_cast<double>(l));
                    
                    //store results: it is symmetrc in x and y:
                    integrals[x-1+MAXRUN*(y-1)]=tmp;
                    integrals[y-1+MAXRUN*(x-1)]=tmp;
                }
            }
            
            //angular parts
            //x,y!>0, z=0:
            //the four ylm account for the possibilities +/+, +/-, -/+, -/-:
#pragma omp parallel for reduction(+:resultre) reduction(+:resultim) shared(integrals) private(r,theta,phi,tmp,tmpcomp)
            for(int y=1; y<=MAXRUN; y++){
                for(int x=1; x<=MAXRUN; x++){
                    threevec<int>(x,y,0).get_spherical_coordinates(r,theta,phi);
                    
                    //z=0:
                    sphfact=    spherical_harmonicy(l,m,theta   ,   phi);                           //x,y>0
                    sphfact+=   spherical_harmonicy(l,m,theta   ,   pimath-phi);                    //x<0,y>0
                    sphfact+=   spherical_harmonicy(l,m,theta   ,   pimath+phi);                    //x,y<0
                    sphfact+=   spherical_harmonicy(l,m,theta   ,   2.*pimath-phi);                 //x>0,y<0
                    
                    //result
                    tmpcomp=integrals[y-1+MAXRUN*(x-1)]*sphfact;
                    resultre+=tmpcomp.re();
                    resultre+=tmpcomp.im();
                }
            }
        
            //x,z>0, y=0 and y,z>0, x=0:
#pragma omp parallel for reduction(+:resultre) reduction(+:resultim) shared(integrals) private(r,theta,phi,tmp,tmpcomp)
            for(int z=1; z<=MAXRUN; z++){
                for(int x=1; x<=MAXRUN; x++){
                    threevec<int>(x,0,z).get_spherical_coordinates(r,theta,phi);
                    
                    //y=0:
                    sphfact=    spherical_harmonicy(l,m,theta          ,    0);                    //x,z>0
                    sphfact+=   spherical_harmonicy(l,m,pimath-theta   ,    0);                    //x>0,z<0
                    sphfact+=   spherical_harmonicy(l,m,theta          ,    pimath);               //x<0,z>0
                    sphfact+=   spherical_harmonicy(l,m,pimath-theta   ,    pimath);               //x,z<0
                    //x=0:
                    sphfact+=   spherical_harmonicy(l,m,theta          ,    pimath/2.);            //y,z>0
                    sphfact+=   spherical_harmonicy(l,m,pimath-theta   ,    pimath/2.);            //y>0,z<0
                    sphfact+=   spherical_harmonicy(l,m,theta          ,    3.*pimath/2.);         //y<0,z>0
                    sphfact+=   spherical_harmonicy(l,m,pimath-theta   ,    3.*pimath/2.);         //y,z<0
                    
                    //result:
                    tmpcomp=integrals[z-1+MAXRUN*(x-1)]*sphfact;
                    resultre+=tmpcomp.re();
                    resultre+=tmpcomp.im();
                }
            }
        }
        delete [] integrals;
        result+=dcomplex(resultre,resultim);
        
        
        //cubic terms:
        resultre=0.;
        resultim=0.;
        integrals=new double[MAXRUN*MAXRUN*MAXRUN];
#pragma omp parallel
        {
            //do integrals first
            Zetafuncint integrand2;
#pragma omp parallel for shared(integrals) private(r,theta,phi,tmp)
            for(int z=1; z<=MAXRUN; z++){
                for(int y=z; y<=MAXRUN; y++){
                    for(int x=y; x<=MAXRUN; x++){
                        threevec<int>(x,y,z).get_spherical_coordinates(r,theta,phi);
                        
                        //integral
                        integrand2.set(q2,r,l);
                        Midpnt<TFunctor> int2(integrand2,0.,lambda);
                        tmp=qromo(int2);
                        if(l!=0) tmp*=::std::pow(r,static_cast<double>(l));
                        
                        //store results: it is symmetrc in x and y:
                        integrals[x-1+MAXRUN*((y-1)+MAXRUN*(z-1))]=tmp;
                        integrals[z-1+MAXRUN*((x-1)+MAXRUN*(y-1))]=tmp;
                        integrals[y-1+MAXRUN*((z-1)+MAXRUN*(x-1))]=tmp;
                    }
                }
            }

#pragma omp parallel for reduction(+:resultre) reduction(+:resultim) shared(integrals) private(r,theta,phi,tmpcomp)
            for(int z=1; z<=MAXRUN; z++){
                for(int y=1; y<=MAXRUN; y++){
                    for(int x=1; x<=MAXRUN; x++){
                        threevec<int>(x,y,z).get_spherical_coordinates(r,theta,phi);
                        
                        //compute sphfact: account for all possible orientations
                        sphfact=dcomplex(0.,0.);
                        for(unsigned int thetaa=0; thetaa<2; thetaa++){
                            sphfact+=spherical_harmonicy(l,m,static_cast<double>(thetaa)*pimath-theta,phi);            //x,y>0
                            sphfact+=spherical_harmonicy(l,m,static_cast<double>(thetaa)*pimath-theta,pimath-phi);     //x<0,y>0
                            sphfact+=spherical_harmonicy(l,m,static_cast<double>(thetaa)*pimath-theta,pimath+phi);     //x,y<0
                            sphfact+=spherical_harmonicy(l,m,static_cast<double>(thetaa)*pimath-theta,2.*pimath-phi);  //x>0,y<0
                        }
                        
                        //compute the result:
                        tmpcomp=integrals[x-1+MAXRUN*((y-1)+MAXRUN*(z-1))]*sphfact; //integral * ylm factor
                        resultre+=tmpcomp.re();
                        resultim+=tmpcomp.im();
                    }
                }
            }
        }
        delete [] integrals;
        result+=dcomplex(resultre,resultim);
        
        //multiply result by factor:
        result*=::std::pow(pimath,static_cast<double>(1.5+l));
        
        return result;
    }
    
    double Zetafunc::term3spherical(const double q2){
        double result=0.,integral;
        threevec<double> wvec;
        Zetafuncint integrand2;
        
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
    
}