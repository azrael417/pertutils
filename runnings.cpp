#include "mathutils.hpp"
#include "pertutils.hpp"

namespace pertutils{
    
    void alpha::set_betavec(){
        
        std::vector<double> betatmp;
        betatmp.push_back(betazero());
        if(looporder>1){
            betatmp.push_back(betaone());
        }
        if(looporder>2){
            betatmp.push_back(betatwo());
        }
        if(looporder>3){
            betatmp.push_back(betathree());
        }
        
        betavec.assign(betatmp);
        betfunc.set(betavec);
    }
    
    alpha::alpha(double MUSTART, double ALPHASTART, double MUMIN, double MUMAX, int MUCOUNT, int NF, int LO, int id) : mustart(MUSTART), alphastart(ALPHASTART), mumin(MUMIN), mumax(MUMAX), mucount(MUCOUNT), numflavours(NF), looporder(LO), intid(id), stepmin(0.0), integrated(false), interpolated(false){
        
        atol=1e-8;
        rtol=atol;
        
        double astart;
        rhomin=::std::log(::std::pow(mumin,2.));
        rhomax=::std::log(::std::pow(mumax,2.));
        if( std::abs(mustart)>atol && std::abs(alphastart)>atol ){
            set_betavec();
            perform_integration(alphastart/pimath,::std::log(::std::pow(mustart,2.)),rhomin,mucount);
            astart=yvec[yvec.size()-1]/pimath;
        }
        else{
            switch(numflavours){
                case 3: mustart=0.99*mc_msbar;
                    break;
                case 4: mustart=0.99*mb_msbar;
                    break;
                case 5: mustart=mz;
                    break;
                default: mustart=2.;
                    break;
            }
            rhostart=::std::log(::std::pow(mustart,2.));
            
            //obtaining starting values by running from mz:
            if(numflavours!=0){
                alphastart=alpha_step_from_mz(rhostart);
            }
            else{
                alphastart=alpha_from_lambda(rhostart);
            }
            astart=alphastart/pimath;
            set_betavec();
            
            //evolve from rhostart to rhomin:
            perform_integration(astart,rhostart,rhomin,mucount);
            astart=yvec[yvec.size()-1]/pimath;
        }
        ::std::cout << "alpha::alpha: starting value alpha( " << mumin << "^2 GeV^2 ) = " << astart*pimath << ::std::endl;
        
        //integration
        perform_integration(astart,rhomin,rhomax,mucount);
        
        //interpolation
        perform_interpolation();
    }
    
    alpha::alpha(double MUMIN, double MUMAX, int MUCOUNT, int NF, int LO, int id) : mumin(MUMIN), mumax(MUMAX), mucount(MUCOUNT), numflavours(NF), looporder(LO), intid(id), stepmin(0.0), integrated(false), interpolated(false){
        
        atol=1e-8;
        rtol=atol;
        
        double astart;
        rhomin=::std::log(::std::pow(mumin,2.));
        rhomax=::std::log(::std::pow(mumax,2.));
        switch(numflavours){
            case 3: mustart=0.99*mc_msbar;
                break;
            case 4: mustart=0.99*mb_msbar;
                break;
            case 5: mustart=mz;
                break;
            default: mustart=2.;
                break;
        }
        rhostart=::std::log(::std::pow(mustart,2.));
        
        //obtaining starting values by running from mz:
        if(numflavours!=0){
            alphastart=alpha_step_from_mz(rhostart);
        }
        else{
            alphastart=alpha_from_lambda(rhostart);
        }
        astart=alphastart/pimath;
        set_betavec();
        
        //evolve from rhostart to rhomin:
        perform_integration(astart,rhostart,rhomin,mucount);
        astart=yvec[yvec.size()-1]/pimath;
        ::std::cout << "alpha::alpha: starting value alpha( " << mumin << "^2 GeV^2 ) = " << astart*pimath << ::std::endl;
        
        //integration:
        perform_integration(astart,rhomin,rhomax,mucount);
        
        //interpolation
        perform_interpolation();
    }
    
    void alpha::perform_integration(double astart, double rmin, double rmax, int stepcount){
        Vector<double> avector(1);
        double rstep=(rmax-rmin)/(static_cast<double>(stepcount));
        avector[0]=astart;
        Output out(stepcount);
        
        if(std::abs(rmin-rmax)<atol){
            xvec.resize(1);
            yvec.resize(1);
            xvec[0]=rmin;
            yvec[0]=astart*pimath;
        }
        
        else{
            if(intid==0){
                Odeint<StepperDopr853<betafunction> > doprintegrator(avector,rmin,rmax,atol,rtol,rstep,stepmin,out,betfunc);
                doprintegrator.integrate();
            }
            else{
                Odeint<StepperBS<betafunction> > bsintegrator(avector,rmin,rmax,atol,rtol,rstep*10.,stepmin,out,betfunc);
                bsintegrator.integrate();
            }
            
            xvec.resize(out.count);
            yvec.resize(out.count);
            for(int i=0; i<out.count; i++){
                xvec[i]=out.xsave[i];
                yvec[i]=out.ysave[0][i]*pimath;
            }
        }
        integrated=true;
    }
    
    void alpha::perform_interpolation(){
        if(integrated && yvec.size()>1){
            if(interpolated) delete interpolation;
            interpolation=new spline_interp(xvec,yvec);
            
            interpolated=true;
        }
    }
    
    void alpha::set_integrator(int id){
        intid=id;
        
        if(intid==0){
            ::std::cout << "alpha::set_integrator: integrator set to Dormand-Prince." << ::std::endl;
        }
        else{
            ::std::cout << "alpha::set_integrator: integrator set to Bulrisch-Stoer." << ::std::endl;
        }
        integrated=false;
    }
    
    void alpha::set_LO(int loopord){
        looporder=loopord;
        integrated=false;
        interpolated=false;
        
        //obtaining starting values by running from mz:
        if(numflavours!=0) alphastart=alpha_step_from_mz(rhomin);
        else alphastart=alpha_from_lambda(rhomin);
        double astart=alphastart/pimath;
        set_betavec();
        
        //integration:
        perform_integration(astart,rhomin,rhomax,mucount);
        
        //interpolation
        perform_interpolation();
    }
    
    void alpha::set_NF(int nf){
        numflavours=nf;
        integrated=false;
        interpolated=false;
        
        //obtaining starting values by running from mz:
        if(numflavours!=0) alphastart=alpha_step_from_mz(rhomin);
        else alphastart=alpha_from_lambda(rhomin);
        double astart=alphastart/pimath;
        set_betavec();
        
        //integration:
        perform_integration(astart,rhomin,rhomax,mucount);
        
        //interpolation
        perform_interpolation();
    }
    
    alpha::alpha(::std::string filename){
        ::std::ifstream input;
        double token;
        
        input.open(filename.c_str());
        if(!input.good()){
            ::std::cerr << "alpha: error, cannot open " << filename << "for reading!" << ::std::endl;
        }
        else{
            xvec.resize(mucount+1);
            yvec.resize(mucount+1);
            input >> mustart;
            input >> alphastart;
            input >> mumin;
            input >> mumax;
            input >> mucount;
            for(int i=0; i<=mucount; i++){
                input >> token;
                xvec[i]=token;
                input >> token;
                yvec[i]=token;
            }
            input.close();
            
            //do interpolation:
            integrated=true;
            perform_interpolation();
        }
    }
    
    //Print stored xvec,yvec into file:
    int alpha::writefile(::std::string filename){
        ::std::ofstream output;
        output.precision(10);
        output.setf(::std::ios::scientific,::std::ios::floatfield);
        
        output.open(filename.c_str());
        if(!output.good()){
            ::std::cerr << "alpha::writefile: error, cannot open " << filename << " for writing!" << ::std::endl;
            return EXIT_FAILURE;
        }
        else{
            if(integrated){
                output << mustart << ::std::endl;
                output << alphastart << ::std::endl;
                output << mumin << ::std::endl;
                output << mumax << ::std::endl;
                output << mucount << ::std::endl;
                for(unsigned int i=0; i<xvec.size(); i++){
                    output << xvec[i] << " " << yvec[i] << ::std::endl;
                }
            }
            else{
                ::std::cerr << "alpha::writefile: error, please do the integration first!" << ::std::endl;
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
        double result;
        double rho=::std::log(::std::pow(mu,2.));
        
        if(interpolated){
            result=interpolation->interp(rho);
        }
        else if(yvec.size()==1){
            if(std::abs(rho-xvec[0])<atol) result=yvec[0];
            else result=-1.;
        }
        else{
            ::std::cerr << "alpha::operator(): please perform interpolation first!" << ::std::endl;
            result=-1.;
        }
        return result;
    }
    
    double alpha::get(const double mu){
        double result;
        double rho=::std::log(::std::pow(mu,2.));
        
        if(interpolated){
            result=interpolation->interp(rho);
        }
        else if(yvec.size()==1){
            if(std::abs(rho-xvec[0])<atol) result=yvec[0];
            else result=-1.;
        }
        else{
            ::std::cerr << "alpha::get: please perform interpolation first!" << ::std::endl;
            result=-1.;
        }
        return result;
    }
    
    alpha::~alpha(){
        xvec.clear();
        yvec.clear();
    }
    
    double alpha::alpha_from_lambda(double rmin){
        int numflavoursbackup=numflavours;
        numflavours=0;
        set_betavec();
        double L=rmin-::std::log(LQCD*LQCD);
        
        double result=1./(betavec[0]*L)*(1.-betavec[1]*::std::log(L)/(::std::pow(betavec[0],2.)*L)+1./(::std::pow(betavec[0]*L,2.))*(::std::pow(betavec[1],2.)/::std::pow(betavec[0],2.)*(::std::pow(::std::log(L),2.)-::std::log(L)-1.)+betavec[2]/betavec[0]));
        
        numflavours=numflavoursbackup;
        set_betavec();
        
        return result*pimath;
    }
    
    double alpha::alpha_step_from_mz(double rmin){
        double rup,rdown;
        double astart=alphas_mz/pimath;
        int nflavourbackup=numflavours;
        double result;
        int nsteps;
        
        ::std::cout << "Running down from Mz" << ::std::endl;
        
        ::std::vector<double> thresholds;
        thresholds.push_back(::std::log(mz*mz));
        if(rmin< ::std::log(mb_msbar*mb_msbar)) thresholds.push_back(::std::log(mb_msbar*mb_msbar));
        if(rmin< ::std::log(mc_msbar*mc_msbar)) thresholds.push_back(::std::log(mc_msbar*mc_msbar));
        thresholds.push_back(0.);
        
        for(unsigned int i=0; i<(thresholds.size()-1); i++){
            
            //update flavour content:
            rup=thresholds[i];
            rdown=max(rmin,thresholds[i+1]);
            numflavours=5-i;
            set_betavec();
            
            //integrate:
            nsteps=static_cast<int>(lround(std::abs(rup-rdown)*1.e2));
            perform_integration(astart,rup,rdown,nsteps);
            
            //update starting values:
            astart=yvec[yvec.size()-1]/pimath;
            ::std::cout << "nf= " << numflavours << ", alpha(" << ::std::exp(rdown/2) << "^2 GeV^2)= " << astart*pimath << ::std::endl;
        }
        
        //reset stuff:
        result=astart*pimath;
        numflavours=nflavourbackup;
        set_betavec();
        
        return result;
    }
    
    //QCD Beta-function at 4-loop (hep-ph/9910332,p.20, 72):
    double alpha::betazero(){
        double nf=static_cast<double>(numflavours);
        return (11.-2.*nf/3.)/4.;
    }
    
    double alpha::betaone(){
        double nf=static_cast<double>(numflavours);
        return (102.-38.*nf/3.)/16.;
    }
    
    double alpha::betatwo(){
        double nf=static_cast<double>(numflavours);
        return (2857./2.-5033.*nf/18.+325.*nf*nf/54.)/64.;
    }
    
    double alpha::betathree(){
        double nf=static_cast<double>(numflavours);
        return (149753./6.+3564.*zeta3-(1078361./162.+6508.*zeta3/27.)*nf+(50065./162.+6472.*zeta3/81.)*nf*nf+1093./729.*nf*nf*nf)/256.;
    }
    
    double alpha::betafunc(double a){
        double result=0.0;
        for(unsigned int i=0; i<betavec.size(); i++){
            result-=betavec[i]*::std::pow(a,static_cast<double>(i+2));
        }
        return result;
    }
    
    double alpha::pdbetafunc(double a){
        double result=0.0;
        for(unsigned int i=0; i<betavec.size(); i++){
            result-=static_cast<double>(i+2)*betavec[i]*::std::pow(a,static_cast<double>(i+1));
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
            result-=static_cast<double>((i+2)*(i+1))*betavec[i]*::std::pow(a,static_cast<double>(i));
        }
        return result;
    }
    
    double alpha::d2betafunc(double a){
        double result=pdbetafunc(a)*dbetafunc(a);
        result+=::std::pow(betafunc(a),2.)*pd2betafunc(a);
        return result;
    }
    
}
