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

alpha::alpha(double MUSTART, double ALPHASTART, double MUMIN, double MUMAX, int MUCOUNT, int NF, int LO, int id) : mustart(MUSTART), alphastart(ALPHASTART), mumin(MUMIN), mumax(MUMAX), mucount(MUCOUNT), numflavours(NF), looporder(LO), intid(id), stepmin(0.0), integrated(false), interpolated(false){
    
    atol=1e-6;
    rtol=atol;
    
    rhomin=log(pow(mumin,2));
    rhomax=log(pow(mumax,2));
    set_betavec();
    perform_integration(alphastart/pimath,log(pow(mustart,2)),rhomin,mucount);
	double astart=yvec[yvec.size()-1]/pimath;
    
    //integration
    perform_integration(astart,rhomin,rhomax,mucount);
    
    //interpolation
    perform_interpolation();
}

alpha::alpha(double MUMIN, double MUMAX, int MUCOUNT, int NF, int LO, int id) : mumin(MUMIN), mumax(MUMAX), mucount(MUCOUNT), numflavours(NF), looporder(LO), intid(id), stepmin(0.0), integrated(false), interpolated(false){
	
    atol=1e-6;
    rtol=atol;

	mustart=mumin;
    rhomin=log(pow(mumin,2));
	rhomax=log(pow(mumax,2));
    
    //obtaining starting values by running from mz:
	alphastart=alpha_step_from_mz(rhomin);
    double astart=alphastart/pimath;
    set_betavec();
    
    //integration:
    perform_integration(astart,rhomin,rhomax,mucount);
    
    //interpolation
    perform_interpolation();
}

void alpha::perform_integration(double astart, double rmin, double rmax, int stepcount){
    std::vector<double> avector(1);
    double rstep=fabs(rmax-rmin)/((double)stepcount);
    avector[0]=astart;
    Output out(stepcount);
    
    if(intid==0){
        Odeint<StepperDopr853<betafunction> > doprintegrator(avector,rmin,rmax,atol,rtol,rstep,stepmin,out,betfunc);
        doprintegrator.integrate();
    }
    else{
        Odeint<StepperBS<betafunction> > bsintegrator(avector,rmin,rmax,atol,rtol,rstep,stepmin,out,betfunc);
        bsintegrator.integrate();
    }
    
    xvec.clear();
    yvec.clear();
    for(int i=0; i<(out.count); i++){
        xvec.push_back((out.xsave[i]));
        yvec.push_back((out.ysave[0][i])*pimath);
    }
    integrated=true;
}

void alpha::perform_interpolation(){
    if(integrated){
        if(interpolated) delete interpolation;
        interpolation=new spline_interp(xvec,yvec);
        
        interpolated=true;
    }
}

void alpha::set_integrator(int id){
    intid=id;
    
    if(intid==0){
        std::cout << "alpha::set_integrator: integrator set to Dormand-Prince." << std::endl;
    }
    else{
        std::cout << "alpha::set_integrator: integrator set to Bulrisch-Stoer." << std::endl;
    }
    integrated=false;
}

void alpha::set_LO(int loopord){
	looporder=loopord;
	integrated=false;
    interpolated=false;
    
    //obtaining starting values by running from mz:
	alphastart=alpha_step_from_mz(rhomin);
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
	alphastart=alpha_step_from_mz(rhomin);
    double astart=alphastart/pimath;
    set_betavec();
    
    //integration:
    perform_integration(astart,rhomin,rhomax,mucount);
    
    //interpolation
    perform_interpolation();
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
        
        //do interpolation:
        integrated=true;
        perform_interpolation();
	}
}

//Print stored xvec,yvec into file:
int alpha::writefile(std::string filename){
	std::ofstream output;
	output.precision(10);
    output.setf(std::ios::scientific,std::ios::floatfield);
	
	output.open(filename.c_str());
	if(!output.good()){
		std::cerr << "alpha::writefile: error, cannot open " << filename << " for writing!" << std::endl;
		return EXIT_FAILURE;
	}
	else{
        if(integrated){
            output << mustart << std::endl;
            output << alphastart << std::endl;
            output << mumin << std::endl;
            output << mumax << std::endl;
            output << mucount << std::endl;
            for(unsigned int i=0; i<xvec.size(); i++){
                output << xvec[i] << " " << yvec[i] << std::endl;
            }
        }
        else{
            std::cerr << "alpha::writefile: error, please do the integration first!" << std::endl;
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
    double rho, result;
    if(interpolated){
        rho=log(pow(mu,2));
        result=interpolation->interp(rho);
    }
    else{
        std::cerr << "alpha::operator(): please perform interpolation first!" << std::endl;
        result=-1.;
    }
	return result;
}

double alpha::get(const double mu){
    double rho, result;
    if(interpolated){
        rho=log(pow(mu,2));
        result=interpolation->interp(rho);
    }
    else{
        std::cerr << "alpha::get: please perform interpolation first!" << std::endl;
        result=-1.;
    }
	return result;
}

alpha::~alpha(){
	xvec.clear();
	yvec.clear();
}

double alpha::alpha_step_from_mz(double rmin){
    double rhost=log(pow(mz,2)), rm;
	double astart=alphas_mz/pimath;
	int nflavourbackup=numflavours;
    double result;
    
    std::vector<double> thresholds;
    if(rmin>=log(mb_msbar*mb_msbar)){
        set_betavec();
        
        //integration:
        perform_integration(astart,rhost,rmin,1000);
        result=yvec[yvec.size()-1];
    }
    else{
        thresholds.push_back(log(mb_msbar*mb_msbar));
        if(rmin<log(mc_msbar*mc_msbar)) thresholds.push_back(log(mc_msbar*mc_msbar));
        thresholds.push_back(0.);
        
        rm=rmin;
        for(unsigned int i=0; i<(thresholds.size()-1); i++){
            
            //update flavour content:
            numflavours=5-i;
            set_betavec();
            
            //integrate:
            perform_integration(astart,rhost,rm,1000);
            
            //update starting values:
            astart=yvec[yvec.size()-1]/pimath;
            rhost=thresholds[i];
            rm=max(rmin,thresholds[i+1]);
        }
        result=astart*pimath;
    }
    numflavours=nflavourbackup;
    set_betavec();
    
    return result;
}

/*double alpha::alpha_step_from_mz(double rmin){
	double rhost=log(pow(mz,2));
	double astart=alphas_mz/pimath;
	int nflavour=numflavours;
	

	
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
}*/

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

/*double alpha::afunc(double astart, double sstart, double sziel, int numsteps){ //liefert a(log(mu^2)) bei gegebenem a(log(mu_0^2)) und log(mu_0^2).                                                                                                  
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
}*/

