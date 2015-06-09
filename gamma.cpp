//
//  gamma.cpp
//  pertutils
//
//  Created by Thorsten Kurth on 23.04.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

#include "mathutils.hpp"
#include "pertutils.hpp"


namespace pertutils{

    Matrix<dcomplex> get_gamma(const gamma_id gid){
        Matrix<dcomplex> result(4,4,0.);
        
        switch(gid){
            case ONE:
                result[0][0]=dcomplex(1.,0.);
                result[1][1]=dcomplex(1.,0.);
                result[2][2]=dcomplex(1.,0.);
                result[3][3]=dcomplex(1.,0.);
                break;
            case GX:
                result[0][3]=dcomplex(0,1.);
                result[1][2]=dcomplex(0,1.);
                result[2][1]=dcomplex(0,-1.);
                result[3][0]=dcomplex(0,-1.);
                break;
            case GY:
                result[0][3]=dcomplex(-1.,0.);
                result[1][2]=dcomplex(1.,0.);
                result[2][1]=dcomplex(1.,0.);
                result[3][0]=dcomplex(-1.,0.);
                break;
            case GZ:
                result[0][2]=dcomplex(0.,1.);
                result[1][3]=dcomplex(0,-1.);
                result[2][0]=dcomplex(0.,-1.);
                result[3][1]=dcomplex(0.,1.);
                break;
            case GT:
                result[0][2]=dcomplex(1.,0.);
                result[1][3]=dcomplex(1.,0.);
                result[2][0]=dcomplex(1.,0.);
                result[3][1]=dcomplex(1.,0.);
                break;
            case G5:
                result[0][0]=dcomplex(1.,0.);
                result[1][1]=dcomplex(1.,0.);
                result[2][2]=dcomplex(-1.,0.);
                result[3][3]=dcomplex(-1.,0.);
                break;
            case CG5:
                result[0][1]=dcomplex(-1.,0.);
                result[1][0]=dcomplex(1.,0.);
                result[2][3]=dcomplex(-1.,0.);
                result[3][2]=dcomplex(1.,0.);
                break;
            case C:
                result[0][1]=dcomplex(-1.,0.);
                result[1][0]=dcomplex(1.,0.);
                result[2][3]=dcomplex(1.,0.);
                result[3][2]=dcomplex(-1.,0.);
            case IG4CG5:
                result[0][3]=dcomplex(0,-1.);
                result[1][2]=dcomplex(0,1.);
                result[2][1]=dcomplex(0,-1.);
                result[3][0]=dcomplex(0,1.);
                break;
            case PARPLUS:
                result[0][0]=dcomplex(0.5,0.);
                result[1][1]=dcomplex(0.5,0.);
                result[2][2]=dcomplex(0.5,0.);
                result[3][3]=dcomplex(0.5,0.);
                
                result[0][2]=dcomplex(0.5,0.);
                result[1][3]=dcomplex(0.5,0.);
                result[2][0]=dcomplex(0.5,0.);
                result[3][1]=dcomplex(0.5,0.);
                break;
            case PARMINUS:
                result[0][0]=dcomplex(0.5,0.);
                result[1][1]=dcomplex(0.5,0.);
                result[2][2]=dcomplex(0.5,0.);
                result[3][3]=dcomplex(0.5,0.);
                
                result[0][2]=dcomplex(-0.5,0.);
                result[1][3]=dcomplex(-0.5,0.);
                result[2][0]=dcomplex(-0.5,0.);
                result[3][1]=dcomplex(-0.5,0.);
                break;
        }
        return result;
    }
    
    Matrix<double> get_gamma(const basak_gamma_id gid){
        Matrix<double> result(4,4,0.);
        
        switch(gid){
            case UTRANS:
                result[0][1]=1./sqrt2;
                result[0][3]=-1./sqrt2;
                result[1][0]=-1./sqrt2;
                result[1][2]=1./sqrt2;
                result[2][1]=1./sqrt2;
                result[2][3]=1./sqrt2;
                result[3][0]=-1./sqrt2;
                result[3][2]=-1./sqrt2;
                break;
            case G1G1:
                result[0][1]=1./(2.*sqrt2);
                result[0][3]=1./(2.*sqrt2);
                result[1][0]=-1./(2.*sqrt2);
                result[1][2]=-1./(2.*sqrt2);
                result[2][1]=1./(2.*sqrt2);
                result[2][3]=1./(2.*sqrt2);
                result[3][0]=-1./(2.*sqrt2);
                result[3][2]=-1./(2.*sqrt2);
                break;
            case G1U1:
                result[0][1]=1./(2.*sqrt2);
                result[0][3]=-1./(2.*sqrt2);
                result[1][0]=-1./(2.*sqrt2);
                result[1][2]=1./(2.*sqrt2);
                result[2][1]=-1./(2.*sqrt2);
                result[2][3]=1./(2.*sqrt2);
                result[3][0]=1./(2.*sqrt2);
                result[3][2]=-1./(2.*sqrt2);
                break;
        }
        return result;
    }
    
    Matrix<double> get_polarization(const polarization_id pid){
        Matrix<double> result(4,4,0.);
        
        switch(pid){
            case UPDOWNAV:
                result[0][0]=0.25;
                result[0][2]=0.25;
                result[1][1]=0.25;
                result[1][3]=0.25;
                result[2][0]=0.25;
                result[2][2]=0.25;
                result[3][1]=0.25;
                result[3][3]=0.25;
                break;
            case SING0:
                result[0][1]=1./(2.*sqrt2);
                result[0][3]=1./(2.*sqrt2);
                result[1][0]=(-1.)/(2.*sqrt2);
                result[1][2]=(-1.)/(2.*sqrt2);
                result[2][1]=1./(2.*sqrt2);
                result[2][3]=1./(2.*sqrt2);
                result[3][0]=(-1.)/(2.*sqrt2);
                result[3][2]=(-1.)/(2.*sqrt2);
                break;
            case TRIP0:
                result[0][1]=-1./(2.*sqrt2);
                result[0][3]=-1./(2.*sqrt2);
                result[1][0]=-1./(2.*sqrt2);
                result[1][2]=-1./(2.*sqrt2);
                result[2][1]=-1./(2.*sqrt2);
                result[2][3]=-1./(2.*sqrt2);
                result[3][0]=-1./(2.*sqrt2);
                result[3][2]=-1./(2.*sqrt2);
                break;
            case TRIPP1:
                result[1][1]=0.5;
                result[1][3]=0.5;
                result[3][1]=0.5;
                result[3][3]=0.5;
                break;
            case TRIPM1:
                result[0][0]=0.5;
                result[0][2]=0.5;
                result[2][0]=0.5;
                result[2][2]=0.5;
                break;
        }
        return result;
    }
    
}
