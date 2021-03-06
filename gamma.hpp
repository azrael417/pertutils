//
//  gamma.hpp
//  pertutils
//
//  Created by Thorsten Kurth on 23.04.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

//this files contains 4D spinor matrix generation:

#ifndef _GAMMA
#define _GAMMA

namespace pertutils{
    
    enum gamma_id {
        ONE,
        GX,
        GY,
        GZ,
        GT,
        G5,
        CG5,
        C,
        IG4CG5,
        PARPLUS,
        PARMINUS,
    };
    
    enum basak_gamma_id {
        UTRANS,
        G1G1,
        G1U1
    };
    
    enum polarization_id{
        UPDOWNAV,
        SING0,
        TRIP0,
        TRIPP1,
        TRIPM1
    };
    
    Matrix<dcomplex> get_gamma(const gamma_id gid);
    Matrix<double> get_gamma(const basak_gamma_id gid);
    Matrix<double> get_polarization(const polarization_id gid);
    
}
#endif
