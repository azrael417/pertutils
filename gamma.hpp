//
//  gamma.hpp
//  pertutils
//
//  Created by Thorsten Kurth on 23.04.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

#ifndef _GAMMA
#define _GAMMA
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
    PARMINUS
};

NRmatrix<dcomplex> get_gamma(const gamma_id gid);

#endif
