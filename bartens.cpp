//
//  bartens.cpp
//  pertutils
//
//  Created by Thorsten Kurth on 24.06.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

#include "mathutils.hpp"
#include "pertutils.hpp"

void BBTensor::join(const BBTensor& rhs, const bool& asym){
    //Compute Kronecker Product of tensors:
    std::vector<TTTensor> newdata;
    for(unsigned int t=0; t<nt; t++){
        TTTensor ltr(kron(data[t],rhs.data[t]));
        newdata.push_back(ltr);
    }
    data=newdata;
    
    //Join vectors of source positions:
    std::vector<quark> newquarks;
    newquarks.reserve(quarks.size()+rhs.quarks.size());
    newquarks.insert(newquarks.end(),quarks.begin(),quarks.end());
    newquarks.insert(newquarks.end(),rhs.quarks.begin(),rhs.quarks.end());
    quarks=newquarks;
    
    if(asym){
        
    }
}