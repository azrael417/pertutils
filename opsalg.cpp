//
//  opsalg.cpp
//  pertutils
//
//  Created by Thorsten Kurth on 23.09.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

#include "mathutils.hpp"
#include "pertutils.hpp"

//constructors
baryon_op::baryon_op(const std::vector<std::string>& names, const std::vector<double>& coeffs, const std::vector<NRvector<unsigned int> >& spins): opnames(names), coefficients(coeffs), spinids(spins){
    //check sanity:
    if(opnames.size()!=coefficients.size() || opnames.size()!=spinids.size()){
        std::cerr << "baryon_op::baryon_op: inconsistent sizes!" << std::endl;
        this->clear();
    }
};

baryon_op::baryon_op(const baryon_op& rhs): opnames(rhs.opnames), coefficients(rhs.coefficients), spinids(rhs.spinids){};

//destructors:
baryon_op::~baryon_op(){
    this->clear();
}

//member functions:
void baryon_op::clear(){
    opnames.clear();
    coefficients.clear();
    spinids.clear();
}

//operators:
//unary:
baryon_op& baryon_op::operator=(const baryon_op& rhs){
    this->clear();
    opnames=rhs.opnames;
    coefficients=rhs.coefficients;
    spinids=rhs.spinids;
    return *this;
}

baryon_op& baryon_op::operator+=(const baryon_op& rhs){
    opnames.insert( opnames.end(), rhs.opnames.begin(), rhs.opnames.end() );
    coefficients.insert( coefficients.end(), rhs.coefficients.begin(), rhs.coefficients.end() );
    spinids.insert( spinids.end(), rhs.spinids.begin(), rhs.spinids.end() );
    return *this;
}

baryon_op& baryon_op::operator-=(const baryon_op& rhs){
    std::vector<double> tmpcoeffs(rhs.coefficients);
    for(unsigned int n=0; n<tmpcoeffs.size(); n++) tmpcoeffs[n]*=-1.;
    opnames.insert( opnames.end(), rhs.opnames.begin(), rhs.opnames.end() );
    coefficients.insert( coefficients.end(), tmpcoeffs.begin(), tmpcoeffs.end() );
    spinids.insert( spinids.end(), rhs.spinids.begin(), rhs.spinids.end() );
    return *this;
}

baryon_op& baryon_op::operator*=(const double& rhs){
    for(unsigned int n=0; n<coefficients.size(); n++) coefficients[n]*=rhs;
    return *this;
}

baryon_op& baryon_op::operator*=(const baryon_op& rhs){
    unsigned int nn=static_cast<unsigned int>(opnames.size());
    unsigned int mm=static_cast<unsigned int>(rhs.opnames.size());
    unsigned int spincount=spinids[0].dim()+rhs.spinids[0].dim();
    std::vector<std::string> tmpnames;
    std::vector<NRvector<unsigned int> > tmpspins;
    std::vector<double> tmpcoeffs;
    
    NRvector<unsigned> tmpsp(spincount);
    for(unsigned int n=0;n<nn; n++) for(unsigned int m=0; m<mm; m++){
        tmpnames.push_back(opnames[n].append(" "+rhs.opnames[m]));
        tmpcoeffs.push_back(coefficients[n]*coefficients[m]);
        tmpspins.push_back(spinids[n].append(spinids[m]));
    }
    return *this;
}

//binary:
baryon_op operator+(const baryon_op& lhs, const baryon_op& rhs){
    baryon_op result(lhs);
    result+=rhs;
    return result;
}

baryon_op operator-(const baryon_op& lhs, const baryon_op& rhs){
    baryon_op result(lhs);
    result-=rhs;
    return result;
}

baryon_op operator*(const baryon_op& lhs, const double& rhs){
    baryon_op result(lhs);
    result*=rhs;
    return result;
}

baryon_op operator*(const double& lhs, const baryon_op& rhs){
    baryon_op result(rhs);
    result*=lhs;
    return result;
}

baryon_op operator*(const baryon_op& lhs, const baryon_op& rhs){
    baryon_op result(lhs);
    result*=rhs;
    return result;
}

std::ostream& operator<<(std::ostream &os,const baryon_op &obj){
    std::string tmp;
    for(unsigned int n=0; n<obj.opnames.size(); n++){
        if(n>0 && obj.coefficients[n]>=0) tmp="+";
        else tmp="";
        os << tmp << obj.coefficients[n] << "*" << obj.opnames[n] << "_(";
        for(unsigned int s=0; s<obj.spinids[n].dim(); s++) os << obj.spinids[n][s];
        os << ") ";
    }
    return os;
}
