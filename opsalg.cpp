//
//  opsalg.cpp
//  pertutils
//
//  Created by Thorsten Kurth on 23.09.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

#include "mathutils.hpp"
#include "pertutils.hpp"

//******************************************************************
//******************************************************************
//START baryon_op
//******************************************************************
//******************************************************************
//constructors
baryon_op::baryon_op(const std::vector<std::string>& names, const std::vector<double>& coeffs, const std::vector<NRvector<unsigned int> >& spins): spinids(spins), opnames(names), coefficients(coeffs){
    //check sanity:
    if(opnames.size()!=coefficients.size() || opnames.size()!=spinids.size()){
        std::cerr << "baryon_op::baryon_op: inconsistent sizes!" << std::endl;
        this->clear();
    }
};

baryon_op::baryon_op(const baryon_op& rhs): spinids(rhs.spinids), opnames(rhs.opnames), coefficients(rhs.coefficients){};

//member functions:
void baryon_op::clear(){
    opnames.clear();
    coefficients.clear();
    spinids.clear();
}

baryon_op& baryon_op::apply(const std::string& op){
    for(unsigned int n=0; n<opnames.size(); n++){
        opnames[n].append(op);
    }
    return *this;
}

baryon_op& baryon_op::bar(){
    for(unsigned int n=0; n<opnames.size(); n++){
        std::string tmpstring=opnames[n];
        double sign=1.;
        size_t pos=-1;
        while((pos=tmpstring.find("Dp",pos+1))!=std::string::npos){
            sign*=-1.;
            tmpstring.replace(pos, 2, "DM");
        }
        pos=-1;
        while((pos=tmpstring.find("Dm",pos+1))!=std::string::npos){
            sign*=-1.;
            tmpstring.replace(pos, 2, "Dp");
        }
        pos=-1;
        while((pos=tmpstring.find("DM",pos+1))!=std::string::npos){
            tmpstring.replace(pos, 2, "Dm");
        }
        pos=-1;
        while((pos=tmpstring.find("Nb",pos+1))!=std::string::npos){
            tmpstring.replace(pos, 2, "nb");
        }
        pos=-1;
        while((pos=tmpstring.find("N",pos+1))!=std::string::npos){
            tmpstring.insert(pos+1, "b");
        }
        pos=-1;
        while((pos=tmpstring.find("nb",pos+1))!=std::string::npos){
            tmpstring.replace(pos, 2, "N");
            //tmpstring.erase(pos+1);
        }
        opnames[n]=tmpstring;
        coefficients[n]*=sign;
    }
    return *this;
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
    std::vector<std::string> tmpnames;
    std::vector<NRvector<unsigned int> > tmpspins;
    std::vector<double> tmpcoeffs;
    
    for(unsigned int n=0; n<nn; n++) for(unsigned int m=0; m<mm; m++){
        tmpnames.push_back(opnames[n]+" "+rhs.opnames[m]);
        tmpcoeffs.push_back(coefficients[n]*coefficients[m]);
        NRvector<unsigned int> tmpvec(spinids[n]);
        tmpspins.push_back(tmpvec.append(spinids[m]));
    }
    opnames=tmpnames;
    coefficients=tmpcoeffs;
    spinids=tmpspins;
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
    os.precision(10);
    std::string tmp;
    for(unsigned int n=0; n<obj.opnames.size(); n++){
        if(obj.coefficients[n]>=0) tmp=" + ";
        else tmp=" - ";
        
        //check for coeff=0 and coeff=+/-1:
        if(fabs(obj.coefficients[n])<1.e-8) continue;
        else if(fabs(fabs(obj.coefficients[n])-1.)<1.e-8) os << tmp << obj.opnames[n] << "_(";
        else os << tmp << fabs(obj.coefficients[n]) << " * " << obj.opnames[n] << "_(";
        for(unsigned int s=0; s<obj.spinids[n].dim(); s++) os << obj.spinids[n][s];
        os << ")";
    }
    return os;
}

//******************************************************************
//******************************************************************
//END baryon_op
//******************************************************************
//******************************************************************

//******************************************************************
//******************************************************************
//START quark_cont
//******************************************************************
//******************************************************************
quark_cont::quark_cont(const std::vector<NRvector<std::string> >& quarkss, const std::vector<NRvector<std::string> >& attributess, const std::vector<NRvector<std::string> >& idcss, const std::vector<std::string>& sym_coefff, const std::vector<double>& num_coefff): quarks(quarkss), attributes(attributess), idcs(idcss), sym_coeff(sym_coefff), num_coeff(num_coefff){
    //check for consistency:
    if(quarks.size()!=attributes.size() || quarks.size()!=idcs.size() || quarks.size()!=sym_coeff.size() || quarks.size()!=num_coeff.size()){
        std::cerr << "quark_cont::quark_cont: internal structure error!" << std::endl;
        this->clear();
    }
}

quark_cont get_op(const unsigned int& opnumber, const std::string& opname, const NRvector<unsigned int>& spins, const double& coeff){
    std::vector<NRvector<std::string> > quarks;
    std::vector<NRvector<std::string> > attributess;
    std::vector<NRvector<std::string> > idcs;
    std::vector<std::string> sym_coeff;
    std::vector<double> num_coeff;
    
    NRvector<std::string> qvec(3);
    NRvector<std::string> avec(3);
    NRvector<std::string> ivec(3);
    
    //indices:
    std::string col[3];
    col[0]="a";
    col[1]="b";
    col[2]="c";
    for(unsigned int s=0; s<3; s++) ivec[s]=col[s]+std::to_string(opnumber)+";"+std::to_string(spins[s+opnumber]);
    idcs.push_back(ivec);
    idcs.push_back(ivec);
    
    //attributes:
    avec[0]="loc";
    avec[1]="loc";
    avec[2]="loc";
    
    //prefactor eps_abc for every single term:
    std::string tmp="eps_(";
    for(unsigned int s=0; s<2; s++) tmp+=col[s]+std::to_string(opnumber)+";";
    tmp+=col[2]+std::to_string(opnumber)+")";
    sym_coeff.push_back(tmp);
    sym_coeff.push_back(tmp);
    
    size_t pos;
    bool bar=false;
    std::string tmpstring=opname;
    if((pos=tmpstring.find("P"))!=std::string::npos){
        if(tmpstring.find("b")==(pos+1)) bar=true;
        
        qvec[0]="u";
        qvec[1]="d";
        qvec[2]="u";
        if(bar){
            for(unsigned int s=0; s<3; s++) qvec[s]+="b";
        }
        quarks.push_back(qvec);
        num_coeff.push_back(sqrt(1./2.)*coeff);
        
        qvec[0]="d";
        qvec[1]="u";
        qvec[2]="u";
        if(bar){
            for(unsigned int s=0; s<3; s++) qvec[s]+="b";
        }
        quarks.push_back(qvec);
        num_coeff.push_back(-sqrt(1./2.)*coeff);
        
        if(bar) tmpstring.erase(tmpstring.begin()+pos,tmpstring.begin()+(pos+2));
        else tmpstring.erase(tmpstring.begin()+pos,tmpstring.begin()+(pos+1));
    }
    else if((pos=tmpstring.find("N"))!=std::string::npos){
        if(tmpstring.find("b")==(pos+1)) bar=true;
        
        qvec[0]="d";
        qvec[1]="u";
        qvec[2]="d";
        if(bar){
            for(unsigned int s=0; s<3; s++) qvec[s]+="b";
        }
        quarks.push_back(qvec);
        num_coeff.push_back(sqrt(1./2.)*coeff);
        
        qvec[0]="u";
        qvec[1]="d";
        qvec[2]="d";
        if(bar){
            for(unsigned int s=0; s<3; s++) qvec[s]+="b";
        }
        quarks.push_back(qvec);
        num_coeff.push_back(-sqrt(1./2.)*coeff);
        
        if(bar) tmpstring.erase(tmpstring.begin()+pos,tmpstring.begin()+(pos+2));
        else tmpstring.erase(tmpstring.begin()+pos,tmpstring.begin()+(pos+1));
    }
    
    if(tmpstring.compare("")!=0){
        avec[2]=tmpstring;
    }
    attributess.push_back(avec);
    attributess.push_back(avec);
    
    return quark_cont(quarks, attributess, idcs, sym_coeff, num_coeff);
}

quark_cont::quark_cont(const baryon_op& barop){
    std::vector<std::string> token;
    for(unsigned int n=0; n<barop.opnames.size(); n++){
        token.clear();
        tokenize(barop.opnames[n], token);
        quark_cont qcont(get_op(0,token[0],barop.spinids[n],barop.coefficients[n]));
        for(unsigned int s=1; s<token.size(); s++){
            qcont*=get_op(s,token[s],barop.spinids[n],barop.coefficients[n]);
        }
        (*this)+=qcont;
    }
}

quark_cont& quark_cont::operator+=(const quark_cont& rhs){
    quarks.insert( quarks.end(), rhs.quarks.begin(), rhs.quarks.end() );
    attributes.insert( attributes.end(), rhs.attributes.begin(), rhs.attributes.end() );
    idcs.insert( idcs.end(), rhs.idcs.begin(), rhs.idcs.end() );
    sym_coeff.insert( sym_coeff.end(), rhs.sym_coeff.begin(), rhs.sym_coeff.end() );
    num_coeff.insert( num_coeff.end(), rhs.num_coeff.begin(), rhs.num_coeff.end() );
    return *this;
}

quark_cont& quark_cont::operator*=(const quark_cont& rhs){
    unsigned int nn=static_cast<unsigned int>(quarks.size());
    unsigned int mm=static_cast<unsigned int>(rhs.quarks.size());
    std::vector<NRvector<std::string> > tmpquarks;
    std::vector<NRvector<std::string> > tmpattributes;
    std::vector<NRvector<std::string> > tmpidcs;
    std::vector<std::string> tmpsym_coeff;
    std::vector<double> tmpnum_coeff;
    
    for(unsigned int n=0; n<nn; n++) for(unsigned int m=0; m<mm; m++){
        //quarks:
        NRvector<std::string> tmpvec(quarks[n]);
        tmpvec.append(rhs.quarks[m]);
        tmpquarks.push_back(tmpvec);
        
        //attributes:
        tmpvec.assign(attributes[n]);
        tmpvec.append(rhs.attributes[m]);
        tmpattributes.push_back(tmpvec);
        
        //indices:
        tmpvec.assign(idcs[n]);
        tmpvec.append(rhs.idcs[m]);
        tmpidcs.push_back(tmpvec);
        
        std::string tmp=sym_coeff[n]+" "+rhs.sym_coeff[m];
        tmpsym_coeff.push_back(tmp);
        tmpnum_coeff.push_back(num_coeff[n]*rhs.num_coeff[m]);
    }
    quarks=tmpquarks;
    attributes=tmpattributes;
    idcs=tmpidcs;
    sym_coeff=tmpsym_coeff;
    num_coeff=tmpnum_coeff;
    
    return *this;
}

std::ostream& operator<<(std::ostream &os,const quark_cont &obj){
    os.precision(10);
    std::string tmp;
    for(unsigned int n=0; n<obj.quarks.size(); n++){
        if(obj.num_coeff[n]>=0) tmp=" + ";
        else tmp=" - ";
        tmp+=obj.sym_coeff[n]+" *";
        
        //check for coeff=0 and coeff=+/-1:
        if(fabs(obj.num_coeff[n])<1.e-8) continue;
        os << tmp;
        if(fabs(fabs(obj.num_coeff[n])-1.)>1.e-8) os << " " << fabs(obj.num_coeff[n]) << " *";
        for(unsigned int s=0; s<obj.quarks[n].dim(); s++){
            os << " " << obj.quarks[n][s] << "_(" << obj.idcs[n][s] << ")";
        }
    }
    return os;
}

void quark_cont::clear(){
    quarks.clear();
    attributes.clear();
    idcs.clear();
    sym_coeff.clear();
    num_coeff.clear();
}
//******************************************************************
//******************************************************************
//END quark_op
//******************************************************************
//******************************************************************
