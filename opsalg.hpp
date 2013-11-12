//
//  opsalg.hpp
//  pertutils
//
//  Created by Thorsten Kurth on 23.09.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

#ifndef _OPSALG
#define _OPSALG

class quark_cont;

class baryon_op{
private:
    std::vector<NRvector<std::string> > spinids;
    std::vector<std::string> opnames;
    std::vector<double> coefficients;
    
public:
    //constructors:
    baryon_op(){};
    baryon_op(const std::vector<std::string>&,const std::vector<double>&,const std::vector<NRvector<std::string> >&);
    baryon_op(const baryon_op&);
    
    //usual members:
    void clear();
    baryon_op& apply(const std::string&);
    baryon_op& bar();
    
    //operators:
    baryon_op& operator=(const baryon_op& rhs);
    baryon_op& operator+=(const baryon_op& rhs);
    baryon_op& operator-=(const baryon_op& rhs);
    baryon_op& operator*=(const double& rhs);
    baryon_op& operator*=(const baryon_op& rhs);
    friend baryon_op operator+(const baryon_op& lhs, const baryon_op& rhs);
    friend baryon_op operator-(const baryon_op& lhs, const baryon_op& rhs);
    friend baryon_op operator*(const baryon_op& lhs, const double& rhs);
    friend baryon_op operator*(const double& lhs, const baryon_op& rhs);
    friend baryon_op operator*(const baryon_op& lhs, const baryon_op& rhs);
    friend std::ostream& operator<<(std::ostream &os,const baryon_op &obj);
    
    //destructors:
    ~baryon_op(){ this->clear(); };
    
    //other friends:
    friend class quark_cont;
};

class quark_cont{
private:
    std::vector<NRvector<std::string> > quarks;
    std::vector<NRvector<std::string> > attributes;
    std::vector<NRvector<std::string> > idcs;
    std::vector<std::string> sym_coeff;
    std::vector<double> num_coeff;
    
    std::vector<NRvector<std::string> > props;
    std::vector<NRvector<std::string> > props_idcs;
    std::vector<double> props_signs;
    
public:
    quark_cont(){};
    quark_cont(const std::vector<NRvector<std::string> >&, const std::vector<NRvector<std::string> >&, const std::vector<NRvector<std::string> >&, const std::vector<std::string>&, const std::vector<double>&);
    quark_cont(const baryon_op&);
    ~quark_cont(){ this->clear(); };
    
    //operators:
    quark_cont& operator+=(const quark_cont& rhs);
    quark_cont& operator*=(const quark_cont& rhs);
    friend std::ostream& operator<<(std::ostream &os,const quark_cont &obj);
    
    //member functions:
    void clear();
    int contract();
    int reorder();
    int print_contractions(std::ostream& os, const std::string mode="human-readable");
};

quark_cont get_op(const unsigned int& opnumber, const std::string& opname, const NRvector<std::string>& spins, const double& coeff);
#endif
