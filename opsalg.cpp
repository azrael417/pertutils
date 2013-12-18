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
baryon_op::baryon_op(const std::vector<std::string>& names, const std::vector<double>& coeffs, const std::vector<NRvector<std::string> >& spins): spinids(spins), opnames(names), coefficients(coeffs){
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
    std::vector<NRvector<std::string> > tmpspins;
    std::vector<double> tmpcoeffs;
    
    for(unsigned int n=0; n<nn; n++) for(unsigned int m=0; m<mm; m++){
        tmpnames.push_back(opnames[n]+" "+rhs.opnames[m]);
        tmpcoeffs.push_back(coefficients[n]*rhs.coefficients[m]);
        NRvector<std::string> tmpvec(spinids[n]);
        tmpspins.push_back(tmpvec.append(rhs.spinids[m]));
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
quark_cont::quark_cont(const std::vector<NRvector<std::string> >& quarkss, const std::vector<NRvector<std::string> >& attributess, const std::vector<NRvector<std::string> >& idcss, const std::vector<std::string>& sym_coefff, const std::vector<double>& num_coefff): quarks(quarkss), attributes(attributess), idcs(idcss), sym_coeff(sym_coefff), num_coeff(num_coefff), isolimit(false), operator_id(-1){
    //check for consistency:
    if(quarks.size()!=attributes.size() || quarks.size()!=idcs.size() || quarks.size()!=sym_coeff.size() || quarks.size()!=num_coeff.size()){
        std::cerr << "quark_cont::quark_cont: internal structure error!" << std::endl;
        this->clear();
    }
}

static quark_cont get_op(const unsigned int& opnumber, const std::string& opname, const NRvector<std::string>& spins){
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
    for(unsigned int s=0; s<3; s++) ivec[s]=col[s]+std::to_string(opnumber)+","+spins[s+3*opnumber];
    idcs.push_back(ivec);
    idcs.push_back(ivec);
    
    //attributes:
    avec[0]="loc";
    avec[1]="loc";
    avec[2]="loc";
    
    //prefactor eps_abc for every single term:
    std::string tmp="eps_(";
    for(unsigned int s=0; s<2; s++) tmp+=col[s]+std::to_string(opnumber)+",";
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
        num_coeff.push_back(sqrt(1./2.));
        
        qvec[0]="d";
        qvec[1]="u";
        qvec[2]="u";
        if(bar){
            for(unsigned int s=0; s<3; s++) qvec[s]+="b";
        }
        quarks.push_back(qvec);
        num_coeff.push_back(-sqrt(1./2.));
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
        num_coeff.push_back(sqrt(1./2.));
        
        qvec[0]="u";
        qvec[1]="d";
        qvec[2]="d";
        if(bar){
            for(unsigned int s=0; s<3; s++) qvec[s]+="b";
        }
        quarks.push_back(qvec);
        num_coeff.push_back(-sqrt(1./2.));
    }
    if(bar) tmpstring.erase(tmpstring.begin()+pos,tmpstring.begin()+(pos+2));
    else tmpstring.erase(tmpstring.begin()+pos,tmpstring.begin()+(pos+1));
    
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
        
        //start with first operator:
        quark_cont qcont(get_op(0,token[0],barop.spinids[n]));
        for(unsigned int s=1; s<token.size(); s++){
            qcont*=get_op(s,token[s],barop.spinids[n]);
        }
        //multiply overall prefactor:
        qcont*=barop.coefficients[n];
        
        (*this)+=qcont;
    }
    isolimit=false;
    operator_id=-1;
}

quark_cont& quark_cont::operator+=(const quark_cont& rhs){
    quarks.insert( quarks.end(), rhs.quarks.begin(), rhs.quarks.end() );
    attributes.insert( attributes.end(), rhs.attributes.begin(), rhs.attributes.end() );
    idcs.insert( idcs.end(), rhs.idcs.begin(), rhs.idcs.end() );
    sym_coeff.insert( sym_coeff.end(), rhs.sym_coeff.begin(), rhs.sym_coeff.end() );
    num_coeff.insert( num_coeff.end(), rhs.num_coeff.begin(), rhs.num_coeff.end() );
    return *this;
}


quark_cont& quark_cont::operator*=(const double& rhs){
    for(unsigned int n=0; n<num_coeff.size(); n++) num_coeff[n]*=rhs;
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
            if(obj.attributes[n][s].compare("loc")==0) os << " " << obj.quarks[n][s] << "_(" << obj.idcs[n][s] << ")";
            else os << " " << obj.quarks[n][s] << obj.attributes[n][s] << "_(" << obj.idcs[n][s] << ")";
        }
        os << std::endl;
    }
    return os;
}

//member functions:
void quark_cont::clear(){
    quarks.clear();
    attributes.clear();
    idcs.clear();
    sym_coeff.clear();
    num_coeff.clear();
    operator_id=-1;
    isolimit=false;
}

void quark_cont::toggle_isospin_limit(){
    isolimit=!isolimit;
}

void quark_cont::set_operator_id(const int& id){
    operator_id=id;
}

static int contract_helper(const NRvector<std::string>& quarks, const NRvector<std::string>& attributes, const NRvector<std::string>& idcs, std::vector<NRvector<std::string> >& resprops, std::vector<NRvector<std::string> >& residcs, std::vector<double>& signs){
    //count quarks and barred-quarks first and compare;
    unsigned int nu=0, nubar=0, nd=0, ndbar=0, ns=0, nsbar=0;
    for(unsigned int s=0; s<quarks.dim(); s++){
        if(quarks[s].find("ub")==0){
            nubar++;
        }
        else if(quarks[s].find("u")==0){
            nu++;
        }
        
        if(quarks[s].find("db")==0){
            ndbar++;
        }
        else if(quarks[s].find("d")==0){
            nd++;
        }
        
        if(quarks[s].find("sb")==0){
            nsbar++;
        }
        else if(quarks[s].find("s")==0){
            ns++;
        }
    }
    
    if(nu!=nubar || nd!=ndbar || ns!=nsbar){
        std::cerr << "contract_helper: error, quark content does not allow for contraction!" << std::endl;
        return EXIT_FAILURE;
    }
    
    //sort vector: innermost: u, outermost s:
    std::vector<unsigned int> table(quarks.dim());
    unsigned int runs=0,rund=ns,runu=ns+nd,runubar=ns+nd+nu,rundbar=ns+nd+nu+nubar,runsbar=ns+nd+nu+nubar+ndbar;
    for(unsigned int s=0; s<quarks.dim(); s++){
        if(quarks[s].find("s")==0 && quarks[s].find("sb")!=0){
            table[s]=runs;
            runs++;
        }
        else if(quarks[s].find("d")==0 && quarks[s].find("db")!=0){
            table[s]=rund;
            rund++;
        }
        else if(quarks[s].find("u")==0 && quarks[s].find("ub")!=0){
            table[s]=runu;
            runu++;
        }
        else if(quarks[s].find("ub")==0){
            table[s]=runubar;
            runubar++;
        }
        else if(quarks[s].find("db")==0){
            table[s]=rundbar;
            rundbar++;
        }
        else if(quarks[s].find("sb")==0){
            table[s]=runsbar;
            runsbar++;
        }
    }
    
    //sort vectors;
    NRvector<std::string> tmpquarks(quarks.dim()), tmpattributes(quarks.dim()), tmpidcs(quarks.dim());
    for(unsigned int s=0; s<tmpquarks.dim(); s++){
        tmpquarks[table[s]]=quarks[s];
        tmpattributes[table[s]]=attributes[s];
        tmpidcs[table[s]]=idcs[s];
    }
    int globpar=sort(table);
    
    //compute number of contractions:
    unsigned int nscont=static_cast<unsigned int>(factrl(static_cast<int>(ns)));
    unsigned int ndcont=static_cast<unsigned int>(factrl(static_cast<int>(nd)));
    unsigned int nucont=static_cast<unsigned int>(factrl(static_cast<int>(nu)));
    unsigned int ncont=nscont*ndcont*nucont;
    
    //initialize vectors:
    signs.resize(ncont);
    resprops.resize(ncont);
    residcs.resize(ncont);
    
    //go through contractions:
    std::vector<unsigned int> svec(ns),svectmp(ns);
    std::vector<unsigned int> dvec(nd),dvectmp(nd);
    std::vector<unsigned int> uvec(nu),uvectmp(nu);
    for(unsigned int s=0; s<ns; s++) svec[s]=s;
    for(unsigned int s=0; s<nd; s++) dvec[s]=s;
    for(unsigned int s=0; s<nu; s++) uvec[s]=s;
    int spar,dpar,upar;
    
    NRvector<std::string> props(ns+nd+nu);
    NRvector<std::string> indices(ns+nd+nu);
    for(unsigned int nsc=0; nsc<nscont; nsc++){
        reverse_copy(svec.begin(),svec.end(),svectmp.begin());
        for(unsigned int s=0; s<ns; s++){
            unsigned int offs=ns+nd+nu+nubar+ndbar+svectmp[s];
            props[s]="Ss("+tmpattributes[s]+";"+tmpattributes[offs]+")";
            indices[s]="("+tmpidcs[s]+";"+tmpidcs[offs]+")";
        }
        svectmp=svec;
        spar=sort(svectmp);
        
        for(unsigned int ndc=0; ndc<ndcont; ndc++){
            reverse_copy(dvec.begin(),dvec.end(),dvectmp.begin());
            for(unsigned int d=0; d<nd; d++){
                unsigned int offd=ns+nd+nu+nubar+dvectmp[d];
                props[ns+d]="Sd("+tmpattributes[ns+d]+";"+tmpattributes[offd]+")";
                indices[ns+d]="("+tmpidcs[ns+d]+";"+tmpidcs[offd]+")";
            }
            dvectmp=dvec;
            dpar=sort(dvectmp);
            
            for(unsigned int nuc=0; nuc<nucont; nuc++){
                reverse_copy(uvec.begin(),uvec.end(),uvectmp.begin());
                for(unsigned int u=0; u<nu; u++){
                    unsigned int offu=ns+nd+nu+uvectmp[u];
                    props[ns+nd+u]="Su("+tmpattributes[ns+nd+u]+";"+tmpattributes[offu]+")";
                    indices[ns+nd+u]="("+tmpidcs[ns+nd+u]+";"+tmpidcs[offu]+")";
                }
                uvectmp=uvec;
                upar=sort(uvectmp);
                
                unsigned int idx=nuc+nucont*(ndc+ndcont*nsc);
                signs[idx]=static_cast<double>(globpar*upar*dpar*spar);
                resprops[idx]=props;
                residcs[idx]=indices;
                
                std::next_permutation(uvec.begin(),uvec.end());
            }
            std::next_permutation(dvec.begin(),dvec.end());
        }
        std::next_permutation(svec.begin(),svec.end());
    }
    
    return EXIT_SUCCESS;
}

int quark_cont::contract(){
    props.clear();
    props_idcs.clear();
    props_signs.clear();
    int error;
    
    for(unsigned int n=0; n<quarks.size(); n++){
        std::vector<NRvector<std::string> > tmpprops;
        std::vector<NRvector<std::string> > tmpidcs;
        std::vector<double> tmpsigns;
        error=contract_helper(quarks[n],attributes[n],idcs[n],tmpprops,tmpidcs,tmpsigns);
        if(error==EXIT_FAILURE){
            std::cerr << "quark_cont::contract: error, contraction could not be performed" << std::endl;
            return EXIT_FAILURE;
        }
        props.insert(props.end(), tmpprops.begin(), tmpprops.end());
        props_idcs.insert(props_idcs.end(), tmpidcs.begin(), tmpidcs.end());
        props_signs.insert(props_signs.end(), tmpsigns.begin(), tmpsigns.end());
    }
    return reorder();
}

//reorder quark propagators, such that sink side is regrouped into baryons:
int quark_cont::reorder(){
    if(props.size()==0 || props_idcs.size()==0 || props_signs.size()==0){
        std::cerr << "quark_cont::reorder: please perform contractions first!" << std::endl;
        return EXIT_FAILURE;
    }
    
    unsigned int numprops=static_cast<unsigned int>(props.size()/quarks.size());
    
    //compute pivoting table:
    std::vector<NRvector<unsigned int> > proppivots;
    unsigned int count=0;
    for(unsigned int n=0; n<quarks.size(); n++){
        for(unsigned int s=0; s<numprops; s++){
            NRvector<unsigned int> pivot(props[s+n*numprops].dim());
            for(unsigned int p=0; p<props[s+n*numprops].dim(); p++){
                std::string tmp=props_idcs[s+n*numprops][p];
                size_t pos=tmp.find_first_of(",");
                tmp=tmp.substr(1,pos);
                if(tmp.find("a")==0) count=0;
                if(tmp.find("b")==0) count=1;
                if(tmp.find("c")==0) count=2;
                tmp.erase(0,1);
                pivot[p]=static_cast<unsigned int>(strtoul(tmp.c_str(),NULL,10)*3+count);
            }
            proppivots.push_back(pivot);
        }
    }
    
    //reorder:
    for(unsigned int n=0; n<quarks.size(); n++){
        for(unsigned int s=0; s<numprops; s++){
            NRvector<std::string> tmpprops(props[s+n*numprops].dim());
            NRvector<std::string> tmpprops_idcs(props[s+n*numprops].dim());
            for(unsigned int p=0; p<props[s+n*numprops].dim(); p++){
                tmpprops[proppivots[s+n*numprops][p]]=props[s+n*numprops][p];
                tmpprops_idcs[proppivots[s+n*numprops][p]]=props_idcs[s+n*numprops][p];
            }
            for(unsigned int p=0; p<props[s+n*numprops].dim(); p++){
                props[s+n*numprops][p]=tmpprops[p];
                props_idcs[s+n*numprops][p]=tmpprops_idcs[p];
            }
        }
    }
    
    return EXIT_SUCCESS;
}

static void get_indices_laph(const std::string& prop, const std::string& idcs, unsigned int& massid, unsigned int& spin1id, unsigned int& spin2id, unsigned int& n1id, unsigned int& n2id, bool isolimit=false){
    unsigned int sub=0;
    if(isolimit) sub=1;

    if(prop.find("Su")==0) massid=0;
    if(prop.find("Sd")==0) massid=1-sub;
    if(prop.find("Ss")==0) massid=2-sub;
    std::vector<std::string> tmptoken1, tmptoken2;
    tokenize(idcs, tmptoken1,";");
    
    //examine single parts:
    //part 1:
    tmptoken2.clear();
    tokenize(tmptoken1[0],tmptoken2,",");
    spin1id=static_cast<unsigned int>(strtoul(tmptoken2[1].c_str(),NULL,10));
    spin1id-=1;
    unsigned int count;
    tmptoken2[0].erase(0,1);
    if(tmptoken2[0].find("a")==0) count=0;
    if(tmptoken2[0].find("b")==0) count=1;
    if(tmptoken2[0].find("c")==0) count=2;
    tmptoken2[0].erase(0,1);
    n1id=static_cast<unsigned int>(strtoul(tmptoken2[0].c_str(),NULL,10)*3+count);
    
    //part 2:
    tmptoken2.clear();
    tokenize(tmptoken1[1],tmptoken2,",");
    tmptoken2[1].erase(tmptoken2[1].end()-1);
    spin2id=static_cast<unsigned int>(strtoul(tmptoken2[1].c_str(),NULL,10));
    spin2id-=1;
    if(tmptoken2[0].find("a")==0) count=0;
    if(tmptoken2[0].find("b")==0) count=1;
    if(tmptoken2[0].find("c")==0) count=2;
    tmptoken2[0].erase(0,1);
    n2id=static_cast<unsigned int>(strtoul(tmptoken2[0].c_str(),NULL,10)*3+count);
}


int quark_cont::get_laph_sinks(std::string mode){
    unsigned int numfacts=props[0].dim();
    unsigned int numprops=static_cast<unsigned int>(props.size()/quarks.size());
    
    if(mode.compare("laph2")==0){
        for(unsigned int n=0; n<num_coeff.size(); n++){
            for(unsigned int p=0; p<numprops; p++){
                NRvector<std::string> wwwstrings(numfacts/3);
                unsigned int count=0;
                for(unsigned int s=0; s<numfacts; s+=3){
                    
                    //get indices:
                    unsigned int massid[3],spin1id[3],spin2id[3],n1id[3],n2id[3];
                    get_indices_laph(props[p+numprops*n][s+0],props_idcs[p+numprops*n][s+0],massid[0],spin1id[0],spin2id[0],n1id[0],n2id[0],isolimit);
                    get_indices_laph(props[p+numprops*n][s+1],props_idcs[p+numprops*n][s+1],massid[1],spin1id[1],spin2id[1],n1id[1],n2id[1],isolimit);
                    get_indices_laph(props[p+numprops*n][s+2],props_idcs[p+numprops*n][s+2],massid[2],spin1id[2],spin2id[2],n1id[2],n2id[2],isolimit);

                    //get attributes:
                    std::string tmp[3]={"","",""};
                    std::string attrstring="";
                    for(unsigned int l=0; l<3; l++){
                        std::string tmpstring=props[p+numprops*n][s+l];
                        tmpstring.erase(0,3);
                        tmpstring.erase(tmpstring.length()-1,1);
                        tmpstring=tmpstring.substr(0,tmpstring.find_first_of(";"));
                        unsigned int id=n1id[l]%3;
                        if(tmpstring.compare("loc")==0) tmp[id]+="w";
                        else if(tmpstring.compare("A")==0) tmp[id]+="wa";
                        else if(tmpstring.compare("Dp")==0) tmp[id]+="wdp";
                        else if(tmpstring.compare("D0")==0) tmp[id]+="wd0";
                        else if(tmpstring.compare("Dm")==0) tmp[id]+="wdm";
                        else{
                            std::cerr << "quark_cont::get_laph_sinks: error, specified attribute " << tmpstring << " not in database!" << std::endl;
                            return EXIT_FAILURE;
                        }
                    }

                    //correct for indices, since only half of the indices is used, the others are summed first:
                    for(unsigned int l=0; l<3; l++){
                        n2id[l]-=3*(s+1);
                        attrstring+=tmp[l];
                    }
                    
                    //check for attributes at prop3:
                    attrstring+="0pt["+std::to_string(massid[0])+"]["+std::to_string(massid[1])+"]["+std::to_string(massid[2])+"][src]["+std::to_string(spin2id[0])+"][n"+std::to_string(n2id[0])+"]["+std::to_string(spin2id[1])+"][n"+std::to_string(n2id[1])+"][tf]["+std::to_string(spin2id[2])+"][n"+std::to_string(n2id[2])+"]["+std::to_string(spin1id[0])+"]["+std::to_string(spin1id[1])+"]["+std::to_string(spin1id[2])+"]";
                    
                    wwwstrings[count]=attrstring;
                    count++;
                }
                laph_sinks.push_back(wwwstrings);
            }
        }
    }
    return EXIT_SUCCESS;
}

int quark_cont::get_laph_sources(std::string mode){
    unsigned int numfacts=props[0].dim();
    unsigned int numprops=static_cast<unsigned int>(props.size()/quarks.size());
    
    if(mode.compare("laph2")==0){
        for(unsigned int n=0; n<num_coeff.size(); n++){
            for(unsigned int p=0; p<numprops; p++){
                NRvector<std::string> vvvstrings(numfacts/3);
                unsigned int count=0;
                for(unsigned int s=0; s<numfacts; s+=3){
                    
                    //get indices:
                    unsigned int massid[3],spin1id[3],spin2id[3],n1id[3],n2id[3];
                    get_indices_laph(props[p+numprops*n][s+0],props_idcs[p+numprops*n][s+0],massid[0],spin1id[0],spin2id[0],n1id[0],n2id[0]);
                    get_indices_laph(props[p+numprops*n][s+1],props_idcs[p+numprops*n][s+1],massid[1],spin1id[1],spin2id[1],n1id[1],n2id[1]);
                    get_indices_laph(props[p+numprops*n][s+2],props_idcs[p+numprops*n][s+2],massid[2],spin1id[2],spin2id[2],n1id[2],n2id[2]);
                    
                    //get attributes:
                    std::string tmp[3]={"","",""};
                    std::string attrstring="";
                    for(unsigned int l=0; l<3; l++){
                        std::string tmpstring=props[p+numprops*n][s+l];
                        tmpstring.erase(0,3);
                        tmpstring.erase(tmpstring.length()-1,1);
                        tmpstring=tmpstring.substr(tmpstring.find_first_of(";")+1,tmpstring.length());
                        unsigned int id=n1id[l]%3;
                        if(tmpstring.compare("loc")==0) tmp[id]+="v";
                        else if(tmpstring.compare("A")==0) tmp[id]+="va";
                        else if(tmpstring.compare("Dp")==0) tmp[id]+="vdp";
                        else if(tmpstring.compare("D0")==0) tmp[id]+="vd0";
                        else if(tmpstring.compare("Dm")==0) tmp[id]+="vdm";
                        else{
                            std::cerr << "quark_cont::get_laph_sinks: error, specified attribute " << tmpstring << " not in database!" << std::endl;
                            return EXIT_FAILURE;
                        }
                    }
                    //sort arrays since permutations are already included in the www fields:
                    sort2(n2id,tmp,3);
                    
                    //correct for indices, since only half of the indices is used, the others are summed first:
                    for(unsigned int l=0; l<3; l++){
                        n2id[l]-=3*(s+1);
                        attrstring+=tmp[l];
                    }
                    
                    //check for attributes at prop3:
                    attrstring+="0pt[n"+std::to_string(n2id[0])+"][n"+std::to_string(n2id[1])+"][ti][n"+std::to_string(n2id[2])+"]";
                    
                    vvvstrings[count]=attrstring;
                    count++;
                }
                laph_sources.push_back(vvvstrings);
            }
        }
    }
    
    return EXIT_SUCCESS;
}

//int quark_cont::print_laph_baryon_source(std::ostream& os, const std::string mode){
//    //********************************************************
//    //********************************************************
//    //set up source:
//    //********************************************************
//    //********************************************************
//    std::string indent;
//    unsigned int numfacts=props[0].dim();
//    
//    if(mode.compare("local")==0){
//        os << "//set up source:" << std::endl;
//        os << "double complex *vvv1;\n";
//        os << "MALLOC(vvv1, nt*LAPH*LAPH*LAPH)\n";
//        os << "double complex (*vvv1pt)[LAPH][LAPH][LAPH]= (double complex (*)[LAPH][LAPH][LAPH])vvv1;\n";
//        for(unsigned int i=0; i<3; i++){
//            os << indent << "unsigned int n" << i << ";\n";
//            os << indent << "for(n" << i << "=0; n" << i << "<LAPH; n" << i << "++){\n";
//            indent+="\t";
//        }
//        
//        os << indent << "unsigned int tf;\n";
//        os << indent << "for (tf=0; tf<lt; tf++){\n";
//        indent+="\t";
//        
//        //compute vvv1pt string from vvv0pt-string:
//        os << indent << "vvv1pt[tf][n0][n1][n2]= vvv0pt[n0][n1][tf][n2];\n";
//        indent.erase(0,1);
//        os << indent << "}\n";
//        for(unsigned int i=0; i<3; i++){
//            indent.erase((numfacts-i-1),1);
//            os << indent << "} //end loop n" << (numfacts-i-1) << "\n";
//        }
//        os << "memcpy( vvv0, vvv1, lt*LAPH*LAPH*LAPH*sizeof(double complex) );\n";
//        os << "comm_collect_char_dir( (char *)vvv1, (char *)vvv0, lt*LAPH*LAPH*LAPH*sizeof(double complex), TUP );\n";
//    }
//    return EXIT_SUCCESS;
//}

//
//int quark_cont::print_laph_sums(std::ostream& os, const std::string& wwwsummed){
//    //********************************************************
//    //********************************************************
//    //Laph sums:
//    //********************************************************
//    //********************************************************
//    std::string indent;
//    unsigned int numfacts=props[0].dim();
//    
//    os << std::endl;
//    os << "/* LAPH sums */\n";
//    os << "stopper_start( &tmr_laphsum );\n";
//    os << "int tf, ti, src;\n";
//    os << "for (ti=srcstart, src=0; ti<nt; ti+=srcinc, src++) for (tf=0; tf<lt; tf++){\n";
//    indent+="\t";
//    os << indent << "int tf1= lt*mynode_dir[TUP]+tf;\n";
//    os << indent << "int tdiff= (tf1-ti+nt)%nt;\n";
//    os << std::endl << indent << "double complex sum= 0.0;\n";
//    for(unsigned int i=0; i<numfacts; i++){
//        os << indent << "unsigned int n" << i << ";\n";
//        os << indent << "for(n" << i << "=0; n" << i << "<LAPH; n" << i << "++){\n";
//        indent+="\t";
//    }
//    
//    os << indent << "sum += " << wwwsummed;
//    for(unsigned int p=0; p<numfacts; p+=3){
//        os << " * vvv1pt[ti][n" << p << "][n" << p+1 << "][n" << p+2 << "]";
//    }
//    os << ";";
//    
//    os << std::endl;
//    for(unsigned int i=0; i<numfacts; i++){
//        indent.erase((numfacts-i-1),1);
//        os << indent << "} //end loop n" << (numfacts-i-1) << "\n";
//    }
//    
//    os << indent << "P->bar[0][0][tdiff]+= sum;\n";
//    indent.erase(0,1);
//    os << indent << "} //end t-loop\n";
//    indent.erase(0,1);
//    os << indent << "stopper_stop( &tmr_laphsum );\n";
//    
//    return EXIT_SUCCESS;
//}

int quark_cont::print_contractions(std::ostream& os, const std::string mode){
    if(props.size()==0 || props_idcs.size()==0 || props_signs.size()==0){
        std::cerr << "quark_cont::print_contractions: please perform contractions first!" << std::endl;
        return EXIT_FAILURE;
    }
    
    unsigned int numprops=static_cast<unsigned int>(props.size()/quarks.size());
    
    if(mode.compare("human-readable")==0){
        //print contractions in human readable form
        std::string tmp="";
        for(unsigned int n=0; n<quarks.size(); n++){
            tmp+=(num_coeff[n]<0. ? " - " : " + ");
            tmp+=std::to_string(fabs(num_coeff[n]))+" * "+sym_coeff[n]+" * [ ";
            for(unsigned int s=0; s<numprops; s++){
                tmp+=(props_signs[s+n*numprops]>0. ? "+ " : "- ");
                for(unsigned int p=0; p<props[s+n*numprops].dim(); p++){
                    tmp+=props[s+n*numprops][p]+"_"+props_idcs[s+n*numprops][p]+" ";
                }
            }
            tmp+="]\n";
        }
        os << tmp << std::endl;
    }
    else if(mode.compare("laph")==0){
        //print code based on laph in order to compute contractions
        unsigned int numfacts=props[0].dim();
        indent="";
        os << std::endl << "//compute sink blocks and diagrams:\n{\n";
        indent+="\t";
        os << indent+"int tf;\n";
        os << indent+"for (tf=0; tf<lt; tf++){" << std::endl;
        indent+="\t";
        os << indent+"COMPLEX sum=0.;\n";
        for(unsigned int i=0; i<(numfacts*2); i++){
            os << indent << "unsigned int n" << i << ";\n";
            os << indent << "for(n" << i << "=0; n" << i << "<LAPH; n" << i << "++){\n";
            indent+="\t";
        }
        unsigned int massid,spin1id,spin2id,n1id,n2id;
        for(unsigned int n=0; n<num_coeff.size(); n++){
            if(n!=0) os << std::endl;
            if(num_coeff[n]>0.) os << indent << "sum += ";
            else os << indent << "sum -= ";
            os << std::to_string(fabs(num_coeff[n])) << " * ";
            for(unsigned int p=0; p<numfacts; p+=3){
                os << "(~vvv[n" << p << "][n" << p+1 << "][tf][n" << p+2 << "])" << " * ";
            }
            for(unsigned int p=numfacts; p<(2*numfacts); p+=3){
                os << "vvv[n" << p << "][n" << p+1 << "][ti][n" << p+2 << "]" << " * ";
            }
            os << "(" << std::endl;
            for(unsigned int p=0; p<numprops; p++){
                os << indent+"\t";
                std::string sign=(props_signs[p+n*numprops]>0. ? " + " : " - ");
                os << sign;
                for(unsigned int s=0; s<numfacts; s++){
                    get_indices_laph(props[p+numprops*n][s],props_idcs[p+numprops*n][s],massid,spin1id,spin2id,n1id,n2id,isolimit);
                    std::string tmp="vMv["+std::to_string(massid)+"][n"+std::to_string(n2id)+"]["+std::to_string(spin2id)+"][tf][n"+std::to_string(n1id)+"]["+std::to_string(spin1id)+"]";
                    if(s==0) os << tmp;
                    else os << " * " << tmp;
                }
                os << std::endl;
            }
            os << indent << ");" << std::endl;
        }
        os << std::endl;
        for(unsigned int i=0; i<(numfacts*2); i++){
            indent.erase((numfacts*2-i-1),1);
            os << indent << "} //end loop n" << (numfacts*2-i-1) << "\n";
        }
        
        //finalize:
        os << std::endl;
        os << indent+"int tf1= lt*mynode_dir[TUP]+tf;\n";
        os << indent+"int tdiff= (tf1-ti+nt)%nt;\n";
        os << indent+"P->bar["+std::to_string(operator_id)+"][tdiff]+= sum;\n";
        indent.erase(0,1);
        os << indent+"} //end loop ti, tf\n";
        indent.erase(0,1);
        os << indent+"}\n";
    }
    else if(mode.compare("laph2")==0){
        if(laph_sinks.size()==0) get_laph_sinks(mode);
        if(laph_sources.size()==0) get_laph_sources(mode);
        
        unsigned int numfacts=props[0].dim();
        indent="";
        
        //********************************************************
        //********************************************************
        //sink blocks and diagrams:
        //********************************************************
        //********************************************************
        os << std::endl << "//compute sink blocks and diagrams:\n{\n";
        indent+="\t";
        os << indent+"int tf, src;\n";
        os << indent+"NODE0_PRINTF(\"Computing contractions for operator %d\",operator_id)\n;";
        os << indent << "for(tf=0; tf<lt; tf++) for(src=0; src<nsrc; src++){" << std::endl;
        indent+="\t";
        //set up temporary storage:
        std::string wwwsummed("www1pt[tf][src]");
        if(operator_id>=0){
            wwwsummed="www1pt[tf][src]["+std::to_string(operator_id)+"]";
        }
        for(unsigned int i=0; i<numfacts; i++){
            os << indent << "unsigned int n" << i << ";\n";
            os << indent << "for(n" << i << "=0; n" << i << "<LAPH; n" << i << "++){\n";
            indent+="\t";
            
            wwwsummed+="[n"+std::to_string(i)+"]";
        }
        
        for(unsigned int n=0; n<num_coeff.size(); n++){
            if(n!=0) os << std::endl;
            if(num_coeff[n]>0.) os << indent << wwwsummed << " += ";
            else os << std::endl << indent << wwwsummed << " -= ";
            os << std::to_string(fabs(num_coeff[n])) << " * (" << std::endl;
            for(unsigned int p=0; p<numprops; p++){
                os << indent+"\t";
                std::string sign=(props_signs[p+n*numprops]>0. ? " + " : " - ");
                os << sign;
                unsigned int count=0;
                for(unsigned int s=0; s<numfacts; s+=3){
                    
                    if(s==0) os << laph_sinks[p+n*numprops][count];
                    else os << " * " << laph_sinks[p+n*numprops][count];
                    count++;
                }
                os << std::endl;
            }
            os << indent << ");" << std::endl;
        }
        os << std::endl;
        for(unsigned int i=0; i<numfacts; i++){
            indent.erase((numfacts-i-1),1);
            os << indent << "} //end loop n" << (numfacts-i-1) << "\n";
        }
        indent.erase(0,1);
        os << indent << "} //end loop tf, nsrc\n";
        indent.erase(0,1);
        os << "}\n";
    }
    
    return EXIT_SUCCESS;
}
//******************************************************************
//******************************************************************
//END quark_op
//******************************************************************
//******************************************************************
