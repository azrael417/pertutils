//
//  bartens.cpp
//  pertutils
//
//  Created by Thorsten Kurth on 24.06.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

#include "mathutils.hpp"
#include "pertutils.hpp"

BBTensor::BBTensor(const std::vector<dcomplex>& array, const std::vector< fourvec<int> >& sourcepositions, const std::vector<std::string>& ordering, const std::vector<std::string>& bartypes, const bool silent) : baryons(bartypes), spos(sourcepositions), numsources(static_cast<unsigned int>(sourcepositions.size())), numbaryons((0)){
    unsigned int colcount=0,spincount=0,sourcecount=0;
    bool fail=false;
    
    std::vector<std::string> order(ordering);
    std::vector<unsigned int> modesizes;
    for(unsigned int i=0; i<static_cast<unsigned int>(order.size()); i++){
        if(order[i].find("baryon")==0){
            numbaryons++;
            modesizes.push_back(4);
        }
        else if(order[i].find("color")==0){
            colcount++;
            modesizes.push_back(3);
        }
        else if(order[i].find("spin")==0){
            spincount++;
            modesizes.push_back(4);
        }
        else if(order[i].find("source")==0){
            sourcecount++;
            modesizes.push_back(numsources);
        }
    }
    if( (colcount!=spincount) || (colcount!=sourcecount) ){
        std::cerr << "BBTensor::BBtensor: error, you did not specify all the spin/color/source combinations for the quark sources!" << std::endl;
        fail=true;
    }
    if(3*numbaryons!=colcount){
        std::cerr << "BBTensor::BBtensor: error, your number of baryons is not equal three times the number of quarks!" << std::endl;
        fail=true;
    }
    if(numbaryons!=baryons.size()){
        std::cerr << "BBTensor::BBtensor: error, you have to specify exactly as many baryons in the spec list as present in the data!" << std::endl;
    }
    if(order.size()!=numbaryons+spincount+colcount+sourcecount){
        std::cerr << "BBTensor::BBtensor: warning, your string list also contains other unsupported entries! These will be removed now!" << std::endl;
        for(unsigned int i=0; i<static_cast<unsigned int>(order.size()); i++){
            if( (order[i].find("baryon")!=0) && (order[i].find("spin")!=0) && (order[i].find("color")!=0) && (order[i].find("source")!=0) ){
                order.erase(order.begin()+i);
            }
        }
    }
    
    //search whether all indices appear:
    fail*=!find(order,numbaryons,"baryon");
    fail*=!find(order,colcount,"color");
    fail*=!find(order,spincount,"spin");
    fail*=!find(order,sourcecount,"source");
    
    //setup:
    if(!fail){
        TTTensor tmp(array,modesizes);
        
        //group indices in order to obtain form: B(B_1,...,B_n|A_1,A_2,...A_n) with A_i=(spin_i,color_i,src_i):
        for(unsigned int b=0; b<numbaryons; b++){
            std::stringstream searchstring;
            for(unsigned int i=b; i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "baryon" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,b);
                    //move element of order also:
                    order.insert(order.begin()+b,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
        }
        for(unsigned int b=0; b<spincount; b++){
            std::stringstream searchstring;
            for(unsigned int i=(numbaryons+3*b+0); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "spin" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,numbaryons+3*b+0);
                    //move element of order also:
                    order.insert(order.begin()+numbaryons+3*b+0,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
            for(unsigned int i=(numbaryons+3*b+1); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "color" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,numbaryons+3*b+1);
                    //move element of order also:
                    order.insert(order.begin()+numbaryons+3*b+1,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
            for(unsigned int i=(numbaryons+3*b+2); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "source" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,numbaryons+3*b+2);
                    //move element of order also:
                    order.insert(order.begin()+numbaryons+3*b+2,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
        }
        data=tmp;
        
        //set quark-id vector:
        quark qrk;
        for(unsigned int b=numbaryons; b<(9*numbaryons); b+=3){
            qrk.spinid=b+0;
            qrk.colorid=b+1;
            qrk.sourceid=b+2;
            quarks.push_back(qrk);
        }
        for(unsigned int b=0; b<numbaryons; b++){
            if(baryons[b].compare("neutron")==0){
                quarks[0+3*b].flavourid=DOWN;
                quarks[1+3*b].flavourid=DOWN;
                quarks[2+3*b].flavourid=UP;
            }
            if(baryons[b].compare("proton")==0){
                quarks[0+3*b].flavourid=UP;
                quarks[1+3*b].flavourid=UP;
                quarks[2+3*b].flavourid=DOWN;
            }
        }
        
        if(!silent){
            std::cout << "The following quark content has been specified:" << std::endl;
            for(unsigned int b=0; b<numbaryons; b++){
                std::cout << baryons[b] << b << ":" << std::endl;
                for(unsigned int i=0; i<3; i++){
                    std::cout << "\t " << return_flavour(quarks[i+3*b]) << i << ":" << std::endl;
                    std::cout << "\t\t spin: " << quarks[i+3*b].spinid << std::endl;
                    std::cout << "\t\t color: " << quarks[i+3*b].colorid << std::endl;
                    std::cout << "\t\t source: " << quarks[i+3*b].sourceid << std::endl;
                }
            }
        }
    }
};


BBTensor::BBTensor(const TTTensor& atens, const std::vector< fourvec<int> >& sourcepositions, const std::vector<std::string>& ordering, const std::vector<std::string>& bartypes, const bool silent) : baryons(bartypes), spos(sourcepositions), numsources(static_cast<unsigned int>(sourcepositions.size())), numbaryons(0){
    unsigned int colcount=0,spincount=0,sourcecount=0;
    bool fail=false;
    
    std::vector<std::string> order(ordering);
    std::vector<unsigned int> modesizes;
    for(unsigned int i=0; i<static_cast<unsigned int>(order.size()); i++){
        if(order[i].find("baryon")==0){
            numbaryons++;
            modesizes.push_back(4);
        }
        else if(order[i].find("color")==0){
            colcount++;
            modesizes.push_back(3);
        }
        else if(order[i].find("spin")==0){
            spincount++;
            modesizes.push_back(4);
        }
        else if(order[i].find("source")==0){
            sourcecount++;
            modesizes.push_back(numsources);
        }
    }
    if( (colcount!=spincount) || (colcount!=sourcecount) ){
        std::cerr << "BBTensor::BBtensor: error, you did not specify all the spin/color/source combinations for the quark sources!" << std::endl;
        fail=true;
    }
    if(3*numbaryons!=colcount){
        std::cerr << "BBTensor::BBtensor: error, your number of baryons is not equal three times the number of quarks!" << std::endl;
        fail=true;
    }
    if(numbaryons!=baryons.size()){
        std::cerr << "BBTensor::BBtensor: error, you have to specify exactly as many baryons in the spec list as present in the data!" << std::endl;
    }
    if(order.size()!=numbaryons+spincount+colcount+sourcecount){
        std::cerr << "BBTensor::BBtensor: warning, your string list also contains other unsupported entries! These will be removed now!" << std::endl;
        for(unsigned int i=0; i<static_cast<unsigned int>(order.size()); i++){
            if( (order[i].find("baryon")!=0) && (order[i].find("spin")!=0) && (order[i].find("color")!=0) && (order[i].find("source")!=0) ){
                order.erase(order.begin()+i);
            }
        }
    }
    
    //search whether all indices appear:
    fail*=!find(order,numbaryons,"baryon");
    fail*=!find(order,colcount,"color");
    fail*=!find(order,spincount,"spin");
    fail*=!find(order,sourcecount,"source");
    
    //setup:
    if(!fail){
        TTTensor tmp(atens);
        
        //group indices in order to obtain form: B(B_1,...,B_n|A_1,A_2,...A_n) with A_i=(spin_i,color_i,src_i):
        for(unsigned int b=0; b<numbaryons; b++){
            std::stringstream searchstring;
            for(unsigned int i=b; i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "baryon" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,b);
                    //move element of order also:
                    order.insert(order.begin()+b,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
        }
        for(unsigned int b=0; b<spincount; b++){
            std::stringstream searchstring;
            for(unsigned int i=(numbaryons+3*b+0); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "spin" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,numbaryons+3*b+0);
                    //move element of order also:
                    order.insert(order.begin()+numbaryons+3*b+0,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
            for(unsigned int i=(numbaryons+3*b+1); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "color" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,numbaryons+3*b+1);
                    //move element of order also:
                    order.insert(order.begin()+numbaryons+3*b+1,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
            for(unsigned int i=(numbaryons+3*b+2); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "source" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,numbaryons+3*b+2);
                    //move element of order also:
                    order.insert(order.begin()+numbaryons+3*b+2,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
        }
        data=tmp;
        
        //set quark-id vector:
        quark qrk;
        for(unsigned int b=numbaryons; b<(9*numbaryons); b+=3){
            qrk.spinid=b+0;
            qrk.colorid=b+1;
            qrk.sourceid=b+2;
            quarks.push_back(qrk);
        }
        for(unsigned int b=0; b<numbaryons; b++){
            if(baryons[b].compare("neutron")==0){
                quarks[0+3*b].flavourid=DOWN;
                quarks[1+3*b].flavourid=DOWN;
                quarks[2+3*b].flavourid=UP;
            }
            if(baryons[b].compare("proton")==0){
                quarks[0+3*b].flavourid=UP;
                quarks[1+3*b].flavourid=UP;
                quarks[2+3*b].flavourid=DOWN;
            }
        }
        
        if(!silent){
            std::cout << "The following quark content has been specified:" << std::endl;
            for(unsigned int b=0; b<numbaryons; b++){
                std::cout << baryons[b] << b << ":" << std::endl;
                for(unsigned int i=0; i<3; i++){
                    std::cout << "\t " << return_flavour(quarks[i+3*b]) << i << ":" << std::endl;
                    std::cout << "\t\t spin: " << quarks[i+3*b].spinid << std::endl;
                    std::cout << "\t\t color: " << quarks[i+3*b].colorid << std::endl;
                    std::cout << "\t\t source: " << quarks[i+3*b].sourceid << std::endl;
                }
            }
        }
    }
};

//API functions:
TTTensor BBTensor::get_data()const{
    return data;
}

fourvec<int> BBTensor::get_spos(unsigned int sourceid)const{
    if(sourceid<numsources){
        return spos[sourceid];
    }
    else return fourvec<int>(-1,-1,-1,-1);
}

std::vector<std::string> BBTensor::get_baryons()const{
    return baryons;
}

unsigned int BBTensor::get_numbaryons()const{
    return numbaryons;
}

void BBTensor::print_info(const bool verbose)const{
    std::cout << "The following quark content has been specified:" << std::endl;
    for(unsigned int b=0; b<numbaryons; b++){
        std::cout << baryons[b] << b << ":" << std::endl;
        for(unsigned int i=0; i<3; i++){
            std::cout << "\t " << return_flavour(quarks[i+3*b]) << i << ":" << std::endl;
            std::cout << "\t\t spin: " << quarks[i+3*b].spinid << std::endl;
            std::cout << "\t\t color: " << quarks[i+3*b].colorid << std::endl;
            std::cout << "\t\t source: " << quarks[i+3*b].sourceid << std::endl;
        }
    }
    std::cout << "The source positions are:" << std::endl;
    for(unsigned int s=0; s<numsources; s++){
        std::cout << "\t\t" << spos[s] << std::endl;
    }
    if(verbose){
        std::cout << "The Data Layout is:" << std::endl;
        data.print_info();
    }
}

//friend functions:
BBTensor extract_sources(const BBTensor& a, const std::vector<bool>& idvec){
    BBTensor result(a);
    
    if(idvec.size()!=a.numsources){
        return a;
    }
    
    unsigned int numsnew=0;
    std::vector< fourvec<int> > srcspos;
    for(unsigned int s=0; s<a.numsources; s++){
        if(idvec[s]){
            numsnew++;
            srcspos.push_back(a.spos[s]);
        }
    }
    
    if(numsnew==0){
        std::cerr << "extract_sources: error, you cannot drop all sources!" << std::endl;
        return a;
    }
    else if(numsnew==a.numsources) return a;
    
    TTTensor tmp(a.data);
    for(unsigned int q=0; q<static_cast<unsigned int>(a.quarks.size()); q++){
        tmp=extract(tmp,a.quarks[q].sourceid,idvec);
    }
    result.numsources=numsnew;
    result.data=tmp;
    result.spos=srcspos;
    
    return result;
}

//dot product:
TTTensor dot(const BBTensor& t1, const BBTensor& t2){
    TTTensor result;
    
    //determine whether the number of internal quarks is the same:
    if(t1.numbaryons!=t2.numbaryons){
        std::cerr << "BBTensor::dot: error, partial contractions not yet implemented!" << std::endl;
        return result;
    }
    unsigned int nbary=static_cast<unsigned int>(t1.baryons.size());
    for(unsigned int b=0; b<nbary; b++){
        if(t1.baryons[b].compare(t2.baryons[b])!=0){
            std::cerr << "Warning, the order of baryons is not the same!" << std::endl;
        }
    }

    //determine the quark flavours which will be contracted
    std::vector<unsigned int> intquarks1,intquarks2;
    unsigned int nqfcont1[6], nqfcont2[6];
    for(unsigned int f=0; f<6; f++){
        nqfcont1[f]=0;
        nqfcont2[f]=0;
    }
    for(unsigned int i=0; i<static_cast<unsigned int>(t1.quarks.size()); i++){
        nqfcont1[t1.quarks[i].flavourid]++;
        intquarks1.push_back(i);
    }
    for(unsigned int i=0; i<static_cast<unsigned int>(t2.quarks.size()); i++){
        intquarks2.push_back(i);
        nqfcont2[t2.quarks[i].flavourid]++;
    }
    if(intquarks1.size()!=intquarks2.size()){
        std::cerr << "BBTensor::dot: integrity error, the number of quarks which are to be contracted does not match!" << std::endl;
        return result;
    }
    for(unsigned int f=0; f<6; f++){
        if(nqfcont1[f]!=nqfcont2[f]){
            std::cerr << "BBTensor::dot: integrity error, the number of quarks with the same flavour and which will be contracted do not match!" << std::endl;
            return result;
        }
    }
    
    //set up array of indices which will be contracted:
    std::vector<unsigned int> idt1, idt2;
    for(int f=0; f<6; f++){
        for(unsigned int i=0; i<static_cast<unsigned int>(intquarks1.size()); i++){
            if(t1.quarks[intquarks1[i]].flavourid==f){
                idt1.push_back(t1.quarks[intquarks1[i]].spinid);
                idt1.push_back(t1.quarks[intquarks1[i]].colorid);
                idt1.push_back(t1.quarks[intquarks1[i]].sourceid);
                continue;
            }
        }
        for(unsigned int i=0; i<static_cast<unsigned int>(intquarks2.size()); i++){
            if(t2.quarks[intquarks2[i]].flavourid==f){
                idt2.push_back(t2.quarks[intquarks2[i]].spinid);
                idt2.push_back(t2.quarks[intquarks2[i]].colorid);
                idt2.push_back(t2.quarks[intquarks2[i]].sourceid);
                continue;
            }
        }
    }
    
    //perform actual contraction:
    if(t1.numsources==t2.numsources){
        result=dot(t1.data,idt1,t2.data,idt2);
    }
    else if(t1.numsources!=t2.numsources){
        //reduction necessary:
        unsigned int nsmin=min(t1.numsources,t2.numsources);
        unsigned int nsmax=max(t1.numsources,t2.numsources);
        
        //draw all possible combinations of nsmin sources from nsmax samples:
        unsigned int combcount=static_cast<unsigned int>(bico(nsmax, nsmin));
        combination comb(nsmax,nsmin);
        std::vector<bool> combi(comb.return_combination());
        
        if(t1.numsources>t2.numsources){
            BBTensor tmpbb(extract_sources(t1,combi));
            
            result=dot(tmpbb.data,idt1,t2.data,idt2);
            for(unsigned int i=1; i<combcount; i++){
                combi=comb.return_combination();
                tmpbb=extract_sources(t1,combi);
                result+=dot(tmpbb.data,idt1,t2.data,idt2);
            }
            result/=static_cast<double>(combcount);
        }
        else{
            BBTensor tmpbb(extract_sources(t2,combi));
            result=dot(t1.data,idt1,tmpbb.data,idt2);
            
            for(unsigned int i=1; i<combcount; i++){
                combi=comb.return_combination();
                tmpbb=extract_sources(t2,combi);
                result+=dot(t1.data,idt1,tmpbb.data,idt2);
            }
            result/=static_cast<double>(combcount);
        }
    }    
    return result;
}

//contract all in- and external quarks:
dcomplex projectdot(const BBTensor& t1, const BBTensor& t2, const TTTensor& poltens){
    TTTensor contraction(dot(t1,t2));
    TTTensor restens(dot(contraction,poltens));
    std::vector<unsigned int> index(1);
    index[0]=0;
    return restens(index);
}


//computes the kronecker product of two baryon blocks and performs an anti-symmetrization on demand:
BBTensor join_barblocks(const BBTensor& lhs, const BBTensor& rhs, const bool antisym){
    if(lhs.numsources!=rhs.numsources){
        std::cerr << "join_barblocks: error, the number of sources has to be the same!" << std::endl;
        return lhs;
    }
    
    std::vector<std::string> baryonslhs(lhs.baryons), baryonsrhs(rhs.baryons);
    std::vector<std::string> resbaryons(baryonslhs);
    resbaryons.insert(resbaryons.end(),baryonsrhs.begin(),baryonsrhs.end());
    unsigned int bclhs=static_cast<unsigned int>(baryonslhs.size());
    unsigned int bcrhs=static_cast<unsigned int>(baryonsrhs.size());
    
    //join quark vector and shift coordinates of second quark:
    std::vector<quark> quarkslhs(lhs.quarks), quarksrhs(rhs.quarks);
    unsigned int qclhs=static_cast<unsigned int>(quarkslhs.size());
    unsigned int qcrhs=static_cast<unsigned int>(quarksrhs.size());
    
    for(unsigned int q=0; q<quarkslhs.size(); q++){
        quarkslhs[q].spinid+=bcrhs;
        quarkslhs[q].colorid+=bcrhs;
        quarkslhs[q].sourceid+=bcrhs;
    }
    for(unsigned int q=0; q<quarksrhs.size(); q++){
        quarksrhs[q].spinid+=(bclhs+3*qclhs);
        quarksrhs[q].colorid+=(bclhs+3*qclhs);
        quarksrhs[q].sourceid+=(bclhs+3*qclhs);
    }
    std::vector<quark> resquarks(quarkslhs);
    resquarks.insert(resquarks.end(),quarksrhs.begin(),quarksrhs.end());
    
    //resulting tensor:
    TTTensor lhstens(lhs.data), rhstens(rhs.data);
    TTTensor restens(kron(lhstens,rhstens));
    
    //move external indices to the beginning of new tensor:
    for(unsigned int b=0; b<bcrhs; b++){
        restens=move_block(restens, b+bclhs+3*qclhs, bclhs+b);
    }
    
    restens.print_info();
    
    //if antisymmetrization is required, perform anti-symmetrization:
    if(antisym){
        TTTensor tmptens;
        for(unsigned int ql=0; ql<qclhs; ql++){
            for(unsigned int qr=0; qr<qcrhs; qr++){
                if(quarkslhs[qr].flavourid==quarksrhs[qr].flavourid){
                    std::cout << "ANTISYM " << ql << " " << qr << std::endl;
                    std::cout << "SPINSWAP " << quarkslhs[ql].spinid << " <-> " << quarksrhs[qr].spinid << std::endl;
                    tmptens=swap(restens,quarkslhs[ql].spinid, quarksrhs[qr].spinid);
                    std::cout << "COLSWAP " << quarkslhs[ql].colorid << " <-> " << quarksrhs[qr].colorid << std::endl;
                    tmptens=swap(tmptens,quarkslhs[ql].colorid, quarksrhs[qr].colorid);
                    std::cout << "SRCSWAP" << quarkslhs[ql].sourceid << " <-> " << quarksrhs[qr].sourceid << std::endl;
                    tmptens=swap(tmptens,quarkslhs[ql].sourceid, quarksrhs[qr].sourceid);
                    std::cout << "DIFF " << std::endl;
                    restens-=tmptens;
                }
            }
        }
    }
    
    //construct output tensor:
    BBTensor result;
    result.data=restens;
    result.quarks=resquarks;
    result.spos=lhs.spos;
    result.baryons=resbaryons;
    result.numbaryons=static_cast<unsigned int>(resbaryons.size());
    for(unsigned int f=0; f<6; f++) result.numquarksperflavour[f]=lhs.numquarksperflavour[f]+rhs.numquarksperflavour[f];
    
    return result;
}

//I/O:
int write_barblock(const std::string filename, const BBTensor& tens){
    std::ofstream output;
    
    output.open(filename.c_str(),std::ios_base::binary);
    if(!output.good()){
        std::cerr << "write_barblock: error while trying to open file " << filename << "!" << std::endl;
        return EXIT_FAILURE;
    }
    
    int result=write_barblock(output,tens);
    output.close();
    
    return result;
}

int write_barblock(std::ofstream& output, const BBTensor& tens){
    double* head;
    //baryon header:
    head=new double[1+tens.baryons.size()];
    head[0]=static_cast<double>(tens.baryons.size());
    for(unsigned int b=0; b<static_cast<unsigned int>(tens.baryons.size()); b++){
        if(tens.baryons[b].compare("neutron")==0) head[1+b]=0.;
        else if(tens.baryons[b].compare("proton")==0) head[1+b]=1.;
    }
    output.write(reinterpret_cast<char*>(head), (1+tens.baryons.size())*sizeof(double));
    delete [] head;
    
    //source header:
    head=new double[1+4*tens.numsources];
    head[0]=static_cast<double>(tens.numsources);
    for(unsigned int i=0; i<tens.numsources; i++){
        for(unsigned int l=0; l<4; l++){
            head[1+l+4*i]=static_cast<double>(tens.spos[i][l]);
        }
    }
    output.write(reinterpret_cast<char*>(head), (1+4*tens.numsources)*sizeof(double));
    delete [] head;
    
    //quarks header:
    head=new double[1+4*tens.quarks.size()];
    head[0]=static_cast<double>(tens.quarks.size());
    for(unsigned int i=0; i<tens.quarks.size(); i++){
        head[1+0+4*i]=static_cast<double>(tens.quarks[i].flavourid);
        head[1+1+4*i]=static_cast<double>(tens.quarks[i].spinid);
        head[1+2+4*i]=static_cast<double>(tens.quarks[i].colorid);
        head[1+3+4*i]=static_cast<double>(tens.quarks[i].sourceid);
    }
    output.write(reinterpret_cast<char*>(head), (1+4*tens.quarks.size())*sizeof(double));
    delete [] head;
    
    //data content:
    write_tensor(output,tens.data);
        
    return EXIT_SUCCESS;
}

int read_barblock(const std::string filename, BBTensor& tens){
    std::ifstream input;
    input.open(filename.c_str(),std::ios_base::binary);
    if(!input.good()){
        std::cerr << "read_barblock: error while trying to open file " << filename << "!" << std::endl;
        return EXIT_FAILURE;
    }
    
    int result=read_barblock(input,tens);
    input.close();
    
    return result;
}

int read_barblock(std::ifstream& input, BBTensor& tens){
    double dummy, *head;
    input.read(reinterpret_cast<char*>(&dummy),sizeof(double));
    std::vector<std::string> baryons(static_cast<unsigned int>(dummy));
    head=new double[baryons.size()];
    input.read(reinterpret_cast<char*>(head),baryons.size()*sizeof(double));
    for(unsigned int b=0; b<baryons.size(); b++){
        if(fabs(head[b])<1.e-9) baryons[b]="neutron";
        else if(fabs(head[b]-1.)<1.e-9) baryons[b]="proton";
    }
    delete [] head;
    tens.baryons=baryons;
    tens.numbaryons=static_cast<unsigned int>(baryons.size());
    
    //source header:
    input.read(reinterpret_cast<char*>(&dummy),sizeof(double));
    tens.numsources=static_cast<unsigned int>(dummy);
    tens.spos.resize(tens.numsources);
    head=new double[4*tens.numsources];
    input.read(reinterpret_cast<char*>(head), (4*tens.numsources)*sizeof(double));
    for(unsigned int i=0; i<tens.numsources; i++){
        for(unsigned int l=0; l<4; l++){
            tens.spos[i][l]=static_cast<int>(head[l+4*i]);
        }
    }
    delete [] head;
    
    //quarks header:
    input.read(reinterpret_cast<char*>(&dummy),sizeof(double));
    tens.quarks.resize(static_cast<unsigned int>(dummy));
    head=new double[4*tens.quarks.size()];
    input.read(reinterpret_cast<char*>(head), (4*tens.quarks.size())*sizeof(double));
    for(unsigned int i=0; i<tens.quarks.size(); i++){
        tens.quarks[i].flavourid=static_cast<flavour>(head[0+4*i]);
        tens.quarks[i].spinid=static_cast<flavour>(head[1+4*i]);
        tens.quarks[i].colorid=static_cast<flavour>(head[2+4*i]);
        tens.quarks[i].sourceid=static_cast<flavour>(head[3+4*i]);
    }
    delete [] head;
    
    //set the numbers of quarks per flavour:
    for(unsigned int f=0; f<6; f++) tens.numquarksperflavour[f]=0;
    for(unsigned int q=0; q<tens.quarks.size(); q++){
        tens.numquarksperflavour[tens.quarks[q].flavourid]++;
    }
    
    //read data:
    read_tensor(input,tens.data);
    
    return EXIT_SUCCESS;
}

