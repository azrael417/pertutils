//
//  bartens.cpp
//  pertutils
//
//  Created by Thorsten Kurth on 24.06.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

#include "mathutils.hpp"
#include "pertutils.hpp"

BBTensor::BBTensor(const std::vector<dcomplex>& array, const std::vector< fourvec<int> >& sourcepositions, const std::vector<std::string>& ordering, const std::vector<std::string>& bartypes, const bool silent) : baryons(bartypes), spos(sourcepositions), numsources(static_cast<unsigned int>(sourcepositions.size())){
    unsigned int colcount=0,spincount=0,barcount=0,sourcecount=0;
    bool fail=false;
    
    std::vector<std::string> order(ordering);
    std::vector<unsigned int> modesizes;
    for(unsigned int i=0; i<static_cast<unsigned int>(order.size()); i++){
        if(order[i].find("baryon")==0){
            barcount++;
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
    if(3*barcount!=colcount){
        std::cerr << "BBTensor::BBtensor: error, your number of baryons is not equal three times the number of quarks!" << std::endl;
        fail=true;
    }
    if(barcount!=baryons.size()){
        std::cerr << "BBTensor::BBtensor: error, you have to specify exactly as many baryons in the spec list as present in the data!" << std::endl;
    }
    if(order.size()!=barcount+spincount+colcount+sourcecount){
        std::cerr << "BBTensor::BBtensor: warning, your string list also contains other unsupported entries! These will be removed now!" << std::endl;
        for(unsigned int i=0; i<static_cast<unsigned int>(order.size()); i++){
            if( (order[i].find("baryon")!=0) && (order[i].find("spin")!=0) && (order[i].find("color")!=0) && (order[i].find("source")!=0) ){
                order.erase(order.begin()+i);
            }
        }
    }
    
    //search whether all indices appear:
    fail*=!find(order,barcount,"baryon");
    fail*=!find(order,colcount,"color");
    fail*=!find(order,spincount,"spin");
    fail*=!find(order,sourcecount,"source");
    
    //setup:
    if(!fail){
        TTTensor tmp(array,modesizes);
        
        //group indices in order to obtain form: B(B_1,...,B_n|A_1,A_2,...A_n) with A_i=(spin_i,color_i,src_i):
        for(unsigned int b=0; b<barcount; b++){
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
            for(unsigned int i=(barcount+3*b+0); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "spin" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,barcount+3*b+0);
                    //move element of order also:
                    order.insert(order.begin()+barcount+3*b+0,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
            for(unsigned int i=(barcount+3*b+1); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "color" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,barcount+3*b+1);
                    //move element of order also:
                    order.insert(order.begin()+barcount+3*b+1,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
            for(unsigned int i=(barcount+3*b+2); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "source" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,barcount+3*b+2);
                    //move element of order also:
                    order.insert(order.begin()+barcount+3*b+2,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
        }
        data=tmp;
        
        //set quark-id vector:
        quark qrk;
        for(unsigned int b=barcount; b<(9*barcount); b+=3){
            qrk.spinid=b+0;
            qrk.colorid=b+1;
            qrk.sourceid=b+2;
            quarks.push_back(qrk);
        }
        for(unsigned int b=0; b<barcount; b++){
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
            for(unsigned int b=0; b<barcount; b++){
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


BBTensor::BBTensor(const TTTensor& atens, const std::vector< fourvec<int> >& sourcepositions, const std::vector<std::string>& ordering, const std::vector<std::string>& bartypes, const bool silent) : baryons(bartypes), spos(sourcepositions), numsources(static_cast<unsigned int>(sourcepositions.size())){
    unsigned int colcount=0,spincount=0,barcount=0,sourcecount=0;
    bool fail=false;
    
    std::vector<std::string> order(ordering);
    std::vector<unsigned int> modesizes;
    for(unsigned int i=0; i<static_cast<unsigned int>(order.size()); i++){
        if(order[i].find("baryon")==0){
            barcount++;
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
    if(3*barcount!=colcount){
        std::cerr << "BBTensor::BBtensor: error, your number of baryons is not equal three times the number of quarks!" << std::endl;
        fail=true;
    }
    if(barcount!=baryons.size()){
        std::cerr << "BBTensor::BBtensor: error, you have to specify exactly as many baryons in the spec list as present in the data!" << std::endl;
    }
    if(order.size()!=barcount+spincount+colcount+sourcecount){
        std::cerr << "BBTensor::BBtensor: warning, your string list also contains other unsupported entries! These will be removed now!" << std::endl;
        for(unsigned int i=0; i<static_cast<unsigned int>(order.size()); i++){
            if( (order[i].find("baryon")!=0) && (order[i].find("spin")!=0) && (order[i].find("color")!=0) && (order[i].find("source")!=0) ){
                order.erase(order.begin()+i);
            }
        }
    }
    
    //search whether all indices appear:
    fail*=!find(order,barcount,"baryon");
    fail*=!find(order,colcount,"color");
    fail*=!find(order,spincount,"spin");
    fail*=!find(order,sourcecount,"source");
    
    //setup:
    if(!fail){
        TTTensor tmp(atens);
        
        //group indices in order to obtain form: B(B_1,...,B_n|A_1,A_2,...A_n) with A_i=(spin_i,color_i,src_i):
        for(unsigned int b=0; b<barcount; b++){
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
            for(unsigned int i=(barcount+3*b+0); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "spin" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,barcount+3*b+0);
                    //move element of order also:
                    order.insert(order.begin()+barcount+3*b+0,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
            for(unsigned int i=(barcount+3*b+1); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "color" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,barcount+3*b+1);
                    //move element of order also:
                    order.insert(order.begin()+barcount+3*b+1,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
            for(unsigned int i=(barcount+3*b+2); i<static_cast<unsigned int>(order.size()); i++){
                searchstring.clear();
                searchstring.str("");
                searchstring << "source" << b << "\0";
                if(order[i].compare(searchstring.str())==0){
                    tmp=move_block(tmp,i,barcount+3*b+2);
                    //move element of order also:
                    order.insert(order.begin()+barcount+3*b+2,order[i]);
                    order.erase(order.begin()+i+1);
                    continue;
                }
            }
        }
        data=tmp;
        
        //set quark-id vector:
        quark qrk;
        for(unsigned int b=barcount; b<(9*barcount); b+=3){
            qrk.spinid=b+0;
            qrk.colorid=b+1;
            qrk.sourceid=b+2;
            quarks.push_back(qrk);
        }
        for(unsigned int b=0; b<barcount; b++){
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
            for(unsigned int b=0; b<barcount; b++){
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


TTTensor BBTensor::extract_data()const{
    return data;
}


TTTensor dot(const BBTensor& t1, const BBTensor& t2){
    TTTensor result;
    
    //determine whether the number of internal quarks is the same:
    if(t1.baryons.size()!=t2.baryons.size()){
        std::cerr << "BBTensor::dot: error, partial contractions not yet implemented!" << std::endl;
        return result;
    }
    if(t1.numsources!=t2.numsources){
        std::cerr << "BBTensor:dot: error, the number of sources used should be the same for both tensors!" << std::endl;
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
    
    //if the number of baryons is less than the number of sources, then we should restrict the tensor accoordingly: if e.g. we want do compute 2 baryon correlation functions and the tensor contains 3 source positions, then we can select two different contractions, where all quarks come from either pos1 and pos2, pos2 and pos3, or pos1 and pos3:
    
    //set up array of indices which will be contracted:
    std::vector<unsigned int> idt1, idt2;
    for(unsigned int f=0; f<6; f++){
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
    result=dot(t1.data,idt1,t2.data,idt2);
    return result;
}


//void BBTensor::join(const BBTensor& rhs, const bool& asym){
//    //Compute Kronecker Product of tensors:
//    std::vector<TTTensor> newdata;
//    unsigned int oldsize=data[0].get_dim();
//    for(unsigned int t=0; t<nt; t++){
//        TTTensor ltr(kron(data[t],rhs.data[t]));
//        newdata.push_back(ltr);
//    }
//    data=newdata;
//    newdata.clear();
//
//    //anti-symmetrize in all old and new quark indices:
//    if(asym){
//        unsigned int index1, index2;
//        for(unsigned int t=0; t<nt; t++){
//            for(unsigned int f=0; f<6; f++){
//                for(unsigned int i1=0; i1<numquarksperflavour[f]; i1++){
//                    if( (quarks[i1].isexternal) || (quarks[i1].flavourid!=f)) continue;
//                    for(unsigned int i2=0; i2<rhs.numquarksperflavour[f]; i2++){
//                        if( !(rhs.quarks[i2].isexternal) && (rhs.quarks[i2].flavourid==f) ){
//                            index1=quarks[i1].spinid;
//                            index2=rhs.quarks[i2].spinid+oldsize;
//                            TTTensor tmp(swap(data[t],index1,index2));
//                            index1=quarks[i1].colorid;
//                            index2=rhs.quarks[i2].colorid+oldsize;
//                            tmp=swap(tmp,index1,index2);
//                            tmp*=-1.;
//                            TTTensor tmp2(data[t]+tmp);
//                            data.erase(data.begin()+t);
//                            data.insert(data.begin()+t,tmp2);
//                        }
//                    }
//                }
//            }
//        }
//    }
//    
//    //Join vectors of source positions:
//    std::vector<quark> newquarks;
//    newquarks.reserve(quarks.size()+rhs.quarks.size());
//    newquarks.insert(newquarks.end(),quarks.begin(),quarks.end());
//    newquarks.insert(newquarks.end(),rhs.quarks.begin(),rhs.quarks.end());
//    //add offset to new quark locations:
//    unsigned int offset=static_cast<unsigned int>(quarks.size());
//    for(unsigned int i=0; i<rhs.quarks.size(); i++){
//        newquarks[offset+i].spinid+=oldsize;
//        newquarks[offset+i].colorid+=oldsize;
//    }
//    quarks=newquarks;
//    for(unsigned int i=0; i<rhs.quarks.size(); i++){
//        numquarksperflavour[rhs.quarks[i].flavourid]++;
//    }
//}
//
////contract all internal quarks:
//
////contract all external quarks. internal quarks have to be contracted before:
//std::vector<dcomplex> project(const BBTensor& t1, const TTTensor& proj){
//    std::vector<dcomplex> result;
//    
//    //collect external quarks:
//    std::vector<unsigned int> extquarks,intquarks;
//    unsigned int nqfcont[6];
//    for(unsigned int f=0; f<6; f++) nqfcont[f]=0;
//    for(unsigned int i=0; i<static_cast<unsigned int>(t1.quarks.size()); i++){
//        if(t1.quarks[i].isexternal) extquarks.push_back(i);
//        else{
//            nqfcont[t1.quarks[i].flavourid]++;
//            intquarks.push_back(i);
//        }
//    }
//    if(intquarks.size()!=0){
//        std::cerr << "BBTensor::dot: please contract the internal quarks first!" << std::endl;
//        return result;
//    }
//    if( extquarks.size()!=proj.get_dim() ){
//        std::cerr << "BBTensor::dot: error, the number of external quarks does not match or the dimension of your projection tensor!" << std::endl;
//        return result;
//    }
//    
//    //set up contraction indices:
//    std::vector<unsigned int> ind;
//    for(unsigned int i=0; i<static_cast<unsigned int>(extquarks.size()); i++){
//        ind.push_back(i);
//    }
//    
//    //perform contractions:
//    std::vector<unsigned int> dum(1);
//    dum[0]=0;
//    for(unsigned int t=0; t<t1.nt; t++){
//        TTTensor res(dot(t1.data[t],ind,proj,ind));
//        result.push_back(res(dum));
//    }
//    return result;
//}

