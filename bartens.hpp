//
//  bartens.h
//  pertutils
//
//  Created by Thorsten Kurth on 24.06.13.
//  Copyright (c) 2013 Bergische Universität Wuppertal. All rights reserved.
//

#ifndef _BARTENS
#define _BARTENS
const unsigned int bartens_blocksize=256*27;
const unsigned int bartens_modesizes[7]={4,4,4,4,3,3,3};

//quark flavours:
enum flavour{
    UP=0,
    DOWN=1,
    STRANGE=2,
    CHARM=3,
    BOTTOM=4,
    TOP=5
};

//struct which identifies a quark inside the baryon:
typedef struct{
    flavour flavourid;
    bool isexternal;
    unsigned int spinid;
    unsigned int colorid;
    unsigned int sourceid;
}quark;


static bool find(const std::vector<std::string>& ordering, const unsigned int& count, const std::string searchstring){
    bool found=true,tmpfound,fail=false;
    for(unsigned int b=0; b<count; b++){
        tmpfound=false;
        for(unsigned int i=0; i<static_cast<unsigned int>(ordering.size()); i++){
            std::stringstream srch;
            srch << searchstring << b << "\0";
            if(ordering[i].compare(srch.str())==0){
                tmpfound=true;
                continue;
            }
        }
        found*=tmpfound;
    }
    if(!found){
        std::cerr << "BBTensor::BBtensor: error, you did not specify all " << searchstring << " indices correctly!" << std::endl;
        fail=true;
    }
    return fail;
}


//baryon tensor class
class BBTensor{
    
private:
    TTTensor data;
    std::vector<quark> quarks;
    unsigned int numsources;
    unsigned int numquarksperflavour[6];
    
public:
    
    //constructors:
    BBTensor(const std::vector<dcomplex>& array, const unsigned int& nsources, const std::vector<std::string>& ordering);
    
    ~BBTensor(){
        quarks.clear();
    }
    
    //operations:
    void join(const BBTensor& rhs, const bool& asym=false);
    
    //friend functions: this dot product exclusively acts on internal indices:
    friend BBTensor dot(const BBTensor& t1, const BBTensor& t2);
    //friend functions: this dot product exclusively acts on external indices:
    friend std::vector<dcomplex> project(const BBTensor& t1, const TTTensor& proj);
};

#endif
