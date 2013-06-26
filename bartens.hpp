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
    BBTensor(const std::vector<dcomplex>& array, const unsigned int& nsources, const std::vector<std::string>& ordering) : numsources(nsources){
        unsigned int colcount=0,spincount=0,barcount=0,sourcecount=0;
        bool fail=false;
        
        std::vector<unsigned int> modesizes;
        for(unsigned int i=0; i<static_cast<unsigned int>(ordering.size()); i++){
            if(ordering[i].compare("baryon")==0){
                barcount++;
                modesizes.push_back(4);
            }
            else if(ordering[i].compare("color")==0){
                colcount++;
                modesizes.push_back(3);
            }
            else if(ordering[i].compare("spin")==0){
                spincount++;
                modesizes.push_back(4);
            }
            else if(ordering[i].compare("source")==0){
                sourcecount++;
                modesizes.push_back(numsources);
            }
        }
        if( (colcount!=spincount) || (colcount!=sourcecount) ){
            std::cerr << "BBTensor::BBtensor: error, you did not specify all the spin/color/source combinations for the quark sources!" << std::endl;
            fail=true;
        }
        if(barcount!=colcount*3){
            std::cerr << "BBTensor::BBtensor: error, your number of baryons is not equal three times the number of quarks!" << std::endl;
            fail=true;
        }
        
        //search whether all indices appear:
        fail*=!find(ordering,barcount,"baryon");
        fail*=!find(ordering,colcount,"color");
        fail*=!find(ordering,spincount,"spin");
        fail*=!find(ordering,sourcecount,"source");
        
        //setup:
        if(!fail){
            TTTensor tmp(array,modesizes);
            
            //group indices in order to obtain form: B(bar1,...,barn|A_1,A_2,...A_n) with A_i=(s_i,c_i,src_i):
            for(unsigned int b=0; b<barcount; b++){
                std::stringstream searchstring;
                for(unsigned int i=0; i<static_cast<unsigned int>(ordering.size()); i++){
                    searchstring.clear();
                    searchstring.str("");
                    searchstring << "baryon" << b << "\0";
                    if(ordering[i].compare(searchstring.str())==0){
                        tmp=move_block(tmp,i,b);
                        //move element of ordering also:
                        continue;
                    }
                }
            }
        }
    };
    
    
   
    
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
