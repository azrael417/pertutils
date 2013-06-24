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
    bool isextern;
    flavour flavourid;
    unsigned int spinid;
    unsigned int colorid;
}quark;

class BBTensor{
private:
    std::vector<TTTensor> data;
    std::vector<quark> quarks;
    unsigned int nt;
public:
    
    BBTensor(const std::vector<TTTensor>& ddata, const std::vector<quark>& qquarks) : data(ddata), quarks(qquarks), nt(static_cast<unsigned int>(data.size())){
        for(unsigned int t=0; t<nt; t++){
            if(data[t].get_dim()!=7){
                std::cerr << "BBTensor: warning, dimensions for t=" << t << " not equal to 7!" << std::endl;
            }
            else{
                for(unsigned int d=0; d<7; d++){
                    if(bartens_modesizes[d]!=data[t].get_nk(d)){
                        std::cerr << "BBTensor: warning, the mode sizes for t=" << t << " of mode " << d << " is not equal to " << bartens_modesizes[d] << "!" << std::endl;
                    }
                }
            }
        }
        if(quarks.size()!=4){
            std::cerr << "BBTensor: error, you did not fully specify which index corresponds to which quark flavour!" << std::endl;
        }
    };
    
    ~BBTensor(){
        data.clear();
        quarks.clear();
    }
    
};

#endif
