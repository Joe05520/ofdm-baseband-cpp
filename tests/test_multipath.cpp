#include "multipath.hpp"
#include <iostream>
#include <stdexcept>
using namespace hw2;
namespace mp=hw2::multipath;
int main(){
    std::size_t checks=0;double worst=0;
    auto check=[&](bool ok){++checks;if(!ok)throw std::runtime_error("multipath check failed");};
    for(auto mode:{Modulation::QPSK,Modulation::QAM16})for(auto step:{4u,8u}) {
        auto bits=random_bits(6*N*bits_per_symbol(mode),1234);auto h=mp::rayleigh(step,9876);auto H=mp::response(h);
        for(auto cp:{0u,8u,12u,16u,24u,32u}){
            auto tx=mp::transmit(bits,mode,cp),rx=mp::filter(tx,h);double residual=0;
            // Direct convolution at samples around a symbol boundary.
            for(std::size_t n=N+cp-3;n<N+cp+30;++n){Complex ref{};
                for(std::size_t l=0;l<h.size();++l)if(n>=l && n-l<tx.size())ref+=h[l]*tx[n-l];
                check(std::abs(rx[n]-ref)<1e-14);}
            for(std::size_t f=0;f<6;++f){auto y=mp::window_fft(rx,f,cp),z=mp::zero_force(y,H);
                for(std::size_t k=0;k<N;++k){auto x=map_symbol(&bits[(f*N+k)*bits_per_symbol(mode)],mode);
                    residual+=std::norm(y[k]-H[k]*x);if(cp>=3*step)worst=std::max(worst,std::abs(z[k]-x));}
                if(cp>=3*step)check(mp::errors(z,bits,mode,f)==0);}
            check(cp>=3*step ? residual<1e-24 : residual>1e-3);
        }
    }
    // Variable-CP path must match the original transmitter exactly at G=16.
    auto bits=random_bits(8*N*4,42);check(mp::transmit(bits,Modulation::QAM16,16)==modulate(bits,Modulation::QAM16));
    check(mp::filter(Vector{1.,2.,3.},Vector{1.,0.,2.})==Vector({1.,2.,5.,4.,6.}));
    check(mp::filter(Vector{1.,2.,3.},Vector{1.})==Vector({1.,2.,3.}));
    check(worst<1e-12);
    std::cout<<"PASS "<<checks<<" checks; max noise-free ZF error "<<worst<<"\n";
}
