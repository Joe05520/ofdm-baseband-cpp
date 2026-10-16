#include "multipath.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace hw2::multipath {
Vector transmit(const Bits& bits, Modulation mode, std::size_t cp) {
    const auto block=N*bits_per_symbol(mode);
    if(cp>N || bits.empty() || bits.size()%block) throw std::invalid_argument("invalid CP or bits");
    Vector tx;tx.reserve(bits.size()/block*(N+cp));
    for(std::size_t f=0;f<bits.size()/block;++f) {
        Vector x(N);
        for(std::size_t k=0;k<N;++k)x[k]=map_symbol(&bits[f*block+k*bits_per_symbol(mode)],mode);
        const auto t=hw1::ifft(x);
        tx.insert(tx.end(),t.end()-cp,t.end());tx.insert(tx.end(),t.begin(),t.end());
    }
    return tx;
}
Vector filter(const Vector& tx,const Vector& h) {
    if(tx.empty() || h.empty())throw std::invalid_argument("empty FIR input");
    Vector rx(tx.size()+h.size()-1);
    for(std::size_t l=0;l<h.size();++l)if(h[l]!=Complex{})
        for(std::size_t n=0;n<tx.size();++n)rx[n+l]+=h[l]*tx[n];
    return rx;
}
Vector rayleigh(std::size_t step,std::uint64_t seed) {
    if(step==0 || 3*step>=N)throw std::invalid_argument("invalid delay profile");
    double p[4],sum=0;
    for(int l=0;l<4;++l){p[l]=std::pow(10.,-3.*l/10.);sum+=p[l];}
    Gaussian rng(seed);Vector h(3*step+1);
    // Normalize expected PDP energy only. Do NOT normalize individual fades.
    for(int l=0;l<4;++l)h[l*step]=rng.sample(p[l]/sum);
    return h;
}
Vector response(const Vector& h) {
    if(h.empty() || h.size()>N)throw std::invalid_argument("invalid channel length");
    Vector p(N);std::copy(h.begin(),h.end(),p.begin());return hw1::fft(p);
}
Vector window_fft(const Vector& rx,std::size_t frame,std::size_t cp) {
    const auto start=frame*(N+cp)+cp;
    if(cp>N || start+N>rx.size())throw std::invalid_argument("FFT window out of bounds");
    return hw1::fft(Vector(rx.begin()+start,rx.begin()+start+N));
}
Vector zero_force(const Vector& y,const Vector& H) {
    if(y.size()!=N || H.size()!=N)throw std::invalid_argument("ZF size mismatch");
    Vector z(N);
    for(std::size_t k=0;k<N;++k){
        if(std::norm(H[k])==0)throw std::runtime_error("exact spectral null cannot be zero forced");
        z[k]=y[k]/H[k];
    }
    return z;
}
std::uint64_t errors(const Vector& y,const Bits& bits,Modulation mode,std::size_t frame) {
    const auto bps=bits_per_symbol(mode);const auto block=N*bps;
    if(y.size()!=N || (frame+1)*block>bits.size())throw std::invalid_argument("decision size mismatch");
    std::uint64_t e=0;
    for(std::size_t k=0;k<N;++k){const auto b=demap_symbol(y[k],mode);
        for(unsigned j=0;j<bps;++j)e+=b[j]!=bits[frame*block+k*bps+j];}
    return e;
}
}
