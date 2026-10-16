#include "ofdm.hpp"
#include <cmath>
#include <algorithm>

#include <limits>
#include <stdexcept>

namespace hw2 {
unsigned bits_per_symbol(Modulation m) { return m == Modulation::QPSK ? 2 : 4; }
std::string name(Modulation m) { return m == Modulation::QPSK ? "QPSK" : "16QAM"; }
Bits random_bits(std::size_t count, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    Bits bits(count);
    std::uint64_t word = 0;
    for (std::size_t i=0; i<count; ++i) {
        if (i%64 == 0) word = rng();
        bits[i] = (word >> (i%64)) & 1;
    }
    return bits;
}
Complex map_symbol(const std::uint8_t* b, Modulation m) {
    if (m == Modulation::QPSK)
        return Complex(1.0-2*b[0], 1.0-2*b[1]) / std::sqrt(2.0);
    constexpr double level[4] = {-3,-1,3,1}; // Gray: 00,01,10,11
    return Complex(level[2*b[0]+b[1]],level[2*b[2]+b[3]])/std::sqrt(10.0);
}
std::array<std::uint8_t,4> demap_symbol(Complex z, Modulation m) {
    if (m == Modulation::QPSK)
        return {std::uint8_t(z.real()<0),std::uint8_t(z.imag()<0),0,0};
    std::array<std::uint8_t,4> b{};
    const double v[2] = {z.real()*std::sqrt(10.0),z.imag()*std::sqrt(10.0)};
    for (unsigned a=0; a<2; ++a) {
        b[2*a] = v[a]>=0;
        b[2*a+1] = std::abs(v[a])<2;
    }
    return b;
}
Vector modulate(const Bits& bits, Modulation m) {
    const auto bps=bits_per_symbol(m), block=unsigned(N)*bps;
    if (bits.empty() || bits.size()%block) throw std::invalid_argument("bit count must be a positive multiple of 64*bits_per_symbol");
    Vector tx(bits.size()/block*L), freq(N);
    for (std::size_t f=0; f<bits.size()/block; ++f) {
        for (std::size_t k=0;k<N;++k) freq[k]=map_symbol(&bits[f*block+k*bps],m);
        const auto time=hw1::ifft(freq);
        for (std::size_t n=0;n<CP;++n) tx[f*L+n]=time[N-CP+n];
        for (std::size_t n=0;n<N;++n) tx[f*L+CP+n]=time[n];
    }
    return tx;
}
Vector add_delay(const Vector& tx,std::size_t delay) {
    if(delay>=L) throw std::invalid_argument("acquisition model requires delay in [0,79]");
    Vector out(delay+tx.size()+L); // Tail permits measuring errors even for a wrong timing estimate.
    std::copy(tx.begin(),tx.end(),out.begin()+delay);
    return out;
}
double Gaussian::uniform_open() {
    // 52 random mantissa bits; midpoint bins avoid log(0), and remain strictly below one.
    return (double(rng_()>>12)+0.5)*0x1p-52;
}
Complex Gaussian::sample(double variance) {
    if (variance<0 || !std::isfinite(variance)) throw std::invalid_argument("invalid noise variance");
    const double radius=std::sqrt(-variance*std::log(uniform_open()));
    const double angle=2*hw1::pi*uniform_open();
    return {radius*std::cos(angle),radius*std::sin(angle)};
}
Vector awgn(const Vector& tx,double variance,std::uint64_t seed,NoiseStats* stats) {
    Gaussian gen(seed); Vector rx(tx.size()); NoiseStats s;
    for(std::size_t i=0;i<tx.size();++i) {
        const auto w=gen.sample(variance); rx[i]=tx[i]+w;
        s.power+=std::norm(w);s.mean_i+=w.real();s.mean_q+=w.imag();s.iq+=w.real()*w.imag();
    }
    if(stats) {
        const double n=double(tx.size());
        if(n==0) throw std::invalid_argument("empty AWGN statistics");
        *stats={s.power/n,s.mean_i/n,s.mean_q/n,s.iq/n};
    }
    return rx;
}
Timing synchronize(const Vector& rx,std::size_t symbols) {
    if(symbols==0 || rx.size()<L || symbols>(rx.size()-(L-1))/L)
        throw std::invalid_argument("insufficient samples for all 80 timing candidates");
    Timing result; double best=-std::numeric_limits<double>::infinity();
    for(std::size_t d=0;d<L;++d) {
        Complex corr{};double energy=0;
        for(std::size_t s=0;s<symbols;++s) for(std::size_t n=0;n<CP;++n) {
            const auto a=rx[d+s*L+n],b=rx[d+s*L+n+N];
            // Lecture 02, p.48: Phi_d = sum r[d+k] conj(r[d+k+N]).
            corr+=a*std::conj(b);energy+=std::norm(b);
        }
        // Gamma_d = |Phi_d|^2 / P_d^2; P_d uses the delayed samples only.
        // This one-sided normalization is not generally bounded above by 1.
        result.metric[d]=energy>0 ? std::norm(corr/energy) : 0;
        if(result.metric[d]>best) {best=result.metric[d];result.cp_start=d;}
    }
    return result;
}
std::vector<Timing> synchronize_sweep(const Vector& rx,const std::vector<std::size_t>& counts) {
    if(counts.empty() || counts.front()==0 ||
       !std::is_sorted(counts.begin(),counts.end()) ||
       std::adjacent_find(counts.begin(),counts.end())!=counts.end() ||
       rx.size()<L || counts.back()>(rx.size()-(L-1))/L)
        throw std::invalid_argument("invalid timing sweep or insufficient samples");
    // Rolling C=CP windows evaluate exactly the lecture sum. Reuse them for
    // all M values; validate this acceleration against the direct function.
    const auto max_m=counts.back();
    const auto last_start=(max_m-1)*L+(L-1);
    const auto product_count=last_start+CP;
    Vector products(product_count);std::vector<double> energies(product_count);
    for(std::size_t i=0;i<product_count;++i) {
        products[i]=rx[i]*std::conj(rx[i+N]);energies[i]=std::norm(rx[i+N]);
    }
    std::array<Complex,L> corr{};std::array<double,L> power{};
    Complex window{};double window_power=0;
    for(std::size_t i=0;i<CP;++i){window+=products[i];window_power+=energies[i];}
    std::vector<Timing> result;std::size_t next=0;
    for(std::size_t start=0;start<=last_start;++start) {
        const auto d=start%L;corr[d]+=window;power[d]+=window_power;
        if(d==L-1 && (start/L+1)==counts[next]) {
            Timing t;double best=-std::numeric_limits<double>::infinity();
            for(std::size_t j=0;j<L;++j) {
                t.metric[j]=power[j]>0 ? std::norm(corr[j]/power[j]) : 0;
                if(t.metric[j]>best){best=t.metric[j];t.cp_start=j;}
            }
            result.push_back(t);++next;
            if(next==counts.size())break;
        }
        if(start<last_start) {
            window+=products[start+CP]-products[start];
            window_power+=energies[start+CP]-energies[start];
        }
    }
    return result;
}
Decode demodulate(const Vector& rx,const Bits& bits,Modulation m,std::size_t estimated_cp) {
    const auto bps=bits_per_symbol(m),block=unsigned(N)*bps;
    if(bits.empty() || bits.size()%block) throw std::invalid_argument("invalid reference bit count");
    const std::size_t frames=bits.size()/block;
    if(estimated_cp>=L || rx.size()<estimated_cp+frames*L) throw std::invalid_argument("truncated receiver stream");
    Decode result;Vector time(N);
    for(std::size_t f=0;f<frames;++f) {
        std::copy_n(rx.begin()+estimated_cp+f*L+CP,N,time.begin());
        const auto freq=hw1::fft(time);
        if(f<3) result.first_symbols.insert(result.first_symbols.end(),freq.begin(),freq.end());
        for(std::size_t k=0;k<N;++k) {
            const auto got=demap_symbol(freq[k],m);
            for(unsigned b=0;b<bps;++b) result.errors+=got[b]!=bits[f*block+k*bps+b];
        }
    }
    result.compared_bits=bits.size();return result;
}
double mean_power(const Vector& x) {
    if(x.empty()) throw std::invalid_argument("empty power input");
    double e=0;for(auto z:x)e+=std::norm(z);return e/double(x.size());
}
}
