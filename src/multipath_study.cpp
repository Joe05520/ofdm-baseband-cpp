#include "multipath.hpp"
#include "sync_config.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <array>
using namespace hw2;
namespace mp=hw2::multipath;
constexpr std::size_t packets=1024,frames=64;
constexpr std::uint64_t bit_root=2026100701000000ULL,channel_root=2026100702000000ULL,noise_root=2026100703000000ULL,timing_root=2026100704000000ULL;
struct Stat {std::uint64_t errors=0;void add(std::uint64_t e){errors+=e;} };
int main(int argc,char**argv){try{
    const std::filesystem::path out=argc>1?argv[1]:"results/multipath";std::filesystem::create_directories(out);
    auto open=[&](const std::string&name){std::ofstream f(out/name);f<<std::setprecision(17);if(!f)throw std::runtime_error("cannot open output");return f;};
    auto ber=open("ber.csv"),scatter=open("constellation.csv"),timing=open("timing_trials.csv"),trace=open("timing_metric.csv");
    ber<<"modulation,profile,max_delay,cp,snr_db,receiver,errors,bits,packets,ber\n";
    scatter<<"modulation,profile,cp,snr_db,symbol,x_i,x_q,y_i,y_q,z_i,z_q,h_i,h_q\n";
    timing<<"modulation,profile,snr_db,packet,true_cp,estimated_cp,delta,fft_window_valid,errors,bits\n";
    trace<<"profile,snr_db,d,gamma,true_cp,estimated_cp\n";
    for(auto mode:{Modulation::QPSK,Modulation::QAM16}){
      const auto bps=bits_per_symbol(mode);const std::uint64_t packet_bits=frames*N*bps;
      std::array<std::array<std::array<std::array<Stat,2>,12>,6>,3> stats{};
      for(std::size_t p=0;p<packets;++p){
        auto bits=random_bits(packet_bits,bit_root+10000*bps+p);Vector X(frames*N),W(frames*N);
        Gaussian noise(noise_root+10000*bps+p);
        for(std::size_t f=0;f<frames;++f){Vector w(N);
            for(std::size_t k=0;k<N;++k){X[f*N+k]=map_symbol(&bits[(f*N+k)*bps],mode);w[k]=noise.sample(1./N);}
            auto wf=hw1::fft(w);std::copy(wf.begin(),wf.end(),W.begin()+f*N);}
        for(int profile=0;profile<3;++profile){
          const std::string label=profile==0?"AWGN":profile==1?"Rayleigh_D12":"Rayleigh_D24";
          auto h=profile==0?Vector{1.}:mp::rayleigh(profile==1?4:8,channel_root+p);auto H=mp::response(h);
          int ci=0;for(auto cp:{0u,8u,12u,16u,24u,32u}){
            if(profile==0 && cp!=16){++ci;continue;}
            auto tx=mp::transmit(bits,mode,cp),rx=mp::filter(tx,h);
            std::uint64_t counts[12][2]{};
            for(std::size_t f=0;f<frames;++f){auto y=mp::window_fft(rx,f,cp),z=mp::zero_force(y,H);
                for(int s=0;s<12;++s){double scale=s==11?0:std::pow(10.,-3.*s/20.);Vector yn(N);
                    for(std::size_t k=0;k<N;++k)yn[k]=y[k]+scale*W[f*N+k];auto zn=mp::zero_force(yn,H);
                    counts[s][0]+=mp::errors(zn,bits,mode,f);counts[s][1]+=mp::errors(yn,bits,mode,f);
                    if(p==0 && f<4 && profile==2 && (cp==16 || cp==24) && (s==8 || s==11))
                      for(std::size_t k=0;k<N;++k){auto x=X[f*N+k];scatter<<name(mode)<<','<<label<<','<<cp<<','<<(s==11?-1:3*s)<<','<<f*N+k<<','<<x.real()<<','<<x.imag()<<','<<yn[k].real()<<','<<yn[k].imag()<<','<<zn[k].real()<<','<<zn[k].imag()<<','<<H[k].real()<<','<<H[k].imag()<<'\n';}
                }
            }
            for(int s=0;s<12;++s)for(int e=0;e<2;++e)stats[profile][ci][s][e].add(counts[s][e]);
            if(cp==16){
              // Unchanged lecture formula at G=C=16,M=64, unknown CP=37.
              // This diagnostic does not feed the controlled CP BER experiment.
              auto delayed=add_delay(rx,37);
              for(int snr:{-1,15}){auto r=snr<0?delayed:awgn(delayed,1./(N*std::pow(10.,snr/10.)),timing_root+10000*bps+p);
                  auto t=synchronize_sweep(r,{64}).front();int delta=int(t.cp_start)-37;const int D=int(h.size()-1);
                  bool valid=delta<=0 && delta>=D-16;
                  auto decoded=demodulate(r,bits,mode,t.cp_start);
                  timing<<name(mode)<<','<<label<<','<<snr<<','<<p<<",37,"<<t.cp_start<<','<<delta<<','<<valid<<','<<decoded.errors<<','<<decoded.compared_bits<<'\n';
                  if(p==0 && mode==Modulation::QPSK)for(std::size_t d=0;d<L;++d)trace<<label<<','<<snr<<','<<d<<','<<t.metric[d]<<",37,"<<t.cp_start<<'\n';
              }
            }
            ++ci;
          }
        }
        if((p+1)%256==0)std::cout<<name(mode)<<" "<<p+1<<"/"<<packets<<std::endl;
      }
      int D[3]={0,12,24};for(int profile=0;profile<3;++profile){std::string label=profile==0?"AWGN":profile==1?"Rayleigh_D12":"Rayleigh_D24";int ci=0;
        for(auto cp:{0u,8u,12u,16u,24u,32u}){if(profile==0 && cp!=16){++ci;continue;}
                for(int s=0;s<12;++s)for(int e=0;e<2;++e){auto a=stats[profile][ci][s][e];
            ber<<name(mode)<<','<<label<<','<<D[profile]<<','<<cp<<','<<(s==11?-1:3*s)<<','<<(e==0?"ZF":"NoEQ")<<','<<a.errors<<','<<packet_bits*packets<<','<<packets<<','<<double(a.errors)/(packet_bits*packets)<<'\n';}
          ++ci;}
      }
    }
    std::cout<<"Study complete\n";
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
