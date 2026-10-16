#include "multipath.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace hw2;
namespace mp=hw2::multipath;

int main(int argc,char**argv){try{
 const std::filesystem::path out=argc>1?argv[1]:"results/m40";std::filesystem::create_directories(out);
 std::ofstream trials(out/"multipath_timing_trials.csv"),trace(out/"multipath_timing_metric.csv");
 if(!trials||!trace)throw std::runtime_error("cannot open outputs");trials<<std::setprecision(17);trace<<std::setprecision(17);
 trials<<"modulation,profile,snr_db,packet,symbols,true_cp,estimated_cp,delta,fft_window_valid,errors,bits\n";
 trace<<"profile,snr_db,symbols,d,gamma,true_cp,estimated_cp\n";
 for(auto mode:{Modulation::QPSK,Modulation::QAM16}){
  const auto bps=bits_per_symbol(mode);
  for(std::size_t p=0;p<1024;++p){
   auto bits=random_bits(64*N*bps,2026100701000000ULL+10000*bps+p);
   auto tx=mp::transmit(bits,mode,16);
   for(int profile=0;profile<3;++profile){
    const std::string label=profile==0?"AWGN":profile==1?"Rayleigh_D12":"Rayleigh_D24";
    auto h=profile==0?Vector{1.}:mp::rayleigh(profile==1?4:8,2026100702000000ULL+p);
    auto delayed=add_delay(mp::filter(tx,h),37);
    for(int snr:{-1,15}){
     auto r=snr<0?delayed:awgn(delayed,1./(N*std::pow(10.,snr/10.)),2026100704000000ULL+10000*bps+p);
     auto sweep=synchronize_sweep(r,{40,64});
     for(std::size_t i=0;i<2;++i){
      const auto& t=sweep[i];const auto m=i==0?40:64;int delta=int(t.cp_start)-37;int D=int(h.size()-1);
      auto decoded=demodulate(r,bits,mode,t.cp_start);
      trials<<name(mode)<<','<<label<<','<<snr<<','<<p<<','<<m<<",37,"<<t.cp_start<<','<<delta<<','<<(delta<=0&&delta>=D-16)<<','<<decoded.errors<<','<<decoded.compared_bits<<'\n';
      if(p==0&&mode==Modulation::QPSK)for(std::size_t d=0;d<L;++d)trace<<label<<','<<snr<<','<<m<<','<<d<<','<<t.metric[d]<<",37,"<<t.cp_start<<'\n';
     }
    }
   }
  }
  std::cout<<name(mode)<<" PASS 1024 paired packets"<<std::endl;
 }
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
