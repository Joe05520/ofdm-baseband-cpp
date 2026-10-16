#include "ofdm.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
namespace fs=std::filesystem;
std::uint64_t mix(std::uint64_t x) {
    x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
    x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);
}
int main(int argc,char** argv){try{
    if(argc<3 || argc>5)throw std::invalid_argument("sync_study output.csv selection|fine [M] [trials=2000]");
    const std::string phase=argv[2];
    if(phase!="selection" && phase!="fine")throw std::invalid_argument("unknown phase");
    const auto chosen=argc>3 ? std::stoull(argv[3]) : 1024;
    const auto trials=argc>4 ? std::stoull(argv[4]) : 2000;
    if(chosen==0 || chosen>1024 || trials==0)throw std::invalid_argument("invalid M or trial count");
    std::vector<std::size_t> counts=phase=="selection" ?
        std::vector<std::size_t>{1,2,4,8,16,32,64,128,256,512,1024} :
        std::vector<std::size_t>{chosen};
    if(phase=="fine") { counts.clear(); for(std::size_t m=32;m<=64;++m)counts.push_back(m); }
    const std::vector<int> snrs{0,3,15};
    fs::create_directories(fs::path(argv[1]).parent_path());
    std::ofstream out(argv[1]);if(!out)throw std::runtime_error("cannot open study CSV");
    out<<std::setprecision(17);
    out<<"phase,modulation,snr_db,trial,symbols,true_cp,estimated_cp,error,circular_error,correct,peak_metric,true_metric,bit_seed,noise_seed\n";
    const std::uint64_t root=2026100601000000ULL;
    const auto begin=std::chrono::steady_clock::now();
    for(auto mode:{hw2::Modulation::QPSK,hw2::Modulation::QAM16}) {
        for(std::size_t trial=0;trial<trials;++trial) {
            const auto label=root+10000000ULL*hw2::bits_per_symbol(mode)+100ULL*trial;
            const auto bit_seed=mix(label);
            const auto delay=(37*trial+13)%hw2::L;
            const auto bits=hw2::random_bits(((phase=="fine" ? 1024 : counts.back())+2)*hw2::N*hw2::bits_per_symbol(mode),bit_seed);
            const auto tx=hw2::add_delay(hw2::modulate(bits,mode),delay);
            for(int snr:snrs) {
                const auto noise_seed=mix(label+1+snr);
                const auto rx=hw2::awgn(tx,1.0/(hw2::N*std::pow(10.0,snr/10.0)),noise_seed);
                const auto timing=hw2::synchronize_sweep(rx,counts);
                for(std::size_t i=0;i<counts.size();++i) {
                    const auto& t=timing[i];const int error=int(t.cp_start)-int(delay);
                    const int circular=(error+int(hw2::L)+int(hw2::L)/2)%int(hw2::L)-int(hw2::L)/2;
                    out<<phase<<','<<hw2::name(mode)<<','<<snr<<','<<trial<<','<<counts[i]<<','<<delay<<','
                       <<t.cp_start<<','<<error<<','<<circular<<','<<(error==0)<<','
                       <<t.metric[t.cp_start]<<','<<t.metric[delay]<<','<<bit_seed<<','<<noise_seed<<'\n';
                }
            }
            if((trial+1)%200==0) {
                out.flush();std::cout<<phase<<' '<<hw2::name(mode)<<" trial="<<trial+1<<'/'<<trials
                    <<" elapsed="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<"s"<<std::endl;
            }
        }
    }
    return 0;
}catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return 1;}}
