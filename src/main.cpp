#include "ofdm.hpp"
#include "sync_config.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace fs=std::filesystem;
std::ofstream output(const fs::path& p) {
    std::ofstream f(p);if(!f)throw std::runtime_error("cannot write "+p.string());
    f<<std::setprecision(17);return f;
}
int main(int argc,char** argv) { try {
    std::size_t used_k=hw2::ACQUISITION_SYMBOLS;
    std::size_t bit_count=1ULL<<25; fs::path out="results";
    if(argc>1 && std::string(argv[1])=="--help") {
        std::cout<<"Usage: ofdm [output_directory] [bits] [M]\nDefault bits=33554432; optional bits is for validation only.\n";return 0;
    }
    if(argc>4)throw std::invalid_argument("Usage: ofdm [output_directory] [bits] [M]");
    if(argc>1)out=argv[1];if(argc>2)bit_count=std::stoull(argv[2]);
    if(argc>3)used_k=std::stoull(argv[3]);
    if(used_k==0 || used_k>1024)throw std::invalid_argument("M must be 1..1024");
    if(bit_count<1024*64*4 || bit_count%(64*4))throw std::invalid_argument("bits must be >=262144 and divisible by 256");
    fs::create_directories(out);
    constexpr std::uint64_t bit_seed=313512050, noise_base=2026101200;
    constexpr std::size_t delay=37;
    auto ber=output(out/"ber.csv");
    ber<<"modulation,snr_db,bits,errors,ber,ofdm_symbols,true_cp,estimated_cp,fft_start,sync_symbols,tx_power,noise_variance,noise_power,measured_snr_db,noise_mean_i,noise_mean_q,noise_iq,bit_seed,noise_seed,seconds\n";
    auto syn=output(out/"synchronization.csv");syn<<"modulation,snr_db,symbols,candidate,metric,selected,true_cp\n";
    const auto bits=hw2::random_bits(bit_count,bit_seed);
    for(auto mode:{hw2::Modulation::QPSK,hw2::Modulation::QAM16}) {
        auto tx=hw2::modulate(bits,mode);const auto power=hw2::mean_power(tx);const auto frames=tx.size()/hw2::L;
        tx=hw2::add_delay(tx,delay);
        for(int snr=0;snr<=24;snr+=3) {
            const auto start=std::chrono::steady_clock::now();
            const auto variance=1.0/(hw2::N*std::pow(10.0,snr/10.0));
            const auto seed=noise_base+(mode==hw2::Modulation::QAM16 ? 1000 : 0)+snr;
            hw2::NoiseStats stats;const auto rx=hw2::awgn(tx,variance,seed,&stats);
            hw2::Timing timing;
            for(std::size_t k:{1,2,4,8,16,32,40,64,128,256,512,1024}) {
                const auto t=hw2::synchronize(rx,k);
                for(std::size_t d=0;d<hw2::L;++d)
                    syn<<hw2::name(mode)<<','<<snr<<','<<k<<','<<d<<','<<t.metric[d]<<','<<(d==t.cp_start)<<','<<delay<<'\n';
                if(k==used_k)timing=t;
            }
            // Only the estimated timing, never the transmitter's delay, enters the receiver.
            const auto decoded=hw2::demodulate(rx,bits,mode,timing.cp_start);
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            ber<<hw2::name(mode)<<','<<snr<<','<<decoded.compared_bits<<','<<decoded.errors<<','<<double(decoded.errors)/decoded.compared_bits<<','<<frames<<','<<delay<<','<<timing.cp_start<<','<<timing.cp_start+hw2::CP<<','<<used_k<<','<<power<<','<<variance<<','<<stats.power<<','<<10*std::log10(power/stats.power)<<','<<stats.mean_i<<','<<stats.mean_q<<','<<stats.iq<<','<<bit_seed<<','<<seed<<','<<seconds<<'\n';ber.flush();
            std::cout<<hw2::name(mode)<<" SNR="<<snr<<" errors="<<decoded.errors<<'/'<<decoded.compared_bits<<" CP="<<timing.cp_start<<" elapsed="<<seconds<<"s"<<std::endl;
            if(snr==3 || snr==15) {
                const auto tag=hw2::name(mode)+"_snr"+std::to_string(snr);
                auto wave=output(out/(tag+"_waveform.csv"));wave<<"sample,absolute_sample,ofdm_symbol,is_cp,tx_i,tx_q,rx_i,rx_q\n";
                for(std::size_t n=0;n<3*hw2::L;++n) {
                    const auto a=tx[delay+n],b=rx[delay+n];
                    wave<<n<<','<<delay+n<<','<<n/hw2::L<<','<<(n%hw2::L<hw2::CP)<<','<<a.real()<<','<<a.imag()<<','<<b.real()<<','<<b.imag()<<'\n';
                }
                auto con=output(out/(tag+"_constellation.csv"));con<<"symbol,ofdm_symbol,subcarrier,tx_i,tx_q,rx_i,rx_q\n";
                for(std::size_t n=0;n<decoded.first_symbols.size();++n) {
                    const auto a=hw2::map_symbol(&bits[n*hw2::bits_per_symbol(mode)],mode),b=decoded.first_symbols[n];
                    con<<n<<','<<n/hw2::N<<','<<n%hw2::N<<','<<a.real()<<','<<a.imag()<<','<<b.real()<<','<<b.imag()<<'\n';
                }
            }
        }
    }
    return 0;
} catch(const std::exception& e) {std::cerr<<"Error: "<<e.what()<<'\n';return 1;} }
