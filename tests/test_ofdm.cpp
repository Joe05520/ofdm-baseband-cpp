#include "ofdm.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace hw2;
int checks=0;
void require(bool ok,const char* label){++checks;if(!ok)throw std::runtime_error(label);}
template<class F>void throws(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}require(caught,"expected rejection");}
int main(){try{
 double max_dft=0,max_roundtrip=0;
 for(std::size_t n:{1,3,16,64,127}){
  Vector x(n);for(std::size_t j=0;j<n;++j)x[j]={std::sin(j*0.73)+0.2*j,std::cos(j*0.91)};
  const auto y=hw1::fft(x),back=hw1::ifft(y);
  for(std::size_t k=0;k<n;++k){Complex ref{};for(std::size_t j=0;j<n;++j)ref+=x[j]*std::polar(1.0,-2*hw1::pi*double(j*k)/n);
   max_dft=std::max(max_dft,std::abs(y[k]-ref));max_roundtrip=std::max(max_roundtrip,std::abs(back[k]-x[k]));}
 }
 require(max_dft<1e-9,"independent direct DFT");require(max_roundtrip<1e-11,"FFT/IFFT reconstruction");
 std::cout<<"DFT max absolute error="<<max_dft<<"; roundtrip="<<max_roundtrip<<'\n';
 for(auto m:{Modulation::QPSK,Modulation::QAM16}){
  const auto bps=bits_per_symbol(m);double energy=0;
  for(unsigned word=0;word<(1u<<bps);++word){std::uint8_t b[4]{};for(unsigned i=0;i<bps;++i)b[i]=(word>>i)&1;
   const auto z=map_symbol(b,m);energy+=std::norm(z);const auto got=demap_symbol(z,m);
   for(unsigned i=0;i<bps;++i)require(got[i]==b[i],"all mapping labels roundtrip");}
  require(std::abs(energy/(1u<<bps)-1)<1e-12,"unit expected constellation energy");
 }
 std::uint8_t zeros[4]={0,0,0,0};
 require(std::abs(map_symbol(zeros,Modulation::QAM16)-Complex(-3,-3)/std::sqrt(10.0))<1e-14,"16QAM specified Gray label");
 constexpr double levels[]={-3,-1,1,3};constexpr unsigned labels[]={0,1,3,2};
 for(int i=0;i<4;++i)for(int q=0;q<4;++q)for(double delta:{-0.999,0.999}){
  auto b=demap_symbol(Complex(levels[i]+delta,levels[q]-delta)/std::sqrt(10.0),Modulation::QAM16);
  require(2*b[0]+b[1]==labels[i] && 2*b[2]+b[3]==labels[q],"nearest neighbor threshold");}
 for(auto m:{Modulation::QPSK,Modulation::QAM16}){
  const auto bits=random_bits(128*N*bits_per_symbol(m),313512050);const auto tx=modulate(bits,m);
  for(std::size_t f=0;f<128;++f)for(std::size_t i=0;i<CP;++i)require(tx[f*L+i]==tx[f*L+i+N],"CP copies last 16 samples");
  for(std::size_t d=0;d<L;++d){const auto rx=add_delay(tx,d);const auto t=synchronize(rx,64);
   require(t.cp_start==d,"all noiseless timing offsets");const auto dec=demodulate(rx,bits,m,t.cp_start);
   require(dec.errors==0 && dec.compared_bits==bits.size(),"noiseless end to end all bits");}
  const auto wrong=demodulate(add_delay(tx,37),bits,m,38);
  require(wrong.errors>bits.size()/4,"one sample wrong timing causes errors");
 }
 Vector blank(1000000);NoiseStats stat;const auto noisy=awgn(blank,0.02,9981,&stat);
 require(std::abs(stat.power/0.02-1)<0.005,"AWGN complex variance");
 require(std::abs(stat.mean_i)<0.0005 && std::abs(stat.mean_q)<0.0005,"AWGN zero mean");
 require(std::abs(stat.iq)<0.00005,"AWGN independent I Q");
 double lag=0;for(std::size_t i=1;i<noisy.size();++i)lag+=(noisy[i]*std::conj(noisy[i-1])).real();
 require(std::abs(lag/(noisy.size()-1))<0.0001,"AWGN lag one correlation");
 std::cout<<"Noise mean="<<stat.mean_i<<","<<stat.mean_q<<" power="<<stat.power<<" IQ="<<stat.iq<<'\n';
 // Formula-specific checks distinguish the lecture metric from the old real metric.
 Vector asymmetric(160);for(std::size_t n=0;n<CP;++n){asymmetric[n]={2,0};asymmetric[n+N]={1,0};}
 require(std::abs(synchronize(asymmetric,1).metric[0]-4)<1e-14,"lecture one-sided normalization can exceed one");
 const auto refbits=random_bits(1026*N*4,313512051);
 const auto reference=awgn(add_delay(modulate(refbits,Modulation::QAM16),37),1.0/N,1313);
 const std::vector<std::size_t> counts{1,2,4,8,16,32,64,128,256,512,1024};
 const auto sweep=synchronize_sweep(reference,counts);
 Vector rotated=reference;
 for(std::size_t i=0;i<rotated.size();++i)rotated[i]*=std::polar(1.0,2*hw1::pi*.2*i/N);
 double rolling_error=0,phase_error=0;
 for(std::size_t i=0;i<counts.size();++i){
  const auto direct=synchronize(reference,counts[i]),phase=synchronize(rotated,counts[i]);
  require(direct.cp_start==sweep[i].cp_start,"rolling sweep and direct sum choose same boundary");
  for(std::size_t d=0;d<L;++d){
   rolling_error=std::max(rolling_error,std::abs(direct.metric[d]-sweep[i].metric[d]));
   phase_error=std::max(phase_error,std::abs(direct.metric[d]-phase.metric[d]));
  }
 }
 require(rolling_error<1e-11,"rolling window matches lecture direct sum");
 require(phase_error<1e-11,"magnitude metric invariant to common lag-N CFO phase");
 require(synchronize(Vector(160),1).metric[0]==0,"zero energy returns finite zero metric");
 std::cout<<"Lecture metric rolling max error="<<rolling_error<<"; CFO phase invariance="<<phase_error<<'\n';
 throws([&]{synchronize_sweep(reference,{});});throws([&]{synchronize_sweep(reference,{4,1});});
 throws([&]{synchronize_sweep(reference,{1,1});});throws([&]{synchronize_sweep(Vector(80),{1});});
 int sync_good=0;
 for(auto m:{Modulation::QPSK,Modulation::QAM16}){
  const auto bits=random_bits(1026*N*bits_per_symbol(m),77123);const auto tx=modulate(bits,m);
  for(std::size_t d=0;d<L;++d){const auto rx=awgn(add_delay(tx,d),1.0/N,81000+d+bits_per_symbol(m)*100);
   const auto t=synchronize(rx,1024);sync_good+=t.cp_start==d;require(t.cp_start==d,"0dB acquisition all 80 offsets each modulation");}}
 throws([]{hw1::fft({});});throws([]{modulate({},Modulation::QPSK);});
 throws([]{modulate(Bits(129),Modulation::QPSK);});throws([]{synchronize(Vector(80),1);});
 throws([]{synchronize(Vector(100),0);});throws([]{add_delay(Vector(80),80);});
 throws([]{demodulate(Vector(79),Bits(128),Modulation::QPSK,0);});
 std::cout<<"Low SNR sync="<<sync_good<<"/160; PASS "<<checks<<" checks\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<" after "<<checks<<" checks\n";return 1;}}
