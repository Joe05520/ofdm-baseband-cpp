function verify_ofdm_figures
% Check measured series in all fourteen report figures.
root=fileparts(fileparts(mfilename('fullpath')));o=figure_folder(root);r=fullfile(root,'results');checks=0;count=0;
for mode=["QPSK","16QAM"]
 for snr=[3,15]
  name=sprintf('%s_snr%d_waveform',mode,snr);f=loadfig(name);
  for m=[64,40]
   p=r;prefix='waveform';if m==40,p=fullfile(r,'m40');prefix='waveform40';end
   w=readtable(fullfile(p,[name '.csv']));v={w.tx_i,w.rx_i,w.tx_q,w.rx_q};
   for k=1:4,check(f,sprintf('%s_%d',prefix,k),w.sample,v{k});end
  end
  close(f);
 end
 f=loadfig(char(mode+"_constellations"));
 for m=[64,40]
  p=r;tag='received_symbols';if m==40,p=fullfile(r,'m40');tag='received_symbols40';end
  obj=findall(f,'Tag',tag);assert(numel(obj)==2);
  for snr=[3,15]
   c=readtable(fullfile(p,sprintf('%s_snr%d_constellation.csv',mode,snr)));matched=false;
   for a=obj',matched=matched||(same(a.XData,c.rx_i)&&same(a.YData,c.rx_q));end
   assert(matched);checks=checks+1;
  end
 end
 close(f);
end
f=loadfig('ber_comparison');
for mode=["QPSK","16QAM"],for m=[64,40]
 p=r;tag=char(mode+"_measured");if m==40,p=fullfile(r,'m40');tag=char(mode+"_measured40");end
 b=readtable(fullfile(p,'ber.csv'),'TextType','string');b=b(b.modulation==mode&b.errors>0,:);check(f,tag,b.snr_db,b.ber);
end,end
ax=findobj(f,'Type','axes');assert(strcmp(ax.YLabel.String,'BER'));close(f);
s=readtable(fullfile(r,'synchronization.csv'),'TextType','string');a=readtable(fullfile(r,'m40','synchronization.csv'),'TextType','string');s=[s;a(a.symbols==40,:)];f=loadfig('synchronization');
for mode=["QPSK","16QAM"],for m=[1,16,32,40,64]
 a=sortrows(s(s.modulation==mode&s.snr_db==0&s.symbols==m,:),'candidate');check(f,sprintf('%s_K%d',mode,m),a.candidate,a.metric);
end,end
close(f);
t=readtable(fullfile(r,'sync_selection_summary.csv'),'TextType','string');fine=readtable(fullfile(r,'sync_fine_summary.csv'),'TextType','string');
f=loadfig('sync_m_selection');
for mode=["QPSK","16QAM"],for snr=[0,3,15]
 a=sortrows(t(t.modulation==mode&t.snr_db==snr,:),'symbols');check(f,sprintf('%s_accuracy_%d',mode,snr),a.symbols,a.correct_rate);check(f,sprintf('%s_rmse_%d',mode,snr),a.symbols,a.rmse);
 a=fine(fine.modulation==mode&fine.snr_db==snr&fine.symbols==40,:);check(f,sprintf('%s_added40_1_%d',mode,snr),40,a.correct_rate);check(f,sprintf('%s_added40_2_%d',mode,snr),40,a.rmse);
end,end
assert(isempty(findall(f,'Type','errorbar')));close(f);
f=loadfig('sync_m_fine');
for mode=["QPSK","16QAM"],for snr=[0,3,15]
 a=sortrows(fine(fine.modulation==mode&fine.snr_db==snr,:),'symbols');assert(height(a)==33);check(f,sprintf('%s_fine1_%d',mode,snr),a.symbols,a.correct_rate);check(f,sprintf('%s_fine2_%d',mode,snr),a.symbols,a.rmse);
end,end
close(f);
b=readtable(fullfile(r,'multipath','ber.csv'),'TextType','string');
for mode=["QPSK","16QAM"]
 f=loadfig(char("multipath/ber_"+mode));
 for profile=["Rayleigh_D12","Rayleigh_D24"],for cp=[0,8,12,16,24,32]
  a=b(b.modulation==mode&b.profile==profile&b.cp==cp&b.receiver=="ZF"&b.snr_db>=0&b.errors>0,:);check(f,sprintf('%s_%s_CP%d',mode,profile,cp),a.snr_db,a.ber);
 end,end
 a=b(b.modulation==mode&b.profile=="AWGN"&b.receiver=="ZF"&b.snr_db>=0&b.errors>0,:);objects=findobj(f,'Tag',char(mode+"_AWGN"));assert(numel(objects)==2);
 for obj=objects',assert(same(obj.XData,a.snr_db)&&same(obj.YData,a.ber));checks=checks+1;end
 close(f);
end
c=readtable(fullfile(r,'multipath','constellation.csv'),'TextType','string');f=loadfig('multipath/equalization_constellations');
for mode=["QPSK","16QAM"],for k=1:3
 cp=24;if k==2,cp=16;end;a=c(c.modulation==mode&c.cp==cp&c.snr_db==-1,:);
 if k==1,z=complex(a.y_i,a.y_q);else,z=complex(a.z_i,a.z_q);end
 check(f,sprintf('%s_scatter%d',mode,k),real(z),imag(z));ideal=unique(complex(a.x_i,a.x_q));check(f,sprintf('%s_ideal%d',mode,k),real(ideal),imag(ideal));
end,end
close(f);f=loadfig('multipath/timing_multipath');
for m=[64,40]
 if m==64,p=fullfile(r,'multipath');tf='timing_metric.csv';af='timing_trials.csv';tt='timing_';ht='hist_';else,p=fullfile(r,'m40');tf='multipath_timing_metric.csv';af='multipath_timing_trials.csv';tt='timing40_';ht='hist40_';end
 t=readtable(fullfile(p,tf),'TextType','string');a=readtable(fullfile(p,af),'TextType','string');
 if m==40,t=t(t.symbols==40,:);a=a(a.symbols==40,:);end
 for profile=["AWGN","Rayleigh_D12","Rayleigh_D24"]
  x=t(t.profile==profile&t.snr_db==-1,:);check(f,char(tt+profile),x.d,x.gamma);
 end
 for profile=["Rayleigh_D12","Rayleigh_D24"]
  x=a(a.modulation=="QPSK"&a.profile==profile&a.snr_db==15,:);check(f,char(ht+profile),0:24,histcounts(x.delta,-.5:1:24.5)/1024);
 end
end
legends=findobj(f,'Type','legend');for l=legends',assert(strcmp(l.Interpreter,'none'));end
close(f);assert(count==14);fprintf('PASS: %d report figures; %d measured series checked.\n',count,checks);
 function f=loadfig(name)
  for ext={'.png','.svg','.fig'},assert(isfile(figure_file(o,[name ext{1}])));end
  f=openfig(figure_file(o,[name '.fig']),'invisible');count=count+1;
 end
 function check(f,tag,x,y)
  obj=findall(f,'Tag',tag);assert(numel(obj)==1,'Missing/duplicate series: %s',tag);
  assert(same(obj.XData,x)&&same(obj.YData,y),'CSV mismatch: %s',tag);checks=checks+1;
 end
end
function ok=same(a,b)
ok=numel(a)==numel(b)&&all(abs(double(a(:))-double(b(:)))<1e-12);
end
