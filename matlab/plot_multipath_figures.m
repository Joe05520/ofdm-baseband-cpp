function plot_multipath_figures
% MATLAB-authoritative Rayleigh multipath / CP extension plots.
% C++ implements the communications system; MATLAB plots exported measurements.
% No Communications Toolbox is required. Run from any working directory.
root=fileparts(fileparts(mfilename('fullpath')));r=fullfile(root,'results','multipath');o=figure_folder(root,'multipath');if ~isfolder(o),mkdir(o);end
B=readtable(fullfile(r,'ber.csv'),'TextType','string');S=readtable(fullfile(r,'constellation.csv'),'TextType','string');T=readtable(fullfile(r,'timing_metric.csv'),'TextType','string');A=readtable(fullfile(r,'timing_trials.csv'),'TextType','string');
assert(height(B)==624);colors=lines(6);names=strings(0,1);
for mode=["QPSK","16QAM"]
f=newfig(10.8,4.6);t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for profile=["Rayleigh_D12","Rayleigh_D24"],ax=nexttile(t);hold(ax,'on');
for j=1:6,cp=[0,8,12,16,24,32];cp=cp(j);a=B(B.modulation==mode & B.profile==profile & B.cp==cp & B.receiver=="ZF" & B.snr_db>=0,:);a=a(a.errors>0,:);semilogy(ax,a.snr_db,a.ber,'o-','Color',colors(j,:),'LineWidth',1.2,'MarkerSize',3,'DisplayName',sprintf('CP = %d',cp),'Tag',sprintf('%s_%s_CP%d',mode,profile,cp));end
a=B(B.modulation==mode & B.profile=="AWGN" & B.receiver=="ZF" & B.snr_db>=0 & B.errors>0,:);semilogy(ax,a.snr_db,a.ber,'k--','LineWidth',1.3,'DisplayName','AWGN CP = 16','Tag',char(mode+"_AWGN"));
ax.YScale='log';ylim(ax,[1e-5,.6]);xlim(ax,[0,30]);xlabel(ax,'SNR (dB)');ylabel(ax,'BER');style(ax);title(ax,replace(profile,'_',' '));legend(ax,'Location','southwest','FontSize',8,'NumColumns',2,'Box','off');
end
title(t,mode+"  perfect timing and CSI  ZF equalizer");savefigs(f,char("ber_"+mode));
end

f=newfig(10.8,6.3);t=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
for mode=["QPSK","16QAM"],for choice=1:3,ax=nexttile(t);hold(ax,'on');cp=24;if choice==2,cp=16;end;a=S(S.modulation==mode & S.cp==cp & S.snr_db==-1,:);
if choice==1,z=complex(a.y_i,a.y_q);label='CP24 without equalization';else,z=complex(a.z_i,a.z_q);label=sprintf('CP%d with ZF',cp);end
plot(ax,real(z),imag(z),'.','MarkerSize',5,'Tag',sprintf('%s_scatter%d',mode,choice));ideal=unique(complex(a.x_i,a.x_q));plot(ax,real(ideal),imag(ideal),'kx','MarkerSize',6,'LineWidth',1.2,'Tag',sprintf('%s_ideal%d',mode,choice));axis(ax,'equal');lim=max(1.3,1.08*max(abs([real(z);imag(z)])));xlim(ax,[-lim,lim]);ylim(ax,[-lim,lim]);xlabel(ax,'I');ylabel(ax,'Q');title(ax,mode+"  "+label,'FontSize',10);style(ax);
end,end
title(t,'D = 24  no noise  realization 0  first 256 symbols');savefigs(f,'equalization_constellations');

T40=readtable(fullfile(root,'results','m40','multipath_timing_metric.csv'),'TextType','string');
A40=readtable(fullfile(root,'results','m40','multipath_timing_trials.csv'),'TextType','string');
f=newfig(10.8,4);t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');ax=nexttile(t);hold(ax,'on');cols=lines(3);profileIndex=0;
for profile=["AWGN","Rayleigh_D12","Rayleigh_D24"]
 profileIndex=profileIndex+1;a=T(T.profile==profile&T.snr_db==-1,:);plot(ax,a.d,a.gamma,'Color',cols(profileIndex,:),'LineWidth',1.2,'DisplayName',char(replace(profile,"_D"," D = ")+" M = 64"),'Tag',char("timing_"+profile));
 a=T40(T40.profile==profile&T40.snr_db==-1&T40.symbols==40,:);plot(ax,a.d,a.gamma,'--','Color',cols(profileIndex,:),'LineWidth',1.2,'DisplayName',char(replace(profile,"_D"," D = ")+" M = 40"),'Tag',char("timing40_"+profile));
end
xline(ax,37,'k--','True CP','HandleVisibility','off');xlabel(ax,'d');ylabel(ax,'\Gamma_d');title(ax,'QPSK realization 0  no noise');legend(ax,'Location','northeast','FontSize',7,'Box','off','Interpreter','none');style(ax);
ax=nexttile(t);hold(ax,'on');profileIndex=0;
for profile=["Rayleigh_D12","Rayleigh_D24"]
 profileIndex=profileIndex+1;a=A(A.modulation=="QPSK"&A.profile==profile&A.snr_db==15,:);n=histcounts(a.delta,-.5:1:24.5);plot(ax,0:24,n/1024,'o-','Color',cols(profileIndex,:),'MarkerSize',3,'DisplayName',char(replace(profile,"_D"," D = ")+" M = 64"),'Tag',char("hist_"+profile));
 a=A40(A40.modulation=="QPSK"&A40.profile==profile&A40.snr_db==15&A40.symbols==40,:);n=histcounts(a.delta,-.5:1:24.5);plot(ax,0:24,n/1024,'--d','Color',cols(profileIndex,:),'MarkerSize',3,'DisplayName',char(replace(profile,"_D"," D = ")+" M = 40"),'Tag',char("hist40_"+profile));
end
xlabel(ax,'Timing error');ylabel(ax,'Fraction');title(ax,'M = 40 / 64, C = G = 16, 15 dB');legend(ax,'Location','northeast','FontSize',7,'Box','off','Interpreter','none');style(ax);
title(t,'Paired M = 40 / 64 timing in multipath');savefigs(f,'timing_multipath');

assert(numel(names)==4);
fprintf('PASS: Generated 4 multipath report figures.\n');
    function savefigs(f,name)
        exportgraphics(f,figure_file(o,[name '.png']),'Resolution',300);print(f,figure_file(o,[name '.svg']),'-dsvg');savefig(f,figure_file(o,[name '.fig']));names(end+1)=string(name);close(f);
    end
end
function f=newfig(w,h)
f=figure('Visible','off','Color','w','Units','inches','Position',[1,1,w,h]);
end
function style(ax)
set(ax,'FontName','Arial','FontSize',9,'LineWidth',.8,'Box','off','XGrid','on','YGrid','on','GridAlpha',.15);
end
