function plot_ofdm_figures(outputDir)
%PLOT_HW2_FIGURES Reproduce all OFDM figures from measured C++ CSV data.
% Author: Jin-Teng Sung.
% Tested with MATLAB R2024a. No Communications Toolbox is needed.
%
% In the project directory:
%   addpath('matlab'); plot_ofdm_figures
% Or save copies elsewhere:
%   plot_ofdm_figures(fullfile(pwd, 'my_figures'))
%
% Outputs: ten 300-dpi PNGs, ten SVGs, ten editable MATLAB FIGs.
% The C++ simulator implements mapping, IFFT, CP, AWGN, synchronization,
% FFT, demapping and BER. This function only reads and plots its results.
% The analytical reference curves below are not simulated BER values.
% Zero-error observations remain zero in CSV and are omitted on the log axis.

root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, outputDir = figure_folder(root); end
if ~isfolder(outputDir), mkdir(outputDir); end
results = fullfile(root, 'results');
B = readtable(fullfile(results, 'ber.csv'), 'TextType', 'string');
B40=readtable(fullfile(results,'m40','ber.csv'),'TextType','string');
F=readtable(fullfile(results,'sync_fine_summary.csv'),'TextType','string');
S = readtable(fullfile(results, 'synchronization.csv'), 'TextType', 'string');
S40=readtable(fullfile(results,'m40','synchronization.csv'),'TextType','string');
S=[S;S40(S40.symbols==40,:)];
selectedM=64;
T=readtable(fullfile(results,'sync_selection_summary.csv'),'TextType','string');
assert(height(B) == 18 && all(B.bits == 33554432), 'Incomplete BER dataset.');
assert(all(abs(B.ber - B.errors ./ B.bits) < 1e-15), 'BER count mismatch.');
assert(all(B.estimated_cp == B.true_cp), 'Dataset contains timing failures.');
blue = [0.14 0.36 0.54]; orange = [0.72 0.36 0.14]; gray = [.34 .36 .40];
modes = ["QPSK", "16QAM"];
figureNames = strings(0,1);

for mode = modes
    for snr = [3,15]
        tag = sprintf('%s_snr%d', mode, snr);
        W = readtable(fullfile(results, [tag '_waveform.csv']));
        C = readtable(fullfile(results, [tag '_constellation.csv']));
        W40=readtable(fullfile(results,'m40',[tag '_waveform.csv']));
        C40=readtable(fullfile(results,'m40',[tag '_constellation.csv']));
        assert(height(W) == 240 && height(C) == 192, 'Wrong sample count.');
        assert(isequal(W.sample, (0:239)'), 'Waveform sample order mismatch.');
        assert(isequal(C.symbol, (0:191)'), 'Constellation symbol order mismatch.');
        f = newFigure(10.8,4.1);
        t = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
        lim = 1.35 * max(abs([W.tx_i;W.tx_q;W.rx_i;W.rx_q]));
        values = {W.tx_i,W.rx_i,W.tx_q,W.rx_q};
        for tile = 1:4
            ax = nexttile(t,tile); styleAxes(ax); hold(ax,'on');
            for block = 0:2
                start = block*80;
                patch(ax,[start-.5 start+15.5 start+15.5 start-.5], ...
                    [-lim -lim lim lim],orange,'FaceAlpha',.13, ...
                    'EdgeColor','none','HandleVisibility','off');
                xline(ax,start+15.5,'--','Color',orange,'LineWidth',.65, ...
                    'HandleVisibility','off');
                if tile <= 2
                    text(ax,start+7.5,.88*lim,'CP','Color',orange, ...
                        'HorizontalAlignment','center','FontSize',8);
                    text(ax,start+47.5,.88*lim,sprintf('Data %d',block+1), ...
                        'HorizontalAlignment','center','FontSize',8);
                end
            end
            color = gray; if mod(tile,2)==0, color=blue; end
            plot(ax,W.sample,values{tile},'Color',color,'LineWidth',.75, ...
                'Tag',sprintf('waveform_%d',tile),'DisplayName','M = 64');
            values40={W40.tx_i,W40.rx_i,W40.tx_q,W40.rx_q};
            plot(ax,W40.sample,values40{tile},'--','Color',orange,'LineWidth',.7,'Tag',sprintf('waveform40_%d',tile),'DisplayName','M = 40');
            if tile==2,legend(ax,'Location','southwest','FontSize',7,'Box','off');end
            xlim(ax,[-.5,239.5]); ylim(ax,[-lim,lim]);
            xticks(ax,[0 16 80 96 160 176 239]); ax.YGrid='on';
            if tile<=2
                if tile==1, title(ax,'Before AWGN'); else, title(ax,'After AWGN'); end
                ax.XTickLabel=[];
            else
                xlabel(ax,'n');
            end
            if tile==1, ylabel(ax,'I'); end
            if tile==3, ylabel(ax,'Q'); end
        end
        title(t,sprintf('%s OFDM waveforms at SNR %d dB',mode,snr), ...
            'FontName','Arial','FontSize',13,'FontWeight','bold');
        name=[tag '_waveform']; saveFigure(f,outputDir,name);
        figureNames(end+1)=string(name); %#ok<AGROW>

    end
    f=newFigure(10.8,4.5);
    t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    for snr=[3,15]
        C=readtable(fullfile(results,sprintf('%s_snr%d_constellation.csv',mode,snr)));
        C40=readtable(fullfile(results,'m40',sprintf('%s_snr%d_constellation.csv',mode,snr)));
        ax=nexttile(t);styleAxes(ax);constellation(ax,C,mode,blue,orange,C40);
        title(ax,sprintf('SNR %d dB',snr));
        if snr==15,legend(ax,'Location','northeast','FontSize',8,'Box','off');end
    end
    title(t,sprintf('%s received constellations',mode),'FontName','Arial', ...
        'FontSize',14,'FontWeight','bold');
    name=char(mode+"_constellations");saveFigure(f,outputDir,name);
    figureNames(end+1)=string(name); %#ok<AGROW>
end

% BER curves. QPSK: Q(sqrt(gamma)). Gray 16-QAM: exact 4-PAM bit BER.
f=newFigure(10.4,5);ax=axes(f);styleAxes(ax);hold(ax,'on');
ax.YScale='log'; colors={blue,orange}; markers={'o','s'};
theoryRows=table();
for i=1:2
    mode=modes(i); color=colors{i}; marker=markers{i};
    b=sortrows(B(B.modulation==mode,:),'snr_db');
    assert(isequal(b.snr_db,(0:3:24)'), 'Missing SNR point.');
    snr=linspace(0,24,241); gamma=10.^(snr/10);
    Q=@(x) 0.5*erfc(x/sqrt(2));
    if mode=="QPSK"
        theory=Q(sqrt(gamma)); bps=2;
    else
        a=sqrt(gamma/5);theory=.75*Q(a)+.5*Q(3*a)-.25*Q(5*a);bps=4;
    end
    plot(ax,snr,theory,'--','Color',color,'LineWidth',1.4, ...
        'DisplayName',char(mode+" exact AWGN theory"),'Tag',char(mode+"_theory"));
    nz=b.errors>0;
    plot(ax,b.snr_db(nz),b.ber(nz),'LineStyle','none','Marker',marker, ...
        'MarkerSize',6,'MarkerFaceColor',color,'Color',color, ...
        'DisplayName',char(mode+" M = 64"),'Tag',char(mode+"_measured"));
    c=B40(B40.modulation==mode,:);nz40=c.errors>0;
    plot(ax,c.snr_db(nz40),c.ber(nz40),'LineStyle','none','Marker',marker,'MarkerSize',9,'Color',color,'LineWidth',1,'DisplayName',char(mode+" M = 40"),'Tag',char(mode+"_measured40"));
    theoryRows=[theoryRows;table(repmat(mode,numel(snr),1),snr',theory', ...
        'VariableNames',{'modulation','snr_db','theory_ber'})]; %#ok<AGROW>
end
xlim(ax,[-.5,24.5]);ylim(ax,[1e-8,1]);xticks(ax,0:3:24);grid(ax,'on');
ax.YMinorGrid='off';xlabel(ax,'SNR (dB)');ylabel(ax,'BER');
title(ax,'OFDM BER with estimated symbol timing');
legend(ax,'Location','southwest','FontSize',8.5,'Box','off');
saveFigure(f,outputDir,'ber_comparison');figureNames(end+1)="ber_comparison";
writetable(theoryRows,fullfile(results,'matlab_theory_reference.csv'));

f=newFigure(10.8,4);t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for mode=modes
    ax=nexttile(t);styleAxes(ax);hold(ax,'on');
    ks=[1 16 32 40 selectedM]; colors={[.60 .60 .60],orange,[.20 .52 .33],[.55 .20 .55],blue};styles={':','--','-.','--','-'};
    for i=1:5
        s=sortrows(S(S.modulation==mode & S.snr_db==0 & S.symbols==ks(i),:),'candidate');
        assert(height(s)==80 && isequal(s.candidate,(0:79)'), 'Timing candidates missing.');
        plot(ax,s.candidate,s.metric,'Color',colors{i},'LineStyle',styles{i}, ...
            'LineWidth',1.3,'DisplayName',sprintf('M = %d',ks(i)), ...
            'Tag',sprintf('%s_K%d',mode,ks(i)));
    end
    xline(ax,37,'-.','Color',[.2 .2 .2],'DisplayName','True CP start = 37');
    xlim(ax,[0 79]);ylim(ax,[0 max(.6,1.05*max(S.metric(S.modulation==mode&S.snr_db==0&ismember(S.symbols,ks))))]);grid(ax,'on');
    xlabel(ax,'d');title(ax,mode);
    if mode=="QPSK",ylabel(ax,'\Gamma_d');
    else,legend(ax,'Location','northeast','FontSize',8,'Box','off');end
end
title(t,'Timing acquisition at SNR 0 dB','FontName','Arial','FontSize',14,'FontWeight','bold');
saveFigure(f,outputDir,'synchronization');figureNames(end+1)="synchronization";
snrs=[0 3 15];studyColors={blue,orange,[.20 .52 .33]};
% Full-range measured accuracy and absolute timing RMSE.
f=newFigure(10.8,7);t=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
for mi=1:2
    mode=modes(mi);
    for metricIndex=1:2
        ax=nexttile(t,mi+2*(metricIndex-1));styleAxes(ax);hold(ax,'on');
        for si=1:3
            r=sortrows(T(T.modulation==mode&T.snr_db==snrs(si),:),'symbols');
            if metricIndex==1
                plot(ax,r.symbols,r.correct_rate,'-o','Color',studyColors{si}, ...
                    'LineWidth',1.1,'MarkerSize',4,'DisplayName',sprintf('%d dB',snrs(si)), ...
                    'Tag',sprintf('%s_accuracy_%d',mode,snrs(si)));
            else
                plot(ax,r.symbols,r.rmse,'-o','Color',studyColors{si},'LineWidth',1.1, ...
                    'MarkerSize',4,'DisplayName',sprintf('%d dB',snrs(si)), ...
                    'Tag',sprintf('%s_rmse_%d',mode,snrs(si)));
            end
        end
        for si=1:3
            r40=F(F.modulation==mode&F.snr_db==snrs(si)&F.symbols==40,:);
            if metricIndex==1,v40=r40.correct_rate;else,v40=r40.rmse;end
            plot(ax,40,v40,'d','MarkerSize',7,'Color',studyColors{si},'LineWidth',1.3,'HandleVisibility','off','Tag',sprintf('%s_added40_%d_%d',mode,metricIndex,snrs(si)));
        end
        xline(ax,40,':','M = 40','HandleVisibility','off');
        ax.XScale='log';ax.XTick=[1 4 16 40 64 256 1024];xlim(ax,[.9 1150]);grid(ax,'on');
        xline(ax,selectedM,'--','M = 64','HandleVisibility','off');
        xlabel(ax,'M');
        if metricIndex==1
            ylim(ax,[0 1.025]);ylabel(ax,'Accuracy');title(ax,mode);
            if mi==2,legend(ax,'Location','southeast','Box','off');end
        else
            ylabel(ax,'RMSE (samples)');ylim(ax,[0 max(T.rmse)*1.06]);
        end
    end
end
title(t,'M sweep: 2,000 trials per condition', ...
    'FontName','Arial','FontSize',12,'FontWeight','bold');
saveFigure(f,outputDir,'sync_m_selection');figureNames(end+1)="sync_m_selection";

F=readtable(fullfile(results,'sync_fine_summary.csv'),'TextType','string');
f=newFigure(10.8,6.8);t=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
for mi=1:2
 for row=1:2
  ax=nexttile(t,mi+2*(row-1));styleAxes(ax);hold(ax,'on');
  for si=1:3
   r=sortrows(F(F.modulation==modes(mi)&F.snr_db==snrs(si),:),'symbols');
   if row==1,v=r.correct_rate;else,v=r.rmse;end
   plot(ax,r.symbols,v,'-o','Color',studyColors{si},'MarkerSize',3,'DisplayName',sprintf('%d dB',snrs(si)),'Tag',sprintf('%s_fine%d_%d',modes(mi),row,snrs(si)));
  end
  xline(ax,40,'--','M = 40','HandleVisibility','off');
  xlabel(ax,'M');xlim(ax,[32 64]);xticks(ax,32:4:64);grid(ax,'on');
  if row==1,ylabel(ax,'Accuracy');ylim(ax,[.988 1.001]);title(ax,modes(mi));else,ylabel(ax,'RMSE (samples)');end
  if mi==2&&row==1,legend(ax,'Location','southeast','Box','off');end
 end
end
title(t,'Integer M sweep: 32 to 64','FontName','Arial','FontSize',13,'FontWeight','bold');
saveFigure(f,outputDir,'sync_m_fine');figureNames(end+1)="sync_m_fine";
assert(numel(figureNames)==10 && numel(unique(figureNames))==10);
fprintf('PASS: Generated 10 report figures in PNG, SVG and FIG formats.\n');
end

function f=newFigure(width,height)
f=figure('Visible','off','Color','white','Units','inches', ...
    'Position',[1 1 width height],'PaperPositionMode','auto', ...
    'InvertHardcopy','off');
end
function styleAxes(ax)
set(ax,'FontName','Arial','FontSize',10,'Box','off','LineWidth',.65, ...
    'GridAlpha',.17,'TickDir','out','TitleFontSizeMultiplier',1.1);
end
function constellation(ax,C,mode,blue,orange,C40)
hold(ax,'on');
scatter(ax,C.rx_i,C.rx_q,18,blue,'filled','MarkerFaceAlpha',.7, ...
    'DisplayName','M = 64','Tag','received_symbols');
scatter(ax,C40.rx_i,C40.rx_q,24,orange,'o','LineWidth',.55,'DisplayName','M = 40','Tag','received_symbols40');
ideal=unique(complex(C.tx_i,C.tx_q));
scatter(ax,real(ideal),imag(ideal),65,orange,'x','LineWidth',1.5, ...
    'DisplayName','Ideal points','Tag','ideal_symbols');
bound=2.6;if mode=="16QAM",bound=2.8;end
assert(all(abs(C.rx_i)<bound & abs(C.rx_q)<bound),'Constellation point clipped.');
xlim(ax,[-bound bound]);ylim(ax,[-bound bound]);axis(ax,'square');grid(ax,'on');
xlabel(ax,'I');ylabel(ax,'Q');
end
function saveFigure(f,outputDir,name)
f.Name=name;
f.UserData=struct('author','Jin-Teng Sung', ...
    'generator','matlab/plot_ofdm_figures.m','data_source','results/*.csv');
drawnow;
exportgraphics(f,figure_file(outputDir,[name '.png']),'Resolution',300,'BackgroundColor','white');
print(f,figure_file(outputDir,[name '.svg']),'-dsvg','-painters');
savefig(f,figure_file(outputDir,[name '.fig']));
close(f);
fprintf('Saved %s\n',name);
end
