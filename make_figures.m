function make_figures(outdir)
%MAKE_FIGURES Single-scenario publication figures 2,2b,3,4,5.
% Fig5: direct epsPY response and sign-conditioned responses of 4 inputs.
if nargin<1||isempty(outdir)
    outdir=fullfile(fileparts(mfilename('fullpath')),'output');
end
figdir=fullfile(outdir,'figures');if ~exist(figdir,'dir'),mkdir(figdir);end
Z=load(fullfile(outdir,'results.mat'),'names','config','conv');
A=load(fullfile(outdir,'analysis_results.mat'), ...
    'SobolCI','SobolCI_Y2','Heat','Dir');
P=Z.config.priors;names=Z.names;cfg=Z.config;
blue=[.16 .40 .68];orange=[.86 .37 .19];teal=[.10 .56 .45];
set(groot,'defaultAxesFontName','Arial', ...
    'defaultAxesFontSize',10,'defaultFigureColor','white');
scenarios={'all priors','epsPY < 0','epsPY > 0'};

% Figure 2: Y1 Sobol indices, 1 x 3 (no environment dimension).
for which=1:2
    if which==1
        data=A.SobolCI;ylabelText='Y1';
        outname='Fig02_sobol_sorted_95CI';label='Fig 2';
    else
        data=A.SobolCI_Y2;ylabelText='Y2';
        outname='Fig02b_Y2_sobol_sorted_95CI';label='Fig 2b';
    end
    fh=figure('Visible','off','Position',[20 20 1620 610]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    for sg=0:2
        nexttile;hold on;
        rr=data(data(:,1)==sg,:);
        assert(size(rr,1)==10,'Expected 10 parameters for each epsPY scenario');
        for j=1:10
            plot_one_bar(j-.18,rr(j,4),rr(j,5),rr(j,6),blue);
            plot_one_bar(j+.18,rr(j,7),rr(j,8),rr(j,9),orange);
        end
        if all(~isfinite(rr(:,4))) && all(~isfinite(rr(:,7)))
            xlim([-.1 1.0]);text(.04,.50,'Undefined: output variance below threshold', ...
                'Units','normalized','FontSize',9,'Color',[.3 .3 .3]);
        else
            axis tight;
        end
        set(gca,'YDir','reverse','YTick',1:10,'YTickLabel',names(rr(:,3)));
        ylim([.45 10.55]);xlabel('Sobol index (95% t CI)');
        title(scenarios{sg+1});xline(0,':','Color',[.4 .4 .4]);
        h1=plot(NaN,NaN,'s','Color',blue,'MarkerFaceColor',blue);
        h2=plot(NaN,NaN,'s','Color',orange,'MarkerFaceColor',orange);
        legend([h1 h2],{'S_1','S_T'},'Location','best','FontSize',8);
        grid on;set(gca,'YGrid','off');
    end
    sgtitle(sprintf('%s | %s, Ptarget=0.80; descending S_T; 95%% CIs from six scrambles', ...
        label,ylabelText));
    export_all(figdir,outname,fh,false);
end

% Figure 3: 3 growth targets x 10 inputs, all 4 response variables.
fh=figure('Visible','off','Position',[20 20 1430 860]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
outputs={'Y1','Y2','Dplus','Dminus'};
for out=1:4
    nexttile;vals=nan(3,10);
    for t=1:3
        for idx=1:10
            rr=A.Heat(A.Heat(:,1)==t & A.Heat(:,2)==out & A.Heat(:,3)==idx,:);
            assert(size(rr,1)==1);
            vals(t,idx)=rr(5);
        end
    end
    hh=imagesc(vals);set(hh,'AlphaData',isfinite(vals));
    set(gca,'Color',[.9 .9 .9],'XTick',1:10,'XTickLabel',names, ...
        'YTick',1:3,'YTickLabel',{'Ptarget=.50','Ptarget=.80','Ptarget=.95'});
    xtickangle(40);clim([0 1]);colorbar;title(outputs{out});
end
sgtitle('Fig 3 | Single-scenario total-effect sensitivity across growth targets and outputs');
export_all(figdir,'Fig03_sensitivity_heatmap',fh,true);

% Figure 4: convergence of three most influential full-prior Y1 inputs.
rr=A.SobolCI(A.SobolCI(:,1)==0 & A.SobolCI(:,2)<=3,:);
chosen=rr(:,3)';palette=[blue;orange;teal];
if strcmp(cfg.mode,'quick'),nvec=[8 16 32 64];
else,nvec=[256 512 1024 2048];end
fh=figure('Visible','off','Position',[20 20 1240 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for out=1:2
    nexttile;hold on;hh=gobjects(1,numel(chosen));
    for j=1:numel(chosen)
        idx=chosen(j);
        v=Z.conv(Z.conv(:,5)==idx,:);
        for np=nvec
            found=v(v(:,2)==np,:);
            if np==nvec(end)
                if ~isempty(found)
                    y=found(end,5+out);
                    if isfinite(y)
                        plot(np,y,'d','Color',palette(j,:), ...
                            'MarkerFaceColor',palette(j,:), ...
                            'LineWidth',1.3,'MarkerSize',7);
                    end
                end
            else
                found=found(found(:,1)<=cfg.nCI,:);
                [mu,lo,hi]=replicate_ci(found(:,5+out));
                if isfinite(mu)
                    errorbar(np,mu,mu-lo,hi-mu,'o', ...
                        'Color',palette(j,:),'MarkerFaceColor',palette(j,:), ...
                        'LineWidth',1.15);
                end
            end
        end
        hh(j)=plot(NaN,NaN,'-o','Color',palette(j,:), ...
            'MarkerFaceColor',palette(j,:));
    end
    set(gca,'XScale','log','XTick',nvec, ...
        'XTickLabel',cellstr(num2str(nvec')));
    xlim([min(nvec)*.85 max(nvec)*1.15]);yline(0,':');
    xlabel('Sobol base sample size N per scramble');
    if out==1,ylabel('S_1');title('First-order index S_1');
    else,ylabel('S_T');title('Total-effect index S_T');end
    legend(hh,names(chosen),'Location','best');grid on;
end
sgtitle('Fig 4 | Sobol convergence: equal-N 95% CIs and independent largest-N estimate (diamonds)');
export_all(figdir,'Fig04_sobol_convergence',fh,false);

% Figure 5: 4x3 panels: 1 direct epsPY + 9 sign-conditioned parameters.
% Solid = epsPY negative, dashed = epsPY positive; magnitude remains
% independently sampled in each subset. Blue=Increase, red=Decrease.
fh=figure('Visible','off','Position',[20 20 1580 1480]);
tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
plotVars=[8 1 2 3 10 4 5 6 7 9];
for tile=1:10
    idx=plotVars(tile);nexttile;hold on;
    if idx==8
        for d=1:2
            rr=A.Dir(A.Dir(:,1)==8 & A.Dir(:,2)==0 & A.Dir(:,4)==d,:);
            [~,order]=sort(rr(:,3));rr=rr(order,:);
            if d==1,cc=blue;else,cc=orange;end
            plot(rr(1:4,7),rr(1:4,8),'-o','Color',cc, ...
                'LineWidth',1.4,'MarkerFaceColor','w');
            plot(rr(5:8,7),rr(5:8,8),'-o','Color',cc, ...
                'LineWidth',1.4,'MarkerFaceColor','w');
        end
        xlim([P.epsPYNegative(1) P.epsPYPositive(2)]);
        title('epsPY | direct response');
    else
        key=names{idx};range=P.(key);
        for sg=1:2
            for d=1:2
                rr=A.Dir(A.Dir(:,1)==idx & A.Dir(:,2)==sg & A.Dir(:,4)==d,:);
                [~,order]=sort(rr(:,3));rr=rr(order,:);
                if d==1,cc=blue;else,cc=orange;end
                if sg==1,style='-o';else,style='--o';end
                plot(rr(:,7),rr(:,8),style,'Color',cc,'LineWidth',1.3, ...
                    'MarkerSize',5,'MarkerFaceColor','w');
            end
        end
        if strcmpi(P.dist.(key),'loguniform'),set(gca,'XScale','log');end
        xlim(range);
        title(sprintf('%s | conditional on epsPY sign',key));
    end
    ylim([-.02 1.02]);ylabel('Direction probability');xlabel(names{idx});
    grid on;
end
nexttile(11);axis off;hold on;
h1=plot(NaN,NaN,'-o','Color',blue,'MarkerFaceColor','w','LineWidth',1.4);
h2=plot(NaN,NaN,'-o','Color',orange,'MarkerFaceColor','w','LineWidth',1.4);
h3=plot(NaN,NaN,'--o','Color',blue,'MarkerFaceColor','w','LineWidth',1.4);
h4=plot(NaN,NaN,'--o','Color',orange,'MarkerFaceColor','w','LineWidth',1.4);
legend([h1 h2 h3 h4],{'Increase | epsPY<0','Decrease | epsPY<0', ...
    'Increase | epsPY>0','Decrease | epsPY>0'}, ...
    'Location','northwest','FontSize',10,'Box','off');
text(.03,.44,'Increase: Y_1 > +0.05 per mil', ...
    'Units','normalized','FontSize',10);
text(.03,.34,'Decrease: Y_1 < -0.05 per mil', ...
    'Units','normalized','FontSize',10);
text(.03,.20,'For all other nine inputs: epsPY SIGN is fixed,', ...
    'Units','normalized','FontSize',9);
text(.03,.13,'but its magnitude and other priors still vary.', ...
    'Units','normalized','FontSize',9);
nexttile(12);axis off;
text(.02,.82,'One model, ten independent inputs', ...
    'Units','normalized','FontWeight','bold','FontSize',10);
text(.02,.66,'Direct epsPY panel: four bins per sign.', ...
    'Units','normalized','FontSize',9);
text(.02,.53,'All other panels: eight bins per parameter.', ...
    'Units','normalized','FontSize',9);
text(.02,.40,'Full active prior ranges on each x axis.', ...
    'Units','normalized','FontSize',9);
text(.02,.26,'Probabilities marginalize all remaining inputs.', ...
    'Units','normalized','FontSize',9);
sgtitle('Fig 5 | Single-scenario direction probabilities: epsPY and sign-conditioned parameter effects');
export_all(figdir,'Fig05_direction_probabilities_epsPY_sign_conditioned',fh,false);
fprintf('Exported Figures 2, 2b, 3, 4, 5 (PNG, PDF and SVG).\n');
end

function plot_one_bar(y,x,low,high,c)
if ~isfinite(x),return;end
w=.135;
patch([0 x x 0],[y-w y-w y+w y+w],c,'EdgeColor','none');
if isfinite(low)&&isfinite(high)
    plot([low high],[y y],'k-','LineWidth',.9);
    plot([low low],[y-.055 y+.055],'k-','LineWidth',.9);
    plot([high high],[y-.055 y+.055],'k-','LineWidth',.9);
end
end

function [mu,low,high]=replicate_ci(v)
v=v(isfinite(v));mu=NaN;low=NaN;high=NaN;
if numel(v)<2,return;end
mu=mean(v);err=tinv(.975,numel(v)-1)*std(v,0)/sqrt(numel(v));
low=mu-err;high=mu+err;
end

function export_all(figdir,name,fh,rasterPDF)
exportgraphics(fh,fullfile(figdir,[name '.png']),'Resolution',300);
if rasterPDF
    exportgraphics(fh,fullfile(figdir,[name '.pdf']), ...
        'ContentType','image','Resolution',600);
else
    exportgraphics(fh,fullfile(figdir,[name '.pdf']),'ContentType','vector');
end
set(fh,'PaperPositionMode','auto');
print(fh,fullfile(figdir,[name '.svg']),'-dsvg');
close(fh);
end
