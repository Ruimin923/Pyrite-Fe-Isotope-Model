function make_figures(outdir)
%MAKE_FIGURES Manuscript plots: Figure 1, 2 (Y1), 2b (Y2), 3, 4, 5.
% Figures are PRIOR-PREDICTIVE. No observed Fe isotope data fitted.
if nargin<1 || isempty(outdir),outdir=fullfile(fileparts(mfilename('fullpath')),'output');end
figdir=fullfile(outdir,'figures');if ~exist(figdir,'dir'),mkdir(figdir);end
Z=load(fullfile(outdir,'results.mat'),'names','config','indices','conv');
A=load(fullfile(outdir,'analysis_results.mat'),...
    'SobolCI','SobolCI_Y2','Heat','Dir','Cases','EnsembleBand','EnsembleBandTargets');
names=Z.names;cfg=Z.config;labels={'SL','YX'};
blue=[.16 .40 .68];orange=[.86 .37 .19];teal=[.10 .56 .45];
set(groot,'defaultAxesFontName','Arial','defaultAxesFontSize',10,'defaultFigureColor','white');

% Remove deprecated figures from earlier runs in the same output folder.
for oldName={'Fig03b_sign_scenarios_heatmap','FigS01_reaction_contributions', ...
        'Fig05_pairwise_interactions','Fig06_direction_probabilities_f_H_epsPY_Psi', ...
        'Fig05_direction_probabilities_f_H_epsPY_Psi','Fig05_direction_probabilities_f_H_epsPY_Psi_kappaOX'}
    for ext={'.png','.pdf','.svg'}
        obsolete=fullfile(figdir,[oldName{1} ext{1}]);
        if isfile(obsolete),delete(obsolete);end
    end
end

% FIG 01: first TWO columns unchanged: reservoirs and one matched case.
% ONLY the THIRD column splits all A-draw pyrite Fe-isotope trajectories
% into negative/positive epsPY subgroups within each SL/YX axis.
% Each sign has instant (solid) and cumulative (dashed) median plus
% pointwise central 95% prior-predictive intervals (not CIs/errors).
f=figure('Visible','off','Position',[20 20 1510 790]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for e=1:2
    p=readmatrix(fullfile(outdir,sprintf('representative_profile_%d.csv',e)));
    d=readmatrix(fullfile(outdir,sprintf('representative_curve_%d.csv',e)));
    cr=A.Cases(A.Cases(:,1)==e,:);quality=cr(end);
    nexttile;plot(p(:,1),p(:,2:4),'LineWidth',1.6);
    xlabel('Dimensionless time');ylabel('Normalized Fe inventory');
    title(sprintf('%s  |  O, L and P (matched case)',labels{e}));
    legend({'O: Fe(III)','L: Fe(II)','P: pyrite'},'Location','best');
    nexttile;hold on;
    plot(d(:,1),d(:,2),'-','Color',blue,'LineWidth',1.8);
    plot(d(:,1),d(:,3),'-','Color',orange,'LineWidth',1.8);
    xlabel('q = P / P_{target}');ylabel('Pyrite \delta^{56}Fe (per mil)','Interpreter','tex');
    title(sprintf('%s  |  direction-matched example (flag %d)',labels{e},quality));
    legend({'Instantaneous new pyrite','Cumulative bulk pyrite'},'Location','best');
    nexttile;hold on;
    signColor=[blue;orange];
    legendHandles=gobjects(4,1);
    legendNames={ 'epsPY < 0: new', 'epsPY < 0: cumulative', ...
        'epsPY > 0: new', 'epsPY > 0: cumulative' };
    counts=zeros(1,2);
    for sg=1:2
        % Extract all A draws for THIS environment and THIS epsPY sign.
        band=A.EnsembleBand(A.EnsembleBand(:,1)==e & ...
            A.EnsembleBand(:,2)==sg,:);
        assert(size(band,1)==cfg.nq,'Missing signed all-A isotope envelope');
        q=band(:,3);counts(sg)=round(band(1,10));
        color=signColor(sg,:);
        % New pyrite: solid median + central 95% predictive band.
        hNew=fill([q;flipud(q)],[band(:,5);flipud(band(:,6))], ...
            color,'FaceAlpha',.14,'EdgeColor','none');
        set(hNew,'HandleVisibility','off');
        % Cumulative pyrite: dashed median + its OWN predictive band.
        hCum=fill([q;flipud(q)],[band(:,8);flipud(band(:,9))], ...
            color,'FaceAlpha',.08,'EdgeColor','none');
        set(hCum,'HandleVisibility','off');
        legendHandles((sg-1)*2+1)=plot(q,band(:,4),'-', ...
            'Color',color,'LineWidth',1.8);
        legendHandles((sg-1)*2+2)=plot(q,band(:,7),'--', ...
            'Color',color,'LineWidth',1.8);
    end
    xlabel('q = P / P_{target}');
    ylabel('Pyrite \delta^{56}Fe (per mil)','Interpreter','tex');
    title({sprintf('%s  |  all A draws by epsPY sign',labels{e}), ...
        sprintf('negative n=%d; positive n=%d',counts(1),counts(2))});
    legend(legendHandles,legendNames,'Location','best','FontSize',7);

end
sgtitle('Fig 1  |  Reservoir evolution, sign-compatible example and prior-predictive envelopes (Ptarget = 0.80)');
emit(figdir,'Fig01_representative_trajectories',f);

% SUPPLEMENTARY FIG S1: compare Ptarget=0.50, 0.80 and 0.95 in one
% 2 x 3 figure (rows SL/YX; columns Ptarget). ALL A draws are included.
% Within EACH panel, epsPY<0 and epsPY>0 are separate; each sign has
% instantaneous/new (solid) and cumulative/bulk (dashed) median and
% its OWN pointwise central 95% PRIOR-PREDICTIVE interval (not CI).
f=figure('Visible','off','Position',[20 20 1660 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for e=1:2
    for t=1:3
        nexttile;hold on;
        hh=gobjects(4,1);
        counts=zeros(1,2);
        for sg=1:2
            dat=A.EnsembleBandTargets( ...
                A.EnsembleBandTargets(:,1)==e & A.EnsembleBandTargets(:,2)==sg & ...
                abs(A.EnsembleBandTargets(:,3)-cfg.targets(t))<1e-12,:);
            assert(size(dat,1)==cfg.nq,'Incomplete Supplementary Figure S1 data');
            qq=dat(:,4);cc=blue;if sg==2,cc=orange;end
            counts(sg)=dat(1,11);
            hfill=fill([qq;flipud(qq)],[dat(:,6);flipud(dat(:,7))],cc, ...
                'FaceAlpha',.14,'EdgeColor','none');
            set(hfill,'HandleVisibility','off');
            hfill=fill([qq;flipud(qq)],[dat(:,9);flipud(dat(:,10))],cc, ...
                'FaceAlpha',.07,'EdgeColor','none');
            set(hfill,'HandleVisibility','off');
            hh((sg-1)*2+1)=plot(qq,dat(:,5),'-','Color',cc,'LineWidth',1.8);
            hh((sg-1)*2+2)=plot(qq,dat(:,8),'--','Color',cc,'LineWidth',1.8);
        end
        title(sprintf('%s | P_{target}=%.2f | n_{-}=%d, n_{+}=%d', ...
            labels{e},cfg.targets(t),counts(1),counts(2)),'Interpreter','tex');
        xlabel('q = P / P_{target}');
        ylabel('Pyrite \delta^{56}Fe (per mil)','Interpreter','tex');
        if e==1 && t==3
            legend(hh,{'epsPY < 0: new','epsPY < 0: cumulative', ...
                'epsPY > 0: new','epsPY > 0: cumulative'}, ...
                'Location','best','FontSize',8);
        end
        xlim([.05 .95]);
    end
end
sgtitle('Supplementary Fig S1 | All A draws: instant vs cumulative pyrite, separated by epsPY sign and Ptarget');
emit(figdir,'FigS01_Ptarget_isotope_envelopes',f);

% FIG 02: S1/ST with 95% t CIs; all 10 inputs ordered by mean ST.
f=figure('Visible','off','Position',[20 20 1620 1070]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
scenarioNames={'all priors','epsPY < 0','epsPY > 0'};
for e=1:2
    for sg=0:2
        nexttile;hold on;
        rows=A.SobolCI(A.SobolCI(:,1)==e&A.SobolCI(:,2)==sg,:);
        assert(size(rows,1)==10,'Incomplete Sobol CI table');
        for r=1:10
            % columns: e sign rank idx S1 lo hi ST lo hi nRep
            plotOneBar(r-.18,rows(r,5),rows(r,6),rows(r,7),blue);
            plotOneBar(r+.18,rows(r,8),rows(r,9),rows(r,10),orange);
        end
        axis tight;
        set(gca,'YDir','reverse','YTick',1:10,'YTickLabel',names(rows(:,4)));
        ylim([.45 10.55]);xlabel('Sobol index (95% t CI)');
        title(sprintf('%s: %s',labels{e},scenarioNames{sg+1}));
        xline(0,':','Color',[.3 .3 .3]);
        h1=plot(NaN,NaN,'s','MarkerFaceColor',blue,'Color',blue);
        h2=plot(NaN,NaN,'s','MarkerFaceColor',orange,'Color',orange);
        legend([h1 h2],{'S_1','S_T'},'Location','best','FontSize',8);
        grid on;set(gca,'YGrid','off');
    end
end
sgtitle('Fig 2  |  Y1, Ptarget=0.80; descending S_T; CIs from independent scrambles');
emit(figdir,'Fig02_sobol_sorted_95CI',f);

% FIG 02b: EXACTLY the Fig 02 design for Y2 at Ptarget=.80.
% Y2 is non-monotonicity, not geometric curvature.  If Var(Y2)~0,
% indices/intervals remain NaN (no fabricated zero-sensitivity bars).
f=figure('Visible','off','Position',[20 20 1620 1070]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for e=1:2
    for sg=0:2
        nexttile;hold on;
        rows=A.SobolCI_Y2(A.SobolCI_Y2(:,1)==e&A.SobolCI_Y2(:,2)==sg,:);
        assert(size(rows,1)==10,'Incomplete Y2 Sobol CI table');
        for r=1:10
            plotOneBar(r-.18,rows(r,5),rows(r,6),rows(r,7),blue);
            plotOneBar(r+.18,rows(r,8),rows(r,9),rows(r,10),orange);
        end
        if all(~isfinite(rows(:,5))) && all(~isfinite(rows(:,8)))
            xlim([-.1 1.0]);
            text(.03,.50,'Undefined: Var(Y2) below threshold', ...
                'Units','normalized','Color',[.35 .35 .35],'FontSize',9);
        else
            axis tight;
        end
        set(gca,'YDir','reverse','YTick',1:10,'YTickLabel',names(rows(:,4)));
        ylim([.45 10.55]);xlabel('Sobol index (95% t CI)');
        title(sprintf('%s: %s',labels{e},scenarioNames{sg+1}));
        xline(0,':','Color',[.3 .3 .3]);
        h1=plot(NaN,NaN,'s','MarkerFaceColor',blue,'Color',blue);
        h2=plot(NaN,NaN,'s','MarkerFaceColor',orange,'Color',orange);
        legend([h1 h2],{'S_1','S_T'},'Location','best','FontSize',8);
        grid on;set(gca,'YGrid','off');
    end
end
sgtitle('Fig 2b  |  Y2, Ptarget=0.80; descending S_T; CIs from independent scrambles');
emit(figdir,'Fig02b_Y2_sobol_sorted_95CI',f);

% FIG 03: all outputs and SL/YX / three Ptarget values, fixed parameter order.
f=figure('Visible','off','Position',[20 20 1440 890]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
outNames={'Y1','Y2','Dplus','Dminus'};
for out=1:4
    nexttile;
    data=nan(6,10);
    for e=1:2
        for t=1:3
            for idx=1:10
                rr=A.Heat(A.Heat(:,1)==e&A.Heat(:,2)==t&A.Heat(:,3)==out&A.Heat(:,4)==idx,:);
                data((e-1)*3+t,idx)=rr(6);
            end
        end
    end
    hh=imagesc(data);set(hh,'AlphaData',isfinite(data));
    set(gca,'Color',[.9 .9 .9],'XTick',1:10,'XTickLabel',names,...
        'YTick',1:6,'YTickLabel',{'SL .50','SL .80','SL .95','YX .50','YX .80','YX .95'});
    xtickangle(40);clim([0 1]);cb=colorbar;cb.Label.String='S_T';
    title(outNames{out});
end
sgtitle('Fig 3  |  Total-effect sensitivity by environment, target and response');
emit(figdir,'Fig03_sensitivity_heatmap',f);

% FIG 04: convergence for top three unconditional Y1 parameters at .80.
f=figure('Visible','off','Position',[20 20 1280 830]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
palette=[blue;orange;teal];
for e=1:2
    rr=A.SobolCI(A.SobolCI(:,1)==e&A.SobolCI(:,2)==0&A.SobolCI(:,3)<=3,:);
    chosen=rr(:,4)';
    for out=1:2  % S1 and ST, respectively
        nexttile;hold on;
        for j=1:numel(chosen)
            idx=chosen(j);base=Z.conv(Z.conv(:,1)==e&Z.conv(:,6)==idx,:);
            if strcmp(cfg.mode,'quick'),nvec=[8 16 32 64];else,nvec=[256 512 1024 2048];end
            for nn=nvec
                selected=base(base(:,3)==nn,:);
                if nn==nvec(end)
                    % Independent largest-N check: one run, NO CI.
                    if isempty(selected),continue;end
                    r=selected(end,6+out);
                    if isfinite(r),plot(nn,r,'d','MarkerSize',7,'LineWidth',1.3, ...
                            'Color',palette(j,:),'MarkerFaceColor',palette(j,:));end
                else
                    selected=selected(selected(:,2)<=cfg.nCI,:);
                    [mu,lo,hi,~]=tinterval(selected(:,6+out));
                    if isfinite(mu)
                        errorbar(nn,mu,mu-lo,hi-mu,'o','MarkerSize',5,...
                            'Color',palette(j,:),'MarkerFaceColor',palette(j,:),'LineWidth',1.1);
                    end
                end
            end
            h(j)=plot(NaN,NaN,'-o','Color',palette(j,:),'MarkerFaceColor',palette(j,:)); %#ok<AGROW>
        end
        set(gca,'XScale','log','XTick',nvec,'XTickLabel',cellstr(num2str(nvec')));
        xlim([min(nvec)*.85 max(nvec)*1.15]);yline(0,':');
        xlabel('Sobol sample N per scramble');
        if out==1,ylabel('S_1');else,ylabel('S_T');end
        if out==1,kind='S1';else,kind='ST';end
        title(sprintf('%s  |  %s',labels{e},kind));
        legend(h,names(chosen),'Location','best','FontSize',8);
    end
end
sgtitle('Fig 4  |  95% t CIs at equal N; diamonds = independent larger N, no CI');
emit(figdir,'Fig04_sobol_convergence',f);

% FIG 05: direction probabilities in a 4 x 3 layout.
% Rows 1-2: SL/YX for epsPY, f and Psi.
% Rows 3-4: beta and kappaOX, with epsPY sign fixed separately.
% For f, Psi, beta and kappaOX, solid lines use epsPY<0 and dashed
% lines use epsPY>0. Blue=Increase (Y1>0.05), orange=Decrease (Y1<-0.05).
f=figure('Visible','off','Position',[20 20 1560 1440]);
tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
panelMap={ ...
    1,8,'eps';  1,1,'cond'; 1,2,'cond'; ...
    2,8,'eps';  2,1,'cond'; 2,2,'cond'; ...
    1,3,'cond'; 2,3,'cond'; 0,0,'legend'; ...
    1,10,'cond';2,10,'cond';0,0,'note'};
for tile=1:size(panelMap,1)
    e=panelMap{tile,1}; idx=panelMap{tile,2}; mode=panelMap{tile,3};
    nexttile(tile); hold on;
    if strcmp(mode,'legend')
        axis off;
        h1=plot(NaN,NaN,'-o','Color',blue,'MarkerFaceColor','w','LineWidth',1.3);
        h2=plot(NaN,NaN,'-o','Color',orange,'MarkerFaceColor','w','LineWidth',1.3);
        h3=plot(NaN,NaN,'--o','Color',blue,'MarkerFaceColor','w','LineWidth',1.3);
        h4=plot(NaN,NaN,'--o','Color',orange,'MarkerFaceColor','w','LineWidth',1.3);
        legend([h1 h2 h3 h4],{'Increase | epsPY < 0','Decrease | epsPY < 0', ...
            'Increase | epsPY > 0','Decrease | epsPY > 0'}, ...
            'Location','northwest','Box','off','FontSize',9);
        text(0.02,0.35,'For f, Psi, beta and kappaOX, each panel is recalculated', ...
            'Units','normalized','FontSize',9);
        text(0.02,0.26,'after fixing the epsPY branch sign.', ...
            'Units','normalized','FontSize',9);
        text(0.02,0.10,'epsPY panels show the direct effect of epsPY itself.', ...
            'Units','normalized','FontSize',9);
        continue;
    elseif strcmp(mode,'note')
        axis off;
        text(0.02,0.82,'Probability definitions', 'Units','normalized', ...
            'FontWeight','bold','FontSize',10);
        text(0.02,0.68,'Increase:  Y_1 > 0.05 per mil', ...
            'Units','normalized','FontSize',9);
        text(0.02,0.58,'Decrease:  Y_1 < -0.05 per mil', ...
            'Units','normalized','FontSize',9);
        text(0.02,0.42,'Other priors are marginalized within each', ...
            'Units','normalized','FontSize',9);
        text(0.02,0.33,'environment and fixed epsPY subset.', ...
            'Units','normalized','FontSize',9);
        continue;
    end

    if strcmp(mode,'eps')
        for d=1:2
            dat=A.Dir(A.Dir(:,1)==e & A.Dir(:,2)==idx & A.Dir(:,3)==0 & A.Dir(:,5)==d,:);
            [~,ord]=sort(dat(:,4));dat=dat(ord,:);
            if d==1,cc=blue;else,cc=orange;end
            plot(dat(1:4,8),dat(1:4,9),'-o','Color',cc,'MarkerFaceColor','w','LineWidth',1.2);
            plot(dat(5:8,8),dat(5:8,9),'-o','Color',cc,'MarkerFaceColor','w','LineWidth',1.2);
        end
        ylabel('Direction probability');xlabel(names{idx});
        title(sprintf('%s | %s',labels{e},names{idx}));
    else
        for sg=1:2
            for d=1:2
                dat=A.Dir(A.Dir(:,1)==e & A.Dir(:,2)==idx & A.Dir(:,3)==sg & A.Dir(:,5)==d,:);
                [~,ord]=sort(dat(:,4));dat=dat(ord,:);
                if d==1,cc=blue;else,cc=orange;end
                if sg==1,ls='-';else,ls='--';end
                plot(dat(:,8),dat(:,9),[ls 'o'],'Color',cc,'MarkerFaceColor','w','LineWidth',1.2);
            end
        end
        ylim([-.02 1.02]);
        if idx==2
            if e==1,dst=cfg.priors.dist.Psi_SL;else,dst=cfg.priors.dist.Psi_YX;end
        elseif idx==3
            if e==1,dst=cfg.priors.dist.beta_SL;else,dst=cfg.priors.dist.beta_YX;end
        elseif idx==10
            if e==1,dst=cfg.priors.dist.kappaOX_SL;else,dst=cfg.priors.dist.kappaOX_YX;end
        else
            dst='uniform';
        end
        if strcmpi(dst,'loguniform')
            v=A.Dir(A.Dir(:,1)==e & A.Dir(:,2)==idx & isfinite(A.Dir(:,8)),8);
            if ~isempty(v) && all(v>0),set(gca,'XScale','log');end
        end
        ylabel('Direction probability');xlabel(names{idx});
        title(sprintf('%s | %s | fixed epsPY sign',labels{e},names{idx}));
    end
    ylim([-.02 1.02]);
    grid on;
end
sgtitle('Fig 5 | Direction probabilities after fixing epsPY endmembers for f, Psi, beta and kappaOX');
emit(figdir,'Fig05_direction_probabilities_4x3_fixed_epsPY',f);

fprintf('Exported Fig 1, 2, 2b, 3, 4, 5 and Supplementary Fig S1 as PNG/PDF/SVG.\n');
end

function plotOneBar(y,v,lower,upper,color)
if ~isfinite(v),return;end
w=.135;
patch([0 v v 0],[y-w y-w y+w y+w],color,'EdgeColor','none');
if isfinite(lower)&&isfinite(upper)
    plot([lower upper],[y y],'k-','LineWidth',.9);
    plot([lower lower],[y-.055 y+.055],'k-','LineWidth',.9);
    plot([upper upper],[y-.055 y+.055],'k-','LineWidth',.9);
end
end

function [mu,lo,hi,n]=tinterval(v)
v=v(isfinite(v));n=numel(v);mu=NaN;lo=NaN;hi=NaN;
if n<2,return;end
mu=mean(v);err=tinv(.975,n-1)*std(v,0)/sqrt(n);lo=mu-err;hi=mu+err;
end

function emit(figdir,base,fh)
% MATLAB R2024b: do NOT use print(...,'-dpdf','-vector') for large
% figures. print uses a fixed paper size and may crop tiledlayouts.
% exportgraphics crops to the figure content, avoiding paper clipping.
% Complex translucent uncertainty bands (Fig 1 and Supplementary S1),
% and the dense heatmap (Fig 3), are exported to PDF as 600-dpi images.
% Forcing them to vector can trigger MATLAB's ContentTypeImageSuggested
% warning and lead to very large PDFs. All other PDFs remain vector.
% SVG is retained as the separate editable vector format in MATLAB R2024b.
exportgraphics(fh,fullfile(figdir,[base '.png']),'Resolution',300);
pdfFile=fullfile(figdir,[base '.pdf']);
rasterPDF={'Fig01_representative_trajectories', ...
           'FigS01_Ptarget_isotope_envelopes', ...
           'Fig03_sensitivity_heatmap'};
if any(strcmp(base,rasterPDF))
    exportgraphics(fh,pdfFile,'ContentType','image','Resolution',600);
else
    exportgraphics(fh,pdfFile,'ContentType','vector');
end
% SVG is a non-page vector format. Match its exported size to the actual
% on-screen figure instead of relying on default paper dimensions.
set(fh,'PaperPositionMode','auto');
print(fh,fullfile(figdir,[base '.svg']),'-dsvg');
close(fh);
end
