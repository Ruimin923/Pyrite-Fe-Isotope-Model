function analyze_results(outdir)
%ANALYZE_RESULTS Publication tables for SL/YX Sobol and direction diagnostics.
% All direction probabilities are PRIOR PREDICTIVE, not fitted posteriors.
if nargin<1 || isempty(outdir),outdir=fullfile(fileparts(mfilename('fullpath')),'output');end
z=load(fullfile(outdir,'results.mat'),'samples','indices','config','names');
S=z.samples;I=z.indices;cfg=z.config;P=cfg.priors;
ciReps=1:cfg.nCI; K=10; targets=[.50 .80 .95];

% Direction probability for all targets and PY sign scenarios.
prob=zeros(2*3*3,8);np=0;
for e=1:2
    ss=S(S(:,1)==e,:);
    for t=1:3
        for sg=0:2
            sub=ss;
            if sg==1,sub=sub(sub(:,11)<0,:);end
            if sg==2,sub=sub(sub(:,11)>0,:);end
            v=sub(:,13+t);np=np+1;
            prob(np,:)=[e t sg numel(v) mean(v) mean(v>.05) mean(v<-.05) mean(abs(v)<=.05)];
        end
    end
end
writematrix(prob,fullfile(outdir,'probabilities.csv'));

% Per-scramble, fixed-N approximate 95% Student t confidence intervals.
% Separate Y1 and Y2 tables, both with the SAME 11-column table layout (10 parameters).
% Replicate-based intervals quantify scramble variability, not geological error.
SobolCI=sobolTable(I,cfg,outdir,1,K);
SobolCI_Y2=sobolTable(I,cfg,outdir,2,K);
writematrix(SobolCI,fullfile(outdir,'sobol_ci_sorted.csv'));
writematrix(SobolCI_Y2,fullfile(outdir,'sobol_ci_sorted_Y2.csv'));
writematrix(SobolCI(SobolCI(:,3)<=3,:),fullfile(outdir,'dominant_parameters.csv'));
writematrix(SobolCI_Y2(SobolCI_Y2(:,3)<=3,:),fullfile(outdir,'dominant_parameters_Y2.csv'));

% Across scenario (environment x target x output) ST/S1 heatmap data.
Heat=zeros(2*3*4*K,6);ih=0;
for e=1:2
    for t=1:3
        for out=1:4
            for idx=1:K
                m=I(:,1)==e&I(:,2)<=cfg.nCI&I(:,4)==t&I(:,5)==out&I(:,6)==idx;
                v=I(m,7:8);v=v(all(isfinite(v),2),:);
                if isempty(v),u=[NaN NaN];else,u=mean(v,1);end
                ih=ih+1;Heat(ih,:)=[e t out idx u];
            end
        end
    end
end
writematrix(Heat,fullfile(outdir,'scenario_sensitivity_heatmap.csv'));

% Direction probabilities for Fig. 5.
% epsPY is shown directly by environment. For f, Psi, beta and kappaOX,
% direction probabilities are recalculated AFTER fixing the epsPY endmember
% sign, i.e., using separate epsPY<0 and epsPY>0 subsets.
% Output columns in directional_responses.csv:
% env, parameterIndex, epsPYgroup(0 all / 1 neg / 2 pos), bin, direction,
% lower, upper, center, probability.
parAll=[8 1 2 3 10]; bins=8;
% epsPY: 2 env x 1 signGroup x 8 bins x 2 directions = 32 rows.
% Other four parameters: 2 env x 4 par x 2 signGroups x 8 bins x 2 directions = 256 rows.
Dir=nan(32+256,9);id=0;
for e=1:2
    ss=S(S(:,1)==e,:);y=ss(:,15);epsCol=11;
    for ik=1:numel(parAll)
        idx=parAll(ik);col=3+idx;
        if idx==8
            aNeg=linspace(P.epsPYNegative(1),P.epsPYNegative(2),5);
            aPos=linspace(P.epsPYPositive(1),P.epsPYPositive(2),5);
            edgesPairs=[aNeg(1:4)' aNeg(2:5)';aPos(1:4)' aPos(2:5)'];
            signGroups=0;
        else
            switch idx
                case 1
                    bounds=P.f;dist=P.dist.f;
                case 2
                    if e==1,bounds=P.Psi_SL;dist=P.dist.Psi_SL;
                    else,bounds=P.Psi_YX;dist=P.dist.Psi_YX;end
                case 3
                    if e==1,bounds=P.beta_SL;dist=P.dist.beta_SL;
                    else,bounds=P.beta_YX;dist=P.dist.beta_YX;end
                case 10
                    if e==1,bounds=P.kappaOX_SL;dist=P.dist.kappaOX_SL;
                    else,bounds=P.kappaOX_YX;dist=P.dist.kappaOX_YX;end
            end
            if bounds(1)==bounds(2)
                edgesPairs=nan(bins,2);edgesPairs(1,:)=bounds;
            else
                if strcmpi(dist,'loguniform')
                    ed=logspace(log10(bounds(1)),log10(bounds(2)),bins+1);
                else
                    ed=linspace(bounds(1),bounds(2),bins+1);
                end
                edgesPairs=[ed(1:end-1)' ed(2:end)'];
            end
            signGroups=[1 2];
        end
        assert(size(edgesPairs,1)==8);
        for sg=signGroups
            if sg==0
                baseMask=true(size(ss,1),1);
            elseif sg==1
                baseMask=ss(:,epsCol)<0;
            else
                baseMask=ss(:,epsCol)>0;
            end
            for b=1:bins
                a=edgesPairs(b,1);bb=edgesPairs(b,2);
                if ~isfinite(a) || ~isfinite(bb)
                    mask=false(size(ss,1),1);
                elseif idx~=8 && a==bb
                    mask=baseMask & abs(ss(:,col)-a)<=1e-12*max(1,abs(a));
                elseif b==bins || (idx==8 && b==4)
                    mask=baseMask & ss(:,col)>=a & ss(:,col)<=bb;
                else
                    mask=baseMask & ss(:,col)>=a & ss(:,col)<bb;
                end
                n=sum(mask);valid=15;if strcmp(cfg.mode,'quick'),valid=2;end
                for d=1:2
                    id=id+1;pp=NaN;
                    if n>=valid
                        if d==1,pp=mean(y(mask)>.05);else,pp=mean(y(mask)<-.05);end
                    end
                    if idx~=8 && exist('dist','var') && strcmpi(dist,'loguniform') && a>0 && bb>0
                        center=sqrt(a*bb);
                    else
                        center=(a+bb)/2;
                    end
                    Dir(id,:)=[e idx sg b d a bb center pp];
                end
            end
        end
    end
end
Dir=Dir(1:id,:);
writematrix(Dir,fullfile(outdir,'directional_responses.csv'));

% Reproducible representative cases selected to match the reported SIGN
% (not fitted isotope values), near median Y1 within eligible first-A cases.
Cases=zeros(2,16);Morph=[];
for e=1:2
    ens=load(fullfile(outdir,sprintf('ensemble_%d_1.mat',e)),'A','Y');
    yy=ens.Y(:,2,1,1); y2=ens.Y(:,2,2,1);
    if e==1,match=yy>.05;else,match=yy<-.05;end
    good=match & y2<.05;
    quality=2;
    if ~any(good),good=match;quality=1;end
    if ~any(good),good=true(size(yy));quality=0;end
    ix=find(good);[~,j]=min(abs(yy(ix)-median(yy(ix))));pick=ix(j);
    x=ens.A(pick,:);
    [y,~,curves,~,profile,cumulativeCurves]=fe_forward(x,1e-8,181);
    Cases(e,:)=[e 1 pick x y(2,1) y(2,2) quality];
    writematrix(profile',fullfile(outdir,sprintf('representative_profile_%d.csv',e)));
    writematrix([linspace(.05,.95,181)' curves'],...
        fullfile(outdir,sprintf('representative_delta_%d.csv',e)));
    writematrix([linspace(.05,.95,181)' cumulativeCurves'],...
        fullfile(outdir,sprintf('representative_cumulative_%d.csv',e)));
    writematrix([linspace(.05,.95,181)' curves(2,:)' cumulativeCurves(2,:)'],...
        fullfile(outdir,sprintf('representative_curve_%d.csv',e)));
    % Retain the morphology diagnostic on the first max 128 A cases.
    for sample=1:min(128,size(ens.A,1))
        [~,~,cc]=fe_forward(ens.A(sample,:),1e-8,181);
        for t=1:3
            Morph(end+1,:)=[e sample t classifyCurve(cc(t,:),.05)]; %#ok<AGROW>
        end
    end
end
writematrix(Cases,fullfile(outdir,'representative_cases.csv'));
writematrix(Morph,fullfile(outdir,'morphology.csv'));

% ALL A draws: sign-separated isotope envelopes for THREE Ptarget windows.
% This does NOT split reservoir evolution, or the matched example, by sign.
% Stored arrays instantA, cumulativeA have shape [N x 3 x nq].
% The cumulative delta is calculated from the P isotope INVENTORY ratio,
% not by averaging instantaneous delta values.
% New table columns:
% env, epsPYsign (1 negative, 2 positive), Ptarget, q,
% new median, new 2.5th, new 97.5th,
% cumulative median, cumulative 2.5th, cumulative 97.5th, subgroup n.
% The legacy table remains the t=2 (Ptarget=0.80) 10-column subset,
% so the first two Figure 1 columns and its third-column layout are intact.
q=linspace(.05,.95,cfg.nq)';
Ptargets=[.50 .80 .95];
EnsembleBandTargets=zeros(2*3*2*cfg.nq,11);
EnsembleBand=zeros(4*cfg.nq,10);
for e=1:2
    for t=1:3
        instantBySign={[],[]}; cumulativeBySign={[],[]};
        for rep=1:numel(cfg.NperRep)
            ens=load(fullfile(outdir,sprintf('ensemble_%d_%d.mat',e,rep)), ...
                'A','instantA','cumulativeA');
            N=size(ens.A,1);
            assert(isequal(size(ens.instantA),[N,3,cfg.nq]) && ...
                isequal(size(ens.cumulativeA),[N,3,cfg.nq]), ...
                'Expected all three Ptarget isotope windows for every A draw');
            neg=ens.A(:,8)<0;pos=ens.A(:,8)>0;
            assert(all(neg|pos),'Unexpected zero epsPY draw');
            instantSlice=reshape(ens.instantA(:,t,:),N,cfg.nq);
            cumulativeSlice=reshape(ens.cumulativeA(:,t,:),N,cfg.nq);
            for sg=1:2
                if sg==1,mask=neg;else,mask=pos;end
                instantBySign{sg}=[instantBySign{sg};instantSlice(mask,:)]; %#ok<AGROW>
                cumulativeBySign{sg}=[cumulativeBySign{sg};cumulativeSlice(mask,:)]; %#ok<AGROW>
            end
        end
        assert(size(instantBySign{1},1)+size(instantBySign{2},1)==sum(cfg.NperRep), ...
            'Missing A draws in the signed Ptarget isotope envelopes');
        for sg=1:2
            instant=instantBySign{sg};cumulativePY=cumulativeBySign{sg};
            assert(isequal(size(instant),size(cumulativePY)) && size(instant,1)>=3, ...
                'Bad sign-specific all-A isotope subgroup');
            qi=prctile(instant,[2.5 50 97.5],1);
            qc=prctile(cumulativePY,[2.5 50 97.5],1);
            rows=((e-1)*6+(t-1)*2+sg-1)*cfg.nq+(1:cfg.nq);
            EnsembleBandTargets(rows,:)=[repmat(e,cfg.nq,1) ...
                repmat(sg,cfg.nq,1) repmat(Ptargets(t),cfg.nq,1) q ...
                qi(2,:)' qi(1,:)' qi(3,:)' ...
                qc(2,:)' qc(1,:)' qc(3,:)' repmat(size(instant,1),cfg.nq,1)];
            if t==2
                rows80=((e-1)*2+sg-1)*cfg.nq+(1:cfg.nq);
                EnsembleBand(rows80,:)=EnsembleBandTargets(rows,[1 2 4:11]);
            end
        end
    end
end
writematrix(EnsembleBand,fullfile(outdir,'ensemble_pyrite_bands.csv'));
writematrix(EnsembleBandTargets,fullfile(outdir,'ensemble_pyrite_bands_targets.csv'));

% The user removed the path-wise attribution SUPPLEMENT figure.  The
% forward solver retains the optional diagnostic but no longer computes
% it or saves it as a manuscript result.
save(fullfile(outdir,'analysis_results.mat'),'prob','SobolCI','SobolCI_Y2',...
    'Heat','Dir','Cases','Morph','EnsembleBand','EnsembleBandTargets','-v7.3');
fprintf('Computed Y1/Y2 CIs, direction probabilities and ensemble 95%% prediction bands.\n');
end

function Result=sobolTable(I,cfg,outdir,out,K)
% All outputs: env, sign, ST rank, parameter, S1, S1low, S1high,
% ST, STlow, SThigh, number of usable CI scrambles.
Result=nan(2*3*K,11);c=0;
for e=1:2
    for sg=0:2
        data=nan(K,8);
        for idx=1:K
            est=nan(cfg.nCI,2);
            for rep=1:cfg.nCI
                if sg==0
                    m=I(:,1)==e&I(:,2)==rep&I(:,4)==2&I(:,5)==out&I(:,6)==idx;
                    rr=I(m,7:8);assert(size(rr,1)==1);est(rep,:)=rr;
                else
                    ens=load(fullfile(outdir,sprintf('ensemble_%d_%d.mat',e,rep)), ...
                        'A','B','Y');
                    if sg==1,mask=ens.A(:,8)<0&ens.B(:,8)<0;
                    else,mask=ens.A(:,8)>0&ens.B(:,8)>0;end
                    a=ens.Y(mask,2,out,1);b=ens.Y(mask,2,out,2);
                    ab=ens.Y(mask,2,out,idx+2);
                    if numel(a)>=3
                        vv=var([a;b],1);
                        if vv>1e-12
                            est(rep,:)=[1-mean((b-ab).^2)/(2*vv) ...
                                mean((a-ab).^2)/(2*vv)];
                        end
                    end
                end
            end
            [m1,l1,h1,n1]=tinterval(est(:,1));
            [mt,lt,ht,n2]=tinterval(est(:,2));assert(n1==n2);
            data(idx,:)=[idx m1 l1 h1 mt lt ht n1];
        end
        key=data(:,5);key(~isfinite(key))=-inf;
        [~,order]=sort(key,'descend');
        for rank=1:K
            c=c+1;Result(c,:)=[e sg rank data(order(rank),:)];
        end
    end
end
end

function [mu,lo,hi,n]=tinterval(v)
v=v(isfinite(v));n=numel(v);
mu=NaN;lo=NaN;hi=NaN;
if n<2,return;end
mu=mean(v);crit=tinv(.975,n-1);
err=crit*std(v,0)/sqrt(n);
lo=mu-err;hi=mu+err;
end

function cl=classifyCurve(y,h)
if max(y)-min(y)<=h,cl=0;return;end
lo=y(1);hi=y(1);direction=0;turns=0;first=0;
for i=2:numel(y)
    v=y(i);
    if direction==0
        lo=min(lo,v);hi=max(hi,v);
        if v-lo>h,direction=1;first=1;hi=v;
        elseif hi-v>h,direction=-1;first=-1;lo=v;end
    elseif direction==1
        hi=max(hi,v);if hi-v>h,turns=turns+1;direction=-1;lo=v;end
    else
        lo=min(lo,v);if v-lo>h,turns=turns+1;direction=1;hi=v;end
    end
end
if turns==0,cl=first;elseif turns==1,cl=2*first;else,cl=3;end
end
