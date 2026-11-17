function analyze_results(outdir)
%ANALYZE_RESULTS Single-scenario Y1/Y2 sensitivity and direction summaries.
% No SL/YX grouping. Conditional epsPY branches fix SIGN, not magnitude.
if nargin<1||isempty(outdir)
    outdir=fullfile(fileparts(mfilename('fullpath')),'output');
end
z=load(fullfile(outdir,'results.mat'),'samples','indices','config');
S=z.samples;I=z.indices;cfg=z.config;P=cfg.priors;K=10;

% Samples format: rep,draw,10 parameters,3x4 outputs (column-major),4 fluxes.
% epsPY = samples(:,10), Y1 at Ptarget=.80 = samples(:,14).
assert(size(S,2)==28 && size(I,2)==7,'This is not a single-scenario ensemble');
prob=zeros(9,8);i=0;
for t=1:3
    for sg=0:2
        mask=true(size(S,1),1);
        if sg==1,mask=S(:,10)<0;elseif sg==2,mask=S(:,10)>0;end
        y=S(mask,12+t);y2=S(mask,15+t);
        i=i+1;
        prob(i,:)=[t sg numel(y) mean(y) mean(y>.05) mean(y<-.05) ...
            mean(abs(y)<=.05) mean(y2)];
    end
end
writematrix(prob,fullfile(outdir,'probabilities.csv'));

SobolCI=sobol_table(I,cfg,outdir,1,K);
SobolCI_Y2=sobol_table(I,cfg,outdir,2,K);
writematrix(SobolCI,fullfile(outdir,'sobol_ci_sorted.csv'));
writematrix(SobolCI_Y2,fullfile(outdir,'sobol_ci_sorted_Y2.csv'));
writematrix(SobolCI(SobolCI(:,2)<=3,:), ...
    fullfile(outdir,'dominant_parameters.csv'));
writematrix(SobolCI_Y2(SobolCI_Y2(:,2)<=3,:), ...
    fullfile(outdir,'dominant_parameters_Y2.csv'));

% target, output, parameter index, mean S1, mean ST
Heat=nan(3*4*K,5);row=0;
for t=1:3
    for out=1:4
        for idx=1:K
            rr=I(I(:,1)<=cfg.nCI & I(:,3)==t & I(:,4)==out & ...
                I(:,5)==idx,6:7);
            rr=rr(all(isfinite(rr),2),:);
            if isempty(rr),v=[NaN NaN];else,v=mean(rr,1);end
            row=row+1;Heat(row,:)=[t out idx v];
        end
    end
end
writematrix(Heat,fullfile(outdir,'scenario_sensitivity_heatmap.csv'));

% Fig. 5: (a) direct epsPY effect, with positive and negative intervals;
% (b-j) all other nine parameters, separately conditional on epsPY sign.
% Columns: parameterIndex, signGroup(0=epsPY direct, 1=neg, 2=pos),
% bin, direction(1=increase,2=decrease), lower, upper, center, probability, n.
Dir=nan(16+9*2*8*2,9);r=0;
shown=[8 1 2 3 10 4 5 6 7 9];bins=8;
for idx=shown
    paramCol=2+idx;
    if idx==8
        negEdge=linspace(P.epsPYNegative(1),P.epsPYNegative(2),5);
        posEdge=linspace(P.epsPYPositive(1),P.epsPYPositive(2),5);
        edges=[negEdge(1:4)' negEdge(2:5)';posEdge(1:4)' posEdge(2:5)'];
        signs=0;
    else
        key=P.names{idx};b=P.(key);distribution=P.dist.(key);
        if strcmpi(distribution,'loguniform')
            ed=logspace(log10(b(1)),log10(b(2)),bins+1);
        else
            ed=linspace(b(1),b(2),bins+1);
        end
        edges=[ed(1:end-1)' ed(2:end)'];
        signs=[1 2];
    end
    for sg=signs
        if sg==0,subset=true(size(S,1),1);
        elseif sg==1,subset=S(:,10)<0;
        else,subset=S(:,10)>0;end
        for b=1:bins
            a=edges(b,1);bb=edges(b,2);
            if idx==8
                % Negative and positive ranges each have their own last bin.
                includeUpper=ismember(b,[4 8]);
            else
                includeUpper=b==bins;
            end
            if includeUpper
                inside=S(:,paramCol)>=a & S(:,paramCol)<=bb;
            else
                inside=S(:,paramCol)>=a & S(:,paramCol)<bb;
            end
            mask=subset & inside;n=sum(mask);
            valid=15;if strcmp(cfg.mode,'quick'),valid=2;end
            if idx~=8 && strcmpi(distribution,'loguniform')
                center=sqrt(a*bb);
            else
                center=(a+bb)/2;
            end
            for d=1:2
                r=r+1;probability=NaN;
                if n>=valid
                    if d==1,probability=mean(S(mask,14)>.05);
                    else,probability=mean(S(mask,14)<-.05);end
                end
                Dir(r,:)=[idx sg b d a bb center probability n];
            end
        end
    end
end
assert(r==size(Dir,1));
writematrix(Dir,fullfile(outdir,'directional_responses.csv'));
save(fullfile(outdir,'analysis_results.mat'),'prob','SobolCI','SobolCI_Y2', ...
    'Heat','Dir','-v7.3');
fprintf('Single-scenario postprocessing: Y1/Y2, 4 Sobol outputs, epsPY-sign direction probabilities.\n');
end

function Result=sobol_table(I,cfg,outdir,out,K)
% scenario (0 full,1 negative,2 positive), rank, parameter index,
% S1, lower, upper, ST, lower, upper, n usable scrambles.
Result=nan(3*K,10);row=0;
for sg=0:2
    dat=nan(K,8);
    for idx=1:K
        est=nan(cfg.nCI,2);
        for rep=1:cfg.nCI
            if sg==0
                rr=I(I(:,1)==rep & I(:,3)==2 & I(:,4)==out & I(:,5)==idx,6:7);
                assert(size(rr,1)==1);
                est(rep,:)=rr;
            else
                E=load(fullfile(outdir,sprintf('ensemble_%d.mat',rep)),'A','B','Y');
                if sg==1,mask=E.A(:,8)<0 & E.B(:,8)<0;
                else,mask=E.A(:,8)>0 & E.B(:,8)>0;end
                a=E.Y(mask,2,out,1);b=E.Y(mask,2,out,2);
                ab=E.Y(mask,2,out,idx+2);
                if numel(a)>=3
                    V=var([a;b],1);
                    if V>1e-12
                        est(rep,:)=[1-mean((b-ab).^2)/(2*V), ...
                            mean((a-ab).^2)/(2*V)];
                    end
                end
            end
        end
        [s1,l1,h1,n1]=replicate_ci(est(:,1));
        [st,lt,ht,n2]=replicate_ci(est(:,2));
        assert(n1==n2);
        dat(idx,:)=[idx s1 l1 h1 st lt ht n1];
    end
    ranking=dat(:,5);ranking(~isfinite(ranking))=-inf;
    [~,order]=sort(ranking,'descend');
    for rank=1:K
        row=row+1;Result(row,:)=[sg rank dat(order(rank),:)];
    end
end
end

function [mu,low,high,n]=replicate_ci(v)
v=v(isfinite(v));n=numel(v);mu=NaN;low=NaN;high=NaN;
if n<2,return;end
mu=mean(v);error=tinv(.975,n-1)*std(v,0)/sqrt(n);
low=mu-error;high=mu+error;
end
