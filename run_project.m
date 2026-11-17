function run_project(outdir,mode)
%RUN_PROJECT Single-scenario, 10-input scrambled Sobol/Jansen forward design.
% No environment loop, no pairwise interaction AB_ij calculations.
if nargin<1||isempty(outdir)
    outdir=fullfile(fileparts(mfilename('fullpath')),'output');
end
if nargin<2||isempty(mode),mode='full';end
if ~exist(outdir,'dir'),mkdir(outdir);end
P=model_parameters();K=numel(P.names);assert(K==10);
cfg.mode=char(mode);
cfg.priors=P;cfg.modelVersion='single_generic_noH_10D_v1';
cfg.targets=[.50 .80 .95];cfg.nq=181;
cfg.seedBase=700;
if strcmp(mode,'quick')
    cfg.NperRep=[32 32 32 64];cfg.nPrecision=4;
else
    cfg.NperRep=[1024 1024 1024 1024 1024 1024 2048];cfg.nPrecision=32;
end
cfg.nCI=numel(cfg.NperRep)-1;
cfg.primaryRuns=sum(cfg.NperRep)*(K+2);
assert(cfg.nCI>=3);
assert(license('test','Statistics_Toolbox')==1, ...
    'The Statistics and Machine Learning Toolbox is required.');
summary=zeros(numel(cfg.NperRep)*3*4,7);
indices=zeros(numel(cfg.NperRep)*3*4*K,7);
diagnostics=zeros(cfg.primaryRuns,9);
samples=zeros(sum(cfg.NperRep),28);
conv=[];si=0;ii=0;di=0;ai=0;
write_snapshot(fullfile(outdir,'parameters_used.csv'),P);
fid=fopen(fullfile(outdir,'execution.log'),'w');
assert(fid>=0);autoClose=onCleanup(@()fclose(fid));
fprintf(fid,'MATLAB %s, mode %s, 1 scenario / 10 parameters\n',version,cfg.mode);
for rep=1:numel(cfg.NperRep)
    N=cfg.NperRep(rep);
    rng(cfg.seedBase+rep,'twister');
    U=net(scramble(sobolset(2*K),'MatousekAffineOwen'),N);
    A=transform_parameters(U(:,1:K),P);
    B=transform_parameters(U(:,K+1:2*K),P);
    Y=zeros(N,3,4,K+2);
    instantA=zeros(N,3,cfg.nq);cumulativeA=zeros(N,3,cfg.nq);
    for block=1:K+2
        X=A;
        if block==2
            X=B;
        elseif block>=3
            idx=block-2;X(:,idx)=B(:,idx);
        end
        for j=1:N
            if block==1
                [y,check,curve,flux,~,cumCurve]=fe_forward(X(j,:),1e-8,cfg.nq);
                instantA(j,:,:)=reshape(curve,1,3,cfg.nq);
                cumulativeA(j,:,:)=reshape(cumCurve,1,3,cfg.nq);
            else
                [y,check,~,flux]=fe_forward(X(j,:),1e-8,cfg.nq);
            end
            Y(j,:,:,block)=reshape(y,1,3,4);
            di=di+1;diagnostics(di,:)=[rep block j check];
            if block==1
                ai=ai+1;
                samples(ai,:)=[rep j X(j,:) reshape(y,1,[]) flux];
            end
        end
        fprintf(fid,'rep=%d/%d block=%d/%d done\n',rep,numel(cfg.NperRep),block,K+2);
    end
    for t=1:3
        for out=1:4
            ya=Y(:,t,out,1);yb=Y(:,t,out,2);vv=var([ya;yb],1);
            si=si+1;summary(si,:)=[rep N t out mean(ya) std(ya) vv];
            for idx=1:K
                [s1,st]=jansen(ya,yb,Y(:,t,out,idx+2));
                ii=ii+1;indices(ii,:)=[rep N t out idx s1 st];
            end
        end
    end
    if strcmp(mode,'quick'),prefix=[8 16 32 64];
    else,prefix=[256 512 1024 2048];end
    for np=prefix(prefix<=N)
        for idx=1:K
            [s1,st]=jansen(Y(1:np,2,1,1),Y(1:np,2,1,2), ...
                Y(1:np,2,1,idx+2));
            conv(end+1,:)=[rep np 2 1 idx s1 st]; %#ok<AGROW>
        end
    end
    save(fullfile(outdir,sprintf('ensemble_%d.mat',rep)), ...
        'A','B','Y','U','instantA','cumulativeA','-v7.3');
    fprintf('Single scenario: rep %d/%d, N=%d complete\n', ...
        rep,numel(cfg.NperRep),N);
end
assert(di==cfg.primaryRuns&&ai==size(samples,1)&&si==size(summary,1) ...
    &&ii==size(indices,1),'Incomplete Sobol design');
precision=zeros(cfg.nPrecision,14);
rng(991,'twister');
Pdraw=transform_parameters(net(scramble(sobolset(K),'MatousekAffineOwen'), ...
    cfg.nPrecision),P);
for j=1:cfg.nPrecision
    [y,c]=fe_forward(Pdraw(j,:),1e-8,181);
    [yt,ct]=fe_forward(Pdraw(j,:),1e-10,181);
    precision(j,:)=[j max(abs(y(:)-yt(:))) c ct];
end
writematrix(summary,fullfile(outdir,'summary.csv'));
writematrix(indices,fullfile(outdir,'sobol.csv'));
writematrix(diagnostics,fullfile(outdir,'conservation.csv'));
writematrix(samples,fullfile(outdir,'samples.csv'));
writematrix(precision,fullfile(outdir,'precision.csv'));
writematrix(conv,fullfile(outdir,'sobol_convergence.csv'));
names=P.names;targets=cfg.targets;outputNames={'Y1','Y2','Dplus','Dminus'};
config=cfg;
save(fullfile(outdir,'results.mat'),'summary','indices','diagnostics', ...
    'samples','precision','names','targets','outputNames','config','conv','-v7.3');
end

function [s1,st]=jansen(a,b,ab)
V=var([a;b],1);
if ~isfinite(V)||V<=1e-12,s1=NaN;st=NaN;return;end
st=mean((a-ab).^2)/(2*V);
s1=1-mean((b-ab).^2)/(2*V);
end

function X=transform_parameters(U,P)
assert(size(U,2)==10);
X=zeros(size(U));
keys={'f','Psi','beta','deltaO','deltaL','epsDIR','epsAIR', ...
      'epsPY','epsOX','kappaOX'};
for idx=[1:7 9:10]
    key=keys{idx};X(:,idx)=draw_prior(U(:,idx),P.(key),P.dist.(key));
end
v=U(:,8);neg=v<0.5;
X(neg,8)=draw_prior(v(neg)/0.5,P.epsPYNegative,'uniform');
X(~neg,8)=draw_prior((v(~neg)-0.5)/0.5,P.epsPYPositive,'uniform');
end

function x=draw_prior(u,b,distribution)
if b(1)==b(2),x=repmat(b(1),size(u));return;end
switch lower(distribution)
    case 'uniform'
        x=b(1)+(b(2)-b(1))*u;
    case 'loguniform'
        x=b(1)*(b(2)/b(1)).^u;
    otherwise
        error('Unsupported distribution %s',distribution);
end
end

function write_snapshot(filename,P)
fid=fopen(filename,'w');assert(fid>=0);
autoClose=onCleanup(@()fclose(fid));
fprintf(fid,'index,parameter,branch,distribution,minimum,maximum,units\n');
for idx=1:10
    key=P.names{idx};
    if idx==8
        for s=1:2
            if s==1,b=P.epsPYNegative;side='negative';
            else,b=P.epsPYPositive;side='positive';end
            fprintf(fid,'%d,epsPY,%s,uniform,%.15g,%.15g,per_mil\n', ...
                idx,side,b(1),b(2));
        end
    else
        b=P.(key);dist=P.dist.(key);unit='dimensionless';
        if any(idx==[4 5 6 7 9]),unit='per_mil';end
        if idx==3,unit='fraction';end
        fprintf(fid,'%d,%s,all,%s,%.15g,%.15g,%s\n', ...
            idx,key,dist,b(1),b(2),unit);
    end
end
end
