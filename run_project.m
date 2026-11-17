function run_project(outdir,mode)
%RUN_PROJECT Independent scrambled Sobol/Jansen ensembles.
% One authoritative forward solver: fe_forward.m.
% No parameter-pair interaction analysis or AB_ij solves are performed.
if nargin<1 || isempty(outdir),outdir=fullfile(fileparts(mfilename('fullpath')),'output');end
if nargin<2 || isempty(mode),mode='full';end
if ~exist(outdir,'dir'),mkdir(outdir);end
priors=model_parameters();
names=priors.names;
targets=[.50 .80 .95]; outputNames={'Y1','Y2','Dplus','Dminus'};
config.mode=char(mode);
% All prior bounds are read from model_parameters.m (never duplicated).
config.priors=priors;
config.priorVersion='ten_inputs_Psi_SL_1to10_YX_0p1to1_v1';
config.modelVersion='noH_normalized_KD_KAH_Psi_KPH_div_R_v1';
config.targets=targets;
if strcmp(mode,'quick')
    config.NperRep=[32 32 32 64];
    config.nPrecisionPerEnv=4;
else
    % Six independent N=1024 scrambled replicates for 95%% t intervals;
    % one independent N=2048 replicate for convergence cross-check.
    config.NperRep=[1024 1024 1024 1024 1024 1024 2048];
    config.nPrecisionPerEnv=32;
end
config.seedBase=700;
config.nCI=numel(config.NperRep)-1;
config.primaryRuns=2*sum(config.NperRep)*12;
config.tol=1e-8;config.nq=181;
assert(config.nCI>=3,'Require >=3 independent scrambles for CI plots.');
summary=zeros(2*numel(config.NperRep)*3*4,8);
indices=zeros(2*numel(config.NperRep)*3*4*10,8);
diagnostics=zeros(config.primaryRuns,10);
samples=zeros(2*sum(config.NperRep),29);
conv=[];si=0;ii=0;di=0;ai=0;
write_parameter_snapshot(fullfile(outdir,'parameters_used.csv'),priors);
logFile=fullfile(outdir,'execution.log');fid=fopen(logFile,'w');
assert(fid>=0,'Could not create execution.log');closeLog=onCleanup(@()fclose(fid));
assert(license('test','Statistics_Toolbox')==1, 'Statistics and Machine Learning Toolbox required.');
fprintf(fid,'MATLAB %s; mode %s; seeds 701..; nCI %d\n',version,config.mode,config.nCI);
for env=1:2
    for rep=1:numel(config.NperRep)
        N=config.NperRep(rep);
        rng(config.seedBase+rep,'twister');
        U=net(scramble(sobolset(20),'MatousekAffineOwen'),N);
        A=transformParameters(U(:,1:10),env,priors);
        B=transformParameters(U(:,11:20),env,priors);
        Y=zeros(N,3,4,12);
        instantA=zeros(N,3,config.nq);cumulativeA=zeros(N,3,config.nq);
        for block=1:12
            X=A;
            if block==2
                X=B;
            elseif block>=3
                X(:,block-2)=B(:,block-2);
            end
            for j=1:N
                if block==1
                    % Store every prior-draw isotope trajectory once, as the
                    % complete A-sample ensemble (not just a few examples).
                    [y,check,curves,flux,~,cumCurves]=fe_forward(X(j,:),config.tol,config.nq);
                    instantA(j,:,:)=reshape(curves,1,3,config.nq);
                    cumulativeA(j,:,:)=reshape(cumCurves,1,3,config.nq);
                else
                    [y,check,~,flux]=fe_forward(X(j,:),config.tol,config.nq);
                end
                Y(j,:,:,block)=reshape(y,1,3,4);
                di=di+1;diagnostics(di,:)=[env rep block j check];
                if block==1
                    ai=ai+1;samples(ai,:)=[env rep j X(j,:) reshape(y,1,[]) flux];
                end
            end
            fprintf(fid,'env=%d rep=%d block=%d completed\n',env,rep,block);
        end
        for t=1:3
            for out=1:4
                ya=Y(:,t,out,1);yb=Y(:,t,out,2);vv=var([ya;yb],1);
                si=si+1;summary(si,:)=[env rep N t out mean(ya) std(ya) vv];
                for idx=1:10
                    [s1,st]=jansen(ya,yb,Y(:,t,out,idx+2));
                    ii=ii+1;indices(ii,:)=[env rep N t out idx s1 st];
                end
            end
        end
        % Prefix-based convergence curves, no extra expensive forward runs.
        if strcmp(mode,'quick'),prefix=[8 16 32 64];else,prefix=[256 512 1024 2048];end
        for np=prefix(prefix<=N)
            for idx=1:10
                ya=Y(1:np,2,1,1);yb=Y(1:np,2,1,2);
                zab=Y(1:np,2,1,idx+2);
                [s1,st]=jansen(ya,yb,zab);
                conv(end+1,:)=[env rep np 2 1 idx s1 st]; %#ok<AGROW>
            end
        end
        save(fullfile(outdir,sprintf('ensemble_%d_%d.mat',env,rep)), ...
            'A','B','Y','U','instantA','cumulativeA','-v7.3');
        fprintf('env %d/%d; rep %d/%d; N=%d done\n',env,2,rep,numel(config.NperRep),N);
        fprintf(fid,'Saved env=%d rep=%d\n',env,rep);
    end
end
assert(di==config.primaryRuns && ai==size(samples,1) && ...
    si==size(summary,1) && ii==size(indices,1));
precision=zeros(2*config.nPrecisionPerEnv,15);kp=0;
for env=1:2
    rng(991,'twister');
    xx=transformParameters(net(scramble(sobolset(10),'MatousekAffineOwen'), ...
        config.nPrecisionPerEnv),env,priors);
    for j=1:config.nPrecisionPerEnv
        [u,c]=fe_forward(xx(j,:),1e-8,181);
        [v,cc]=fe_forward(xx(j,:),1e-10,181);
        kp=kp+1;precision(kp,:)=[env j max(abs(u(:)-v(:))) c cc];
    end
end
writematrix(summary,fullfile(outdir,'summary.csv'));
writematrix(indices,fullfile(outdir,'sobol.csv'));
writematrix(diagnostics,fullfile(outdir,'conservation.csv'));
writematrix(samples,fullfile(outdir,'samples.csv'));
writematrix(precision,fullfile(outdir,'precision.csv'));
writematrix(conv,fullfile(outdir,'sobol_convergence.csv'));
save(fullfile(outdir,'results.mat'),'summary','indices','diagnostics', ...
    'samples','precision','names','targets','outputNames','config','conv','-v7.3');
fprintf(fid,'Completed primary=%d; no parameter-interaction simulations\n',config.primaryRuns);
end

function [s1,st]=jansen(ya,yb,z)
vv=var([ya;yb],1);
if vv<=1e-12,s1=NaN;st=NaN;else
    st=mean((ya-z).^2)/(2*vv);
    s1=1-mean((yb-z).^2)/(2*vv);
end
end

function X=transformParameters(U,env,P)
% Transform a 10D scrambled Sobol draw using the single editable prior file.
% Independent input order: f Psi beta deltaO deltaL epsDIR epsAIR
%                          epsPY epsOX kappaOX. H IS NOT AN INPUT.
assert(size(U,2)==10 && any(env==[1 2]));
X=U;
X(:,1)=drawPrior(U(:,1),P.f,P.dist.f);
if env==1
    X(:,2)=drawPrior(U(:,2),P.Psi_SL,P.dist.Psi_SL);
    X(:,3)=drawPrior(U(:,3),P.beta_SL,P.dist.beta_SL);
    X(:,10)=drawPrior(U(:,10),P.kappaOX_SL,P.dist.kappaOX_SL);
else
    X(:,2)=drawPrior(U(:,2),P.Psi_YX,P.dist.Psi_YX);
    X(:,3)=drawPrior(U(:,3),P.beta_YX,P.dist.beta_YX);
    X(:,10)=drawPrior(U(:,10),P.kappaOX_YX,P.dist.kappaOX_YX);
end
X(:,4)=drawPrior(U(:,4),P.deltaO,P.dist.deltaO);
X(:,5)=drawPrior(U(:,5),P.deltaL,P.dist.deltaL);
X(:,6)=drawPrior(U(:,6),P.epsDIR,P.dist.epsDIR);
X(:,7)=drawPrior(U(:,7),P.epsAIR,P.dist.epsAIR);
X(:,9)=drawPrior(U(:,9),P.epsOX,P.dist.epsOX);
% Equal-probability, disjoint negative/positive pyrite fractionation.
neg=U(:,8)<0.5;
X(neg,8)=drawPrior(U(neg,8)/0.5,P.epsPYNegative,'uniform');
X(~neg,8)=drawPrior((U(~neg,8)-0.5)/0.5,P.epsPYPositive,'uniform');
end

function y=drawPrior(u,b,dist)
if b(1)==b(2),y=repmat(b(1),size(u));return;end
switch lower(dist)
    case 'uniform'
        y=b(1)+(b(2)-b(1))*u;
    case 'loguniform'
        assert(b(1)>0 && b(2)>0,'Positive prior bounds required for loguniform');
        y=b(1)*(b(2)/b(1)).^u;
    otherwise
        error('Unknown distribution: %s',dist);
end
end

function write_parameter_snapshot(filename,P)
% Output audit trail, not an editable configuration file.
fid=fopen(filename,'w');assert(fid>=0,'Cannot write parameter snapshot');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'environment,index,parameter,branch,distribution,minimum,maximum,units\n');
for env=1:2
    if env==1
        tag='SL';Psi=P.Psi_SL;beta=P.beta_SL;k=P.kappaOX_SL;
        dPsi=P.dist.Psi_SL;dBeta=P.dist.beta_SL;dK=P.dist.kappaOX_SL;
    else
        tag='YX';Psi=P.Psi_YX;beta=P.beta_YX;k=P.kappaOX_YX;
        dPsi=P.dist.Psi_YX;dBeta=P.dist.beta_YX;dK=P.dist.kappaOX_YX;
    end
    bounds={P.f,Psi,beta,P.deltaO,P.deltaL,P.epsDIR,P.epsAIR, ...
            P.epsPYNegative,P.epsOX,k};
    dists={P.dist.f,dPsi,dBeta,P.dist.deltaO,P.dist.deltaL, ...
           P.dist.epsDIR,P.dist.epsAIR,'uniform',P.dist.epsOX,dK};
    for idx=1:10
        if idx==8
            fprintf(fid,'%s,%d,epsPY,negative,uniform,%.15g,%.15g,per_mil\n', ...
                tag,idx,P.epsPYNegative(1),P.epsPYNegative(2));
            fprintf(fid,'%s,%d,epsPY,positive,uniform,%.15g,%.15g,per_mil\n', ...
                tag,idx,P.epsPYPositive(1),P.epsPYPositive(2));
        else
            units='dimensionless';
            if any(idx==[4 5 6 7 9]),units='per_mil';end
            if idx==3,units='fraction';end
            fprintf(fid,'%s,%d,%s,all,%s,%.15g,%.15g,%s\n', ...
                tag,idx,P.names{idx},dists{idx},bounds{idx}(1),bounds{idx}(2),units);
        end
    end
end
end
