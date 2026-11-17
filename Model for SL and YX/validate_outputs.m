function validate_outputs(outdir)
%VALIDATE_OUTPUTS Numerical and artifact checks (not geological validation).
if nargin<1 || isempty(outdir)
    outdir=fullfile(fileparts(mfilename('fullpath')),'output');
end
Z=load(fullfile(outdir,'results.mat'),'config','indices','samples');
cfg=Z.config;P=model_parameters();K=numel(P.names);
assert(K==10,'This release expects 10 independent inputs.');
assert(isfield(cfg,'priors') && isequaln(cfg.priors,P), ...
    'Current priors differ from saved ensemble: rerun full or quick.');
assert(isfield(cfg,'modelVersion') && ...
    strcmp(cfg.modelVersion,'noH_normalized_KD_KAH_Psi_KPH_div_R_v1'), ...
    'Old model detected; do not reuse outputs from H-dependent kinetics.');
assert(isfile(fullfile(outdir,'parameters_used.csv')));
D=readmatrix(fullfile(outdir,'conservation.csv'));
assert(size(D,1)==cfg.primaryRuns && size(D,2)==10, ...
    'Unexpected number of simulations or conservation columns.');
assert(cfg.primaryRuns==2*sum(cfg.NperRep)*(K+2), 'Incorrect 10D Sobol design size');
assert(all(abs(D(:,10)-.95)<1e-8), 'Pyrite P=0.95 event not reached');
assert(max(abs(D(:,5:7)),[],'all')<2e-9, 'Fe and isotope conservation failed');
assert(max(D(:,8))<2e-7, 'Matrix exponential / ODE disagreement');
assert(isequal(cfg.targets,[.50 .80 .95]),'Unexpected Ptarget windows');
assert(isequal(size(Z.samples),[2*sum(cfg.NperRep),29]),'Wrong sample table width');
assert(size(Z.indices,1)==2*numel(cfg.NperRep)*3*4*K, ...
    'Wrong Sobol index table length');
for e=1:2
    if e==1
        pb=P.Psi_SL; bb=P.beta_SL; kk=P.kappaOX_SL;
    else
        pb=P.Psi_YX; bb=P.beta_YX; kk=P.kappaOX_YX;
    end
    bounds=[P.f;pb;bb;P.deltaO;P.deltaL;P.epsDIR;P.epsAIR; ...
            [NaN NaN];P.epsOX;kk];
    for rep=1:numel(cfg.NperRep)
        ens=load(fullfile(outdir,sprintf('ensemble_%d_%d.mat',e,rep)), ...
            'A','B','Y','instantA','cumulativeA','U');
        N=cfg.NperRep(rep);
        assert(isequal(size(ens.A),[N 10]) && ...
               isequal(size(ens.B),[N 10]) && ...
               isequal(size(ens.U),[N 20]),'Sobol input dimension must be 10+10');
        assert(isequal(size(ens.Y),[N 3 4 12]),'Expected A,B and 10 AB_i outputs');
        assert(isequal(size(ens.instantA),[N 3 cfg.nq]) && ...
               isequal(size(ens.cumulativeA),[N 3 cfg.nq]), ...
               'Missing instantaneous or cumulative growth trajectories');
        assert(all(isfinite(ens.Y),'all') && ...
               all(isfinite(ens.instantA),'all') && ...
               all(isfinite(ens.cumulativeA),'all'),'Nonfinite ensemble values');
        for mat={ens.A,ens.B}
            X=mat{1};
            for j=[1:7 9:10]
                lb=bounds(j,1);ub=bounds(j,2);sc=max([1 abs(lb) abs(ub)]);
                assert(all(X(:,j)>=lb-1e-12*sc & X(:,j)<=ub+1e-12*sc), ...
                    'Parameter %d out of %s prior bounds',j,ternary(e==1,'SL','YX'));
            end
            ep=X(:,8);neg=ep>=P.epsPYNegative(1)&ep<=P.epsPYNegative(2);
            pos=ep>=P.epsPYPositive(1)&ep<=P.epsPYPositive(2);
            assert(all(neg|pos),'Invalid epsPY mixture draw');
            assert(sum(neg)==N/2, 'Scrambled Sobol epsPY sign split is not 50:50');
        end
    end
end
% SL and YX use the same scrambled coordinates but DIFFERENT Psi priors.
for rep=1:numel(cfg.NperRep)
    a=load(fullfile(outdir,sprintf('ensemble_1_%d.mat',rep)),'A');
    b=load(fullfile(outdir,sprintf('ensemble_2_%d.mat',rep)),'A');
    assert(isequal(a.A(:,[1 4:9]),b.A(:,[1 4:9])), ...
        'Identical common-prior Sobol coordinates should give matching values');
    assert(all(a.A(:,2)>=P.Psi_SL(1) & a.A(:,2)<=P.Psi_SL(2)) && ...
           all(b.A(:,2)>=P.Psi_YX(1) & b.A(:,2)<=P.Psi_YX(2)), ...
        'Environment-specific Psi bounds not applied');
end
Q=readmatrix(fullfile(outdir,'precision.csv'));
assert(size(Q,1)==2*cfg.nPrecisionPerEnv,'Wrong precision sample count');
assert(max(Q(:,3))<2e-4,'Tighter ODE tolerance changed results too much');
ci=readmatrix(fullfile(outdir,'sobol_ci_sorted.csv'));
ciY2=readmatrix(fullfile(outdir,'sobol_ci_sorted_Y2.csv'));
assert(isequal(size(ci),[2*3*K,11]) && isequal(size(ciY2),[2*3*K,11]), ...
    'Incomplete Sobol confidence interval tables');
assert(all(ci(:,3)>=1 & ci(:,3)<=K) && ...
       all(ciY2(:,3)>=1 & ciY2(:,3)<=K),'Invalid Sobol ranks');
heat=readmatrix(fullfile(outdir,'scenario_sensitivity_heatmap.csv'));
assert(size(heat,1)==2*3*4*K,'Wrong heatmap dimensions');
dir=readmatrix(fullfile(outdir,'directional_responses.csv'));
assert(isequal(size(dir),[288,9]),'Expected fixed-endmember direction table for Fig. 5');
% epsPY panels: sg=0 only; 2 env x 8 bins x 2 directions.
rr=dir(dir(:,2)==8,:);
assert(all(rr(:,3)==0) && size(rr,1)==32,'epsPY direction rows malformed');
for e=1:2
    for idx=[1 2 3 10]
        rr=dir(dir(:,1)==e & dir(:,2)==idx,:);
        assert(size(rr,1)==32,'Missing conditioned direction bins');
        assert(all(ismember(unique(rr(:,3))',[1 2])),'Expected neg/pos epsPY subsets');
        if idx==1, bnd=P.f;
        elseif idx==2, if e==1,bnd=P.Psi_SL;else,bnd=P.Psi_YX;end
        elseif idx==3, if e==1,bnd=P.beta_SL;else,bnd=P.beta_YX;end
        else, if e==1,bnd=P.kappaOX_SL;else,bnd=P.kappaOX_YX;end
        end
        v=rr(isfinite(rr(:,8)),8);
        assert(all(v>=bnd(1)-1e-12 & v<=bnd(2)+1e-12), ...
            'Direction bins fall outside active prior');
    end
end
assert(isequal(size(readmatrix(fullfile(outdir,'representative_cases.csv'))),[2 16]), ...
    'Incorrect representative-case file format');
for e=1:2
    profile=readmatrix(fullfile(outdir,sprintf('representative_profile_%d.csv',e)));
    curve=readmatrix(fullfile(outdir,sprintf('representative_curve_%d.csv',e)));
    cases=readmatrix(fullfile(outdir,'representative_cases.csv'));
    par=cases(cases(:,1)==e,4:13);
    assert(size(profile,2)==8 && isequal(size(curve),[cfg.nq 3]));
    assert(max(abs(profile(:,5)+profile(:,6)-profile(:,2)))<1e-10, ...
        'DIR+AIR must equal O (no H in normalized equations)');
    assert(max(abs(profile(:,7)-par(2)*profile(:,3)))<1e-9, ...
        'PY flux must equal Psi*L (no H)');
    assert(max(abs(profile(:,8)-par(10)*profile(:,3)))<1e-9, ...
        'OX flux must equal kappaOX*L');
end
bands=readmatrix(fullfile(outdir,'ensemble_pyrite_bands.csv'));
allBands=readmatrix(fullfile(outdir,'ensemble_pyrite_bands_targets.csv'));
assert(isequal(size(bands),[4*cfg.nq,10]) && ...
       isequal(size(allBands),[2*2*3*cfg.nq,11]),'Incomplete prior predictive envelopes');
assert(all(isfinite(allBands),'all') && all(isfinite(bands),'all'), ...
    'Nonfinite predicted isotope envelopes');
assert(all(allBands(:,6)<=allBands(:,5) & allBands(:,5)<=allBands(:,7)) && ...
       all(allBands(:,9)<=allBands(:,8) & allBands(:,8)<=allBands(:,10)), ...
       'Median falls outside its isotope prediction interval');
for e=1:2
    for t=1:3
        count=0;
        for sg=1:2
            rr=allBands(allBands(:,1)==e & allBands(:,2)==sg & ...
                abs(allBands(:,3)-cfg.targets(t))<1e-12,:);
            assert(size(rr,1)==cfg.nq);
            count=count+rr(1,11);
            if t==2
                old=bands(bands(:,1)==e&bands(:,2)==sg,:);
                assert(max(abs(old-rr(:,[1 2 4:11])),[],'all')<1e-9);
            end
        end
        assert(count==sum(cfg.NperRep),'All A draws must contribute to each window');
    end
end
% Validate saved Y endpoint extraction and numerical stability in each environment.
grid=[];
for e=1:2
    ens=load(fullfile(outdir,sprintf('ensemble_%d_1.mat',e)),'A','Y');
    for j=1:min(4,size(ens.A,1))
        base=fe_forward(ens.A(j,:),1e-8,181);
        tight=fe_forward(ens.A(j,:),1e-10,361);
        ref=reshape(ens.Y(j,:,:,1),3,4);
        grid(end+1,:)=[e j max(abs(base-ref),[],'all') ...
                       max(abs(base(:,1:2)-tight(:,1:2)),[],'all')]; %#ok<AGROW>
    end
end
writematrix(grid,fullfile(outdir,'grid_validation.csv'));
assert(max(grid(:,3))<2e-5 && max(grid(:,4))<2e-3, ...
    'Growth-grid resolution sensitivity failed');
% Jansen estimator orientation on a known additive 10D function.
rng(891,'twister');AA=rand(2000,K);BB=rand(2000,K);
YA=sum(AA,2);YB=sum(BB,2);V=var([YA;YB],1);
test=zeros(K,4);
for k=1:K
    AB=AA;AB(:,k)=BB(:,k);zc=sum(AB,2);
    test(k,:)=[k 1-mean((YB-zc).^2)/(2*V) ...
                mean((YA-zc).^2)/(2*V) 1/K];
end
writematrix(test,fullfile(outdir,'estimator_orientation.csv'));
assert(max(abs(test(:,2:3)-1/K),[],'all')<.09,'Jansen orientation check failed');
figdir=fullfile(outdir,'figures');
expected={'Fig01_representative_trajectories','Fig02_sobol_sorted_95CI', ...
    'Fig02b_Y2_sobol_sorted_95CI','Fig03_sensitivity_heatmap', ...
    'Fig04_sobol_convergence', ...
    'Fig05_direction_probabilities_4x3_fixed_epsPY', ...
    'FigS01_Ptarget_isotope_envelopes'};
for j=1:numel(expected)
    for ext={'.png','.pdf','.svg'}
        assert(isfile(fullfile(figdir,[expected{j} ext{1}])), ...
            'Missing figure %s%s',expected{j},ext{1});
    end
end
fid=fopen(fullfile(outdir,'VALIDATION.txt'),'w');
assert(fid>=0,'Cannot write validation report');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'PASS: %s mode, %d independent inputs, H removed\n',cfg.mode,K);
fprintf(fid,'Primary model solves: %d\n',cfg.primaryRuns);
fprintf(fid,'PASS: normalized fluxes JD=fO, JA=(1-f)O, JP=PsiL, JX=kappaL\n');
fprintf(fid,'PASS: Psi SL [%.6g, %.6g], YX [%.6g, %.6g]\n', ...
    P.Psi_SL(1),P.Psi_SL(2),P.Psi_YX(1),P.Psi_YX(2));
fprintf(fid,'PASS: beta SL [%.6g, %.6g], YX [%.6g, %.6g]\n', ...
    P.beta_SL(1),P.beta_SL(2),P.beta_YX(1),P.beta_YX(2));
fprintf(fid,'PASS: kappaOX SL [%.6g, %.6g], YX [%.6g, %.6g]\n', ...
    P.kappaOX_SL(1),P.kappaOX_SL(2),P.kappaOX_YX(1),P.kappaOX_YX(2));
fprintf(fid,'PASS: P=0.95, total Fe/54Fe/56Fe conserved, matrix solution\n');
fprintf(fid,'PASS: three Ptarget windows; all-A 95%% prior prediction bands\n');
fprintf(fid,'PASS: 10D Sobol/Jansen, Y1/Y2, probability and convergence outputs\n');
fprintf(fid,'PASS: six main figures and Supplementary S1 exported in 3 formats\n');
fprintf(fid,'MATLAB version: %s\n',version);
fprintf(fid,'CAUTION: Numerical consistency does not validate a natural mechanism.\n');
if strcmp(cfg.mode,'quick')
    fprintf(fid,'NOT FOR PUBLICATION: quick mode is a smoke test.\n');
end
fprintf('VALIDATION PASSED (%s).\n',cfg.mode);
end

function out=ternary(test,a,b)
if test,out=a;else,out=b;end
end
