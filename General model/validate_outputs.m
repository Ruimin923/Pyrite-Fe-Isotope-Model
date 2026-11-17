function validate_outputs(outdir)
%VALIDATE_OUTPUTS Independent structural and conservation checks.
% Passing validation does not demonstrate geological mechanism validity.
if nargin<1||isempty(outdir)
    outdir=fullfile(fileparts(mfilename('fullpath')),'output');
end
Z=load(fullfile(outdir,'results.mat'),'config','indices','samples','conv');
C=Z.config;P=model_parameters();K=numel(P.names);
assert(K==10 && isfield(C,'priors') && isequaln(C.priors,P), ...
    'Active 10D priors differ from saved ensemble');
assert(strcmp(C.modelVersion,'single_generic_noH_10D_v1'), ...
    'Wrong ensemble version; do not reuse old SL/YX simulations');
assert(isequal(C.targets,[.50 .80 .95]));
assert(C.primaryRuns==sum(C.NperRep)*(K+2));
assert(isfile(fullfile(outdir,'parameters_used.csv')));
assert(isequal(size(Z.samples),[sum(C.NperRep),28]));
assert(isequal(size(Z.indices),[numel(C.NperRep)*3*4*K,7]));
D=readmatrix(fullfile(outdir,'conservation.csv'));
assert(isequal(size(D),[C.primaryRuns,9]));
assert(max(abs(D(:,9)-.95))<1e-8,'P=0.95 event did not terminate');
assert(max(abs(D(:,4:6)),[],'all')<2e-9, ...
    'Total Fe / isotope inventories not conserved');
assert(max(D(:,7))<2e-7,'Total Fe differs from matrix-exponential solution');
limits=[P.f;P.Psi;P.beta;P.deltaO;P.deltaL;P.epsDIR;P.epsAIR; ...
    [NaN NaN];P.epsOX;P.kappaOX];
for rep=1:numel(C.NperRep)
    E=load(fullfile(outdir,sprintf('ensemble_%d.mat',rep)), ...
        'A','B','U','Y','instantA','cumulativeA');
    N=C.NperRep(rep);
    assert(isequal(size(E.A),[N,K])&&isequal(size(E.B),[N,K]) ...
        &&isequal(size(E.U),[N,2*K]));
    assert(isequal(size(E.Y),[N,3,4,K+2]));
    assert(isequal(size(E.instantA),[N,3,C.nq]) && ...
        isequal(size(E.cumulativeA),[N,3,C.nq]));
    assert(all(isfinite(E.Y),'all') && ...
        all(isfinite(E.instantA),'all') && all(isfinite(E.cumulativeA),'all'));
    for matrix={E.A,E.B}
        X=matrix{1};
        for j=[1:7 9:10]
            b=limits(j,:);sc=max([1 abs(b)]);
            assert(all(X(:,j)>=b(1)-1e-12*sc & ...
                X(:,j)<=b(2)+1e-12*sc),'Input %d out of prior bounds',j);
        end
        ev=X(:,8);negative=ev>=P.epsPYNegative(1) & ...
            ev<=P.epsPYNegative(2);
        positive=ev>=P.epsPYPositive(1) & ...
            ev<=P.epsPYPositive(2);
        assert(all(negative|positive),'Invalid epsPY sign-mixture sample');
        assert(sum(negative)==N/2,'Expected equal epsPY sign counts');
    end
end
Q=readmatrix(fullfile(outdir,'precision.csv'));
assert(isequal(size(Q),[C.nPrecision,14]));
assert(max(Q(:,2))<2e-4,'Tighter ODE tolerance changed output too much');
CI=readmatrix(fullfile(outdir,'sobol_ci_sorted.csv'));
CI2=readmatrix(fullfile(outdir,'sobol_ci_sorted_Y2.csv'));
assert(isequal(size(CI),[3*K,10]) && isequal(size(CI2),[3*K,10]));
assert(all(ismember(CI(:,1),0:2)) && ...
    all(CI(:,2)>=1 & CI(:,2)<=K));
H=readmatrix(fullfile(outdir,'scenario_sensitivity_heatmap.csv'));
assert(isequal(size(H),[3*4*K,5]));
DR=readmatrix(fullfile(outdir,'directional_responses.csv'));
assert(isequal(size(DR),[304,9]));
assert(all(DR(:,9)>=0));
neg=DR(DR(:,1)==8,:);
assert(size(neg,1)==16 && all(neg(:,2)==0), ...
    'Expected eight direct epsPY bins with two direction probabilities');
for idx=[1:7 9:10]
    dat=DR(DR(:,1)==idx,:);
    assert(size(dat,1)==32 && ...
        isequal(sort(unique(dat(:,2))),[1;2]), ...
        'Expected 9 x eight bins x two signs x two directions');
    bounds=P.(P.names{idx});
    assert(all(dat(:,7)>=bounds(1)-1e-12 & ...
        dat(:,7)<=bounds(2)+1e-12), ...
        'Bin centers must lie within active prior for param %d',idx);
end
for rep=1:numel(C.NperRep)
    E=load(fullfile(outdir,sprintf('ensemble_%d.mat',rep)),'A','Y');
    for j=1:min(2,size(E.A,1))
        [newY,~,~]=fe_forward(E.A(j,:),1e-8,C.nq);
        stored=reshape(E.Y(j,:,:,1),3,4);
        assert(max(abs(stored-newY),[],'all')<2e-5, ...
            'Stored Y not consistent with fe_forward');
    end
end
figdir=fullfile(outdir,'figures');
expected={'Fig02_sobol_sorted_95CI','Fig02b_Y2_sobol_sorted_95CI', ...
    'Fig03_sensitivity_heatmap','Fig04_sobol_convergence', ...
    'Fig05_direction_probabilities_epsPY_sign_conditioned'};
for j=1:numel(expected)
    for extension={'.png','.pdf','.svg'}
        path=fullfile(figdir,[expected{j} extension{1}]);
        assert(isfile(path),'Missing figure %s',path);
    end
end
fid=fopen(fullfile(outdir,'VALIDATION.txt'),'w');assert(fid>=0);
clean=onCleanup(@()fclose(fid));
fprintf(fid,'PASS: single generic scenario, 10 independent inputs, no SL/YX or H input\n');
fprintf(fid,'PASS: primary model evaluations = %d\n',C.primaryRuns);
fprintf(fid,'PASS: Psi [%.5g,%.5g], beta [%.5g,%.5g], kappaOX [%.5g,%.5g]\n', ...
    P.Psi(1),P.Psi(2),P.beta(1),P.beta(2),P.kappaOX(1),P.kappaOX(2));
fprintf(fid,'PASS: P=.95 event, iron / isotope conservation, matrix consistency\n');
fprintf(fid,'PASS: 10D Sobol, Y1/Y2 and binary direction metrics\n');
fprintf(fid,'PASS: epsPY sign-conditioned direction probabilities and ten-input S1/ST\n');
fprintf(fid,'PASS: five figures (2,2b,3,4,5) exported as PNG/PDF/SVG\n');
fprintf(fid,'MATLAB: %s\n',version);
fprintf(fid,'CAUTION: numerical tests do not validate geology.\n');
if strcmp(C.mode,'quick')
    fprintf(fid,'NOT FOR PUBLICATION: quick mode is only a smoke test.\n');
end
fprintf('VALIDATION PASSED (%s).\n',C.mode);
end
