function P = model_parameters()
%MODEL_PARAMETERS Edit ALL 10 input ranges and distributions here.
% A SINGLE generic iron-cycling scenario is simulated. NO SL/YX labels.
% Rate normalization absorbs H: JD=f*O, JA=(1-f)*O, JP=Psi*L,
% JX=kappaOX*L. Independent input order is fixed by P.names below.
%
% The default Psi/beta/kappaOX priors cover the UNION of the last SL/YX
% bounds. This is a NEW independent shared prior, NOT a mixture of SL/YX.
% Change any bounds or .dist entry here; rerun run_all('full').
% =================== USER-EDITABLE PARAMETER BLOCK ===================
P.f             = [0.00 1.00];
P.Psi           = [0.10 10.00];    % Combined old SL/YX range, loguniform
P.beta          = [1e-5 0.40];    % Combined old SL/YX range, loguniform
P.deltaO        = [-0.30 0.50];
P.deltaL        = [-1.00 1.00];
P.epsDIR        = [-3.00 -2.00];
P.epsAIR        = [-2.50 -0.50];
P.epsPYNegative = [-3.00 -0.50];
P.epsPYPositive = [0.50 3.00];
P.epsPYPositiveWeight = 0.50;
P.epsOX        = [0.40 3.00];
P.kappaOX      = [0.01 1.00];

% Distribution labels: 'uniform' or 'loguniform'.
P.dist.f        = 'uniform';
P.dist.Psi      = 'loguniform';
P.dist.beta     = 'loguniform';
P.dist.deltaO   = 'uniform';
P.dist.deltaL   = 'uniform';
P.dist.epsDIR   = 'uniform';
P.dist.epsAIR   = 'uniform';
P.dist.epsOX    = 'uniform';
P.dist.kappaOX  = 'loguniform';
% epsPY is a 50:50 mixture of uniform negative and positive intervals.
% ================== END USER-EDITABLE PARAMETER BLOCK ================

P.names={'f','Psi','beta','deltaO','deltaL','epsDIR', ...
         'epsAIR','epsPY','epsOX','kappaOX'};
required={'f','Psi','beta','deltaO','deltaL','epsDIR','epsAIR', ...
          'epsPYNegative','epsPYPositive','epsOX','kappaOX'};
for j=1:numel(required)
    key=required{j};b=P.(key);
    assert(isnumeric(b)&&isequal(size(b),[1,2])&&all(isfinite(b))&&b(1)<=b(2), ...
        'Invalid lower/upper limits for %s',key);
end
assert(P.f(1)>=0 && P.f(2)<=1 && P.f(1)<P.f(2),'f must vary in [0,1]');
assert(P.Psi(1)>0 && P.Psi(2)>P.Psi(1),'Psi must vary over positive values');
assert(P.beta(1)>0 && P.beta(2)<1 && P.beta(2)>P.beta(1), ...
    'beta must vary strictly inside (0,1)');
assert(P.kappaOX(1)>=0 && P.kappaOX(1)<P.kappaOX(2), ...
    'kappaOX must be nonnegative with a nonzero range');
assert(P.epsPYNegative(1)<P.epsPYNegative(2)&&P.epsPYNegative(2)<0 ...
    &&P.epsPYPositive(1)>0&&P.epsPYPositive(1)<P.epsPYPositive(2), ...
    'epsPY must have disjoint negative and positive branches');
assert(P.epsPYPositiveWeight==0.5, ...
    'Current scrambled Sobol branch design requires a 50:50 epsPY mixture');
keys=fieldnames(P.dist);
for j=1:numel(keys)
    key=keys{j};distribution=P.dist.(key);b=P.(key);
    assert(any(strcmpi(distribution,{'uniform','loguniform'})), ...
        'Unknown distribution for %s',key);
    if strcmpi(distribution,'loguniform')
        assert(b(1)>0 && b(2)>0,'loguniform bounds must be positive: %s',key);
    end
end
end
