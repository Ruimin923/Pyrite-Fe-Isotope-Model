function P = model_parameters()
%MODEL_PARAMETERS Single editable source of parameter bounds and priors.
% The normalized model has 10 INDEPENDENT inputs; H is NOT sampled:
% f=KD/(KD+KA*H), Psi=KP*H/(KD+KA*H), kappaOX=KX/(KD+KA*H).
% Rates per dimensionless time: JD=f*O, JA=(1-f)*O,
% JP=Psi*L, JX=kappaOX*L.
%
% To change any prior, edit the bounds and its distribution below, then
% rerun run_all('quick') / run_all('full').  REUSE is rejected when priors
% differ.  'loguniform' requires positive lower/upper bounds.
%
% ================== EDIT PARAMETER PRIORS HERE ==================
% COMMON to both environments
P.f             = [ 0.00,  1.00 ];   % DIR fraction of total reduction
P.deltaO        = [-0.30,  0.50 ];   % initial Fe(III) delta56Fe, per mil
P.deltaL        = [-1.00,  1.00 ];   % initial Fe(II) delta56Fe, per mil
P.epsDIR        = [-3.00, -2.00 ];   % DIR fractionation, per mil
P.epsAIR        = [-2.50, -0.50 ];   % AIR fractionation, per mil
P.epsPYNegative = [-3.00, -0.50 ];   % pyrite: negative branch, per mil
P.epsPYPositive = [ 0.50,  3.00 ];   % pyrite: positive branch, per mil
P.epsOX         = [ 0.40,  3.00 ];   % oxidation fractionation, per mil
P.epsPYPositiveWeight = 0.50;       % positive branch probability

% SONGLIN (SL): Psi, initial Fe(II) fraction, oxidation coefficient
P.Psi_SL        = [ 1.00, 10.00 ];   % KP*H/(KD+KA*H)
P.beta_SL       = [ 1e-5,  0.02 ];   % initial L inventory fraction
P.kappaOX_SL    = [ 0.01,  1.00 ];   % KX/(KD+KA*H)

% YOUXI (YX)
P.Psi_YX        = [ 0.10,  1.00 ];
P.beta_YX       = [ 0.05,  0.40 ];
P.kappaOX_YX    = [ 0.10,  1.00 ];

% Distributions can also be modified here. Allowed: uniform, loguniform.
P.dist.f           = 'uniform';
P.dist.Psi_SL      = 'loguniform';
P.dist.Psi_YX      = 'loguniform';
P.dist.beta_SL     = 'loguniform';
P.dist.beta_YX     = 'uniform';
P.dist.deltaO      = 'uniform';
P.dist.deltaL      = 'uniform';
P.dist.epsDIR      = 'uniform';
P.dist.epsAIR      = 'uniform';
P.dist.epsOX       = 'uniform';
P.dist.kappaOX_SL  = 'loguniform';
P.dist.kappaOX_YX  = 'loguniform';
% epsPY branches are each uniform, with the mixture weight above.
% ================== END OF EDITABLE PRIORS ======================

% Stable solver/ensemble input order; do not edit without updating code.
P.names = {'f','Psi','beta','deltaO','deltaL','epsDIR', ...
           'epsAIR','epsPY','epsOX','kappaOX'};

fields={'f','Psi_SL','Psi_YX','beta_SL','beta_YX','deltaO','deltaL', ...
        'epsDIR','epsAIR','epsPYNegative','epsPYPositive','epsOX', ...
        'kappaOX_SL','kappaOX_YX'};
for j=1:numel(fields)
    b=P.(fields{j});
    assert(isnumeric(b) && isequal(size(b),[1 2]) && all(isfinite(b)) ...
        && b(1)<=b(2), 'Invalid bounds in %s: [lower, upper] required',fields{j});
end
assert(P.f(1)>=0 && P.f(2)<=1,'f must lie within [0,1].');
for field={'beta_SL','beta_YX'}
    b=P.(field{1});assert(b(1)>0 && b(2)<1, 'beta must lie inside (0,1)');
end
for field={'Psi_SL','Psi_YX'}
    b=P.(field{1});assert(b(1)>0 && b(2)>0,'Psi must be positive');
end
for field={'kappaOX_SL','kappaOX_YX'}
    b=P.(field{1});assert(b(1)>=0,'kappaOX cannot be negative');
end
assert(P.epsPYNegative(2)<0 && P.epsPYPositive(1)>0, ...
    'epsPY branches must have strictly negative and positive values');
assert(P.epsPYPositiveWeight>0 && P.epsPYPositiveWeight<1, ...
    'epsPY positive branch probability must be in (0,1)');
assert(P.epsPYPositiveWeight==0.5, ...
    'This scrambled Sobol implementation expects a 50:50 epsPY mixture');
for j=1:numel(fieldnames(P.dist))
    keys=fieldnames(P.dist);key=keys{j};val=P.dist.(key);b=P.(key);
    assert(any(strcmpi(val,{'uniform','loguniform'})), ...
        'Unknown distribution for %s',key);
    if strcmpi(val,'loguniform') && b(1)~=b(2)
        assert(b(1)>0,'loguniform prior %s needs positive lower bound',key);
    end
end
% Figure 5 bins need nonzero width for the shown variables.
assert(P.f(1)<P.f(2) && P.Psi_SL(1)<P.Psi_SL(2) && ...
    P.Psi_YX(1)<P.Psi_YX(2) && ...
    P.epsPYNegative(1)<P.epsPYNegative(2) && ...
    P.epsPYPositive(1)<P.epsPYPositive(2), ...
    'f, both Psi priors, and epsPY branches must have nonzero ranges');
end
