function run_all(mode)
%RUN_ALL Single generic 10-input Fe-isotope model and sensitivity figures.
% run_all('quick') : shortened end-to-end smoke test, output_quick/
% run_all('full')  : full scrambled Sobol ensemble, output/
% run_all('reuse') : refresh tables/figures using existing compatible output/
% No SL/YX split or pairwise-interaction modules.
if nargin<1 || isempty(mode),mode='quick';end
mode=lower(char(mode));
assert(any(strcmp(mode,{'quick','full','reuse'})), ...
    'Mode must be quick, full or reuse');
root=fileparts(mfilename('fullpath'));addpath(root);
if strcmp(mode,'quick'),outdir=fullfile(root,'output_quick');
else,outdir=fullfile(root,'output');end
if ~exist(outdir,'dir'),mkdir(outdir);end
if strcmp(mode,'reuse')
    required={fullfile(outdir,'results.mat'), ...
        fullfile(outdir,'precision.csv')};
    assert(all(cellfun(@isfile,required)), ...
        'No completed full run: execute run_all(''full'') first');
    d=load(fullfile(outdir,'results.mat'),'config');cfg=d.config;
    P=model_parameters();
    assert(isfield(cfg,'priors') && isequaln(cfg.priors,P) && ...
        isfield(cfg,'modelVersion') && ...
        strcmp(cfg.modelVersion,'single_generic_noH_10D_v1'), ...
        'Ensemble does not match the active single-scenario model/priors');
    for rep=1:numel(cfg.NperRep)
        f=fullfile(outdir,sprintf('ensemble_%d.mat',rep));
        assert(isfile(f),'Missing ensemble file %s',f);
        info=whos('-file',f);names={info.name};
        assert(all(ismember({'A','B','Y'},names)), ...
            'Incomplete ensemble; rerun the full model');
    end
    fprintf('Reusing existing single-scenario ensemble (%d model evaluations).\n', ...
        cfg.primaryRuns);
else
    fprintf('Computing SINGLE-SCENARIO %s ensemble.\n',mode);
    run_project(outdir,mode);
end
analyze_results(outdir);
make_figures(outdir);
validate_outputs(outdir);
fprintf('Finished. See %s\n',fullfile(outdir,'figures'));
end
