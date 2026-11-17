function run_all(mode)
%RUN_ALL Reproducible SL/YX Fe-isotope model and publication figures.
% run_all('quick') : abbreviated end-to-end test, output_quick/ (NOT publication data)
% run_all('full')  : independent scrambled replicates and full analysis
% run_all('reuse') : regenerate figures/tables from a completed full run
% Legacy compatibility: run_all(true)='full'; run_all(false)='reuse'.
% The optional attribution supplemental FIGURE is removed from this version.
if nargin<1 || isempty(mode),mode='reuse';end
if islogical(mode) || (isnumeric(mode) && isscalar(mode))
    if mode,mode='full';else,mode='reuse';end
end
mode=lower(char(mode));
assert(any(strcmp(mode,{'quick','full','reuse'})), 'Mode must be quick, full or reuse.');
root=fileparts(mfilename('fullpath'));addpath(root);
if strcmp(mode,'quick'),outdir=fullfile(root,'output_quick');else,outdir=fullfile(root,'output');end
if ~exist(outdir,'dir'),mkdir(outdir);end
if ~strcmp(mode,'reuse')
    fprintf('Preparing %s ensemble in %s\n',mode,outdir);
    run_project(outdir,mode);
else
    required={fullfile(outdir,'results.mat'),fullfile(outdir,'precision.csv')};
    assert(all(cellfun(@isfile,required)), 'Missing full ensemble: run run_all(''full'') first.');
    d=load(fullfile(outdir,'results.mat'),'config');config=d.config;
    % Compare the EXACT prior struct saved in the ensemble with the current
    % editable model_parameters.m, not a hard-coded priorVersion string.
    activePriors=model_parameters();
    assert(isfield(config,'priors') && isequaln(config.priors,activePriors) && ...
        isfield(config,'modelVersion') && ...
        strcmp(config.modelVersion,'noH_normalized_KD_KAH_Psi_KPH_div_R_v1'), ...
        ['Saved results were made with different priors/kinetics. ', ...
         'After changing model_parameters.m, run run_all(''full'') again.']);
    for e=1:2
        for rep=1:numel(config.NperRep)
            file=fullfile(outdir,sprintf('ensemble_%d_%d.mat',e,rep));
            assert(isfile(file),'Missing ensemble file; rerun run_all(''full'').');
            info=whos('-file',file);fn={info.name};
            assert(all(ismember({'instantA','cumulativeA'},fn)), ...
                'Old ensemble lacks cumulative/new pyrite curves; rerun run_all(''full'').');
            sInfo=info(strcmp(fn,'instantA'));
            assert(numel(sInfo.size)==3 && sInfo.size(2)==3, ...
                'Old ensemble lacks three Ptarget windows; rerun run_all(''full'').');
        end
    end
    fprintf('Reusing %d primary runs from completed ensemble.\n',config.primaryRuns);
end
% Remove obsolete interaction tables and plot artifacts from earlier releases.
obsolete={'pair_interaction_replicates.csv','pair_interaction_ci.csv'};
for j=1:numel(obsolete)
    filename=fullfile(outdir,obsolete{j});
    if isfile(filename),delete(filename);end
end
analyze_results(outdir);
make_figures(outdir);
validate_outputs(outdir);
fprintf('Finished. Figures: %s\n',fullfile(outdir,'figures'));
end
