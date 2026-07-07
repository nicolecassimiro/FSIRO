function run_all_liver(basepath)
% RUN_ALL_LIVER  Batch-run all Liver patients (triangular + trapezoidal) and
% record DVH metrics vs TROTS. Results under <basepath>/resultados/Liver/.
%
%   run_all_liver              % basepath = current folder
%   run_all_liver('/path/to/implementacao')

    if nargin < 1 || isempty(basepath), basepath = pwd; end
    cfg.site     = 'Liver';
    cfg.datadir  = fullfile(basepath, 'Liver');
    cfg.prefix   = 'Liver';
    cfg.npat     = 10;
    cfg.protocol = @protocol_liver;
    cfg.Vlevels  = [15 20 23 25 28];   % liver 700cc<15; kidneys V20/V23/V28; heart V25
    cfg.outdir   = fullfile(basepath, 'resultados', 'Liver');
    run_all_site(cfg);
end