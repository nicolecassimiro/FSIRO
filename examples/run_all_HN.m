function run_all_HN(basepath)
% RUN_ALL_HN  Batch-run all Head-and-Neck patients (triangular + trapezoidal)
% and record DVH metrics vs TROTS. Results under <basepath>/resultados/HN/.
%
%   run_all_HN                 % basepath = current folder
%   run_all_HN('/path/to/implementacao')

    if nargin < 1 || isempty(basepath), basepath = pwd; end
    cfg.site     = 'HN';
    cfg.datadir  = fullfile(basepath, 'Head-and-Neck');
    cfg.prefix   = 'Head-and-Neck';
    cfg.npat     = 15;
    cfg.protocol = @protocol_HN;
    cfg.Vlevels  = 50;                 % larynx V50 (edema)
    cfg.outdir   = fullfile(basepath, 'resultados', 'HN');
    run_all_site(cfg);
end