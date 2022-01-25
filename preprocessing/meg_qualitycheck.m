function meg_qualitycheck(cfg)

try
    cd(cfg.dataset)
    
    ft_qualitycheck(cfg)
    
    % copy pdf to outputfolder
    pdflist = dir('*.pdf');
    [~, out] = fileparts(cfg.dataset);
    outfile = fullfile(cfg.qualityoutputdir, [out '.pdf']);
    copyfile(pdflist(1).name, outfile)
catch ME
    disp(getReport(ME))
    fid = fopen('~/qualitycheck_errorlog.txt', 'at');
    fprintf(fid,'%s\n%s\n%s\n\n\n', datestr(now), cfg.dataset, getReport(ME));
    fclose('all');
end
