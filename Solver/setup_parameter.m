function [opts, note_str] = setup_parameter(rec_mode_id, weight_mode_id, post_mode_id, quant_mode_id, isShowPSNR)

rec_mode        = {'TV', 'TVNL1', 'TVNL2', 'TVNL3', 'DTV', 'GSR'} ;
weight_mode     = {'No', 'ENOW', 'HENOW', 'EdgeCS', 'Candes', 'PC', 'PC_ENOW', 'PC_HENOW'};
post_mode       = {'No', 'BM3D', 'MH', 'NLM'};
quant_mode      = {'No', 'SQ'};
filter_mode     = {'No', 'NLM', 'BM3D', 'WNNM', 'MH' , 'Wienner', 'Gauss'};

note_str        = [rec_mode{rec_mode_id} '_w' weight_mode{weight_mode_id} '_p' post_mode{post_mode_id} '_q' quant_mode{quant_mode_id}];

%% For all TV
opts.lambda     = [0.5];        opts.nuy       = 0.05;
opts.mu         = 1;            opts.tol         = 0.2e-5;
opts.nInner     = [4];          opts.nOuter     = [30];
opts.nbrLoop    = 5;            opts.isShowPSNR = isShowPSNR;
opts.startNLR   = 5;
%% For reconstruction mode
switch rec_mode{rec_mode_id}
    case 'GSR'
        opts.rec_mode   = 6;            opts.tol        = 0.005;
        opts.nOuter     = 13;           opts.dtv_mode   = 0;
        
        %% For cartoon
        opts.car.lambda = 2;            opts.car.nuy    = 0.1;
        opts.car.mu     = 10;           opts.car.nInner = 4; % do we use this nInner parameter??? no I guess
        opts.car.nLoop_init = 20;       opts.car.nLoop  = 8;
        
        % for recovery
        opts.car.rec_mode = 2;
        
        % For regularization
        opts.car.gamma  = 100;      opts.car.h_reg    = 0.04;
        
        % For weighting scheme
        opts.car.weight_mode = weight_mode_id;
        opts.car.h_weight    = 0.04;
        opts.car.tau         = 5;
        
        % For filter mode
        opts.car.filter_mode = 3;
        switch opts.car.filter_mode
            case 1 % no
                
            case 2 % NLM
                opts.car.hFilter  = 0.04;
                opts.car.sigma    = opts.car.hFilter;
                
            case 3 % BM3D
                opts.car.sigma    = 15;         opts.tex.profile    = 'np';
                
            case 4 % WNNM
                opts.car.sigma    = 15;
                
            case 5 % MH
                opts.car.ratio    = 4;          opts.car.test_lambda= 1;
                opts.car.w        = 3;
                
            case 6 % Wienner
                
            case 7 % Gauss
                opts.car.sigma = 1;
        end;
        
        %% For texture
        opts.tex.lambda = 2;            opts.tex.nuy     = 0.1;
        opts.tex.mu     = 1;            opts.tex.nInner  = [8];
        opts.tex.nbrLoop_init = 4;      opts.tex.nLoop   = 4;
        
        % For recovery
        opts.tex.rec_mode = 1;
        
        % For regularization
        opts.tex.gamma  = 0;            opts.tex.h_reg    = 0;
        
        % For weighting scheme
        opts.tex.weight_mode = 1;       opts.tex.h_weight = 0;
        opts.tex.tau        = 0;
        
        % For filter
        opts.tex.filter_mode = 3;
        switch opts.tex.filter_mode
            case 1 % no
                
            case 2 % NLM
                opts.tex.hFilter  = 0.04;
                opts.tex.sigma    = opts.tex.hFilter;
                
            case 3 % BM3D
                opts.tex.sigma    = 5;         opts.tex.profile    = 'np';
                
            case 4 % WNNM
                opts.tex.sigma    = 5;
                
            case 5 % MH
                opts.tex.ratio    = 4;          opts.tex.test_lambda= 1;
                opts.tex.w        = 3;
                
            case 6 % Wienner
                
                
            case 7 % Gauss
                opts.tex.sigma  = 1;
        end;
        
        note_str =  [rec_mode{rec_mode_id} '_CarRec' rec_mode{opts.car.rec_mode} '_w' weight_mode{opts.car.weight_mode} ...
            '_filter' filter_mode{opts.car.filter_mode} '_texFilter' filter_mode{opts.tex.filter_mode} '_nOuter' num2str(opts.nOuter)];
        
        
        % for GSR
        par.PatchSize   = 8;
        par.SlidingDis  = 4;
        par.Factor      = 240;
        par.ArrayNo     = 60;
        par.SearchWin   = 20;
        %par.org         = image1;
        par.IterNum     = 50;
        par.mu          = 2.5e-3;
        par.lambda      = 0.082;
        par.tol         = 1e-5; 
        opts.par        = par; 
        
    case 'TV'
        opts.lambda     = [1];        opts.nuy       = 0.05;
        opts.mu         = 1;            opts.tol         = 0.2e-5;
        opts.nInner     = [4];          opts.nOuter     = [30];
        opts.nbrLoop    = 1;            opts.isShowPSNR = isShowPSNR;
        opts.rec_mode   = 1;
    case 'TVNL1'
        opts.lambda     = [2];          opts.nuy       = 0.1;
        opts.mu         = 10;           opts.tol         = 0.1e-5;
        opts.nInner     = [4];          opts.nOuter     = [30];
        opts.nbrLoop    = 1;            opts.isShowPSNR = isShowPSNR;
        
        opts.rec_mode   = 2;
        opts.gamma      = 100;          opts.sigma       = 0.03;
        %opts.filter_mode= 2;            opts.sigma       = 0.03;
        %         opts.filter_mode= 3;            opts.hReg       = 10;
    case 'TVNL2'
        opts.lambda     = [2];          opts.nuy       = 0.05;
        opts.mu         = 10;           opts.tol         = 0.2e-5;
        opts.nInner     = [8];          opts.nOuter     = [30];
        opts.nbrLoop    = 1;            opts.isShowPSNR = isShowPSNR;
        opts.rec_mode   = 3;            opts.gamma      = 10;
        
        %opts.filter_mode= 2;            opts.sigma       = 0.03;
        %         opts.filter_mode= 3;            opts.hReg       = 10;
        
    case 'TVNL3'
        opts.lambda     = [2];          opts.nuy       = 0.1;
        opts.mu         = 10;           opts.tol         = 0.1e-5;
        opts.nInner     = [8];          opts.nOuter     = [25];
        opts.nbrLoop    = 1;            opts.isShowPSNR = isShowPSNR;
        opts.rec_mode   = 4;            opts.gamma      = 100;
        opts.filter_mode= 2;            opts.sigma       = 0.03;
        
    case 'DTV'
        
        opts.rec_mode   = 5;            opts.tol        = 0.005;
        opts.nOuter     = 15;           opts.dtv_mode   = 0;
        
        %% For cartoon
        opts.car.lambda = 2;            opts.car.nuy    = 0.1;
        opts.car.mu     = 10;           opts.car.nInner = 4; % do we use this nInner parameter??? no I guess
        opts.car.nLoop_init = 20;       opts.car.nLoop  = 8;
        
        % for recovery
        opts.car.rec_mode = 1;
        
        % For regularization
        opts.car.gamma  = 100;      opts.car.h_reg    = 0.04;
        
        % For weighting scheme
        opts.car.weight_mode = weight_mode_id;
        opts.car.h_weight    = 0.04;
        opts.car.tau         = 5;
        
        % For filter mode
        opts.car.filter_mode = 3;
        switch opts.car.filter_mode
            case 1 % no
                
            case 2 % NLM
                opts.car.hFilter  = 0.04;
                opts.car.sigma    = opts.car.hFilter;
                
            case 3 % BM3D
                opts.car.sigma    = 15;         opts.tex.profile    = 'np';
                
            case 4 % WNNM
                opts.car.sigma    = 15;
                
            case 5 % MH
                opts.car.ratio    = 4;          opts.car.test_lambda= 1;
                opts.car.w        = 3;
                
            case 6 % Wienner
                
            case 7 % Gauss
                opts.car.sigma = 1;
        end;
        
        %% For texture
        opts.tex.lambda = 2;            opts.tex.nuy     = 0.1;
        opts.tex.mu     = 1;            opts.tex.nInner  = [8];
        opts.tex.nbrLoop_init = 4;      opts.tex.nLoop   = 4;
        
        % For recovery
        opts.tex.rec_mode = 1;
        
        % For regularization
        opts.tex.gamma  = 0;            opts.tex.h_reg    = 0;
        
        % For weighting scheme
        opts.tex.weight_mode = 1;       opts.tex.h_weight = 0;
        opts.tex.tau        = 0;
        
        % For filter
        opts.tex.filter_mode = 3;
        switch opts.tex.filter_mode
            case 1 % no
                
            case 2 % NLM
                opts.tex.hFilter  = 0.04;
                opts.tex.sigma    = opts.tex.hFilter;
                
            case 3 % BM3D
                opts.tex.sigma    = 5;         opts.tex.profile    = 'np';
                
            case 4 % WNNM
                opts.tex.sigma    = 5;
                
            case 5 % MH
                opts.tex.ratio    = 4;          opts.tex.test_lambda= 1;
                opts.tex.w        = 3;
                
            case 6 % Wienner
                
                
            case 7 % Gauss
                opts.tex.sigma  = 1;
        end;
        
        note_str =  [rec_mode{rec_mode_id} '_CarRec' rec_mode{opts.car.rec_mode} '_w' weight_mode{opts.car.weight_mode} ...
            '_filter' filter_mode{opts.car.filter_mode} '_texFilter' filter_mode{opts.tex.filter_mode} '_nOuter' num2str(opts.nOuter)];
        
end;

%% For weighting
opts.iterWeight = 15;
switch weight_mode{weight_mode_id}
    case 'No'
        opts.weight_mode = 1;
        
    case 'ENOW'
        opts.weight_mode = 2;
        opts.tau1        = 5;
        opts.hWeight     = 0.03;
        
    case 'HENOW'
        opts.weight_mode = 3;
        opts.tau1        = 5;
        opts.hWeight     = 0.03;
        
    case 'EdgeCS'
        opts.weight_mode = 4;
    case 'Candes'
        opts.weight_mode = 5;
end;

%% For  Post processing
opts.nInnerPost = 4;            opts.nOuterPost = 10;
opts.nLoopPost  = 4;
switch post_mode{post_mode_id}
    case 'No'
        opts.post_mode      = 1;
    case 'BM3D'
        opts.post_mode      = 2;
        opts.sigma          = 10;           opts.profile    = 'np';
    case 'MH'
        opts.post_mode      = 3;
        opts.ratio      	= 4;             opts.test_lambda= 1;
    case 'NLM'
        opts.post_mode      = 4;
        opts.sigma          = 0.04;
        
end;
opts.quant_mode         =  quant_mode_id;

