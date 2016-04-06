function Opts = GSR_ParSet(Phi, x_MH, block_size,image)
% theta = 2, beta = 7, lambda = 12;

[row, col] = size(image);
Opts = [];
Opts.Phi = Phi;
Opts.block_size = block_size;
Opts.row = row;
Opts.col = col;
if ~isfield(Opts,'initial')
    Opts.initial = double(x_MH);
end

Opts.org = image;

if ~isfield(Opts,'IterNum')
    Opts.IterNum = 50;
end

if ~isfield(Opts,'mu')
    Opts.mu = 2.5e-3;
end

if ~isfield(Opts,'lambda')
    Opts.lambda = 0.082;
end

if ~isfield(Opts,'Inloop')
    Opts.Inloop = 200;
end