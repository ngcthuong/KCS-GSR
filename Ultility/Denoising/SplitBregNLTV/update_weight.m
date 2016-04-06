function wopts=update_weight(Im0,h0,nwin,nbloc,NbBestNeigh)
if (~exist('NbBestNeigh','var'))
NbBestNeigh = 10; % >=6
end    
[M,N]=size(Im0);
NbNearNeigh = 4; % 0 or 4
display_messages = 0;
a = 2;
m=2*nwin+1;
w=2*nbloc+1;
NbNeigh=NbNearNeigh+NbBestNeigh;
G = fspecial('gaussian',[m m],a);
 VecGeneralParameters = [ display_messages; M; N; m; w; NbNearNeigh; NbBestNeigh; h0*h0;];
[W,Y] = compute_nl_weights_mex(single(Im0),single(G),single(VecGeneralParameters));
wopts.NbNeigh=NbNeigh;
wopts.m=m;
wopts.w=w;
wopts.W=W;
wopts.Y=Y;
