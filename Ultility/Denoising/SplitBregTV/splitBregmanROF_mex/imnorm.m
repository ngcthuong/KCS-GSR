function  val= imnorm(u,normtype)
if nargin<2
    normtype = 'L2';
end
normtype=lower(normtype);

u=u(:);
len=numel(u);
val=1;
switch normtype
    case 'l1'
        val=sum(abs(u))./len;
    case 'l2'
        val=sqrt(sum(abs(u).^2))/len;
    otherwise
        disp('unknown');
end

