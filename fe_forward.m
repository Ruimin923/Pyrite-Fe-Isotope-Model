function [Y,check,curves,cumulative,profile,cumulativeCurves] = fe_forward(x,tol,nq)
%FE_FORWARD Single, closed three-reservoir iron/isotope forward solver.
% INPUT ORDER (10):
% x = [f Psi beta deltaO deltaL epsDIR epsAIR epsPY epsOX kappaOX]
%
% Definitions (all normalized by KD+KA*H, with H absorbed):
% f=KD/(KD+KA*H), Psi=KP*H/(KD+KA*H),
% kappaOX=KX/(KD+KA*H); tau=(KD+KA*H)*time.
% JD=f*O; JA=(1-f)*O; JP=Psi*L; JX=kappaOX*L.
%
% Model mass is approximated by 54Fe+56Fe. Delta, epsilon units: per mil.
% Ptargets=[.50,.80,.95] sampled at q=P/Ptarget in [.05,.95].
% Y columns: Y1=endpoint delta change; Y2=total variation-|Y1|;
%            Dplus=(Y1>.05); Dminus=(Y1<-.05).
% cumulativeCurves: bulk P inventory isotope ratio, NOT time average.
if nargin<2 || isempty(tol),tol=1e-8;end
if nargin<3 || isempty(nq),nq=181;end
assert(isvector(x) && numel(x)==10,'Expected 10 independent inputs');
x=x(:)';
f=x(1);Psi=x(2);beta=x(3);k=x(10);
assert(f>=0 && f<=1 && Psi>0 && k>=0 && beta>0 && beta<1, ...
    'Invalid Fe rate or initial inventory parameter');
r=1;p=Psi;Rref=15.698;
T0=[1-beta;beta;0];
R0=Rref*(1+x(4:5)'/1000);
a0=T0(1:2)./(1+R0);
y0=[a0;0; T0(1:2)-a0;0;zeros(4,1)];
M=[-r k 0;r -p-k 0;0 p 0];
slowRate=(r+p+k-sqrt((r+p+k)^2-4*r*p))/2;
assert(slowRate>0,'Nonpositive total-Fe decay rate');
tmax=12/slowRate;
opt=odeset('RelTol',tol,'AbsTol',tol*1e-4,'Events',@stopAt95);
sol=ode15s(@rhs,[0 tmax],y0,opt);
while isempty(sol.xe)
    tmax=2*tmax;sol=ode15s(@rhs,[0 tmax],y0,opt);
end
q=linspace(.05,.95,nq);targets=[.50 .80 .95];
Y=zeros(3,4);curves=zeros(3,nq);
if nargout>=6,cumulativeCurves=zeros(3,nq);else,cumulativeCurves=[];end
l1=-slowRate;l2=-(r+p+k)+slowRate;
cs=(-p*beta-l2)/(l1-l2);
for t=1:3
    goal=targets(t)*q;
    lo=zeros(size(q));hi=sol.x(end)*ones(size(q));
    for iter=1:48
        mid=(lo+hi)/2;
        pNow=1-cs*exp(l1*mid)-(1-cs)*exp(l2*mid);
        below=pNow<goal;
        lo(below)=mid(below);hi(~below)=mid(~below);
    end
    ts=(lo+hi)/2;v=deval(sol,ts);
    delta=1000*((v(5,:)./v(2,:))*exp(x(8)/1000)/Rref-1);
    curves(t,:)=delta;
    if nargout>=6
        cumulativeCurves(t,:)=1000*((v(6,:)./v(3,:))/Rref-1);
    end
    y1=delta(end)-delta(1);
    tv=sum(abs(diff(delta)));
    Y(t,:)=[y1 max(0,tv-abs(y1)) double(y1>.05) double(y1<-.05)];
end
v=sol.y;T=v(1:3,:)+v(4:6,:);
exact=expm(M*sol.x(end))*T0;
check=[max(abs(sum(T,1)-1)), ...
       max(abs(sum(v(1:3,:),1)-sum(y0(1:3)))), ...
       max(abs(sum(v(4:6,:),1)-sum(y0(4:6)))), ...
       max(abs(T(:,end)-exact)),sol.x(end),T(3,end)];
cumulative=v(7:10,end)';
if nargout>=5
    times=linspace(0,sol.x(end),301);
    vv=deval(sol,times);TT=vv(1:3,:)+vv(4:6,:);
    profile=[times;TT;f*TT(1,:);(1-f)*TT(1,:);Psi*TT(2,:);k*TT(2,:)];
else
    profile=[];
end
    function dy=rhs(~,y)
        Tot=y(1:3)+y(4:6);
        J=[f*Tot(1);(1-f)*Tot(1);Psi*Tot(2);k*Tot(2)];
        src=[1 1 2 2];dst=[2 2 3 1];eps=x(6:9);
        dy=zeros(10,1);
        for ii=1:4
            s=src(ii);d=dst(ii);
            Rprod=y(s+3)/y(s)*exp(eps(ii)/1000);
            j54=J(ii)/(1+Rprod);j56=J(ii)-j54;
            dy(s)=dy(s)-j54;dy(d)=dy(d)+j54;
            dy(s+3)=dy(s+3)-j56;dy(d+3)=dy(d+3)+j56;
        end
        dy(7:10)=J;
    end
    function [value,isterminal,direction]=stopAt95(~,y)
        value=y(3)+y(6)-.95;
        isterminal=1;direction=1;
    end
end
