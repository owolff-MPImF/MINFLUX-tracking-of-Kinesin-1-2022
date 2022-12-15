function output=calculateSCE(singleTrace,minPhot,maxPhot,photsmW)

xM=singleTrace.x(1:end-1);yM=singleTrace.y(1:end-1); %FPGA derived position estimate (x,y)
dt=singleTrace.dt;t=[0;cumsum(dt(1:end-2))]; %time in measurement
nx=singleTrace.nx(2:end,:);NX=sum(nx,2); %photons by exposure for x localizations and total number of photons per localization
ny=singleTrace.ny(2:end,:);NY=sum(ny,2); %photons by exposure for y localizations and total number of photons per localization
Lx=2.*singleTrace.Lx(2:end); %MINFLUX L for x localization
Ly=2.*singleTrace.Ly(2:end); %MINFLUX L for y localization

bx=(nx(:,1)+nx(:,3)-2.*nx(:,2));by=(ny(:,1)+ny(:,3)-2.*ny(:,2)); %curvature of photon counts
%determine width for SCE
Ntarget=300;Ncurr=mean([NX;NY]);
SWW=ceil(Ntarget/Ncurr);
if all(NY==0,'all') %1D measurements
    filterN=and( movmedian(sum(nx,2),photsmW)>minPhot , movmedian(sum(nx,2),photsmW)<maxPhot ); %filter by number of photons
    filterT=Lx>0; %filter by actual localizations (L>0)
    bx(~and( filterN , filterT ))=NaN;
    curvx=movmean(bx,SWW,'omitnan');curvy=curvx; %calculate cliding curvature
    filterL=abs(nx(:,3)-nx(:,1))<(4.*curvx); %filter by positive curvature and correction <+-L
    by=bx;
    Ly=Lx;
    ny=nx;
    NY=NX;
else %2D measurements
    filterN=and( and( movmedian(sum(nx,2),photsmW)>minPhot , movmedian(sum(nx,2),photsmW)<maxPhot ) , and( movmedian(sum(ny,2),photsmW)>minPhot , movmedian(sum(ny,2),photsmW)<maxPhot ) ); %filter by number of photons
    filterT=and( Lx>0 , Ly>0 ); %filter by actual localizations (L>0)
    bx(~and( filterN , filterT ))=NaN;
    by(~and( filterN , filterT ))=NaN;
    curvx=movmean(bx,SWW,'omitnan');curvy=movmean(by,SWW,'omitnan'); %calculate cliding curvature
    filterL=and( abs(nx(:,3)-nx(:,1))<(4.*curvx) , abs(ny(:,3)-ny(:,1))<(4.*curvy) ); %filter by positive curvature and correction <+-L
end

vld=and( filterL , and( filterN , filterT ) ); %total valid localizations
% sliding curvature estimator
SCEx=Lx./4.*(nx(:,3)-nx(:,1))./curvx+xM;
SCEy=Ly./4.*(ny(:,3)-ny(:,1))./curvy+yM;

output.t=t(vld);
output.x=SCEx(vld);
output.y=SCEy(vld);
output.nx=nx(vld,:);
output.ny=ny(vld,:);

end