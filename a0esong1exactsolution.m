%a0esong1exactsolution


a=1
bb=5
nu=0.3
E=1e5
P=1

coord4
sdConn4


x=coord4(:,1)
y=coord4(:,2)
r=sqrt(x.^2+y.^2)
urexact=a*r.*(1-nu+bb^2*(1+nu)./(r.^2))./(bb^2-a^2)
sigmarexact=(1-bb.^2./(r.^2))./((bb/a)^2-1)
sigmathetaexact=(1+bb^2./r.^2)./((bb/a)^2-1)
cell2mat(sdConn4)
elementconn=cell(length(sdConn4),1);
for i=1:length(sdConn4)
    elementconn{i}=sdConn4{i}(:,1)';
end
figure
showsolution(coord4,elementconn,sigmarexact)
title(['e2sigmarexact'])
exportgraphics(gcf, ['e2sigmarexact.pdf'], 'ContentType', 'auto', 'BackgroundColor', 'white', 'Resolution', 300);




figure
showsolution(coord4,elementconn,sigmathetaexact)
title(['e2sigmathetaexact'])
exportgraphics(gcf, ['e2sigmathetaexact.pdf'], 'ContentType', 'auto', 'BackgroundColor', 'white', 'Resolution', 300);




dy=d4(2:2:end)
dx=d4(1:2:end)
ursbfem=sqrt(dx.^2+dy.^2)
size(d4)
size(dx)


figure
showsolution(coord4,elementconn,urexact)
title(['e2urexact'])
exportgraphics(gcf, ['e2urexact.pdf'], 'ContentType', 'auto', 'BackgroundColor', 'white', 'Resolution', 300);



size(coord4)
size(elementconn)
size(ursbfem)
figure
showsolution(coord4,elementconn,ursbfem)
title(['e2ursbfem'])
exportgraphics(gcf, ['e2ursbfem.pdf'], 'ContentType', 'auto', 'BackgroundColor', 'white', 'Resolution', 300);



urerror=urexact-ursbfem
figure
showsolution(coord4,elementconn,urerror)
title(['e2urerror'])
exportgraphics(gcf, ['e2urerror.pdf'], 'ContentType', 'auto', 'BackgroundColor', 'white', 'Resolution', 300);


sdSln4
strnNode4

e1sigmarsbfermsolutionx=[zeros(length(coord4),1)];
for i=1:length(sdSln4)
    i;
    aconn=[];
    aconn=sdSln4{i}.conn(:,1);
    anode=[];
    anode=sdSln4{i}.node;
    aconnnode=[];
    aconnnode=anode(aconn);
    atempstressx=mat.D *strnNode4{i};
    atempstressxx=atempstressx(1,:);
    for j=1:length(aconnnode)
        e1sigmarsbfermsolutionx(aconnnode(j))=atempstressxx(j);
    end
end
e1sigmarsbfermsolutionx;

size(e1sigmarsbfermsolutionx);


e1sigmarsbfermsolutiony=[zeros(length(coord4),1)];
for i=1:length(sdSln4)
    i;
    aconn=[];
    aconn=sdSln4{i}.conn(:,1);
    anode=[];
    anode=sdSln4{i}.node;
    aconnnode=[];
    aconnnode=anode(aconn);
    atempstressy=mat.D *strnNode4{i};
    atempstressyy=atempstressy(2,:);
    for j=1:length(aconnnode)
        e1sigmarsbfermsolutiony(aconnnode(j))=atempstressyy(j);
    end
end
e1sigmarsbfermsolutiony;

size(e1sigmarsbfermsolutiony);


e1gammyxy=[zeros(length(coord4),1)];
for i=1:length(sdSln4)
    i;
    aconn=[];
    aconn=sdSln4{i}.conn(:,1);
    anode=[];
    anode=sdSln4{i}.node;
    aconnnode=[];
    aconnnode=anode(aconn);
    atempstressy3=mat.D *strnNode4{i};
    atempstressyy3=atempstressy3(3,:);
    for j=1:length(aconnnode)
        e1gammyxy(aconnnode(j))=atempstressyy3(j);
    end
end
e1gammyxy;

size(e1gammyxy);



atheta=atan2(coord4(:,2),coord4(:,1))



nnnnn=coord4./vecnorm(coord4,2,2)



AAA=repmat([ 0 0 1],length(nnnnn),1)
BBB=[nnnnn  zeros(length(nnnnn),1)]

 

CCC = cross(AAA,BBB)

TAO=CCC(:,1:2)





sigmarsbfem=abs(e1sigmarsbfermsolutionx.*nnnnn(:,1)+ e1sigmarsbfermsolutiony .*nnnnn(:,2))

sigmathetasbfem=e1sigmarsbfermsolutionx.*TAO(:,1)+ e1sigmarsbfermsolutiony .*TAO(:,2)


costa=cos(atheta)
sinta=sin(atheta)
cos2ta=costa.^2
sin2ta=sinta.^2
sincosta=sinta.*costa



%
% Tmatraix=[ cos(atheta)         sin(atheta);
%           -sin(atheta)         cos(atheta)             ]
% size(Tmatraix)
%
% sigmaxy=[e1sigmarsbfermsolutionx e1gammyxy;
%          e1gammyxy               e1sigmarsbfermsolutiony]
% size(sigmaxy)
% sr=Tmatraix.*sigmaxy.*Tmatraix'







% 
% sigmarsbfem=e1sigmarsbfermsolutionx.*cos2ta+e1sigmarsbfermsolutiony.*sin2ta+e1gammyxy.*sincosta
% size(sigmarsbfem)
% size(coord4)
% size(elementconn)



figure
showsolution(coord4,elementconn,sigmarsbfem)
title(['e2sigmarsbfem'])
exportgraphics(gcf, ['e2sigmarsbfem.pdf'], 'ContentType', 'auto', 'BackgroundColor', 'white', 'Resolution', 300);

figure
showsolution(coord4,elementconn,sigmarexact-sigmarsbfem)
title(['e2sigmarerr'])
exportgraphics(gcf, ['e2sigmarerr.pdf'], 'ContentType', 'auto', 'BackgroundColor', 'white', 'Resolution', 300);





sigmathetasbfem=e1sigmarsbfermsolutionx.*sin2ta+e1sigmarsbfermsolutiony.*cos2ta-e1gammyxy.*sincosta
figure
showsolution(coord4,elementconn,sigmathetasbfem)
title(['e2sigmathetasbfem'])
exportgraphics(gcf, ['e2sigmathetasbfem.pdf'], 'ContentType', 'auto', 'BackgroundColor', 'white', 'Resolution', 300);



figure
showsolution(coord4,elementconn,sigmathetaexact-sigmathetasbfem)
title(['e2sigmathetaerr'])
exportgraphics(gcf, ['e2sigmathetaerr.pdf'], 'ContentType', 'auto', 'BackgroundColor', 'white', 'Resolution', 300);


