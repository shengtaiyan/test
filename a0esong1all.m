
close all
clc
clear
dbstop error;
tic
addpath(genpath('distmesh'))


fd=@(p) max(ddiff(dcircle(p,0,0,5),dcircle(p,0,0,1)),-min(p(:,1),p(:,2)));
% [p0,t0]=distmesh2d(fd,@huniform,0.63,[0,0;5,5],[0,1;0,5;5,0;1,0;0,1])

 [p0,t0]=distmesh2d(fd,@huniform,1,[0,0;5,5],[0,1;0,5;5,0;1,0;0,1])


tricoord0x0      =p0(:,1)
tricoord0y0      =p0(:,2)
tricoord00       =[tricoord0x0 tricoord0y0]
trisdcoon0       =[t0(:,1) t0(:,2) t0(:,3)]
figure
plot(tricoord0x0,tricoord0y0,'r*')

for i=1:length(p0)
    str = i
    text(p0(i,1)+0.1,p0(i,2)+0.1,num2str(str))

end
axis equal
title('tri-point')
tx   =   tricoord0x0(trisdcoon0);
ty   =   tricoord0y0(trisdcoon0);
figure
for itri=1:size(trisdcoon0,1)
    x   =   tricoord0x0(trisdcoon0(itri,:));
    y   =   tricoord0y0(trisdcoon0(itri,:));
    patch(x,y,'white')
    axis equal
    hold on
end
title('tri-mesh')
mat.D                       = IsoElasMtrx(1, 0.3);
mat.den                     = 2
p                           = tricoord00
t                           = trisdcoon0
[ coord0,  sdConn0, sdsc0 ] = triToSBFEMeshWithPlot( p, t )
figure
opt        = struct('LineSpec','-k', 'sdsc0',sdsc0);
PlotSBFEMesh(coord0, sdConn0, opt);
for i=1:length(coord0)
     
    str = i
    text(coord0(i,1)+0.1,coord0(i,2)+0.1,num2str(str))
end
title('MESH0');
coord0
sdConn0
sdsc0

sdConn0{6}=[ 43    15
             15    34
             34    32
             32    38 
             38    43]
figure
PlotSBFEMesh(coord0, sdConn0, opt);
for i=1:size(sdConn0,1)
    sdsc0(i,:)
    sdsc0(i,1)
    sdsc0(i,2)
    str = i
    text(sdsc0(i,1),sdsc0(i,2),num2str(str))
end
text(coord0(61,1)+0.1,coord0(61,2)+0.1,'A')
text(coord0(71,1)+0.1,coord0(71,2)+0.1,'B')
 
title('MESH-point');
figure
for i=1:length(coord0)
     
    str = i
    text(coord0(i,1)+0.1,coord0(i,2)+0.1,num2str(str))
end
 
 axis([0,6,0,6])
title('point1');

 

the1=atan2(coord0(38,2),coord0(38,1))

the1x=cos(the1)
the1y=sin(the1)
coord0(38,:)=[the1x the1y]

coord0(43,:)

the2=atan2(coord0(43,2),coord0(43,1))

the2x=cos(the2)
the2y=sin(the2)
coord0(43,:)=[the2x the2y]
figure

for i=1:size(sdConn0,1)
    sdsc0(i,:)
    sdsc0(i,1)
    sdsc0(i,2)
    str = i
    text(sdsc0(i,1),sdsc0(i,2),num2str(str))
end
 axis([0,6,0,6])
 
title('MESH-point2');


figure
opt        = struct('LineSpec','-k', 'sdsc0',sdsc0);
PlotSBFEMesh(coord0, sdConn0, opt);
for i=1:length(coord0)
     
    str = i;
    text(coord0(i,1)+0.1,coord0(i,2)+0.1,num2str(str))
end
 axis([0,5,0,5])
title('MESH2');


figure
opt        = struct('LineSpec','-k', 'sdsc0',sdsc0);
PlotSBFEMesh(coord0, sdConn0, opt);
 
 axis([0,5,0,5])
title('MESH2');


 

%a0esong1step1




% displacement boudary condition 位移边界条件
eps        =    1d-5;
dispbcx    =    find(abs(coord0(:,2))<eps );
dispbcy    =    find(abs(coord0(:,1))<eps );
BC_Dispx   =    [dispbcx     2*  ones(size(dispbcx,1),1)  zeros(size(dispbcx,1),1)];

BC_Dispy   =    [dispbcy     ones(size(dispbcy,1),1)  zeros(size(dispbcy,1),1)];
BC_Disp0   =    [BC_Dispx; BC_Dispy];
% BC_Disp0(19,:)=[71 1 0]
dispbc     =    [dispbcx;   dispbcy];
figure
plot(coord0(dispbc,1),coord0(dispbc,2),'b*' )
axis([0,6,0,6])
title('bc-disp')


%find point for force and pressure
pointnumber   =    size(coord0,1);
point         =    zeros(pointnumber,1);
for ipointnumber=1:pointnumber
    if(coord0(ipointnumber,1)^2+coord0(ipointnumber,2)^2-1<eps)
        point(ipointnumber)=ipointnumber;
        ipointnumber=ipointnumber+1;
    end
end
point;
fpoint=find(point>eps);
figure
plot(coord0(fpoint,1),coord0(fpoint,2),'r+')
title('f-point')


%find edge
edge=[];
sizesdconn=length(sdConn0);
for isizesdconn=1:sizesdconn
    isizesdconn;
    sdConn0{isizesdconn};
    sizesdconn1=size(sdConn0{isizesdconn},1);
    for  jsizesdconn1=1:sizesdconn1
        lia=intersect(fpoint,sdConn0{isizesdconn}(jsizesdconn1,:),'stable');
        if size(lia,1)==2
            edge=[edge;sdConn0{isizesdconn}(jsizesdconn1,:)];
        end
    end
end
edge;
%calculate angle
ledge=[];
edgealpha=[];
for i=1:length(edge)
    a=edge(i,1);
    b=edge(i,2);
    coordax=coord0(a,1);
    coorday=coord0(a,2);
    coordbx=coord0(b,1);
    coordby=coord0(b,2);
    ledge(i)=sqrt((coordbx-coordax)^2+(coordby-coorday)^2);
    edgealpha(i)=atan((coordby-coorday)/(coordax-coordbx));

end
ledge;
edgealpha;

angle=[];
fpoint;
coord0(fpoint,:);
x         =coord0(edge,1);
y         =coord0(edge,2);
angle     =atan(y./x);
edge1     =edge(:,1);
edge2     =edge(:,2);
x1        =coord0(edge1,1);
y1        =coord0(edge1,2);
angle1    =atand(y1./x1);
x2        =coord0(edge2,1);
y2        =coord0(edge2,2);
angle2    =atand(y2./x2);
trac      =[cosd(angle1) sind(angle1)  cosd(angle2) sind(angle2) ]';
 
ndn = 2;
NDof = ndn*size(coord0,1);
F0 = zeros(NDof,1);
F0 = addSurfTraction(coord0, edge, trac, F0);



% solve
%ltx {\bf Solution of S-elements and assemblage of global stiffness}
%ltx {\bf and mass matrices}
[sdSln, K0, M0] = SBFEMAssembly(coord0, sdConn0, sdsc0, mat);
%ltx {\bf Static solution of nodal displacements and forces}
[d0, F0] = SolverStatics(K0, BC_Disp0, F0);


%ltx {\bf Plot deformed mesh}
figure
opt = struct('MagnFct', 0.1, 'Undeformed','--k');
PlotDeformedMesh(d0, coord0, sdConn0, opt)
disp('Nodal displacements')
for ii = 1:length(coord0)
    fprintf('%5d %25.15e %25.15d\n',ii, d0(2*ii-1:2*ii))
end
sdSln;
sdStrnMode         =  SElementStrainMode2NodeEle( sdSln );%ltx integration constants
sdIntgConst        =  SElementIntgConst( d0, sdSln );
%ltx displacements and strains at specified raidal coord0inate
d0;
sdsc0

size(sdsc0,1);
nsd=size(sdsc0,1);
 
 


for isd=1:nsd
    xi = 1
    %ltx radial coord0inate
    [nodexy{isd},dsp{isd},strnNode{isd},GPxy{isd},strnEle{isd}]=SElementInDispStrain(xi,sdSln{isd},sdStrnMode{isd},sdIntgConst{isd});
end


for isd=1:nsd
    ximid = 0.5
    %ltx radial coord0inate
    [nodexymid{isd},dspmid{isd},strnNodemid{isd},GPxymid{isd},strnElemid{isd}]=SElementInDispStrain(ximid,sdSln{isd},sdStrnMode{isd},sdIntgConst{isd});
end
strnElemid;


for isd=1:nsd
    xisdsc = 0
    %ltx radial coord0inate
    [nodexysdsc{isd},dspsdsc{isd},strnNodesdsc{isd},GPxysdsc{isd},strnElesdsc{isd}]=SElementInDispStrain(xisdsc,sdSln{isd},sdStrnMode{isd},sdIntgConst{isd});
end
strnElesdsc;
strnNodesdsc;
 
sdscstrain=[];
for i=1:nsd
sdscstrain(i,:)=strnNodesdsc{i}(:,1)';
end

 

magnFct            = 0.1;
Umax               = max(abs(d0));%ltx maximum displacement
Lmax               = max(max(coord0)-min(coord0));
%ltx maximum dimension of domain
fct                = magnFct*Lmax/Umax;
%ltx factor to magnify the displacement ltx augment nodal coord0inates
dispmentforplot    = fct*(reshape(d0,2,[]))';
strnEle{isd};
 
figure
plot(coord0(:,1),coord0(:,2),'b+')
axis equal
title('point1')
% stress recovery

[recoverystress]=spr(sdSln,strnEle,strnNode,mat,nsd,coord0,d0,GPxy,strnElemid);
recoverystress;

size(recoverystress);

 
[sdscstressre]=sprsdscstress(recoverystress,sdSln,coord0,nsd,mat,strnNodesdsc);
 
 
%error estimation use recovery stress

% [eta,u,e,sdSln]=errorestimation(sdIntgConst,sdStrnMode,strnNodea,strnElea,nsd,recoverystressa,mat,sdSln)
elementnumber=nsd;



  
 [eta,matupoly,absepoly,sdSln]=errorestimation(sdIntgConst,sdStrnMode,strnNode,strnEle,elementnumber,recoverystress,mat,sdSln,sdscstressre,strnElemid,sdscstrain);

 


% refine stragety 找到需要细化的polygon，然后判断是几边形，然后细化，在计算，在判断，直到达到细化要求
uall=sum(matupoly);
eall=sum(absepoly);




etaoverline       =sqrt(eall^2/(uall^2+eall^2) );
etamax            =0.001
if etaoverline > etamax%全局判断
    em            =etamax*sqrt((uall^2+eall^2)/nsd)     ;     %单元平均允许误差
    ksi           =cell(size(strnEle,2),1);
    for i=1:nsd
        ksi{i}    = absepoly(i)/ em;
    end
    mark=[];
    for i=1:nsd
        if ksi{i} > 1
            mark(i)=i;
        end
    end
    mark
    size(mark);


    nonzeromark=nonzeros(mark);



    maxnodenumber={};


    for i=nonzeromark'
        i;
        disp(i);
        maxnodenumber{i}=max(sdSln{i}.node);
    end

    maxnodenumber;

    maxmaxnodenumber= max(cell2mat(maxnodenumber));

    elementsize=[];

    for i=nonzeromark'
        i;
        elementsize(i)=length(sdSln{i}.node);
    end
    elementsize;


    addnodenumberstar(1)=maxmaxnodenumber+1;


    length(nonzeromark);

    for i=2:length(nonzeromark)

        addnodenumberstar(i)=addnodenumberstar( i-1)+2*elementsize(nonzeromark(i-1));
    end
    addnodenumberstar;



    [mm,nn]= find(elementsize);


    [m,n]= size(elementsize);
    addnodenumberstarnew=zeros(m,n);
    addnodenumberstarnew(nn)=addnodenumberstar;




    addnewnodecoord={};
    for i=nonzeromark'
        i;
        [addnewnodecoord{i}]=addnodecoord(sdSln,nsd,mark,mat,i,coord0,sdConn0,sdsc0,addnodenumberstarnew);
    end
    addnewnodecoord;
    addnewnodecoordm=[];
    for i=nonzeromark'
        i;
        Atemp=[];
        Atemp=addnewnodecoord{i};
        addnewnodecoordm=[addnewnodecoordm;Atemp];
    end
    addnewnodecoordm;

    coord0;
    tempdis=[];
    for i=1:length(addnewnodecoordm)
        i;
        icoordx=addnewnodecoordm(i,2);
        icoordy=addnewnodecoordm(i,3);
        for j=i+1:length(addnewnodecoordm)
            j;
            jcoordx=addnewnodecoordm(j,2);
            jcoordy=addnewnodecoordm(j,3);
            tempdis=sqrt((jcoordx-icoordx)^2+(jcoordy-icoordy)^2);
            if tempdis<1e-5
                addnewnodecoordm(j,1)=addnewnodecoordm(i,1);
            end
        end
    end
    addnewnodecoordm;
    addnewnodecoordmsort=sortrows(addnewnodecoordm,1)
    tempd=[];
    for i=1:length(addnewnodecoordmsort)
        i;
        inode=addnewnodecoordmsort(i,1);
        for j=i+1:length(addnewnodecoordmsort)
            j;
            jnode=addnewnodecoordmsort(j,1);
            tempd=jnode-inode;
            if tempd<0.5
                addnewnodecoordmsort(j,:)=[0 0 0];
            end
        end
    end

    addnewnodecoordmsort;

    size(addnewnodecoordmsort);
    addnewnodecoordmsort(all(addnewnodecoordmsort==0,2),:)=[];
    size(addnewnodecoordmsort);
    addnewnodecoordmsort;
    addnewnodecoordmsortu=addnewnodecoordmsort;

    coord0;
    length(coord0);

    addnewnodecoordmsorturenode=[addnewnodecoordmsortu  (length(coord0)+1:length(coord0)+length(addnewnodecoordmsortu))'];

    addnewnodecoordmsorturenode;


    addcoord=[addnewnodecoordmsorturenode(:,end-2:end)];

    addnewnodecoord;
    addnewnodepoly={};
    addnewnodecoordrefine=[];
    for i=nonzeromark'
        i
        addnewnodecoordrefine=addnewnodecoord{i};
        for j=1:length(addnewnodecoordrefine)
            j;
            jcoordx=addnewnodecoordrefine(j,2);
            jcoordy=addnewnodecoordrefine(j,3);
            tmepdisten=[];
            for k=1:length(addcoord)
                kcoordx=addcoord(k,1);
                kcoordy=addcoord(k,2);
                tmepdisten=sqrt((kcoordx-jcoordx)^2+(kcoordy-jcoordy)^2);
                if tmepdisten<1e-5

                    addnewnodecoordrefine(j,1)=addcoord(k,3);
                    addnewnodepoly{i}=addnewnodecoordrefine;

                end
            end
        end
    end


    addnewnodepoly;
    sdConn0;
    coord0;
    sdsc0;
    for i=nonzeromark'

        i;
        [coord0,sdConn0,sdsc0]=addnewpoly(addnewnodepoly,sdSln,nsd,mark,mat,i,maxmaxnodenumber,coord0,sdConn0,sdsc0);


    end
 
    coord0=[coord0;addcoord(:,1:2)];



    %%找到没有细化的单元，给他细化了的边上增加一个节点
    %%注意 有的单元会有多个边被细化，需要增加几个点
    for i=1:length(sdConn0)
        i;
        ismember(nonzeromark,i);
        if ismember(nonzeromark,i)==0
            tempelement =sdConn0{i};
            sdnode      =tempelement(:,1);
            tempcoord0   =coord0(sdnode,:);
            for j= length(tempelement):-1:1;
                j;
                tempnode       = tempelement(j,:);
                tempnodecoord0 = [coord0(tempnode(1),1)+coord0(tempnode(2),1)   coord0(tempnode(1),2)+coord0(tempnode(2),2) ] /2 ;
                for k= 1:length(coord0)
                    addnode        =k;
                    addnodecoord0   =coord0(k,:);
                    tempdisplment  =sqrt((coord0(k,1)-tempnodecoord0(1))^2+(coord0(k,2)-tempnodecoord0(2))^2);

                    if  tempdisplment<1e-5
                        tempelement(length(tempelement)+1,:)=[0 0];
                        tempelement(j+2:end,:)=tempelement(j+1:end-1,:);
                        tempelement(j,:)      =[tempnode(1) k];
                        tempelement(j+1,:)    =[k           tempnode(2)];
                    end
                end
            end
            tempelement;
            sdConn0{i}=tempelement;
        end
    end




elseif etaoverline<=etamax
    disp('convergence')
end


sdConn1=sdConn0;
sdsc1=sdsc0;
coord1=coord0;

figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc1);
for i=1:size(sdConn1,1)
    sdsc1(i,:);
    sdsc1(i,1);
    sdsc1(i,2);
    str = i;
    text(sdsc1(i,1),sdsc1(i,2),num2str(str));
end
PlotSBFEMesh(coord1, sdConn1, opt);
title('MESH1');

figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc1);

PlotSBFEMesh(coord1, sdConn1, opt);
title('MESH1-1');


figure
for i=1:length(sdConn1)
    sdsc1(i,:);
    sdsc1(i,1);
    sdsc1(i,2);
    str = i;
    text(sdsc1(i,1),sdsc1(i,2),num2str(str));
end
axis([0,6,0,6])

title('MESH11');


figure
plot(coord1(:,1),coord1(:,2),'b+')
axis equal
title('point1')



figure
opt        = struct('LineSpec','-k', 'sdsc1',sdsc1);
PlotSBFEMesh(coord1, sdConn1, opt);
for i=1:length(coord1)

    str = i
    text(coord1(i,1)+0.1,coord1(i,2)+0.1,num2str(str))
end
axis([0,5,0,5])
title('point1');



%a0esong1step2


 

% displacement boudary condition 位移边界条件
eps        =    1d-5;
dispbc1x1    =    find(abs(coord1(:,2))<eps );
dispbc1y1    =    find(abs(coord1(:,1))<eps );
BC_Dispx1   =    [dispbc1x1     2*  ones(size(dispbc1x1,1),1)  zeros(size(dispbc1x1,1),1)];
BC_Dispy1   =    [dispbc1y1     ones(size(dispbc1y1,1),1)  zeros(size(dispbc1y1,1),1)];
BC_Disp1   =    [BC_Dispx1;  BC_Dispy1];
dispbc1     =    [dispbc1x1;   dispbc1y1];
figure
plot(coord1(dispbc1,1),coord1(dispbc1,2),'b*' )
axis([0,6,0,6])
title('bc-disp1')



%find point1 for force and pressure
point1number1   =    size(coord1,1);
point1         =    zeros(point1number1,1);
for i=1:point1number1
    if(coord1(i,1)^2+coord1(i,2)^2-1<eps)
        point1(i)=i;
        i=i+1;
    end
end
point1;
fpoint1=find(point1>eps);
figure
plot(coord1(fpoint1,1),coord1(fpoint1,2),'r+')
title('f-point1')


%find edge1
edge1=[];
sizesdconn1=length(sdConn1);
for i=1:sizesdconn1
    i;
    sdConn1{i};
    sizesdconn11=size(sdConn1{i},1)
    for  jsizesdconn11=1:sizesdconn11
        lia=intersect(fpoint1,sdConn1{i}(jsizesdconn11,:),'stable');
        if size(lia,1)==2
            edge1=[edge1;sdConn1{i}(jsizesdconn11,:)];
        end
    end
end
edge1;
%calculate angle1



angle1=[];
fpoint1;
coord1(fpoint1,:)
x         =coord1(edge1,1);
y         =coord1(edge1,2);
angle1     =atan(y./x);
edge11     =edge1(:,1);
edge12     =edge1(:,2);
x1        =coord1(edge11,1);
y1        =coord1(edge11,2);
angle11    =atan(y1./x1);
x2        =coord1(edge12,1);
y2        =coord1(edge12,2);
angle12    =atan(y2./x2);
trac1      =[cos(angle11)   sin(angle11)   cos(angle12)   sin(angle12)]';
%add force or pressure
ndn1 = 2;
NDof1 = ndn1*size(coord1,1);
F1 = zeros(NDof1,1);
F1 = addSurfTraction(coord1, edge1, trac1, F1);



% solve
%ltx {\bf Solution of S-elements and assemblage of global stiffness}
%ltx {\bf and mass matrices}
[sdSln1, K1, M1] = SBFEMAssembly(coord1, sdConn1, sdsc1, mat);
%ltx {\bf Static solution of nodal displacements and forces}
[d1, F1] = SolverStatics(K1, BC_Disp1, F1);
d1;
 

%ltx {\bf Plot deformed mesh}
figure
opt = struct('MagnFct', 0.1, 'Undeformed','--k');
PlotDeformedMesh(d1, coord1, sdConn1, opt)
disp('Nodal displacements')
for ii = 1:length(coord1)
    fprintf('%5d %25.15e %25.15d\n',ii, d1(2*ii-1:2*ii))
end
sdSln1;
sdStrnMode1         =  SElementStrainMode2NodeEle( sdSln1 );%ltx integration constants
sdIntgConst1        =  SElementIntgConst( d1, sdSln1 );
%ltx displacements and strains at specified raidal coord1inate
d1;
sdsc1

size(sdsc1,1)
nsd1=size(sdsc1,1)
for isd=1:nsd1
    xi = 1;    %ltx radial coord1inate
    [nodexy1{isd},dsp1{isd},strnNode1{isd},GPxy1{isd},strnEle1{isd}]=SElementInDispStrain(xi,sdSln1{isd},sdStrnMode1{isd},sdIntgConst1{isd});
end

for isd=1:nsd1
    ximid = 0.5;    %ltx radial coord1inate
    [nodexy1mid{isd},dsp1mid{isd},strnNode1mid{isd},GPxy1mid{isd},strnEle1mid{isd}]=SElementInDispStrain(ximid,sdSln1{isd},sdStrnMode1{isd},sdIntgConst1{isd});
end



for isd=1:nsd1
    xisdsc = 0;    %ltx radial coord1inate
    [nodexy1sdsc{isd},dsp1sdsc{isd},strnNode1sdsc{isd},GPxy1sdsc{isd},strnEle1sdsc{isd}]=SElementInDispStrain(xisdsc,sdSln1{isd},sdStrnMode1{isd},sdIntgConst1{isd});
end




sdscstrain1=[];
for i=1:nsd1
    sdscstrain1(i,:)=strnNode1sdsc{i}(:,1)';
end

sdscstrain1;
 

magnfct11            = 0.1;
Umax1               = max(abs(d1));%ltx maximum displacem1ent
Lmax1               = max(max(coord1)-min(coord1));
%ltx maximum dimension of domain
fct1                = magnfct11*Lmax1/Umax1;
%ltx factor to magnify the displacem1ent ltx augment nodal coord1inates
dispmentforplot1    = fct1*(reshape(d1,2,[]))';


% stress recovery
 
 

[recoverystress1]  =spr(sdSln1,strnEle1,strnNode1,mat,nsd1,coord1,d1,GPxy1,strnEle1mid);

[sdscstressre1]=sprsdscstress(recoverystress1,sdSln1,coord1,nsd1,mat,strnNode1sdsc);



% error estimation use recovery stress

% [eta1,u1,e1,sdSln1]  =errorestimation(sdIntgConst1,sdStrnMode1,strnNode1a,strnEle1a,nsd1,recoverystress1a,mat,sdSln1);
elementnumber1=nsd1;
% [eta1,u1,e1,sdSln1]=errorestimation(sdIntgConst1,sdStrnMode1,strnNode1a,strnEle1a,elementnumber1,recoverystress1a,mat,sdSln1)




[eta1,matupoly1,absepoly1,sdSln1]=errorestimation(sdIntgConst1,sdStrnMode1,strnNode1,strnEle1,elementnumber1,recoverystress1,mat,sdSln1,sdscstressre1,strnEle1mid,sdscstrain1);



% refine stragety 找到需要细化的polygon，然后判断是几边形，然后细化，在计算，在判断，直到达到细化要求

uall1=sum(matupoly1);
eall1=sum(absepoly1);



eta1overline1                =sqrt(eall1^2/(uall1^2+eall1^2) );
eta1max=0.001
if eta1overline1 > eta1max%全局判断
    em1            =eta1max*sqrt((uall1^2+eall1^2)/nsd1)     ;     %单元平均允许误差
    ksi1           =cell(size(strnNode1,2),1);
    for i=1:nsd1
        ksi1{i}    = absepoly1(i)/ em1;
    end
    mark1=[];
    for i=1:nsd1
        if ksi1{i} > 1
            mark1(i)=i;
        end
    end
    mark1;
    size(mark1);


    nonzeromark1=nonzeros(mark1);



    maxnodenumber1={};



    for i=nonzeromark1'
        i;
        disp(i);
        maxnodenumber1{i}=max(sdSln1{i}.node);
    end

    maxnodenumber1;

    maxmaxnodenumber1= max(cell2mat(maxnodenumber1));

    elem1entsize1=[]

    for i=nonzeromark1'
        i;
        elem1entsize1(i)=length(sdSln1{i}.node);
    end
    elem1entsize1;

 

    length(elem1entsize1);
    addnodenumberstar1=[];


    addnodenumberstar1(1)=maxmaxnodenumber1+1;



    for j=2:length(nonzeromark1)
        j;
        addnodenumberstar1(j)=addnodenumberstar1(j-1)+2*elem1entsize1(nonzeromark1(j-1));
    end

  

    addnodenumberstar1y(logical(mark1))=addnodenumberstar1;



    addnewnodecoord1={};
    nonzeromark1;
    for i=nonzeromark1'
        [addnewnodecoord1{i}]=addnodecoord(sdSln1,nsd1,mark1,mat,i,coord1,sdConn1,sdsc1 ,addnodenumberstar1y);
    end
    addnewnodecoord1;
    addnewnodecoord1m1=[];
    for i=nonzeromark1'
        i;
        Atem1p=[];
        Atem1p=addnewnodecoord1{i};
        addnewnodecoord1m1=[addnewnodecoord1m1;Atem1p];
    end
    addnewnodecoord1m1;

    coord1;
    tem1pdis=[];
    for i=1:length(addnewnodecoord1m1)
        i;
        icoordx=addnewnodecoord1m1(i,2);
        icoordy=addnewnodecoord1m1(i,3);
        for j=i+1:length(addnewnodecoord1m1)
            j;
            jcoordx=addnewnodecoord1m1(j,2);
            jcoordy=addnewnodecoord1m1(j,3);
            tem1pdis=sqrt((jcoordx-icoordx)^2+(jcoordy-icoordy)^2);
            if tem1pdis<1e-5
                addnewnodecoord1m1(j,1)=addnewnodecoord1m1(i,1);
            end
        end
    end
    addnewnodecoord1m1
    addnewnodecoord1m1sort1=sortrows(addnewnodecoord1m1,1)
    tem1pd=[];
    for i=1:length(addnewnodecoord1m1sort1)
        i;
        inode=addnewnodecoord1m1sort1(i,1);
        for j=i+1:length(addnewnodecoord1m1sort1)
            j;
            jnode=addnewnodecoord1m1sort1(j,1);
            tem1pd=jnode-inode;
            if tem1pd<0.5
                addnewnodecoord1m1sort1(j,:)=[0 0 0];
            end
        end
    end

    addnewnodecoord1m1sort1;

    size(addnewnodecoord1m1sort1);
    addnewnodecoord1m1sort1(all(addnewnodecoord1m1sort1==0,2),:)=[];
    size(addnewnodecoord1m1sort1);
    addnewnodecoord1m1sort1;
    addnewnodecoord1m1sort1u=addnewnodecoord1m1sort1;

    coord1;
    length(coord1);

    addnewnodecoord1m1sort1urenode1=[addnewnodecoord1m1sort1u  (length(coord1)+1:length(coord1)+length(addnewnodecoord1m1sort1u))'];

    addnewnodecoord1m1sort1urenode1;


    addcoord1=[addnewnodecoord1m1sort1urenode1(:,end-2:end)];

    addnewnodecoord1;
    addnewnodepoly1={};
    addnewnodecoord1refine1=[];
    for i=nonzeromark1'
        i;
        addnewnodecoord1refine1=addnewnodecoord1{i};
        for j=1:length(addnewnodecoord1refine1)
            j;
            jcoordx=addnewnodecoord1refine1(j,2);
            jcoordy=addnewnodecoord1refine1(j,3);
            tmepdisten=[];
            for k=1:length(addcoord1)
                kcoordx=addcoord1(k,1);
                kcoordy=addcoord1(k,2);
                tmepdisten=sqrt((kcoordx-jcoordx)^2+(kcoordy-jcoordy)^2);
                if tmepdisten<1e-5

                    addnewnodecoord1refine1(j,1)=addcoord1(k,3);
                    addnewnodepoly1{i}=addnewnodecoord1refine1;

                end
            end
        end
    end


    addnewnodepoly1;
    sdConn1;
    coord1;
    sdsc1;
    for i=nonzeromark1'

        i;
        [coord1,sdConn1,sdsc1 ]=addnewpoly(addnewnodepoly1,sdSln1,nsd1,mark1,mat,i,maxmaxnodenumber1,coord1,sdConn1,sdsc1 );


    end

 

    coord1=[coord1;addcoord1(:,1:2)];



    %%找到没有细化的单元，给他细化了的边上增加一个节点
    %%注意 有的单元会有多个边被细化，需要增加几个点
    for i=1:length(sdConn1)
        i;
        ismember(nonzeromark1,i);
        if ismember(nonzeromark1,i)==0
            tem1pelem1ent =sdConn1{i};
            sdnode        =tem1pelem1ent(:,1);
            tem1pcoord1   =coord1(sdnode,:);
            for j= length(tem1pelem1ent):-1:1;
                j;
                tem1pnode       = tem1pelem1ent(j,:);
                tem1pnodecoord1 = [coord1(tem1pnode(1),1)+coord1(tem1pnode(2),1)   coord1(tem1pnode(1),2)+coord1(tem1pnode(2),2) ] /2 ;
                for k= 1:length(coord1)
                    addnode        =k;
                    addnodecoord1   =coord1(k,:);
                    tem1pdisplment  =sqrt((coord1(k,1)-tem1pnodecoord1(1))^2+(coord1(k,2)-tem1pnodecoord1(2))^2);

                    if  tem1pdisplment<1e-5
                        tem1pelem1ent(length(tem1pelem1ent)+1,:)=[0 0];
                        tem1pelem1ent(j+2:end,:)=tem1pelem1ent(j+1:end-1,:);
                        tem1pelem1ent(j,:)      =[tem1pnode(1) k];
                        tem1pelem1ent(j+1,:)    =[k           tem1pnode(2)];
                    end
                end
            end
            tem1pelem1ent;
            sdConn1{i}=tem1pelem1ent;
        end
    end




elseif eta1overline1<=eta1max
    disp('convergence')
end


sdConn2=sdConn1;
sdsc2=sdsc1;
coord2=coord1;

figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc2);
for i=1:size(sdConn2,1)
    sdsc2(i,:);
    sdsc2(i,1);
    sdsc2(i,2);
    str = i;
    text(sdsc2(i,1),sdsc2(i,2),num2str(str));
end
PlotSBFEMesh(coord2, sdConn2, opt);
title('MESH2');


figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc2);

PlotSBFEMesh(coord2, sdConn2, opt);
title('MESH2-2');


 


figure
plot(coord2(:,1),coord2(:,2),'b+')
axis equal
title('point2')

figure
opt        = struct('LineSpec','-k', 'sdsc2',sdsc2);
PlotSBFEMesh(coord2, sdConn2, opt);
for i=1:length(coord2)

    str = i;
    text(coord2(i,1)+0.1,coord2(i,2)+0.1,num2str(str));
end
axis([0,5,0,5])
title('point2');


%a0esong1step3




% displacement boudary condition 位移边界条件
eps        =    1d-5
dispbc2x2    =    find(abs(coord2(:,2))<eps );
dispbc2y2    =    find(abs(coord2(:,1))<eps );
BC_Dispx2   =    [dispbc2x2     2*  ones(size(dispbc2x2,1),1)  zeros(size(dispbc2x2,1),1)];
BC_Dispy2   =    [dispbc2y2     ones(size(dispbc2y2,1),1)  zeros(size(dispbc2y2,1),1)];
BC_Disp2   =    [BC_Dispx2;  BC_Dispy2];
dispbc2     =    [dispbc2x2;   dispbc2y2];
figure
plot(coord2(dispbc2,1),coord2(dispbc2,2),'b*' )
axis([0,6,0,6])
title('bc-disp2')



%find point2 for force and pressure
point2number2   =    size(coord2,1);
point2         =    zeros(point2number2,1);
for i=1:point2number2
    if(coord2(i,1)^2+coord2(i,2)^2-1<eps)
        point2(i)=i;
        i=i+1;
    end
end
point2;
fpoint2=find(point2>eps);
figure
plot(coord2(fpoint2,1),coord2(fpoint2,2),'r+')
title('f-point2')


%find edge2
edge2=[];
sizesdconn2=length(sdConn2);
for i=1:sizesdconn2
    i;
    sdConn2{i};
    sizesdconn22=size(sdConn2{i},1);
    for  j=1:sizesdconn22
        lia=intersect(fpoint2,sdConn2{i}(j,:),'stable');
        if size(lia,1)==2
            edge2=[edge2;sdConn2{i}(j,:)];
        end
    end
end
edge2;
%calculate angle2



angle2=[];
fpoint2;
coord2(fpoint2,:);
x         =coord2(edge2,1);
y         =coord2(edge2,2);
angle2     =atan(y./x);
edge21     =edge2(:,1);
edge22     =edge2(:,2);
x1        =coord2(edge21,1);
y1        =coord2(edge21,2);
angle21    =atan(y1./x1);
x2        =coord2(edge22,1);
y2        =coord2(edge22,2);
angle22    =atan(y2./x2);
trac2      =[cos(angle21)   sin(angle21)   cos(angle22)   sin(angle22)]';
%add force or pressure
ndn2 = 2;
NDof2 = ndn2*size(coord2,1);
F2 = zeros(NDof2,1);
F2 = addSurfTraction(coord2, edge2, trac2, F2);



% solve
%ltx {\bf Solution of S-elements and assemblage of global stiffness}
%ltx {\bf and mass matrices}
[sdSln2, K2, M2] = SBFEMAssembly(coord2, sdConn2, sdsc2, mat);
%ltx {\bf Static solution of nodal displacements and forces}
[d2, F2] = SolverStatics(K2, BC_Disp2, F2);
d2;
 


%ltx {\bf Plot deformed mesh}
figure
opt = struct('MagnFct', 0.1, 'Undeformed','--k');
PlotDeformedMesh(d2, coord2, sdConn2, opt)
disp('Nodal displacements')
for ii = 1:length(coord2)
    fprintf('%5d %25.15e %25.15d\n',ii, d2(2*ii-1:2*ii))
end
sdSln2;
sdStrnMode2         =  SElementStrainMode2NodeEle( sdSln2 );%ltx integration constants
sdIntgConst2        =  SElementIntgConst( d2, sdSln2 );
%ltx displacements and strains at specified raidal coord2inate
d2;
sdsc2;

size(sdsc2,1)
nsd2=size(sdsc2,1)
for isd=1:nsd2
    xi = 1;    %ltx radial coord2inate
    [nodexy2{isd},dsp2{isd},strnNode2{isd},GPxy2{isd},strnEle2{isd}]=SElementInDispStrain(xi,sdSln2{isd},sdStrnMode2{isd},sdIntgConst2{isd});
end


for isd=1:nsd2
    ximid = 0.5;    %ltx radial coord2inate
    [nodexy2mid{isd},dsp2mid{isd},strnNode2mid{isd},GPxy2mid{isd},strnEle2mid{isd}]=SElementInDispStrain(ximid,sdSln2{isd},sdStrnMode2{isd},sdIntgConst2{isd});
end


for isd=1:nsd2
    xisdsc = 0;    %ltx radial coord2inate
    [nodexy2sdsc{isd},dsp2sdsc{isd},strnNode2sdsc{isd},GPxy2sdsc{isd},strnEle2sdsc{isd}]=SElementInDispStrain(xisdsc,sdSln2{isd},sdStrnMode2{isd},sdIntgConst2{isd});
end



sdscstrain2=[];
for i=1:nsd2
    sdscstrain2(i,:)=strnNode2sdsc{i}(:,1)';
end

 

magnfct22            = 0.1;
Umax2               = max(abs(d2));%ltx maximum displaceM2ent
Lmax2               = max(max(coord2)-min(coord2));
%ltx maximum dimension of domain
fct2                = magnfct22*Lmax2/Umax2;
%ltx factor to magnify the displaceM2ent ltx augment nodal coord2inates
dispmentforplot2    = fct2*(reshape(d2,2,[]))';


% [recoverystress2]  =spr(sdSln2,strnEle2,strnNode2,mat,nsd2,coord2,d2,GPxy2)


[recoverystress2]  =spr(sdSln2,strnEle2,strnNode2,mat,nsd2,coord2,d2,GPxy2,strnEle2mid);

[sdscstressre2]=sprsdscstress(recoverystress2,sdSln2,coord2,nsd2,mat,strnNode2sdsc);

 


% error estimation use recovery stress
elementnumber2=nsd2;


[eta2,matupoly2,absepoly2,sdSln2]=errorestimation(sdIntgConst2,sdStrnMode2,strnNode2,strnEle2,elementnumber2,recoverystress2,mat,sdSln2,sdscstressre2,strnEle2mid,sdscstrain2);




uall2=sum(matupoly2);
eall2=sum(absepoly2);
 

% refine stragety 找到需要细化的polygon，然后判断是几边形，然后细化，在计算，在判断，直到达到细化要求

eta2overline2                =sqrt(eall2^2/(uall2^2+eall2^2) );
eta2max=0.001
if eta2overline2 > eta2max%全局判断
    eM2            =eta2max*sqrt((uall2^2+eall2^2)/nsd2)     ;     %单元平均允许误差
    ksi2           =cell(size(strnNode2,2),1);
    for i=1:nsd2
        ksi2{i}    = absepoly2(i)/ eM2;
    end
    marK2=[];
    for i=1:nsd2
        if ksi2{i} > 1
            marK2(i)=i;
        end
    end
    marK2;
    size(marK2);


    nonzeromarK2=nonzeros(marK2);
 
    maxnodenumber2={};
 

    for i=nonzeromarK2'
        i;
        disp(i);
        maxnodenumber2{i}=max(sdSln2{i}.node);
    end

    maxnodenumber2;

    maxmaxnodenumber2= max(cell2mat(maxnodenumber2));

    eleM2entsize2=[];

    for i=nonzeromarK2'
        i;
        eleM2entsize2(i)=length(sdSln2{i}.node);
    end
  

    length(eleM2entsize2);
    addnodenumberstar2=[];


    addnodenumberstar2(1)=maxmaxnodenumber2+1;



    for j=2:length(nonzeromarK2)
        j;
        addnodenumberstar2(j)=addnodenumberstar2(j-1)+2*eleM2entsize2(nonzeromarK2(j-1));
    end

  
    addnodenumberstar2y(logical(marK2))=addnodenumberstar2;



    addnewnodecoord2={};
    nonzeromarK2;
    for i=nonzeromarK2'
        [addnewnodecoord2{i}]=addnodecoord(sdSln2,nsd2,marK2,mat,i,coord2,sdConn2,sdsc2,addnodenumberstar2y);
    end
    addnewnodecoord2;
    addnewnodecoord2M2=[];
    for i=nonzeromarK2'
        i;
        AteM2p=[];
        AteM2p=addnewnodecoord2{i};
        addnewnodecoord2M2=[addnewnodecoord2M2;AteM2p];
    end
    addnewnodecoord2M2;

    coord2;
    teM2pdis=[];
    for i=1:length(addnewnodecoord2M2)
        i;
        icoordx=addnewnodecoord2M2(i,2);
        icoordy=addnewnodecoord2M2(i,3);
        for j=i+1:length(addnewnodecoord2M2)
            j;
            jcoordx=addnewnodecoord2M2(j,2);
            jcoordy=addnewnodecoord2M2(j,3);
            teM2pdis=sqrt((jcoordx-icoordx)^2+(jcoordy-icoordy)^2);
            if teM2pdis<1e-5
                addnewnodecoord2M2(j,1)=addnewnodecoord2M2(i,1);
            end
        end
    end
    addnewnodecoord2M2;
    addnewnodecoord2M2sort2=sortrows(addnewnodecoord2M2,1);
    teM2pd=[];
    for i=1:length(addnewnodecoord2M2sort2)
        i;
        inode=addnewnodecoord2M2sort2(i,1);
        for j=i+1:length(addnewnodecoord2M2sort2)
            j;
            jnode=addnewnodecoord2M2sort2(j,1);
            teM2pd=jnode-inode;
            if teM2pd<0.5
                addnewnodecoord2M2sort2(j,:)=[0 0 0];
            end
        end
    end

    addnewnodecoord2M2sort2;

    size(addnewnodecoord2M2sort2);
    addnewnodecoord2M2sort2(all(addnewnodecoord2M2sort2==0,2),:)=[];
    size(addnewnodecoord2M2sort2);
    addnewnodecoord2M2sort2;
    addnewnodecoord2M2sort2u=addnewnodecoord2M2sort2;

    coord2;
    length(coord2);

    addnewnodecoord2M2sort2urenode2=[addnewnodecoord2M2sort2u  (length(coord2)+1:length(coord2)+length(addnewnodecoord2M2sort2u))'];

    addnewnodecoord2M2sort2urenode2;


    addcoord2=[addnewnodecoord2M2sort2urenode2(:,end-2:end)];

    addnewnodecoord2;
    addnewnodepoly2={};
    addnewnodecoord2refine2=[];
    for i=nonzeromarK2'
        i
        addnewnodecoord2refine2=addnewnodecoord2{i};
        for j=1:length(addnewnodecoord2refine2)
            j;
            jcoordx=addnewnodecoord2refine2(j,2);
            jcoordy=addnewnodecoord2refine2(j,3);
            tmepdisten=[];
            for k=1:length(addcoord2)
                kcoordx=addcoord2(k,1);
                kcoordy=addcoord2(k,2);
                tmepdisten=sqrt((kcoordx-jcoordx)^2+(kcoordy-jcoordy)^2);
                if tmepdisten<1e-5

                    addnewnodecoord2refine2(j,1)=addcoord2(k,3);
                    addnewnodepoly2{i}=addnewnodecoord2refine2;

                end
            end
        end
    end


    addnewnodepoly2;
    sdConn2;
    coord2;
    sdsc2;
    for i=nonzeromarK2'

        i;
        [coord2,sdConn2,sdsc2]=addnewpoly(addnewnodepoly2,sdSln2,nsd2,marK2,mat,i,maxmaxnodenumber2,coord2,sdConn2,sdsc2);


    end

 

    coord2=[coord2;addcoord2(:,1:2)];



    %%找到没有细化的单元，给他细化了的边上增加一个节点
    %%注意 有的单元会有多个边被细化，需要增加几个点
    for i=1:length(sdConn2)
        i;
        ismember(nonzeromarK2,i);
        if ismember(nonzeromarK2,i)==0
            teM2peleM2ent =sdConn2{i};
            sdnode        =teM2peleM2ent(:,1);
            teM2pcoord2   =coord2(sdnode,:);
            for j= length(teM2peleM2ent):-1:1;
                j;
                teM2pnode       = teM2peleM2ent(j,:);
                teM2pnodecoord2 = [coord2(teM2pnode(1),1)+coord2(teM2pnode(2),1)   coord2(teM2pnode(1),2)+coord2(teM2pnode(2),2) ] /2 ;
                for k= 1:length(coord2)
                    addnode        =k;
                    addnodecoord2   =coord2(k,:);
                    teM2pdisplment  =sqrt((coord2(k,1)-teM2pnodecoord2(1))^2+(coord2(k,2)-teM2pnodecoord2(2))^2);

                    if  teM2pdisplment<1e-5
                        teM2peleM2ent(length(teM2peleM2ent)+1,:)=[0 0];
                        teM2peleM2ent(j+2:end,:)=teM2peleM2ent(j+1:end-1,:);
                        teM2peleM2ent(j,:)      =[teM2pnode(1) k];
                        teM2peleM2ent(j+1,:)    =[k           teM2pnode(2)];
                    end
                end
            end
            teM2peleM2ent;
            sdConn2{i}=teM2peleM2ent;
        end
    end




elseif eta2overline2<=eta2max
    disp('convergence')
end


sdConn3=sdConn2;
sdsc3=sdsc2;
coord3=coord2;

 


cell2mat(sdConn3);

 


figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc3);

PlotSBFEMesh(coord3, sdConn3, opt);
title('MESH3-2');



figure
plot(coord3(:,1),coord3(:,2),'b+')
axis equal
title('point3')

%a0esong1step4

 

% displacement boudary condition 位移边界条件
eps        =    1d-5
dispbc3x3    =    find(abs(coord3(:,2))<eps );
dispbc3y3    =    find(abs(coord3(:,1))<eps );
BC_Dispx3   =    [dispbc3x3     2*  ones(size(dispbc3x3,1),1)  zeros(size(dispbc3x3,1),1)];
BC_Dispy3   =    [dispbc3y3     ones(size(dispbc3y3,1),1)  zeros(size(dispbc3y3,1),1)];
BC_Disp3   =    [BC_Dispx3;  BC_Dispy3];
dispbc3     =    [dispbc3x3;   dispbc3y3];
figure
plot(coord3(dispbc3,1),coord3(dispbc3,2),'b*' )
axis([0,6,0,6])
title('bc-disp3')


 
%find point3 for force and pressure
point3number3   =    size(coord3,1);
point3         =    zeros(point3number3,1);
for i=1:point3number3
    if(coord3(i,1)^2+coord3(i,2)^2-1<eps)
        point3(i)=i;
        i=i+1;
    end
end
point3;
fpoint3=find(point3>eps);
figure
plot(coord3(fpoint3,1),coord3(fpoint3,2),'r+')
title('f-point3')


%find edge3
edge3=[];
sizesdconn3=length(sdConn3);
for i=1:sizesdconn3
    i;
    sdConn3{i};
    sizesdconn33=size(sdConn3{i},1);
    for  j=1:sizesdconn33
        lia=intersect(fpoint3,sdConn3{i}(j,:),'stable');
        if size(lia,1)==2
            edge3=[edge3;sdConn3{i}(j,:)];
        end
    end
end
edge3;
%calculate angle3



angle3=[];
fpoint3;
coord3(fpoint3,:);
x         =coord3(edge3,1);
y         =coord3(edge3,2);
angle3     =atan(y./x);
edge31     =edge3(:,1);
edge32     =edge3(:,2);
x1        =coord3(edge31,1);
y1        =coord3(edge31,2);
angle31    =atan(y1./x1);
x2        =coord3(edge32,1);
y2        =coord3(edge32,2);
angle32    =atan(y2./x2);
trac3      =[cos(angle31)   sin(angle31)   cos(angle32)   sin(angle32)]';
%add force or pressure
ndn3 = 2;
NDof3 = ndn3*size(coord3,1);
F3 = zeros(NDof3,1);
F3 = addSurfTraction(coord3, edge3, trac3, F3);



% solve
%ltx {\bf Solution of S-elements and assemblage of global stiffness}
%ltx {\bf and mass matrices}
[sdSln3, K3, M3] = SBFEMAssembly(coord3, sdConn3, sdsc3, mat);
%ltx {\bf Static solution of nodal displacements and forces}
[d3, F3] = SolverStatics(K3, BC_Disp3, F3);
d3;

 


%ltx {\bf Plot deformed mesh}
figure
opt = struct('MagnFct', 0.1, 'Undeformed','--k');
PlotDeformedMesh(d3, coord3, sdConn3, opt)
disp('Nodal displacements')
for ii = 1:length(coord3)
    fprintf('%5d %25.15e %25.15d\n',ii, d3(2*ii-1:2*ii))
end
sdSln3;
sdStrnMode3         =  SElementStrainMode2NodeEle( sdSln3 );%ltx integration constants
sdIntgConst3        =  SElementIntgConst( d3, sdSln3 );
%ltx displacements and strains at specified raidal coord3inate
d3;
sdsc3;

size(sdsc3,1);
nsd3=size(sdsc3,1);
for isd=1:nsd3
    xi = 1;    %ltx radial coord3inate
    [nodexy3{isd},dsp3{isd},strnNode3{isd},GPxy3{isd},strnEle3{isd}]=SElementInDispStrain(xi,sdSln3{isd},sdStrnMode3{isd},sdIntgConst3{isd});
end

for isd=1:nsd3
    ximid = 0.5;    %ltx radial coord3inate
    [nodexy3mid{isd},dsp3mid{isd},strnNode3mid{isd},GPxy3mid{isd},strnEle3mid{isd}]=SElementInDispStrain(ximid,sdSln3{isd},sdStrnMode3{isd},sdIntgConst3{isd});
end



for isd=1:nsd3
    xisdsc = 0;    %ltx radial coord3inate
    [nodexy3sdsc{isd},dsp3sdsc{isd},strnNode3sdsc{isd},GPxy3sdsc{isd},strnEle3sdsc{isd}]=SElementInDispStrain(xisdsc,sdSln3{isd},sdStrnMode3{isd},sdIntgConst3{isd});
end

sdscstrain3=[];
for i=1:nsd3
    sdscstrain3(i,:)=strnNode3sdsc{i}(:,1)';
end

 

 

magnfct33            = 0.1;
Umax3               = max(abs(d3));%ltx maximum displaceM3ent
Lmax3               = max(max(coord3)-min(coord3));
%ltx maximum dimension of domain
fct3                = magnfct33*Lmax3/Umax3;
%ltx factor to magnify the displaceM3ent ltx augment nodal coord3inates
dispmentforplot3    = fct3*(reshape(d3,2,[]))';


 

[recoverystress3]  =spr(sdSln3,strnEle3,strnNode3,mat,nsd3,coord3,d3,GPxy3,strnEle3mid);

[sdscstressre3]=sprsdscstress(recoverystress3,sdSln3,coord3,nsd3,mat,strnNode3sdsc);

  

recoverystress3;
% error estimation use recovery stress
elementnumber3=nsd3;

 
  
 
[eta3,matupoly3,absepoly3,sdSln3]=errorestimation(sdIntgConst3,sdStrnMode3,strnNode3,strnEle3,elementnumber3,recoverystress3,mat,sdSln3,sdscstressre3,strnEle3mid,sdscstrain3);





% refine stragety 找到需要细化的polygon，然后判断是几边形，然后细化，在计算，在判断，直到达到细化要求

 uall3=sum(matupoly3);
 eall3=sum(absepoly3);
 

eta3overline3                =sqrt(eall3^2/(uall3^2+eall3^2) );
eta3max=0.001;
if eta3overline3 > eta3max%全局判断
    eM3            =eta3max*sqrt((uall3^2+eall3^2)/nsd3)     ;     %单元平均允许误差
    ksi3           =cell(size(strnNode3,2),1);
    for i=1:nsd3
        ksi3{i}    = absepoly3(i)/ eM3;
    end
    marK3=[];
    for i=1:nsd3
        if ksi3{i} > 1
            marK3(i)=i;
        end
    end
    marK3
    size(marK3);


    nonzeromarK3=nonzeros(marK3);



    maxnodenumber3={};



    for i=nonzeromarK3'
        i;
        disp(i);
        maxnodenumber3{i}=max(sdSln3{i}.node);
    end

    maxnodenumber3;

    maxmaxnodenumber3= max(cell2mat(maxnodenumber3));

    eleM3entsize3=[];

    for i=nonzeromarK3'
        i;
        eleM3entsize3(i)=length(sdSln3{i}.node);
    end
    eleM3entsize3;



    % eleM3entsize3(eleM3entsize3==0)=[];


    length(eleM3entsize3);
    addnodenumberstar3=[];


    addnodenumberstar3(1)=maxmaxnodenumber3+1;



    for j=2:length(nonzeromarK3)
        j;
        addnodenumberstar3(j)=addnodenumberstar3(j-1)+2*eleM3entsize3(nonzeromarK3(j-1));
    end
    addnodenumberstar3;

 
    addnodenumberstar3y(logical(marK3))=addnodenumberstar3;



    addnewnodecoord3={};
    nonzeromarK3;
    for i=nonzeromarK3'
        [addnewnodecoord3{i}]=addnodecoord(sdSln3,nsd3,marK3,mat,i,coord3,sdConn3,sdsc3,addnodenumberstar3y);
    end
    addnewnodecoord3;
    addnewnodecoord3M3=[];
    for i=nonzeromarK3'
        i;
        AteM3p=[];
        AteM3p=addnewnodecoord3{i};
        addnewnodecoord3M3=[addnewnodecoord3M3;AteM3p];
    end
    addnewnodecoord3M3;

    coord3;
    teM3pdis=[];
    for i=1:length(addnewnodecoord3M3)
        i;
        icoordx=addnewnodecoord3M3(i,2);
        icoordy=addnewnodecoord3M3(i,3);
        for j=i+1:length(addnewnodecoord3M3)
            j;
            jcoordx=addnewnodecoord3M3(j,2);
            jcoordy=addnewnodecoord3M3(j,3);
            teM3pdis=sqrt((jcoordx-icoordx)^2+(jcoordy-icoordy)^2);
            if teM3pdis<1e-5
                addnewnodecoord3M3(j,1)=addnewnodecoord3M3(i,1);
            end
        end
    end
    addnewnodecoord3M3
    addnewnodecoord3M3sort3=sortrows(addnewnodecoord3M3,1)
    teM3pd=[];
    for i=1:length(addnewnodecoord3M3sort3)
        i;
        inode=addnewnodecoord3M3sort3(i,1);
        for j=i+1:length(addnewnodecoord3M3sort3)
            j;
            jnode=addnewnodecoord3M3sort3(j,1);
            teM3pd=jnode-inode;
            if teM3pd<0.5
                addnewnodecoord3M3sort3(j,:)=[0 0 0];
            end
        end
    end

    addnewnodecoord3M3sort3;

    size(addnewnodecoord3M3sort3);
    addnewnodecoord3M3sort3(all(addnewnodecoord3M3sort3==0,2),:)=[];
    size(addnewnodecoord3M3sort3);
    addnewnodecoord3M3sort3;
    addnewnodecoord3M3sort3u=addnewnodecoord3M3sort3;

    coord3;
    length(coord3);

    addnewnodecoord3M3sort3urenode3=[addnewnodecoord3M3sort3u  (length(coord3)+1:length(coord3)+length(addnewnodecoord3M3sort3u))'];

    addnewnodecoord3M3sort3urenode3;


    addcoord3=[addnewnodecoord3M3sort3urenode3(:,end-2:end)];

    addnewnodecoord3;
    addnewnodepoly3={};
    addnewnodecoord3refine3=[];
    for i=nonzeromarK3'
        i;
        addnewnodecoord3refine3=addnewnodecoord3{i};
        for j=1:length(addnewnodecoord3refine3)
            j;
            jcoordx=addnewnodecoord3refine3(j,2);
            jcoordy=addnewnodecoord3refine3(j,3);
            tmepdisten=[];
            for k=1:length(addcoord3)
                kcoordx=addcoord3(k,1);
                kcoordy=addcoord3(k,2);
                tmepdisten=sqrt((kcoordx-jcoordx)^2+(kcoordy-jcoordy)^2);
                if tmepdisten<1e-5

                    addnewnodecoord3refine3(j,1)=addcoord3(k,3);
                    addnewnodepoly3{i}=addnewnodecoord3refine3;

                end
            end
        end
    end


    addnewnodepoly3;
    sdConn3;
    coord3;
    sdsc3;
    for i=nonzeromarK3'

        i;
        [coord3,sdConn3,sdsc3]=addnewpoly(addnewnodepoly3,sdSln3,nsd3,marK3,mat,i,maxmaxnodenumber3,coord3,sdConn3,sdsc3);


    end

 

    coord3=[coord3;addcoord3(:,1:2)];



    %%找到没有细化的单元，给他细化了的边上增加一个节点
    %%注意 有的单元会有多个边被细化，需要增加几个点
    for i=1:length(sdConn3)
        i;
        ismember(nonzeromarK3,i);
        if ismember(nonzeromarK3,i)==0
            teM3peleM3ent =sdConn3{i};
            sdnode        =teM3peleM3ent(:,1);
            teM3pcoord3   =coord3(sdnode,:);
            for j= length(teM3peleM3ent):-1:1;
                j;
                teM3pnode       = teM3peleM3ent(j,:);
                teM3pnodecoord3 = [coord3(teM3pnode(1),1)+coord3(teM3pnode(2),1)   coord3(teM3pnode(1),2)+coord3(teM3pnode(2),2) ] /2 ;
                for k= 1:length(coord3)
                    addnode        =k;
                    addnodecoord3   =coord3(k,:);
                    teM3pdisplment  =sqrt((coord3(k,1)-teM3pnodecoord3(1))^2+(coord3(k,2)-teM3pnodecoord3(2))^2);

                    if  teM3pdisplment<1e-5
                        teM3peleM3ent(length(teM3peleM3ent)+1,:)=[0 0];
                        teM3peleM3ent(j+2:end,:)=teM3peleM3ent(j+1:end-1,:);
                        teM3peleM3ent(j,:)      =[teM3pnode(1) k];
                        teM3peleM3ent(j+1,:)    =[k           teM3pnode(2)];
                    end
                end
            end
            teM3peleM3ent;
            sdConn3{i}=teM3peleM3ent;
        end
    end




elseif eta3overline3<=eta3max
    disp('convergence')
end


sdConn4=sdConn3;
sdsc4=sdsc3;
coord4=coord3;

figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc4);
for i=1:size(sdConn4,1)
    sdsc4(i,:);
    sdsc4(i,1);
    sdsc4(i,2);
    str = i;
    text(sdsc4(i,1),sdsc4(i,2),num2str(str));
end
PlotSBFEMesh(coord4, sdConn4, opt);
title('MESH4');


figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc4);

PlotSBFEMesh(coord4, sdConn4, opt);
title('MESH4-2');



figure
plot(coord4(:,1),coord4(:,2),'b+')
axis equal
title('point4')


 
%a0esong1step5

 

% displacement boudary condition 位移边界条件
eps        =    1d-5
dispbc4x4    =    find(abs(coord4(:,2))<eps );
dispbc4y4    =    find(abs(coord4(:,1))<eps );
BC_Dispx4   =    [dispbc4x4     2*  ones(size(dispbc4x4,1),1)  zeros(size(dispbc4x4,1),1)];
BC_Dispy4   =    [dispbc4y4     ones(size(dispbc4y4,1),1)  zeros(size(dispbc4y4,1),1)];
BC_Disp4   =    [BC_Dispx4;  BC_Dispy4];
dispbc4     =    [dispbc4x4;   dispbc4y4];
figure
plot(coord4(dispbc4,1),coord4(dispbc4,2),'b*' );
axis([0,6,0,6]);
title('bc-disp4')


 
%find point4 for force and pressure
point4number4   =    size(coord4,1);
point4         =    zeros(point4number4,1);
for i=1:point4number4
    if(coord4(i,1)^2+coord4(i,2)^2-1<eps)
        point4(i)=i;
        i=i+1;
    end
end
point4;
fpoint4=find(point4>eps);
figure
plot(coord4(fpoint4,1),coord4(fpoint4,2),'r+')
title('f-point4')


%find edge4
edge4=[];
sizesdconn4=length(sdConn4)
for i=1:sizesdconn4
    i;
    sdConn4{i};
    sizesdconn43=size(sdConn4{i},1);
    for  j=1:sizesdconn43
        lia=intersect(fpoint4,sdConn4{i}(j,:),'stable');
        if size(lia,1)==2
            edge4=[edge4;sdConn4{i}(j,:)];
        end
    end
end
edge4;
%calculate angle4



angle4=[];
fpoint4;
coord4(fpoint4,:);
x         =coord4(edge4,1);
y         =coord4(edge4,2);
angle4     =atan(y./x);
edge41     =edge4(:,1);
edge42     =edge4(:,2);
x1        =coord4(edge41,1);
y1        =coord4(edge41,2);
angle41    =atan(y1./x1);
x2        =coord4(edge42,1);
y2        =coord4(edge42,2);
angle42    =atan(y2./x2);
trac4      =[cos(angle41)   sin(angle41)   cos(angle42)   sin(angle42)]';
%add force or pressure
ndn4 = 2;
NDof4 = ndn4*size(coord4,1);
F4 = zeros(NDof4,1);
F4 = addSurfTraction(coord4, edge4, trac4, F4);



% solve
%ltx {\bf Solution of S-elements and assemblage of global stiffness}
%ltx {\bf and mass matrices}
[sdSln4, K4, M4] = SBFEMAssembly(coord4, sdConn4, sdsc4, mat);
%ltx {\bf Static solution of nodal displacements and forces}
[d4, F4] = SolverStatics(K4, BC_Disp4, F4);
d4;
 


 

%ltx {\bf Plot deformed mesh}
figure
opt = struct('MagnFct', 0.1, 'Undeformed','--k');
PlotDeformedMesh(d4, coord4, sdConn4, opt)
disp('Nodal displacements')
for ii = 1:length(coord4)
    fprintf('%5d %25.15e %25.15d\n',ii, d4(2*ii-1:2*ii))
end
sdSln4;
sdStrnMode4         =  SElementStrainMode2NodeEle( sdSln4 );%ltx integration constants
sdIntgConst4        =  SElementIntgConst( d4, sdSln4 );
%ltx displacements and strains at specified raidal coord4inate
d4;
sdsc4;

size(sdsc4,1);
nsd4=size(sdsc4,1);
for isd=1:nsd4
    xi = 1;    %ltx radial coord4inate
    [nodexy4{isd},dsp4{isd},strnNode4{isd},GPxy4{isd},strnEle4{isd}]=SElementInDispStrain(xi,sdSln4{isd},sdStrnMode4{isd},sdIntgConst4{isd});
end


for isd=1:nsd4
    ximid = 0.5;    %ltx radial coord4inate
    [nodexy4mid{isd},dsp4mid{isd},strnNode4mid{isd},GPxy4mid{isd},strnEle4mid{isd}]=SElementInDispStrain(ximid,sdSln4{isd},sdStrnMode4{isd},sdIntgConst4{isd});
end


for isd=1:nsd4
    xisdsc = 0;    %ltx radial coord4inate
    [nodexy4sdsc{isd},dsp4sdsc{isd},strnNode4sdsc{isd},GPxy4sdsc{isd},strnEle4sdsc{isd}]=SElementInDispStrain(xisdsc,sdSln4{isd},sdStrnMode4{isd},sdIntgConst4{isd});
end


sdscstrain4=[];
for i=1:nsd4
    sdscstrain4(i,:)=strnNode4sdsc{i}(:,1)';
end

 

 

magnfct44           = 0.1;
Umax4               = max(abs(d4));%ltx maximum displaceM4ent
Lmax4               = max(max(coord4)-min(coord4));
%ltx maximum dimension of domain
fct4                = magnfct44*Lmax4/Umax4;
%ltx factor to magnify the displaceM4ent ltx augment nodal coord4inates
dispmentforplot4    = fct4*(reshape(d4,2,[]))';


% stress recovery
 

[recoverystress4]  =spr(sdSln4,strnEle4,strnNode4,mat,nsd4,coord4,d4,GPxy4,strnEle4mid);

[sdscstressre4]=sprsdscstress(recoverystress4,sdSln4,coord4,nsd4,mat,strnNode4sdsc);


recoverystress4;
% error estimation use recovery stress

elementnumber4=nsd4;

 
 
[eta4,matupoly4,absepoly4,sdSln4]=errorestimation(sdIntgConst4,sdStrnMode4,strnNode4,strnEle4,elementnumber4,recoverystress4,mat,sdSln4,sdscstressre4,strnEle4mid,sdscstrain4);


% refine stragety 鎵惧埌闇€瑕佺粏鍖栫殑polygon锛岀劧鍚庡垽鏂槸鍑犺竟褰紝鐒跺悗缁嗗寲锛屽湪璁＄畻锛屽湪鍒ゆ柇锛岀洿鍒拌揪鍒扮粏鍖栬姹�
 uall4=sum(matupoly4);
 eall4=sum(absepoly4);
 
 
 
eta4overline4                =sqrt(eall4^2/(uall4^2+eall4^2) );
eta4max=0.001;
if eta4overline4 > eta4max%鍏ㄥ眬鍒ゆ柇
    eM4            =eta4max*sqrt((uall4^2+eall4^2)/nsd4)     ;     %鍗曞厓骞冲潎鍏佽璇樊
    ksi4           =cell(size(strnNode4,2),1);
    for i=1:nsd4
        ksi4{i}    = absepoly4(i)/ eM4;
    end
    marK4=[];
    for i=1:nsd4
        if ksi4{i} > 1
            marK4(i)=i;
        end
    end
    marK4
    size(marK4);


    nonzeromarK4=nonzeros(marK4);



    maxnodenumber4={};



    for i=nonzeromarK4'
        i;
        disp(i);
        maxnodenumber4{i}=max(sdSln4{i}.node);
    end

    maxnodenumber4;

    maxmaxnodenumber4= max(cell2mat(maxnodenumber4));

    eleM4entsize4=[];

    for i=nonzeromarK4'
        i;
        eleM4entsize4(i)=length(sdSln4{i}.node);
    end
    eleM4entsize4;



    % eleM4entsize4(eleM4entsize4==0)=[];


    length(eleM4entsize4);
    addnodenumberstar4=[];


    addnodenumberstar4(1)=maxmaxnodenumber4+1;



    for j=2:length(nonzeromarK4)
        j;
        addnodenumberstar4(j)=addnodenumberstar4(j-1)+2*eleM4entsize4(nonzeromarK4(j-1));
    end

    addnodenumberstar4;

 

    addnodenumberstar4y(logical(marK4))=addnodenumberstar4;



    addnewnodecoord4={};
    nonzeromarK4
    for i=nonzeromarK4'
        [addnewnodecoord4{i}]=addnodecoord(sdSln4,nsd4,marK4,mat,i,coord4,sdConn4,sdsc4,addnodenumberstar4y);
    end
    addnewnodecoord4;
    addnewnodecoord4M4=[];
    for i=nonzeromarK4'
        i;
        AteM4p=[];
        AteM4p=addnewnodecoord4{i};
        addnewnodecoord4M4=[addnewnodecoord4M4;AteM4p];
    end
    addnewnodecoord4M4;

    coord4;
    teM4pdis=[];
    for i=1:length(addnewnodecoord4M4)
        i;
        icoordx=addnewnodecoord4M4(i,2);
        icoordy=addnewnodecoord4M4(i,3);
        for j=i+1:length(addnewnodecoord4M4)
            j;
            jcoordx=addnewnodecoord4M4(j,2);
            jcoordy=addnewnodecoord4M4(j,3);
            teM4pdis=sqrt((jcoordx-icoordx)^2+(jcoordy-icoordy)^2);
            if teM4pdis<1e-5
                addnewnodecoord4M4(j,1)=addnewnodecoord4M4(i,1);
            end
        end
    end
    addnewnodecoord4M4;
    addnewnodecoord4M4sort4=sortrows(addnewnodecoord4M4,1);
    teM4pd=[];
    for i=1:length(addnewnodecoord4M4sort4)
        i;
        inode=addnewnodecoord4M4sort4(i,1);
        for j=i+1:length(addnewnodecoord4M4sort4)
            j;
            jnode=addnewnodecoord4M4sort4(j,1);
            teM4pd=jnode-inode;
            if teM4pd<0.5
                addnewnodecoord4M4sort4(j,:)=[0 0 0];
            end
        end
    end

    addnewnodecoord4M4sort4;

    size(addnewnodecoord4M4sort4);
    addnewnodecoord4M4sort4(all(addnewnodecoord4M4sort4==0,2),:)=[];
    size(addnewnodecoord4M4sort4);
    addnewnodecoord4M4sort4;
    addnewnodecoord4M4sort4u=addnewnodecoord4M4sort4;

    coord4;
    length(coord4);

    addnewnodecoord4M4sort4urenode4=[addnewnodecoord4M4sort4u  (length(coord4)+1:length(coord4)+length(addnewnodecoord4M4sort4u))'];

    addnewnodecoord4M4sort4urenode4;


    addcoord4=[addnewnodecoord4M4sort4urenode4(:,end-2:end)];

    addnewnodecoord4;
    addnewnodepoly4={};
    addnewnodecoord4refine4=[];
    for i=nonzeromarK4'
        i;
        addnewnodecoord4refine4=addnewnodecoord4{i};
        for j=1:length(addnewnodecoord4refine4)
            j;
            jcoordx=addnewnodecoord4refine4(j,2);
            jcoordy=addnewnodecoord4refine4(j,3);
            tmepdisten=[];
            for k=1:length(addcoord4)
                kcoordx=addcoord4(k,1);
                kcoordy=addcoord4(k,2);
                tmepdisten=sqrt((kcoordx-jcoordx)^2+(kcoordy-jcoordy)^2);
                if tmepdisten<1e-5

                    addnewnodecoord4refine4(j,1)=addcoord4(k,3);
                    addnewnodepoly4{i}=addnewnodecoord4refine4;

                end
            end
        end
    end


    addnewnodepoly4;
    sdConn4;
    coord4;
    sdsc4;
    for i=nonzeromarK4'

        i;
        [coord4,sdConn4,sdsc4]=addnewpoly(addnewnodepoly4,sdSln4,nsd4,marK4,mat,i,maxmaxnodenumber4,coord4,sdConn4,sdsc4);


    end



    coord4;
    sdConn4
    sdsc4;


    coord4=[coord4;addcoord4(:,1:2)];
 
    for i=1:length(sdConn4)
        i;
        ismember(nonzeromarK4,i);
        if ismember(nonzeromarK4,i)==0
            teM4peleM4ent =sdConn4{i};
            sdnode        =teM4peleM4ent(:,1);
            teM4pcoord4   =coord4(sdnode,:);
            for j= length(teM4peleM4ent):-1:1;
                j;
                teM4pnode       = teM4peleM4ent(j,:);
                teM4pnodecoord4 = [coord4(teM4pnode(1),1)+coord4(teM4pnode(2),1)   coord4(teM4pnode(1),2)+coord4(teM4pnode(2),2) ] /2 ;
                for k= 1:length(coord4)
                    addnode        =k;
                    addnodecoord4   =coord4(k,:);
                    teM4pdisplment  =sqrt((coord4(k,1)-teM4pnodecoord4(1))^2+(coord4(k,2)-teM4pnodecoord4(2))^2);

                    if  teM4pdisplment<1e-5
                        teM4peleM4ent(length(teM4peleM4ent)+1,:)=[0 0];
                        teM4peleM4ent(j+2:end,:)=teM4peleM4ent(j+1:end-1,:);
                        teM4peleM4ent(j,:)      =[teM4pnode(1) k];
                        teM4peleM4ent(j+1,:)    =[k           teM4pnode(2)];
                    end
                end
            end
            teM4peleM4ent;
            sdConn4{i}=teM4peleM4ent;
        end
    end




elseif eta4overline4<=eta4max
    disp('convergence')
end


sdConn5=sdConn4;
sdsc5=sdsc4;
coord5=coord4;

 


figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc5);

PlotSBFEMesh(coord5, sdConn5, opt);
title('MESH5-2');



figure
plot(coord5(:,1),coord5(:,2),'b+')
axis equal
title('point5')
 


toc
disp(['运行时间：',num2str(toc)])

%exportgraphics(gcf, ['e11' '.tiff'] ,'Resolution', 300);







