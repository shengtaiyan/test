%a0esong1step1
% displacement boudary condition 位移边界条件
eps        =    1d-5;
dispbcx    =    find(abs(coord0(:,2))<eps );
dispbcy    =    find(abs(coord0(:,1))<eps );
BC_Dispx   =    [dispbcx     2*  ones(size(dispbcx,1),1)  zeros(size(dispbcx,1),1)];
BC_Dispy   =    [dispbcy         ones(size(dispbcy,1),1)  zeros(size(dispbcy,1),1)];
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
    if(abs(sqrt(coord0(ipointnumber,1)^2+coord0(ipointnumber,2)^2)-1)<eps)
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
x         =coord0(edge,1)
y         =coord0(edge,2)
angle     =atan(y./x)
edge1     =edge(:,1)
edge2     =edge(:,2)
x1        =coord0(edge1,1)
y1        =coord0(edge1,2)
angle1    =atand(y1./x1)
x2        =coord0(edge2,1)
y2        =coord0(edge2,2)
angle2    =atand(y2./x2)


% % 计算圆心(假设圆心在原点)
% center = [0, 0]; 
% 
% % 计算法向量(指向圆心方向)
% norms = zeros(length(edge), 2)
% for i = 1:length(edge)
%     midPoint = mean(coord0(edge(i,:),:));
%     norms(i,:) = (center - midPoint)/norm(center - midPoint)
% end
% 
% p = -1.0; % 压力值
% trac = zeros(4, length(edge));
% for i = 1:length(edge)
%     trac(:,i) = [norms(i,1); norms(i,2); norms(i,1); norms(i,2)] * p
% end

% figure
% PlotEdges(coord0, edge); % 可视化确认荷载方向
% quiver(coord0(edge(:,1),1), coord0(edge(:,1),2), norms(:,1), norms(:,2))
% 

 
 
trac = [-1, 0, -0.9239, -0.3827;    % 单元1（节点1-2）
        -0.9239, -0.3827, -0.7071, -0.7071;    % 单元2（节点2-3）
        -0.7071, -0.7071, -0.3827, -0.9239;    % 单元3（节点3-4）
        -0.3827, -0.9239, 0, -1]





% trac=[0.5*ledge(1)*sin(edgealpha(1))   0.5*ledge(1)*cos(edgealpha(1)) 0.5*ledge(1)*sin(edgealpha(1)) 0.5*ledge(1)*cos(edgealpha(1)) ;
%       0.5*ledge(2)*sin(edgealpha(2))   0.5*ledge(2)*cos(edgealpha(2)) 0.5*ledge(2)*sin(edgealpha(2)) 0.5*ledge(2)*cos(edgealpha(2)) ;
%       0.5*ledge(3)*sin(edgealpha(3))   0.5*ledge(3)*cos(edgealpha(3)) 0.5*ledge(3)*sin(edgealpha(3)) 0.5*ledge(3)*cos(edgealpha(3)) ;
%       ]'
%add force or pressure
ndn = 2;
NDof = ndn*size(coord0,1);
F0 = zeros(NDof,1);
F0 = addSurfTraction(coord0, edge, trac, F0)


size(F0)



size(F0)
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
    xi = 1;
    %ltx radial coord0inate
    [nodexy{isd},dsp{isd},strnNode{isd},GPxy{isd},strnEle{isd}]=SElementInDispStrain(xi,sdSln{isd},sdStrnMode{isd},sdIntgConst{isd});
end
for isd=1:nsd
    ximid = 0.5;
    %ltx radial coord0inate
    [nodexymid{isd},dspmid{isd},strnNodemid{isd},GPxymid{isd},strnElemid{isd}]=SElementInDispStrain(ximid,sdSln{isd},sdStrnMode{isd},sdIntgConst{isd});
end
strnElemid;
for isd=1:nsd
    xisdsc =  1e-15;
    %ltx radial coord0inate
    [nodexysdsc{isd},dspsdsc{isd},strnNodesdsc{isd},GPxysdsc{isd},strnElesdsc{isd}]=SElementInDispStrain(xisdsc,sdSln{isd},sdStrnMode{isd},sdIntgConst{isd});
end
strnElesdsc;
strnNodesdsc;

sdscstrain=[];
for i=1:nsd
    sdscstrain(i,:)=strnNodesdsc{i}(:,1)';
end

sdscstrain;
% % C = {1, 2, 3, 4, 5};  % 示例 cell 数组
% target = 71          % 要查找的数字
%
% % 示例 cell 数组，每个元素是一个矩阵
% C = sdConn0
%
% % 初始化变量
% found = false;
% matrix_index = 0;
% row_index = 0;
% col_index = 0;
% % 遍历 cell 数组中的每个矩阵
% for i = 1:length(C)
%     % 检查当前矩阵中是否存在目标数字
%     [row, col] = find(C{i} == target);
%
%     % 如果找到目标数字
%     if ~isempty(row)
%         found = true;
%         matrix_index = i;
%         row_index = row(1);  % 取第一个匹配的位置
%         col_index = col(1);  % 取第一个匹配的位置
%         break;
%     end
% end
% % 输出结果
% if found
%     fprintf('数字 %d 找到了！\n', target);
%     fprintf('它位于 cell 数组的第 %d 个矩阵中，位置是 (%d, %d)。\n', ...
%         matrix_index, row_index, col_index);
% else
%     fprintf('数字 %d 没有找到。\n', target);
% end

% strnNode{8}
% Astrain=strnNode{8}(:,end)
% Astress=mat.D * Astrain
% strnNode{28}
% Bstrain=strnNode{28}(:,end)
% Bstress=mat.D * Bstrain

magnFct            = 0.1;
Umax               = max(abs(d0));%ltx maximum displacement
Lmax               = max(max(coord0)-min(coord0));
%ltx maximum dimension of domain
fct                = magnFct*Lmax/Umax;
%ltx factor to magnify the displacement ltx augment nodal coord0inates
dispmentforplot    = fct*(reshape(d0,2,[]))';
strnEle{isd};
%
% r=1
% b=5
% a=1
% p=1
% sigmarA=p*(1-b^2/r^2)/((b/a)^2-1)
%
% sigmatheA=p*(1+b^2/r^2)/((b/a)^2-1)
%
% nu=0.3
% E=1
% ur=0
% ur=p*a*a*r*(1-nu+b^2*(1+nu)/r^2)/(E*(b^2-a^2))
%
% r=5
% b=5
% a=1
% p=1
% sigmarB=p*(1-b^2/r^2)/((b/a)^2-1)
%
% sigmatheB=p*(1+b^2/r^2)/((b/a)^2-1)
%
% nu=0.3
% E=1
% ur=0
% ur=p*a*a*r*(1-nu+b^2*(1+nu)/r^2)/(E*(b^2-a^2))
figure
plot(coord0(:,1),coord0(:,2),'b+')
axis equal
title('point1')
% stress recovery

[recoverystress]=spr(sdSln,strnEle,strnNode,mat,nsd,coord0,d0,GPxy,strnElemid);
recoverystress;
size(recoverystress);
[sdscstressre]=sprsdscstress(recoverystress,sdSln,coord0,nsd,mat,strnNodesdsc);
% [gspointtrisector]=sprgspointtrisector(recoverystress,sdSln,coord0,nsd,mat,strnElemid,sdscstress)
% gspointtrisector

%error estimation use recovery stress

% [eta,u,e,sdSln]=errorestimation(sdIntgConst,sdStrnMode,strnNodea,strnElea,nsd,recoverystressa,mat,sdSln)
elementnumber=nsd;
% [eta,u,e,sdSln]=errorestimation(sdIntgConst,sdStrnMode,strnNodea,strnElea,elementnumber,recoverystressa,mat,sdSln)

[eta,matupoly,absepoly,sdSln]=errorestimation(sdIntgConst,sdStrnMode,strnNode,strnEle,elementnumber,recoverystress,mat,sdSln,sdscstressre,strnElemid,sdscstrain);

% refine stragety 找到需要细化的polygon，然后判断是几边形，然后细化，在计算，在判断，直到达到细化要求
uall=sum(matupoly);
eall=sum(absepoly);

% reen=sqrt(eall)/sqrt(uall)
reen=eall/uall
% reen=sqrt(0.5*eall)/sqrt(0.5*uall)

gammy=0.6
etamax=max(absepoly)

% etaoverline       =sqrt(eall^2/(uall^2+eall^2) );
% etamax            =0.1
% if etaoverline > etamax%全局判断
%     em            =etamax*sqrt((uall^2+eall^2)/nsd); %单元平均允许误差
%     ksi           =cell(size(strnEle,2),1);
%     for i=1:nsd
%         ksi{i}    = absepoly(i)/em;
%     end
%     mark=[];
%     for i=1:nsd
%         if ksi{i} > 1
%             mark(i)=i;
%         end
%     end
mark=[]
for i=1:nsd
    if absepoly(i) >= etamax*gammy
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
coord0;
sdConn0
sdsc0;
% for j=1:length(addnewnodepoly)
%     tempcoordadd= addnewnodepoly{j}
%     tempcoordaddcoord=tempcoordadd(:,2:end)
%     coord0=[coord0;tempcoordaddcoord]
% end
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

% elseif etaoverline<=etamax
%     disp('convergence')
% end
sdConn1=sdConn0;
sdsc1=sdsc0;
coord1=coord0;
% figure
% opt        = struct('LineSpec','-k', 'sdSC',sdsc1);
% for i=1:size(sdConn1,1)
%     sdsc1(i,:);
%     sdsc1(i,1);
%     sdsc1(i,2);
%     str = i;
%     text(sdsc1(i,1),sdsc1(i,2),num2str(str));
% end
% PlotSBFEMesh(coord1, sdConn1, opt);
% title('MESH1');

figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc1);
PlotSBFEMesh(coord1, sdConn1, opt);
title('MESH1');

figure
plot(coord1(:,1),coord1(:,2),'b+')
axis equal
title('point1b+')

% figure
% opt        = struct('LineSpec','-k', 'sdsc1',sdsc1);
% PlotSBFEMesh(coord1, sdConn1, opt);
% for i=1:length(coord1)
%     str = i
%     text(coord1(i,1)+0.1,coord1(i,2)+0.1,num2str(str))
% end
% axis([0,5,0,5])
% title('point1');