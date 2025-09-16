%a0esong1step4
% displacement boudary condition 位移边界条件
eps          =    1d-5
dispbc3x3    =    find(abs(coord3(:,2))<eps );
dispbc3y3    =    find(abs(coord3(:,1))<eps );
BC_Dispx3    =    [dispbc3x3     2*  ones(size(dispbc3x3,1),1)  zeros(size(dispbc3x3,1),1)];
BC_Dispy3    =    [dispbc3y3         ones(size(dispbc3y3,1),1)  zeros(size(dispbc3y3,1),1)];
BC_Disp3     =    [BC_Dispx3;  BC_Dispy3];
dispbc3      =    [dispbc3x3;   dispbc3y3];
figure
plot(coord3(dispbc3,1),coord3(dispbc3,2),'b*' )
axis([0,6,0,6])
title('bc-disp3')
%find point3 for force and pressure
point3number3   =    size(coord3,1);
point3         =    zeros(point3number3,1);
for i=1:point3number3
    if(abs(sqrt(coord3(i,1)^2+coord3(i,2)^2)-1)<eps)
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
x          =coord3(edge3,1);
y          =coord3(edge3,2);
angle3     =atan(y./x);
edge31     =edge3(:,1);
edge32     =edge3(:,2);
x1         =coord3(edge31,1);
y1         =coord3(edge31,2);
angle31    =atan(y1./x1);
x2         =coord3(edge32,1);
y2         =coord3(edge32,2);
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
d3(121);
d3(122);
d3(141);
d3(142);
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
    xisdsc = 1e-15;    %ltx radial coord3inate
    [nodexy3sdsc{isd},dsp3sdsc{isd},strnNode3sdsc{isd},GPxy3sdsc{isd},strnEle3sdsc{isd}]=SElementInDispStrain(xisdsc,sdSln3{isd},sdStrnMode3{isd},sdIntgConst3{isd});
end
sdscstrain3=[];
for i=1:nsd3
    sdscstrain3(i,:)=strnNode3sdsc{i}(:,1)';
end
sdscstrain3;

% % C = {1, 2, 3, 4, 5};  % 示例 cell 数组
% target = 71        % 要查找的数字
%
% % 示例 cell 数组，每个元素是一个矩阵
% C = sdConn3
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

%
%   strnNode3{1707}
% Astrain=strnNode3{1707}(:,1)
%
% Astress=mat.D * Astrain
%  strnNode3{1118}
% Bstrain=strnNode3{1118}(:,3)
%
% Bstress=mat.D * Bstrain
magnfct33            = 0.1;
Umax3               = max(abs(d3));%ltx maximum displaceM3ent
Lmax3               = max(max(coord3)-min(coord3));
%ltx maximum dimension of domain
fct3                = magnfct33*Lmax3/Umax3;
%ltx factor to magnify the displaceM3ent ltx augment nodal coord3inates
dispmentforplot3    = fct3*(reshape(d3,2,[]))';
[recoverystress3]  =spr(sdSln3,strnEle3,strnNode3,mat,nsd3,coord3,d3,GPxy3,strnEle3mid);
[sdscstressre3]=sprsdscstress(recoverystress3,sdSln3,coord3,nsd3,mat,strnNode3sdsc);
% [recoverystress3,sdSln3]      =sbfemSpr(strnEle3,strnNode3,sdSln3,mat,nsd3,coord3,d3)
recoverystress3;
% error estimation use recovery stress
elementnumber3=nsd3;
[eta3,matupoly3,absepoly3,sdSln3]=errorestimation(sdIntgConst3,sdStrnMode3,strnNode3,strnEle3,elementnumber3,recoverystress3,mat,sdSln3,sdscstressre3,strnEle3mid,sdscstrain3);
% refine stragety 找到需要细化的polygon，然后判断是几边形，然后细化，在计算，在判断，直到达到细化要求
uall3=sum(matupoly3);
eall3=sum(absepoly3);
% reen3=sqrt(eall3)/sqrt(uall3)
% reen3=sqrt(0.5*eall3)/sqrt(0.5*uall3)
reen3=eall3/uall3
etamax3=max(absepoly3);
% eta3overline3                =sqrt(eall3^2/(uall3^2+eall3^2) );
% eta3max=0.1
% if eta3overline3 > eta3max%全局判断
%     eM3            =eta3max*sqrt((uall3^2+eall3^2)/nsd3)     ;     %单元平均允许误差
%     ksi3           =cell(size(strnNode3,2),1);
%     for i=1:nsd3
%         ksi3{i}    = absepoly3(i)/ eM3;
%     end
%     marK3=[];
%     for i=1:nsd3
%         if ksi3{i} > 1
%             marK3(i)=i;
%         end
%     end
marK3=[]
for i=1:nsd3
    if absepoly3(i) > etamax3*gammy
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
marK3;
% x=[2 8 4 9 3 ]
% y=[0 0 1 0 0 0 0 1 1 0 0 1 0 0 1]
% y(logical(y))=x

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
addnewnodecoord3M3;
addnewnodecoord3M3sort3=sortrows(addnewnodecoord3M3,1);
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
coord3;
sdConn3
sdsc3;
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

% elseif eta3overline3<=eta3max
%     disp('convergence')
% end
sdConn4=sdConn3;
sdsc4=sdsc3;
coord4=coord3;

% figure
% opt        = struct('LineSpec','-k', 'sdSC',sdsc4);
% for i=1:size(sdConn4,1)
%     sdsc4(i,:);
%     sdsc4(i,1);
%     sdsc4(i,2);
%     str = i;
%     text(sdsc4(i,1),sdsc4(i,2),num2str(str));
% end
% PlotSBFEMesh(coord4, sdConn4, opt);
% title('MESH4');


figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc4);
PlotSBFEMesh(coord4, sdConn4, opt);
title('MESH4');
% figure
% plot(coord4(:,1),coord4(:,2),'b+')
% axis equal
% title('point4')



