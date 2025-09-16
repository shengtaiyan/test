%a0esong1step3
% displacement boudary condition 位移边界条件
eps          =    1d-5
dispbc2x2    =    find(abs(coord2(:,2))<eps );
dispbc2y2    =    find(abs(coord2(:,1))<eps );
BC_Dispx2    =    [dispbc2x2     2*  ones(size(dispbc2x2,1),1)  zeros(size(dispbc2x2,1),1)];
BC_Dispy2    =    [dispbc2y2         ones(size(dispbc2y2,1),1)  zeros(size(dispbc2y2,1),1)];
BC_Disp2     =    [BC_Dispx2;  BC_Dispy2];
dispbc2      =    [dispbc2x2;   dispbc2y2];
figure
plot(coord2(dispbc2,1),coord2(dispbc2,2),'b*' )
axis([0,6,0,6])
title('bc-disp2')
%find point2 for force and pressure
point2number2   =    size(coord2,1);
point2         =    zeros(point2number2,1);
for i=1:point2number2
    if(abs(sqrt(coord2(i,1)^2+coord2(i,2)^2)-1)<eps)
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
x          =coord2(edge2,1);
y          =coord2(edge2,2);
angle2     =atan(y./x);
edge21     =edge2(:,1);
edge22     =edge2(:,2);
x1         =coord2(edge21,1);
y1         =coord2(edge21,2);
angle21    =atan(y1./x1);
x2         =coord2(edge22,1);
y2         =coord2(edge22,2);
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
d2(121);
d2(122);
d2(141);
d2(142);
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
    xisdsc = 1e-15;    %ltx radial coord2inate
    [nodexy2sdsc{isd},dsp2sdsc{isd},strnNode2sdsc{isd},GPxy2sdsc{isd},strnEle2sdsc{isd}]=SElementInDispStrain(xisdsc,sdSln2{isd},sdStrnMode2{isd},sdIntgConst2{isd});
end
sdscstrain2=[];
for i=1:nsd2
    sdscstrain2(i,:)=strnNode2sdsc{i}(:,1)';
end
sdscstrain2;
% % C = {1, 2, 3, 4, 5};  % 示例 cell 数组
% target = 71         % 要查找的数字
%
% % 示例 cell 数组，每个元素是一个矩阵
% C = sdConn2
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
% strnNode2{523}
% Astrain=strnNode2{523}(:,1)
% Astress=mat.D * Astrain
% strnNode2{1118}
% Bstrain=strnNode2{1118}(:,3)
%
% Bstress=mat.D * Bstrain
magnfct22           = 0.1;
Umax2               = max(abs(d2));%ltx maximum displaceM2ent
Lmax2               = max(max(coord2)-min(coord2));
%ltx maximum dimension of domain
fct2                = magnfct22*Lmax2/Umax2;
%ltx factor to magnify the displaceM2ent ltx augment nodal coord2inates
dispmentforplot2    = fct2*(reshape(d2,2,[]))';
% [recoverystress2]  =spr(sdSln2,strnEle2,strnNode2,mat,nsd2,coord2,d2,GPxy2)
[recoverystress2]  =spr(sdSln2,strnEle2,strnNode2,mat,nsd2,coord2,d2,GPxy2,strnEle2mid);
[sdscstressre2]=sprsdscstress(recoverystress2,sdSln2,coord2,nsd2,mat,strnNode2sdsc);
% [recoverystress2,sdSln2]      =sbfemSpr(strnEle2,strnNode2,sdSln2,mat,nsd2,coord2,d2)
% error estimation use recovery stress
elementnumber2=nsd2;
[eta2,matupoly2,absepoly2,sdSln2]=errorestimation(sdIntgConst2,sdStrnMode2,strnNode2,strnEle2,elementnumber2,recoverystress2,mat,sdSln2,sdscstressre2,strnEle2mid,sdscstrain2);
uall2=sum(matupoly2);
eall2=sum(absepoly2);
% reen2=sqrt(eall2)/sqrt(uall2)
% reen2=sqrt(0.5*eall2)/sqrt(0.5*uall2)
reen2=eall2/uall2
etamax2=max(absepoly2)
% % refine stragety 找到需要细化的polygon，然后判断是几边形，然后细化，在计算，在判断，直到达到细化要求
% eta2overline2                =sqrt(eall2^2/(uall2^2+eall2^2) );
% eta2max=0.1
% if eta2overline2 > eta2max%全局判断
%     eM2            =eta2max*sqrt((uall2^2+eall2^2)/51)     ;     %单元平均允许误差
%     ksi2           =cell(size(strnNode2,2),1);
%     for i=1:nsd2
%         ksi2{i}    = absepoly2(i)/ eM2;
%     end
%     marK2=[];
%     for i=1:nsd2
%         if ksi2{i} > 1
%             marK2(i)=i;
%         end
%     end
marK2=[]
for i=1:nsd2
    if absepoly2(i) > etamax2*gammy
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
eleM2entsize2;
% eleM2entsize2(eleM2entsize2==0)=[]
length(eleM2entsize2);
addnodenumberstar2=[];
addnodenumberstar2(1)=maxmaxnodenumber2+1;
for j=2:length(nonzeromarK2)
    j;
    addnodenumberstar2(j)=addnodenumberstar2(j-1)+2*eleM2entsize2(nonzeromarK2(j-1));
end

addnodenumberstar2;
marK2;
% x=[2 8 4 9 3 ]
% y=[0 0 1 0 0 0 0 1 1 0 0 1 0 0 1]
% y(logical(y))=x

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
coord2;
sdConn2;
sdsc2;
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

% elseif eta2overline2<=eta2max
%     disp('convergence')
% end


sdConn3=sdConn2;
sdsc3=sdsc2;
coord3=coord2;

% figure
% opt        = struct('LineSpec','-k', 'sdSC',sdsc2);
% for i=1:size(sdConn3,1)
%     sdsc3(i,:)
%     sdsc3(i,1)
%     sdsc3(i,2)
%     str = i
%     text(sdsc3(i,1),sdsc3(i,2),num2str(str))
% end
% PlotSBFEMesh(coord3, sdConn3, opt);
% title('MESH3');


cell2mat(sdConn3);
figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc3);
PlotSBFEMesh(coord3, sdConn3, opt);
title('MESH3');



figure
plot(coord3(:,1),coord3(:,2),'b+')
axis equal
title('point3')



