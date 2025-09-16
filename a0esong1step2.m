%a0esong1step2
% displacement boudary condition 位移边界条件
eps          =    1d-5;
dispbc1x1    =    find(abs(coord1(:,2))<eps );
dispbc1y1    =    find(abs(coord1(:,1))<eps );
BC_Dispx1    =    [dispbc1x1     2*  ones(size(dispbc1x1,1),1)  zeros(size(dispbc1x1,1),1)];
BC_Dispy1    =    [dispbc1y1         ones(size(dispbc1y1,1),1)  zeros(size(dispbc1y1,1),1)];
BC_Disp1     =    [BC_Dispx1;  BC_Dispy1];
dispbc1      =    [dispbc1x1;   dispbc1y1];
figure
plot(coord1(dispbc1,1),coord1(dispbc1,2),'b*' )
axis([0,6,0,6])
title('bc-disp1')
%find point1 for force and pressure
point1number1   =    size(coord1,1);
point1         =    zeros(point1number1,1);
for i=1:point1number1
    if(abs(sqrt(coord1(i,1)^2+coord1(i,2)^2)-1)<eps)
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
x          =coord1(edge1,1);
y          =coord1(edge1,2);
angle1     =atan(y./x);
edge11     =edge1(:,1);
edge12     =edge1(:,2);
x1         =coord1(edge11,1);
y1         =coord1(edge11,2);
angle11    =atan(y1./x1);
x2         =coord1(edge12,1);
y2         =coord1(edge12,2);
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
d1(121);
d1(122);
d1(141);
d1(142);
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
    xisdsc = 1e-15;    %ltx radial coord1inate
    [nodexy1sdsc{isd},dsp1sdsc{isd},strnNode1sdsc{isd},GPxy1sdsc{isd},strnEle1sdsc{isd}]=SElementInDispStrain(xisdsc,sdSln1{isd},sdStrnMode1{isd},sdIntgConst1{isd});
end
sdscstrain1=[];
for i=1:nsd1
    sdscstrain1(i,:)=strnNode1sdsc{i}(:,1)';
end
% % C = {1, 2, 3, 4, 5};  % 示例 cell 数组
% target = 61         % 要查找的数字
%
% % 示例 cell 数组，每个元素是一个矩阵
% C = sdConn1
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
%   strnNode1{67}
% Astrain=strnNode1{67}(:,1)
%
% Astress=mat.D * Astrain
%  strnNode1{186}
% Bstrain=strnNode1{186}(:,1)
%
% Bstress=mat.D * Bstrain
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
% reen1=sqrt(eall1)/sqrt(uall1)
% reen1=sqrt(0.5*eall1)/sqrt(0.5*uall1)
reen1=eall1/uall1
etamax1=max(absepoly1)
% eta1overline1                =sqrt(eall1^2/(uall1^2+eall1^2) );
% eta1max=0.1
% if eta1overline1 > eta1max%全局判断
%     em1            =eta1max*sqrt((uall1^2+eall1^2)/28)     ;     %单元平均允许误差
%     ksi1           =cell(size(strnNode1,2),1);
%     for i=1:nsd1
%         ksi1{i}    = absepoly1(i)/ em1;
%     end
%     mark1=[];
%     for i=1:nsd1
%         if ksi1{i} > 1
%             mark1(i)=i;
%         end
%     end

mark1=[]
for i=1:nsd1
    if absepoly1(i) > etamax1*gammy
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
% elem1entsize1(elem1entsize1==0)=[];
length(elem1entsize1);
addnodenumberstar1=[];
addnodenumberstar1(1)=maxmaxnodenumber1+1;
for j=2:length(nonzeromark1)
    j;
    addnodenumberstar1(j)=addnodenumberstar1(j-1)+2*elem1entsize1(nonzeromark1(j-1));
end
addnodenumberstar1;
mark1;
% x=[2 8 4 9 3 ]
% y=[0 0 1 0 0 0 0 1 1 0 0 1 0 0 1]
% y(logical(y))=x
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
coord1;
sdConn1
sdsc1 ;
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
% elseif eta1overline1<=eta1max
%     disp('convergence')
% end
sdConn2=sdConn1;
sdsc2=sdsc1;
coord2=coord1;
% figure
% opt        = struct('LineSpec','-k', 'sdSC',sdsc2);
% for i=1:size(sdConn2,1)
%     sdsc2(i,:);
%     sdsc2(i,1);
%     sdsc2(i,2);
%     str = i;
%     text(sdsc2(i,1),sdsc2(i,2),num2str(str));
% end
% PlotSBFEMesh(coord2, sdConn2, opt);
% title('MESH2');

figure
opt        = struct('LineSpec','-k', 'sdSC',sdsc2);
PlotSBFEMesh(coord2, sdConn2, opt);
title('MESH2');

figure
plot(coord2(:,1),coord2(:,2),'b+')
axis equal
title('point2')

% figure
% opt        = struct('LineSpec','-k', 'sdsc2',sdsc2);
% PlotSBFEMesh(coord2, sdConn2, opt);
% for i=1:length(coord2)
%     str = i;
%     text(coord2(i,1)+0.1,coord2(i,2)+0.1,num2str(str));
% end
% axis([0,5,0,5])
% title('point2');


