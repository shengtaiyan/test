%a0esong1step6
% displacement boudary condition 位移边界条件
eps          =    1d-5
dispbc5x5    =    find(abs(coord5(:,2))<eps );
dispbc5y5    =    find(abs(coord5(:,1))<eps );
BC_Dispx5    =    [dispbc5x5     2*  ones(size(dispbc5x5,1),1)  zeros(size(dispbc5x5,1),1)];
BC_Dispy5    =    [dispbc5y5         ones(size(dispbc5y5,1),1)  zeros(size(dispbc5y5,1),1)];
BC_Disp5     =    [BC_Dispx5;  BC_Dispy5];
dispbc5      =    [dispbc5x5;   dispbc5y5];
figure
plot(coord5(dispbc5,1),coord5(dispbc5,2),'b*' );
axis([0,6,0,6]);
title('bc-disp5')
%find point5 for force and pressure
point5number5   =    size(coord5,1);
point5         =    zeros(point5number5,1);
for i=1:point5number5
    if(abs(sqrt(coord5(i,1)^2+coord5(i,2)^2)-1)<eps)
        point5(i)=i;
        i=i+1;
    end
end
point5;
fpoint5=find(point5>eps);
figure
plot(coord5(fpoint5,1),coord5(fpoint5,2),'r+')
title('f-point5')
%find edge5
edge5=[];
sizesdConn5=length(sdConn5)
for i=1:sizesdConn5
    i;
    sdConn5{i};
    sizesdConn53=size(sdConn5{i},1);
    for  j=1:sizesdConn53
        lia=intersect(fpoint5,sdConn5{i}(j,:),'stable');
        if size(lia,1)==2
            edge5=[edge5;sdConn5{i}(j,:)];
        end
    end
end
edge5;
%calculate angle5
angle5=[];
fpoint5;
coord5(fpoint5,:);
x           =coord5(edge5,1);
y           =coord5(edge5,2);
angle5      =atan(y./x);
edge51      =edge5(:,1);
edge52      =edge5(:,2);
x1          =coord5(edge51,1);
y1          =coord5(edge51,2);
angle51     =atan(y1./x1);
x2          =coord5(edge52,1);
y2          =coord5(edge52,2);
angle52     =atan(y2./x2);
trac5       =[cos(angle51)   sin(angle51)   cos(angle52)   sin(angle52)]';
%add force or pressure
ndn5 = 2;
NDof5 = ndn5*size(coord5,1);
F5 = zeros(NDof5,1);
F5 = addSurfTraction(coord5, edge5, trac5, F5);
% solve
%ltx {\bf Solution of S-elements and assemblage of global stiffness}
%ltx {\bf and mass matrices}
[sdSln5, K5, M5] = SBFEMAssembly(coord5, sdConn5, sdsc5, mat);
%ltx {\bf Static solution of nodal displacements and forces}
[d5, F5] = SolverStatics(K5, BC_Disp5, F5);
d5;
d5(121);
d5(122);
d5(141);
d5(142);
%ltx {\bf Plot deformed mesh}
figure
opt = struct('MagnFct', 0.1, 'Undeformed','--k');
PlotDeformedMesh(d5, coord5, sdConn5, opt)
disp('Nodal displacements')
for ii = 1:length(coord5)
    fprintf('%5d %25.15e %25.15d\n',ii, d5(2*ii-1:2*ii))
end
sdSln5;
sdStrnMode5         =  SElementStrainMode2NodeEle( sdSln5 );%ltx integration constants
sdIntgConst5        =  SElementIntgConst( d5, sdSln5 );
%ltx displacements and strains at specified raidal coord5inate
d5;
sdsc5;
size(sdsc5,1);
nsd5=size(sdsc5,1);
for isd=1:nsd5
    xi = 1;    %ltx radial coord5inate
    [nodexy5{isd},dsp5{isd},strnNode5{isd},GPxy5{isd},strnEle5{isd}]=SElementInDispStrain(xi,sdSln5{isd},sdStrnMode5{isd},sdIntgConst5{isd});
end
for isd=1:nsd5
    ximid = 0.5;    %ltx radial coord5inate
    [nodexy5mid{isd},dsp5mid{isd},strnNode5mid{isd},GPxy5mid{isd},strnEle5mid{isd}]=SElementInDispStrain(ximid,sdSln5{isd},sdStrnMode5{isd},sdIntgConst5{isd});
end
for isd=1:nsd5
    xisdsc = 1e-15;    %ltx radial coord5inate
    [nodexy5sdsc{isd},dsp5sdsc{isd},strnNode5sdsc{isd},GPxy5sdsc{isd},strnEle5sdsc{isd}]=SElementInDispStrain(xisdsc,sdSln5{isd},sdStrnMode5{isd},sdIntgConst5{isd});
end
sdscstrain5=[];
for i=1:nsd5
    sdscstrain5(i,:)=strnNode5sdsc{i}(:,1)';
end
sdscstrain5;
%%
% % C = {1, 2, 3, 4, 5};  % 示例 cell 数组
% target = 61        % 要查找的数字
%
% % 示例 cell 数组，每个元素是一个矩阵
% C = sdConn5
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
%%
%   strnNode5{1707}
% Astrain=strnNode5{1707}(:,1)
% % Astress=mat.D * Astrain
% %  strnNode5{1118}
% Bstrain=strnNode5{1118}(:,3)
%
% Bstress=mat.D * Bstrain
magnfct5           = 0.1;
Umax5               = max(abs(d5));%ltx maximum displaceM5ent
Lmax5               = max(max(coord5)-min(coord5));
%ltx maximum dimension of domain
fct5                = magnfct5*Lmax5/Umax5;
%ltx factor to magnify the displaceM5ent ltx augment nodal coord5inates
dispmentforplot5    = fct5*(reshape(d5,2,[]))';
% stress recovery
[recoverystress5]  =spr(sdSln5,strnEle5,strnNode5,mat,nsd5,coord5,d5,GPxy5,strnEle5mid);
[sdscstressre5]=sprsdscstress(recoverystress5,sdSln5,coord5,nsd5,mat,strnNode5sdsc);
recoverystress5;
% error estimation use recovery stress
elementnumber5=nsd5;
[eta5,matupoly5,absepoly5,sdSln5]=errorestimation(sdIntgConst5,sdStrnMode5,strnNode5,strnEle5,elementnumber5,recoverystress5,mat,sdSln5,sdscstressre5,strnEle5mid,sdscstrain5);
% refine stragety 鎵惧埌闇€瑕佺粏鍖栫殑polygon锛岀劧鍚庡垽鏂槸鍑犺竟褰紝鐒跺悗缁嗗寲锛屽湪璁＄畻锛屽湪鍒ゆ柇锛岀洿鍒拌揪鍒扮粏鍖栬姹�
uall5=sum(matupoly5);
eall5=sum(absepoly5);
% reen5=sqrt(eall5)/sqrt(uall5)
%reen5=sqrt(0.5*eall5)/sqrt(0.5*uall5)
reen5=eall5/uall5
etamax5=max(absepoly5)

% %%不做细化
% % eta5overline4                =sqrt(eall5^2/(uall5^2+eall5^2) );
% % eta5max=0.1
% % if eta5overline4 > eta5max%鍏ㄥ眬鍒ゆ柇
% %     eM5            =eta5max*sqrt((uall5^2+eall5^2)/nsd5)     ;     %鍗曞厓骞冲潎鍏佽璇樊
% %     ksi4           =cell(size(strnNode5,2),1);
% %     for i=1:nsd5
% %         ksi4{i}    = absepoly5(i)/ eM5;
% %     end
% %     marK5=[];
% %     for i=1:nsd5
% %         if ksi4{i} > 1
% %             marK5(i)=i;
% %         end
% %     end
% marK5=[]
% for i=1:nsd5
%     if absepoly5(i) > etamax5*gammy
%         marK5(i)=i;
%     end
% end
% marK5
% size(marK5);
% nonzeromarK5=nonzeros(marK5);
% maxnodenumber4={};
% for i=nonzeromarK5'
%     i;
%     disp(i);
%     maxnodenumber4{i}=max(sdSln5{i}.node);
% end
% 
% maxnodenumber4;
% maxmaxnodenumber4= max(cell2mat(maxnodenumber4));
% eleM5entsize4=[];
% for i=nonzeromarK5'
%     i;
%     eleM5entsize4(i)=length(sdSln5{i}.node);
% end
% eleM5entsize4;
% % eleM5entsize4(eleM5entsize4==0)=[];
% length(eleM5entsize4);
% addnodenumberstar4=[];
% addnodenumberstar4(1)=maxmaxnodenumber4+1;
% for j=2:length(nonzeromarK5)
%     j;
%     addnodenumberstar4(j)=addnodenumberstar4(j-1)+2*eleM5entsize4(nonzeromarK5(j-1));
% end
% addnodenumberstar4;
% marK5;
% % x=[2 8 4 9 3 ]
% % y=[0 0 1 0 0 0 0 1 1 0 0 1 0 0 1]
% % y(logical(y))=x
% addnodenumberstar4y(logical(marK5))=addnodenumberstar4;
% addnewnodecoord5={};
% nonzeromarK5
% for i=nonzeromarK5'
%     [addnewnodecoord5{i}]=addnodecoord(sdSln5,nsd5,marK5,mat,i,coord5,sdConn5,sdsc5,addnodenumberstar4y);
% end
% addnewnodecoord5;
% addnewnodecoord5M5=[];
% for i=nonzeromarK5'
%     i;
%     AteM5p=[];
%     AteM5p=addnewnodecoord5{i};
%     addnewnodecoord5M5=[addnewnodecoord5M5;AteM5p];
% end
% addnewnodecoord5M5;
% coord5;
% teM5pdis=[];
% for i=1:length(addnewnodecoord5M5)
%     i;
%     icoordx=addnewnodecoord5M5(i,2);
%     icoordy=addnewnodecoord5M5(i,3);
%     for j=i+1:length(addnewnodecoord5M5)
%         j;
%         jcoordx=addnewnodecoord5M5(j,2);
%         jcoordy=addnewnodecoord5M5(j,3);
%         teM5pdis=sqrt((jcoordx-icoordx)^2+(jcoordy-icoordy)^2);
%         if teM5pdis<1e-5
%             addnewnodecoord5M5(j,1)=addnewnodecoord5M5(i,1);
%         end
%     end
% end
% addnewnodecoord5M5;
% addnewnodecoord5M5sort4=sortrows(addnewnodecoord5M5,1);
% teM5pd=[];
% for i=1:length(addnewnodecoord5M5sort4)
%     i;
%     inode=addnewnodecoord5M5sort4(i,1);
%     for j=i+1:length(addnewnodecoord5M5sort4)
%         j;
%         jnode=addnewnodecoord5M5sort4(j,1);
%         teM5pd=jnode-inode;
%         if teM5pd<0.5
%             addnewnodecoord5M5sort4(j,:)=[0 0 0];
%         end
%     end
% end
% addnewnodecoord5M5sort4;
% size(addnewnodecoord5M5sort4);
% addnewnodecoord5M5sort4(all(addnewnodecoord5M5sort4==0,2),:)=[];
% size(addnewnodecoord5M5sort4);
% addnewnodecoord5M5sort4;
% addnewnodecoord5M5sort4u=addnewnodecoord5M5sort4;
% coord5;
% length(coord5);
% addnewnodecoord5M5sort4urenode4=[addnewnodecoord5M5sort4u  (length(coord5)+1:length(coord5)+length(addnewnodecoord5M5sort4u))'];
% addnewnodecoord5M5sort4urenode4;
% addcoord5=[addnewnodecoord5M5sort4urenode4(:,end-2:end)];
% addnewnodecoord5;
% addnewnodepoly4={};
% addnewnodecoord5refine4=[];
% for i=nonzeromarK5'
%     i;
%     addnewnodecoord5refine4=addnewnodecoord5{i};
%     for j=1:length(addnewnodecoord5refine4)
%         j;
%         jcoordx=addnewnodecoord5refine4(j,2);
%         jcoordy=addnewnodecoord5refine4(j,3);
%         tmepdisten=[];
%         for k=1:length(addcoord5)
%             kcoordx=addcoord5(k,1);
%             kcoordy=addcoord5(k,2);
%             tmepdisten=sqrt((kcoordx-jcoordx)^2+(kcoordy-jcoordy)^2);
%             if tmepdisten<1e-5
%                 addnewnodecoord5refine4(j,1)=addcoord5(k,3);
%                 addnewnodepoly4{i}=addnewnodecoord5refine4;
% 
%             end
%         end
%     end
% end
% addnewnodepoly4;
% sdConn5;
% coord5;
% sdsc5;
% for i=nonzeromarK5'
%     i;
%     [coord5,sdConn5,sdsc5]=addnewpoly(addnewnodepoly4,sdSln5,nsd5,marK5,mat,i,maxmaxnodenumber4,coord5,sdConn5,sdsc5);
% end
% coord5;
% sdConn5
% sdsc5;
% coord5=[coord5;addcoord5(:,1:2)];
% %%鎵惧埌娌℃湁缁嗗寲鐨勫崟鍏冿紝缁欎粬缁嗗寲浜嗙殑杈逛笂澧炲姞涓€涓妭鐐�
% %%娉ㄦ剰 鏈夌殑鍗曞厓浼氭湁澶氫釜杈硅缁嗗寲锛岄渶瑕佸鍔犲嚑涓偣
% for i=1:length(sdConn5)
%     i;
%     ismember(nonzeromarK5,i);
%     if ismember(nonzeromarK5,i)==0
%         teM5peleM5ent =sdConn5{i};
%         sdnode        =teM5peleM5ent(:,1);
%         teM5pcoord5   =coord5(sdnode,:);
%         for j= length(teM5peleM5ent):-1:1;
%             j;
%             teM5pnode       = teM5peleM5ent(j,:);
%             teM5pnodecoord5 = [coord5(teM5pnode(1),1)+coord5(teM5pnode(2),1)   coord5(teM5pnode(1),2)+coord5(teM5pnode(2),2) ] /2 ;
%             for k= 1:length(coord5)
%                 addnode        =k;
%                 addnodecoord5   =coord5(k,:);
%                 teM5pdisplment  =sqrt((coord5(k,1)-teM5pnodecoord5(1))^2+(coord5(k,2)-teM5pnodecoord5(2))^2);
% 
%                 if  teM5pdisplment<1e-5
%                     teM5peleM5ent(length(teM5peleM5ent)+1,:)=[0 0];
%                     teM5peleM5ent(j+2:end,:)=teM5peleM5ent(j+1:end-1,:);
%                     teM5peleM5ent(j,:)      =[teM5pnode(1) k];
%                     teM5peleM5ent(j+1,:)    =[k           teM5pnode(2)];
%                 end
%             end
%         end
%         teM5peleM5ent;
%         sdConn5{i}=teM5peleM5ent;
%     end
% end
% % elseif eta5overline4<=eta5max
% %     disp('convergence')
% % end
% sdConn5=sdConn5;
% sdsc5=sdsc5;
% coord5=coord5;
% % figure
% % opt        = struct('LineSpec','-k', 'sdSC',sdsc5);
% % for i=1:size(sdConn5,1)
% %     sdsc5(i,:)
% %     sdsc5(i,1)
% %     sdsc5(i,2)
% %     str = i
% %     text(sdsc5(i,1),sdsc5(i,2),num2str(str))
% % end
% % PlotSBFEMesh(coord5, sdConn5, opt);
% % title('MESH5');
% figure
% opt        = struct('LineSpec','-k', 'sdSC',sdsc5);
% PlotSBFEMesh(coord5, sdConn5, opt);
% title('MESH5');
% % figure
% % plot(coord5(:,1),coord5(:,2),'b+')
% % axis equal
% % title('point5')



