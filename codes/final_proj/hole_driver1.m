clear
clc
run('meshhh.m')
%用户输入参数，可在gmsh里改边长，孔半径，网格密度
v = 0.3;
E = 1e9;
L = 2;
T_x = 1e4;%可改
R = 0.5;
%精确解换成笛卡尔坐标下表达
r = @(x,y) sqrt((x+1)^2 + (y+1)^2);
th = @(x,y) atan2((y+1),(x+1));

sigma_rr = @(x,y) T_x/2*(1-R^2/r(x,y)^2) + T_x/2*(1-4*R^2/r(x,y)^2+3*R^4/r(x,y)^4)*cos(2*th(x,y));
sigma_thth = @(x,y) T_x/2*(1+R^2/r(x,y)^2) - T_x/2*(1+3*R^4/r(x,y)^4)*cos(2*th(x,y));
sigma_rth = @(x,y) -T_x/2*(1+2*R^2/r(x,y)^2 - 3*R^4/r(x,y)^4)*sin(2*th(x,y));

%下面这组画图出来位移红温部分更符合预期，但是应力应变数量级不对
sigma_xx = @(x,y) sigma_rr(x,y)*cos(-th(x,y))^2 + sigma_thth(x,y)*sin(-th(x,y))^2 + 2*sigma_rth(x,y)*sin(-th(x,y))*cos(-th(x,y));
sigma_yy = @(x,y) sigma_rr(x,y)*sin(-th(x,y))^2 + sigma_thth(x,y)*cos(-th(x,y))^2 - 2*sigma_rth(x,y)*sin(-th(x,y))*cos(-th(x,y));
sigma_xy = @(x,y) -sigma_rr(x,y)*sin(-th(x,y))*cos(-th(x,y)) + sigma_thth(x,y)*sin(-th(x,y))*cos(-th(x,y)) + sigma_rth(x,y)*(cos(-th(x,y))^2-sin(-th(x,y))^2);

D1 = (E/(1-v*v)).*[1, v, 0;
                  v, 1, 0;
                  0, 0, 0.5-0.5*v]; 
%用下面这组表达，精确解最大应力30000没问题
%sigma_xx = @(x,y)  T_x/2 * (1 - R^2/(x^2 + y^2)) + T_x/2 * (1 - 4*R^2/(x^2 + y^2) + ...
  %  3*R^4/(x^2 + y^2)^2) * (x^2 - y^2)/(x^2 + y^2);
%sigma_yy = @(x,y)  T_x/2 * (1 + R^2/(x^2 + y^2)) - T_x/2 * (1 + 3*R^4/(x^2 + y^2)^2) * (x^2 - y^2)/(x^2 + y^2);
%sigma_xy = @(x,y)  -T_x/2 * (1 + 2*R^2/(x^2 + y^2) - 3*R^4/(x^2 + y^2)^2) * (2*x*y/(x^2 + y^2));

%后面用一阶高斯求载荷向量F
n_int_BC = 5;
n_en_BC = 2; %局部节点数
[xi_BC, weight_BC] = Gauss(n_int_BC, -1, 1);
%quadrature rule
n_int_xi = 3;
n_int_eta = 3;
n_int = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);
%mesh generation
n_en = 4; %四边形单元
%直接用run出来的mesh，n-np也是
n_el = size(msh.QUADS,1); %总单元数
n_np = msh.nbNod;         %总节点数

coor(:,1) = msh.POS(:,1) + 1;
coor(:,2) = msh.POS(:,2) + 1;

x_coor = coor(:,1);
y_coor = coor(:,2);

%IEN array
IEN_tri = zeros(1,1);
IEN = msh.QUADS(:,1:4);
for ee = 1:size(IEN,1)
    IEN_tri(ee*2-1,1) = IEN(ee,1);
    IEN_tri(ee*2-1,2) = IEN(ee,2);
    IEN_tri(ee*2-1,3) = IEN(ee,3);
    IEN_tri(ee*2,1)   = IEN(ee,1);
    IEN_tri(ee*2,2)   = IEN(ee,3);
    IEN_tri(ee*2,3)   = IEN(ee,4);
end

%ID array
ID = -1 .* ones(n_np,2);
for ii = 1:size(msh.LINES,1)
    if msh.LINES(ii,3) == 11
        ID(msh.LINES(ii,1),2) = 0;
        ID(msh.LINES(ii,2),2) = 0;
    elseif msh.LINES(ii,3) == 10
        ID(msh.LINES(ii,1),1) = 0;
        ID(msh.LINES(ii,2),1) = 0;
    else
    end       
end

count = 0;
for i = 1:n_np  %总节点数，
    for j = 1:2 %两个自由度
        if ID(i,j) == -1
            count = count + 1;
            ID(i,j) = count;
        end
    end
end
n_eq = count;

%将四边形元素的节点连接信息转换为两个三角形元素的节点连接信息
%理顺一下序号
for i = 1 : n_el/2
    a1 = IEN(i,1); a2 = IEN(i,2); 
    a3 = IEN(i,3); a4 = IEN(i,4); 
    IEN(i,1) = a4; 
    IEN(i,2) = a3; 
    IEN(i,3) = a2; 
    IEN(i,4) = a1; 
    
    j = i + n_el/2;
    b1 = IEN(j,1); b2 = IEN(j,2); 
    b3 = IEN(j,3); b4 = IEN(j,4); 
    IEN(j,1) = b3; 
    IEN(j,2) = b4; 
    IEN(j,3) = b1; 
    IEN(j,4) = b2; 
    a = IEN(i, 1);
end


% LM array多一个自由度，
LM = zeros(n_el, 2*n_en);
% IEN 是一个 n_el x n_en 的矩阵，包含元素的局部节点编号
% ID 是一个 n_np x 2 的矩阵，包含节点的全局编号
for i = 1:n_el
    for j = 1:2:2:2*n_en  % 步长为 2，一次处理两个自由度
        node_index = IEN(i, ceil(j/2));
        dof_index = rem(j-1, 2) + 1;  % 计算自由度索引（1 或 2）
        LM(i, j) = ID(node_index, dof_index);
    end
end

%n_v数组，定边界法向量
n_v = msh.LINES;
for ii = 1: size(n_v,1)
    if n_v(ii,3) == 9
        n_v(ii,4) = 0; n_v(ii,5) = 1;%9为top
    elseif n_v(ii,3) == 8
        n_v(ii,4) = 1; n_v(ii,5) = 0;%8为right
    elseif n_v(ii,3) == 10
        n_v(ii,4) = -1; n_v(ii,5) = 0;%10为left
    elseif n_v(ii,3) == 11
        n_v(ii,4) = 0; n_v(ii,5) = -1;%11为bottom
    elseif n_v(ii,3) == 13           %13为圆弧处，法向量单位化
        coor1 = [coor(n_v(ii,1),1),coor(n_v(ii,1),2)]; 
        coor2 = [coor(n_v(ii,2),1),coor(n_v(ii,2),2)];
        midpoint = (coor1 + coor2)/2;
        n_v(ii,4) = -midpoint(1)/sqrt(midpoint(1)^2 + midpoint(2)^2); 
        n_v(ii,5) = -midpoint(2)/sqrt(midpoint(1)^2 + midpoint(2)^2);
    end
end


%loop over element to assembly the matrix and vector
K = zeros(n_eq,n_eq);
F = zeros(n_eq,1);
part = zeros(n_np,3);
for ee = 1:n_el
    x_ele = x_coor(IEN(ee, :));
    y_ele = y_coor(IEN(ee, :));

    k_ele = zeros(n_en*2,n_en*2);
    f_ele = zeros(n_en*2,1);

    for ll = 1:2
        n_x_l = 0.0; 
        n_y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        xi_n  = [-1,1,-1,1]; 
        eta_n = [-1,-1,1,1];
        for aa = 1:n_en 
            n_x_l = n_x_l + x_ele(aa) * Quad(aa, xi_n(ll), eta_n(ll));
            n_y_l = n_y_l + y_ele(aa) * Quad(aa, xi_n(ll), eta_n(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi_n(ll), eta_n(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ_n = dx_dxi * dy_deta - dx_deta * dy_dxi;

        for aa = 1:n_en
            number = IEN(ee,aa);
            part(number,3) = part(number,3) + 1;

            Na = Quad(aa, xi_n(ll), eta_n(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi_n(ll), eta_n(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ_n;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ_n;
            
            part(number,1) = part(number,1) + Na_x;
            part(number,2) = part(number,2) + Na_y;
        end
    end

    for ll = 1:n_int
        x_l = 0.0; y_l = 0.0; 
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1:n_en 
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end 

        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        reveJ = 1/detJ;
        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) *reveJ;
            Na_y = (Na_eta * dx_dxi - Na_xi * dx_deta) *reveJ;
                      
            for bb = 1 : n_en
                Nb = Quad(bb, xi(ll), eta(ll));
                [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi)  *reveJ;
                Nb_y = (Nb_eta * dx_dxi - Nb_xi * dx_deta)  *reveJ;

                B_a = [Na_x, 0; 
                       0, Na_y; 
                       Na_y, Na_x]; 
                B_b = [Nb_x, 0; 
                       0, Nb_y; 
                       Nb_y, Nb_x];
                a = aa*2;
                b = bb*2;
                BDB = B_a' * D1 * B_b;
                k_ele(a-1:a , b-1:b) = k_ele(a-1:a , b-1:b) + (weight(ll) * detJ) .* BDB;

            end % end of bb loop
        end % end of aa loop
    end % end of quadrature loop

    for ii = 1:2
        for aa = 1 : n_en
            pp = 2*(aa-1) + ii;
            PP = LM(ee, pp);
            if PP > 0
                for jj = 1:2
                    for bb = 1 : n_en
                        qq = 2*(bb-1) + jj;
                        QQ = LM(ee, qq);
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(pp, qq);
                        else
                            % modify F with the boundary data
                            % here we do nothing because the boundary data g is zero or
                            % homogeneous
                        end % end of QQ 
                    end % end of bb loop
                end % end of jj loop处理两个自由度
            end % end of PP c 
        end % end of aa loop
    end % end of ii loop遍历所有单元组装起来
end

%一阶高斯求载荷向量F 
for ee = 1:size(n_v,1)
    f_BC = zeros(2*n_en_BC,1);
    x_BC = x_coor(n_v(ee,1:2));
    y_BC = y_coor(n_v(ee,1:2));
    L = sqrt((x_BC(1) - x_BC(2))^2 + (y_BC(1) - y_BC(2))^2);
    for ll = 1:n_int_BC
        x_l_BC = 0.0;
        y_l_BC = 0.0;
        for aa = 1:n_en_BC
            x_l_BC = x_l_BC + x_BC(aa) * PolyShape(1, aa, xi_BC(ll), 0);
            y_l_BC = y_l_BC + y_BC(aa) * PolyShape(1, aa, xi_BC(ll), 0);   

        end
 
        si_xx = sigma_xx(x_l_BC,y_l_BC); 
        si_yy = sigma_yy(x_l_BC,y_l_BC); 
        si_xy = sigma_xy(x_l_BC,y_l_BC);

        h_v = [si_xx, si_xy; 
               si_xy, si_yy] * [n_v(ee,4);n_v(ee,5)];
        h_x = h_v(1); 
        h_y = h_v(2);

        for bb = 1:n_en_BC
            b = bb*2;
            f_BC(b-1) = f_BC(b-1) + weight_BC(ll) * PolyShape(1, bb, xi_BC(ll), 0) * h_x * 0.5 * L ;
            f_BC(b)   = f_BC(b)   + weight_BC(ll) * PolyShape(1, bb, xi_BC(ll), 0) * h_y * 0.5 * L ;
        end
    end

    for dof = 1:2 
        PP_1 = ID(n_v(ee,1), dof);%获取节点的全局索引
        if PP_1 > 0%如果全局索引大于0，将边界元素的力向量f_BC的贡献添加到全局力向量
        F(PP_1) = F(PP_1) + f_BC(dof);
        end
        PP_2 = ID(n_v(ee,2), dof);
        if PP_2 > 0
        F(PP_2) = F(PP_2) + f_BC(2+dof);
        end
    end
end


dn = K \ F;
%储存全局节点的xy坐标
disp_x = zeros(n_np, 1);
disp_y = zeros(n_np, 1);
for ii = 1 : n_np
    index_x = ID(ii,1);
    index_y = ID(ii,2);
    if index_x > 0
        disp_x(ii) = dn(index_x);
    else
    end
    if index_y > 0
        disp_y(ii) = dn(index_y);
    else
    end
end


% 下面计算应变和应力，
strain = zeros(n_el, 3);
for ee = 1:n_el
    x_ele = x_coor(IEN(ee, 1:4));
    y_ele = y_coor(IEN(ee, 1:4));
    for ll = 1:n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1:4
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        reveJ = detJ;
        strains = zeros(3, 1);
        for aa = 1:4
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) * reveJ;
            Na_y = (Na_eta * dx_dxi - Na_xi * dx_deta) * reveJ;

        B_a = [Na_x, 0; 
               0, Na_y; 
               Na_y, Na_x]; 
            node_index = IEN(ee, aa);
            node_disp_x = disp_x(node_index);
            node_disp_y = disp_y(node_index);
            strains = strains + B_a(:, 1) * node_disp_x + B_a(:, 2) * node_disp_y;
        end
        for cc = 1:3
            strain(ee, cc) = strain(ee, cc) + weight(ll) * strains(cc) * detJ;
        end
    end
end

stress = zeros(n_el, 3);
for ii = 1:n_el
    stress(ii, :) = D1 * strain(ii, :)';
end

% 初始化节点应变数组，每个节点有3个应变分量（xx, yy, xy）
node_strain = zeros(n_np, 3);
ele_area = zeros(n_el, 1);
for ee = 1:n_el
    % 获取四边形的四个顶点坐标
    x_ele = x_coor(IEN(ee, 1:4));
    y_ele = y_coor(IEN(ee, 1:4));

    % 计算两个三角形的面积
    area_triangle1 = 0.5 * abs((x_ele(1) - x_ele(3)) * (y_ele(2) - y_ele(1)) - (x_ele(2) - x_ele(1)) * (y_ele(3) - y_ele(1)));
    area_triangle2 = 0.5 * abs((x_ele(3) - x_ele(1)) * (y_ele(4) - y_ele(3)) - (x_ele(4) - x_ele(3)) * (y_ele(1) - y_ele(3)));

    % 两个三角形面积之和
    ele_area(ee) = area_triangle1 + area_triangle2;
end

for ee = 1:n_el
    for aa = 1:4
        node_index = IEN(ee, aa);
        node_strain(node_index, :) = node_strain(node_index, :) + strain(ee, :) * ele_area(ee);
    end
end


for nn = 1:n_np
    total_area = 0;
    for ee = 1:n_el
        if any(IEN(ee, :) == nn)
            total_area = total_area + ele_area(ee);
        end
    end
    if total_area > 0
        node_strain(nn, :) = node_strain(nn, :) / total_area;
    end
end

% 初始化节点应力数组，每个节点有3个应力分量（xx, yy, xy）
stress_num = zeros(n_np, 3);
for nn = 1:n_np
    stress_num(nn, :) = D1 * node_strain(nn, :)';
end


% %下面是画图代码
% 1应力sigma_xx数值解
figure;  
hold on
trisurf(IEN_tri, x_coor, y_coor, stress_num(:, 1));
axis equal;
colormap jet
shading interp
title('stress (\sigma_{xx})num');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar

% 2应力sigma_yy数值解
figure;
hold on
trisurf(IEN_tri, x_coor, y_coor, stress_num(:, 2));
axis equal;
colormap jet
shading interp
title('stress (\sigma_{yy})num');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar

% 3应力sigma_xy数值解
figure;
hold on
trisurf(IEN_tri, x_coor, y_coor, stress_num(:, 3));
axis equal;
colormap jet
shading interp
title('stress (\sigma_{xy})num');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar

%下面是精确解,遍历每个节点，每个节点有三个应力分量
stress_exact = zeros(n_np,3);
for aa = 1 : n_np
    stress_exact(aa,1) = sigma_xx(x_coor(aa),y_coor(aa));
    stress_exact(aa,2) = sigma_yy(x_coor(aa),y_coor(aa));
    stress_exact(aa,3) = sigma_xy(x_coor(aa),y_coor(aa));
end


% 4应力sigma-xx的精确解
figure;
hold on
trisurf(IEN_tri, x_coor, y_coor, stress_exact(:, 1));
axis equal;
colormap jet
shading interp
title('stress (\sigma_{xx})exact');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar

% 5应力sigma-yy的精确解
figure; 
hold on
trisurf(IEN_tri, x_coor, y_coor, stress_exact(:, 2));
axis equal;
colormap jet
shading interp
title('stress (\sigma_{yy})exact');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar

% 6应力sigma-xy的精确解
figure; 
hold on
trisurf(IEN_tri, x_coor, y_coor, stress_exact(:, 3));
axis equal;
colormap jet
shading interp
title('stress (\sigma_{xy})exact');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar

% 7绘制应变 epsilon_xx 数值解
figure;
hold on;
trisurf(IEN_tri, x_coor, y_coor, node_strain(:, 1));
axis equal;
colormap jet;
shading interp;
title('Strain (\epsilon_{xx}) Numerical');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar;

% 8绘制应变 epsilon_yy 数值解
figure;
hold on;
trisurf(IEN_tri, x_coor, y_coor, node_strain(:, 2));
axis equal;
colormap jet;
shading interp;
title('Strain (\epsilon_{yy}) Numerical');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar;

% 9绘制应变 epsilon_xy 数值解
figure;
hold on;
trisurf(IEN_tri, x_coor, y_coor, node_strain(:, 3));
axis equal;
colormap jet;
shading interp;
title('Strain (\epsilon_{xy}) Numerical');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar;

% 10 x方向数值解位移图
figure;
hold on;
trisurf(IEN_tri, x_coor, y_coor, disp_x);
shading interp;
axis equal;
colormap jet;
title('x - direction displacement');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar

% 11 y方向数值解位移图
figure;
hold on;
trisurf(IEN_tri, x_coor, y_coor, disp_y);
shading interp;
axis equal;
colormap jet;
title('y - direction displacement');
xlabel('x - coordinate');
ylabel('y - coordinate');
colorbar