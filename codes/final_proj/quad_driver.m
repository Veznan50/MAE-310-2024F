clear ; clc;
%复原了proj-02里四边单元的传热的代码，final-proj当时不该在proj-02里直接改，
% 重命名为final那个里面改得乱七八糟的，这里改看得清楚点
%主要是ID和LM改成两个自由度下的
kappa = 1.0; % conductivity
E = 1E9;
v = 0.3;
T_x = 1E4;
R = 0.5;
store = zeros(4,3);
% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

exact_xx = @(x,y) -2*y*(1-y);
exact_xy = @(x,y) (1-2*x)*(1-2*y);
exact_yy = @(x,y) -2*x*(1-x);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term
%用平面应力推导出的不带孔矩形板块的解，u和v都是上面这个exact = @(x,y) x*(1-x)*y*(1-y)
f1 = @(x,y) (2*E*y*(y - 1))/(v^2 - 1) - (E*(v/2 - 1/2)*((x - 1)*(y - 1) + x*y + ...
    2*x*(x - 1) + x*(y - 1) + y*(x - 1)))/(v^2 - 1) + (E*v*((x - 1)*(y - 1) + ...
    x*y + x*(y - 1) + y*(x - 1)))/(v^2 - 1);%x方向分量
f2 = @(x,y) (2*E*x*(x - 1))/(v^2 - 1) - (E*(v/2 - 1/2)*((x - 1)*(y - 1) + x*y + ...
    x*(y - 1) + y*(x - 1) + 2*y*(y - 1)))/(v^2 - 1) + (E*v*((x - 1)*(y - 1) + ...
    x*y + x*(y - 1) + y*(x - 1)))/(v^2 - 1);%y方向分量
%传热问题从10间隔10到100没问题
%这里从10间隔10到100运行内存不足，到70都不行，
%从20间隔20到100也不行，20到80可以，但是只剩四个点了
for iii = 20:20:80
% quadrature rule不用动
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation不动
n_en   = 4;               % number of nodes in an element
n_el_x = iii;               % number of elements in x-dir
n_el_y = iii;               % number of elements in y-dir
n_el   = n_el_x * n_el_y; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates不动
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array不动，四边形的
IEN = zeros(n_el, n_en);
for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    ee = (ey-1) * n_el_x + ex; 
    IEN(ee, 1) = (ey-1) * n_np_x + ex;
    IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
    IEN(ee, 3) =  ey    * n_np_x + ex + 1;
    IEN(ee, 4) =  ey    * n_np_x + ex;
  end
end

% ID array两个自由度
ID = zeros(n_np,2);
counter = 0;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index,1) = counter;  
    counter = counter + 1;
    ID(index,2) = counter;
  end
end
n_eq = counter;
% end of ID edit

% LM array多一个自由度，
LM = zeros(n_el, 2*n_en);
% IEN 是一个 n_el x n_en 的矩阵，包含元素的局部节点编号
% ID 是一个 n_np x 2 的矩阵，包含节点的全局编号
for ii = 1:n_el
    for jj = 1:2:2:2*n_en  % 步长为 2，一次处理两个自由度
        node_index = IEN(ii, ceil(jj/2));
        dof_index = rem(jj-1, 2) + 1;  % 计算自由度索引（1 或 2）
        LM(ii, jj) = ID(node_index, dof_index);
    end
end

% allocate the stiffness matrix and load vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);
%平面应力矩阵
D1 = (E/(1-v^2)).*[1, v,     0;
                   v, 1,     0;
                   0, 0, (1-v)/2]; 
% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, :) );
  y_ele = y_coor( IEN(ee, :) );
  
  k_ele = zeros(n_en*2, n_en*2); % element stiffness matrix44变88数据量4倍
  f_ele = zeros(n_en*2, 1);    % element load vector数据量翻倍
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));

      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end                                       
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;% 一样的雅克比
    
    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

      % f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
      a = 2*aa;%xy方向不同的f分别编号
      f_ele(a-1) = f_ele(a-1) + weight(ll) * detJ * f1(x_l, y_l) * Na;
      f_ele(a) = f_ele(a) + weight(ll) * detJ * f2(x_l, y_l) * Na;

   
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        % 笔记10-76里的Ba和Bb
        B_a = [Na_x, 0; 
               0, Na_y; 
               Na_y, Na_x]; 
        B_b = [Nb_x, 0; 
               0, Nb_y; 
               Nb_y, Nb_x];
        BDB = B_a' * D1 * B_b;
        b = 2*bb;
        k_ele(a-1:a,b-1:b) = k_ele(a-1:a,b-1:b) + (weight(ll) * detJ) .* BDB;
       % k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop
%开始四个循环，这里从ii开始，aa，jj，bb结束
   for ii = 1:2
        for aa = 1 : n_en
            pp = 2*(aa-1) + ii;
            PP = LM(ee, pp);
            if PP > 0
                F(PP) = F(PP) + f_ele(pp);
                for jj = 1:2
                    for bb = 1 : n_en
                        qq = 2*(bb-1) + jj;
                        QQ = LM(ee, qq);
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(pp, qq);
                        else
                         % modify F with the boundary data
                         % here we do nothing because the boundary data g
                         % is zero or homogeneous
                         
                        end % end of QQ 
                    end % end of bb loop
                end % end of jj loop处理两个自由度
            end % end of PP c 
        end % end of aa loop
    end % end of ii loop遍历所有单元组装起来
end


% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp_x = zeros(n_np, 1);
disp_y = zeros(n_np, 1);%多储存一个y方向位移

for ii = 1 : n_np
    index_x = ID(ii,1);
    index_y = ID(ii,2);
     if index_x > 0
        disp_x(ii) = dn(index_x);
        if index_y > 0
            disp_y(ii) = dn(index_y);
        else
        end
    else
    % modify disp with the g data. Here it does nothing because g is zero
    end
end


%解题目HW6中2-b问的
e_0 = 0.0; %||e||0，先这样记着，最后再统一开根号
e_1 = 0.0; %||e||1
u_2 = 0.0; %||u||2

for ee = 1 : n_el
    x_ele = x_coor(IEN(ee, :));%坐标
    y_ele = y_coor(IEN(ee, :));
    dx_ele = disp_x(IEN(ee, :));%位移
    dy_ele = disp_y(IEN(ee, :));
    %计算坐标位移，直接用老师的
    for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0; 
        dy_dxi = 0.0; dy_deta = 0.0;

        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll)); 
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] =  Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        reveJ=1/detJ;
        %计算各个数值解
        uh_l   = 0.0; 
        uh_x_l = 0.0; 
        uh_y_l = 0.0;

        for aa = 1 : n_en 
            %计算数值解的位移
            uh_l = uh_l + dx_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            %计算数值解位移在x和y方向上的导数用于u_2中
            uh_x_l = uh_x_l + dx_ele(aa) * (Na_xi * dy_deta - Na_eta * dy_dxi) * reveJ;
            uh_y_l = uh_y_l + dy_ele(aa) * (Na_eta * dx_dxi - Na_xi * dx_deta) * reveJ;
        end
        %计算精确解u
        u_l    = exact   (x_l, y_l);
        u_x_l  = exact_x (x_l, y_l); 
        u_y_l  = exact_y (x_l, y_l);
        u_xx_l = exact_xx(x_l, y_l); 
        u_xy_l = exact_xy(x_l, y_l); 
        u_yy_l = exact_yy(x_l, y_l);


        %算误差的平方，表达式在bb提交的作业中有写明
        %三个差值
        m1 = uh_l - u_l;
        m2 = uh_x_l - u_x_l;
        m3 = uh_y_l - u_y_l;
        e_0 = e_0 + weight(ll) * (m1)^2;
        e_1 = e_1 + weight(ll) * ((m1)^2 + (m2)^2 + (m3)^2);
        u_2 = u_2 + weight(ll) * (u_l^2 + u_x_l^2 + u_xx_l^2 + 2*u_xy_l^2 + u_y_l^2 + u_yy_l^2);
    end
end

ch2 = sqrt(e_0/u_2);%开根号
ch1 = sqrt(e_1/u_2);
store(iii/20,:) = [log(hx), log(ch2), log(ch1)];
end

plot(store(:,1),store(:,2),'-o','LineWidth',3);
hold on
plot(store(:,1),store(:,3),'-x','LineWidth',3);
%第一个线对应e_0，第二个对应e_1
% 计算 L2 范数误差曲线的斜率
coeff_L2 = polyfit(store(:,1), store(:,2), 1);
slope_L2 = coeff_L2(1); 

% 计算 H1 范数误差曲线的斜率
coeff_H1 = polyfit(store(:,1), store(:,3), 1);
slope_H1 = coeff_H1(1); 

% 输出斜率
fprintf('L2 范数误差曲线的斜率: %f\n', slope_L2);
fprintf('H1 范数误差曲线的斜率: %f\n', slope_H1);




