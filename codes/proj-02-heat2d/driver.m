clear ; clc;

kappa = 1.0; % conductivity
%另存了一个算四边形收敛速度的代码文件，这里专门改三角单元
% exact solution定义四边形精确解
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);%对x的偏导数
exact_y = @(x,y) x*(1-x)*(1-2*y);%对y偏导

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% Trianglerature rule
n_int_xi  = 3;% 定义Gauss积分的xi方向的点数
n_int_eta = 3;% 定义Gauss积分的eta方向的点数
n_int     = n_int_xi * n_int_eta;% 计算总的积分点数
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);
% 调用Gauss2D函数生成Gauss积分点和权重，这里应该不用管吧

% mesh generation网格划分，改三角
n_en   = 3;               % number of nodes in an element，三角改3？
n_el_x = 60;               % number of elements in x-dir，元素数，这咋改，也是两倍吗。
n_el_y = 60;               % number of elements in y-dir
%这改完，后面节点数要改，多除2，网格尺寸也要改
n_el   = 2*n_el_x * n_el_y; % total number of elements三角要两倍单元数

%计算x和y方向的节点数以及总节点数。节点数不变，但是单元数变了
n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points总节点数不变

%创建两个数组来存储所有节点的x和y坐标。x、y的总量为nnp
x_coor = zeros(n_np, 1);
y_coor = x_coor;

%计算x和y方向的网格尺寸。
hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates使用双重循环生成所有节点的坐标，应该不用改
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index，全局坐标索引
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array三角元素要改IEN，一个三角元里只有3个节点
IEN = zeros(n_el, n_en);%IEN=（元素编号，元素内节点编号=3）
for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    ee = (ey-1) * n_el_x + ex; % element index
    IEN(ee*2-1, 1) = (ey-1) * n_np_x + ex;
    IEN(ee*2-1, 2) =  ey    * n_np_x + ex + 1;
    IEN(ee*2-1, 3) =  ey    * n_np_x + ex ;
    IEN(ee*2, 1)   = (ey-1) * n_np_x + ex;
    IEN(ee*2, 2)   = (ey-1) * n_np_x + ex + 1;
    IEN(ee*2, 3)   =  ey    * n_np_x + ex + 1;
  end
end

% ID array
ID = zeros(n_np,1);
counter = 0;%% 初始化计数器
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index) = counter;  
  end
end

n_eq = counter;

LM = ID(IEN);%合并成LM

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(n_en, n_en); % element stiffness matrix
  f_ele = zeros(n_en, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Triangle(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Triangle(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Triangle_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;%雅克比行列式，后面要用
    reveJ=1/detJ;

    for aa = 1 : n_en
      Na = Triangle(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Triangle_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      
      f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
      
      for bb = 1 : n_en
        Nb = Triangle(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Triangle_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        
        k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end % end of bb loop
    end % end of aa loop
  end % end of Trianglerature loop
 
  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      
      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data，改边界条件g
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    disp(ii) = dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end
end

% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "disp", "n_el_x", "n_el_y");


 % e_0 第一个误差,dx_dxi 要改成detJ,在103行
    tops = 0.0;     %分子
    bottoms = 0.0;   %分母
    for ee = 1 : n_el
        for qua = 1 : n_int   
           
            top_1 = 0.0;%数值
            top_2 = 0.0;%精确
            bottom = 0.0;
            for aa = 1 : n_en
               
                top_1 = top_1 + disp(IEN(ee, aa)) * Triangle(aa, xi, eta);
                top_2 = top_2 + x_coor(IEN(ee, aa)) * Triangle(aa, xi, eta);
                bottom = bottom + x_coor(IEN(ee, aa)) * Triangle(aa, xi, eta);
            end
            top = (top_1 - exact(top_2))^2 * detJ * weight(qua);
            tops = tops + top;
            bottom = exact(bottom)^2 * detJ * weight(qua);
            bottoms =bottoms + bottom;
        end
    end
    e_0(ii) = sqrt(tops/bottoms);
   

    % e_1 第二个误差
    tops = 0.0;
    bottoms = 0.0;
    for ee = 1 : n_el
        for qua = 1 : n_int
            
            top_1 = 0.0;
            top_2 = 0.0;
            bottom = 0.0;
            for aa = 1 : n_en
                
                top_1 = top_1 + disp(IEN(ee, aa)) * Triangle_grad(aa, xi, eta);
                top_2 = top_2 + x_coor(IEN(ee, aa)) * Triangle(aa, xi, eta);
                bottom = bottom + x_coor(IEN(ee, aa)) * Triangle(aa, xi, eta);
            end
            
            top = (top_1 * reveJ - exact_x(top_2))^2 * detJ * weight(qua);
            tops = tops + top;
            bottom = (exact_x(bottom))^2 * detJ * weight(qua);
            bottoms =bottoms + bottom;
        end
    end
    e_1(ii) = sqrt(tops/bottoms);

%计算x和y方向的网格尺寸。
hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

%hx = 1./ (2: 2: 16);
plot(log(hx), log(e_0), '-o','LineWidth',3)

hold on

plot(log(hx), log(e_1), '-x','LineWidth',3)

slope_e_0 = (log(e_0(end))-log(e_0(1)))/(log(hx(end))-log(hx(1)));

slope_e_1 = (log(e_1(end))-log(e_1(1)))/(log(hx(end))-log(hx(1)));












% EOF