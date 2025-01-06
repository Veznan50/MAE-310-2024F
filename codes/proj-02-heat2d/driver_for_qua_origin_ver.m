clear ; clc;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

exact_xx = @(x,y) -2*y*(1-y);
exact_xy = @(x,y) (1-2*x)*(1-2*y);
exact_yy = @(x,y) -2*x*(1-x);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term
for iii = 10:10:100
% quadrature rule

n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_en   = 4;               % number of nodes in an element
n_el_x = iii;               % number of elements in x-dir
n_el_y = iii;               % number of elements in y-dir
n_el   = n_el_x * n_el_y; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array
IEN = zeros(n_el, n_en);
for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    ee = (ey-1) * n_el_x + ex; % element index
    IEN(ee, 1) = (ey-1) * n_np_x + ex;
    IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
    IEN(ee, 3) =  ey    * n_np_x + ex + 1;
    IEN(ee, 4) =  ey    * n_np_x + ex;
  end
end

% ID array
ID = zeros(n_np,1);
counter = 0;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index) = counter;  
  end
end

n_eq = counter;

LM = ID(IEN);

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
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    
    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      
      f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
      
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        
        k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop
 
  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      
      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
% 狄利克雷边界条件dirichlet boundary condition，先把边角节点定住，
%加网格后再找规律把各个边线的边界条件完善
dirichlet_location=zeros(n_np,2);
for bc = 1 : n_np
    if x_coor(bc) == 0 
        dirichlet_location(bc,1)=bc;
        dirichlet_location(bc,2)=bc;
    end
end
dirichlet_g=zeros(n_np,2);
for bc=1 : n_np
    if ismember(bc,dirichlet_location)
    dirichlet_g(bc,1)=0;
    dirichlet_g(bc,2)=0;
    end
end
% 设置节点1和5的第二个自由度上的位移为0
dirichlet_location(1, 2) = 1;
dirichlet_location(5, 2) = 1;
dirichlet_g(1, 2) = 0;
dirichlet_g(5, 2) = 0;

% 设置节点3和6的第一个自由度上的位移为0
dirichlet_location(3, 1) = 1;
dirichlet_location(6, 1) = 1;
dirichlet_g(3, 1) = 0;
dirichlet_g(6, 1) = 0;

% 诺伊曼边界条件neumann boundary condition
neumann_h=zeros(n_np,2);
for bc=1 : n_np
    if x_coor(bc)==1
        neumann_h(bc,1)=1;
        neumann_h(bc,2)=0;
    end
end
%设置外力加载为10kPa
neumann_h(1, 1) = 10;  % 节点1的第一个自由度
neumann_h(2, 1) = 10;  % 节点2的第一个自由度

% 设置节点1和5的第一个自由度上的应力为0
neumann_h(1, 1) = 0;
neumann_h(5, 1) = 0;

% 设置节点3和6的第二个自由度上的应力为0
neumann_h(3, 2) = 0;
neumann_h(6, 2) = 0;

% 设置节点6，4和5的应力为0
neumann_h(4, 1) = 0;
neumann_h(4, 2) = 0;
neumann_h(6, 1) = 0;
neumann_h(6, 2) = 0;
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


%解题目HW6中2-b问的
e_0 = 0.0; %||e||0，先这样记着，最后再统一开根号
e_1 = 0.0; %||e||1
u_2 = 0.0; %||u||2

for ee = 1 : n_el
    x_ele = x_coor(IEN(ee, :));
    y_ele = y_coor(IEN(ee, :));
    d_ele = disp(IEN(ee, :));
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
        for aa = 1 : n_en %this loop is to calculate the value of uh_l
            %计算数值解的位移
            uh_l = uh_l + d_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            %计算数值解位移在x和y方向上的导数用于u_2中
            uh_x_l = uh_x_l + d_ele(aa) * (Na_xi * dy_deta - Na_eta * dy_dxi) * reveJ;
            uh_y_l = uh_y_l + d_ele(aa) * (Na_eta * dx_dxi - Na_xi * dx_deta) * reveJ;
        end
        %计算精确解u
        u_l    = exact   (x_l, y_l);
        u_x_l  = exact_x (x_l, y_l); 
        u_y_l  = exact_y (x_l, y_l);
        u_xx_l = exact_xx(x_l, y_l); 
        u_xy_l = exact_xy(x_l, y_l); 
        u_yy_l = exact_yy(x_l, y_l);
        %算误差的平方，表达式在bb提交的作业中有写明
        e_0 = e_0 + weight(ll) * (uh_l - u_l)^2;
        e_1 = e_1 + weight(ll) * ((uh_l - u_l)^2 + (uh_x_l - u_x_l)^2 + (uh_y_l - u_y_l)^2);
        u_2 = u_2 + weight(ll) * (u_l^2 + u_x_l^2 + u_xx_l^2 + 2*u_xy_l^2 + u_y_l^2 + u_yy_l^2);
    end
end

ch2 = sqrt(e_0/u_2);%开根号
ch1 = sqrt(e_1/u_2);
store(iii/10,:) = [log(hx), log(ch2), log(ch1)];
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




