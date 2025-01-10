clear ; clc;
run('meshhh.m')
kappa = 1.0; % conductivity
E = 1E9;
v = 0.3;
T_x = 1E4;
R = 0.5;
x=1;y=1;

%预先设置误差存储
ch22 = [];
ch11 = [];
length = [];
ch2 = 0;
ch1 = 0;

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

exact_xx = @(x,y) -2*y*(1-y);
exact_xy = @(x,y) (1-2*x)*(1-2*y);
exact_yy = @(x,y) -2*x*(1-x);

%f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term
%用平面应力推导出的不带孔矩形板块的解，u和v都是上面这个exact = @(x,y) x*(1-x)*y*(1-y)
f1 = @(x,y) (2*E*y*(y - 1))/(v^2 - 1) - (E*(v/2 - 1/2)*((x - 1)*(y - 1) + x*y + ...
    2*x*(x - 1) + x*(y - 1) + y*(x - 1)))/(v^2 - 1) + (E*v*((x - 1)*(y - 1) + ...
    x*y + x*(y - 1) + y*(x - 1)))/(v^2 - 1);%x方向分量
f2 = @(x,y) (2*E*x*(x - 1))/(v^2 - 1) - (E*(v/2 - 1/2)*((x - 1)*(y - 1) + x*y + ...
    x*(y - 1) + y*(x - 1) + 2*y*(y - 1)))/(v^2 - 1) + (E*v*((x - 1)*(y - 1) + ...
    x*y + x*(y - 1) + y*(x - 1)))/(v^2 - 1);%y方向分量
%用平面应变推导出的
f11 = @(x,y)  (E*v*((x - 1)*(y - 1) + x*y + x*(y - 1) + ...
    y*(x - 1)))/((2*v - 1)*(v + 1)) - (E*(v - 1/2)*((x - 1)*(y - 1) + x*y + ...
    2*x*(x - 1) + x*(y - 1) + y*(x - 1)))/((2*v - 1)*(v + 1)) - (2*E*y*(v - 1)*(y - 1))/((2*v - 1)*(v + 1));
f22 = @(x,y)  (E*v*((x - 1)*(y - 1) + x*y + x*(y - 1) + ...
    y*(x - 1)))/((2*v - 1)*(v + 1)) - (E*(v - 1/2)*((x - 1)*(y - 1) + x*y + ...
    x*(y - 1) + y*(x - 1) + 2*y*(y - 1)))/((2*v - 1)*(v + 1)) - (2*E*x*(v - 1)*(x - 1))/((2*v - 1)*(v + 1));

%下面是带孔板的精确解
%sigma_xy表达式
sigma_xy = -T_x/2 * (1 + 2*R^2/(x^2 + y^2) - 3*R^4/(x^2 + y^2)^2) * (2*x*y/(x^2 + y^2));

%对sigma_xy关于x求偏导得到f1
f111 = @(x,y)  -T_x * ( (2*y^3)/(x^2 + y^2)^2 - (4*R^2*x*y)/(x^2 + y^2)^3 + (6*R^4*x*y)/(x^2 + y^2)^4 );

%对sigma_xy关于y求偏导得到f2
f222 = @(x,y)  -T_x * ( (2*x^3)/(x^2 + y^2)^2 - (4*R^2*x*y)/(x^2 + y^2)^3 + (6*R^4*x*y)/(x^2 + y^2)^4 );

%sigma_xx表达式
sigma_xx = T_x/2 * (1 - R^2/(x^2 + y^2)) + T_x/2 * (1 - 4*R^2/(x^2 + y^2) + 3*R^4/(x^2 + y^2)^2) * (x^2 - y^2)/(x^2 + y^2);

%对sigma_xx关于x求偏导得到f3
f333 = @(x,y)  T_x * ( R^2/(x^2 + y^2)^2 + (4*R^2*x)/(x^2 + y^2)^3 - (6*R^4*x)/(x^2 + y^2)^4 + (4*x*y^2)/(x^2 + y^2)^2 - (8*R^2*x*y^2)/(x^2 + y^2)^3 + (12*R^4*x*y^2)/(x^2 + y^2)^4 );

%sigma_yy表达式
sigma_yy = T_x/2 * (1 + R^2/(x^2 + y^2)) - T_x/2 * (1 + 3*R^4/(x^2 + y^2)^2) * (x^2 - y^2)/(x^2 + y^2);

%对sigma_yy关于y求偏导得到f4
f444 = @(x,y)  T_x * ( -R^2/(x^2 + y^2)^2 + (4*R^2*y)/(x^2 + y^2)^3 - (6*R^4*y)/(x^2 + y^2)^4 - (4*x^2*y)/(x^2 + y^2)^2 + (8*R^2*x^2*y)/(x^2 + y^2)^3 - (12*R^4*x^2*y)/(x^2 + y^2)^4 );

f1111 = @(x,y)  -(f333) -(f222);%x方向的力
f2222 = @(x,y)  -(f111) - (f444);%y方向的力
%观察可知f1111+f2222=0，这和物理意义相符，说明计算结果正确

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation最大值设置为90了，间隔20，不然说我运行内存不足
for iii = 10:20:90
n_en   = 4;               % number of nodes in an element
n_el_x = iii;               % number of elements in x-dir
n_el_y = n_el_x;               % number of elements in y-dir
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

% ID array多加一个维度，注意计数器的表达式
ID = zeros(n_np,2);
counter = 1; % 修改初始值为 1
for ny = 2 : n_np_y - 1
    for nx = 2 : n_np_x - 1
        index = (ny-1)*n_np_x + nx;
        ID(index,1) = counter;
        ID(index,2) = counter + 1; % 直接赋值给第二列
        counter = counter + 2; % 然后 counter 增加 2，准备下一对 ID 的赋值
    end
end
n_eq = counter - 1; % counter 最后会超过实际的 n_eq，所以减 1 以得到正确的 n_eq

%合成LM
LM1 = zeros(size(IEN));
LM2 = zeros(size(IEN));
for BB = 1:n_en
    LM1(:, BB) = ID(IEN(:, BB), 1);
    LM2(:, BB) = ID(IEN(:, BB), 2);
end

%下面开始组装矩阵和载荷向量
%平面应力矩阵
D1 = (E/(1-v^2)).*[1, v,     0;
                   v, 1,     0;
                   0, 0, (1-v)/2]; 
%平面应变矩阵
D2 = (E/((v+1)*(2*v-1))).*[(v-1), -v, 0;
                           -v, (v-1), 0;
                            0,    0,  v-0.5];    
% allocate the stiffness matrix and load vector
K = zeros(n_eq,n_eq);
F = zeros(n_eq,1);

Dof = 2;%degree of freedom
% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(n_en*2, n_en*2); % element stiffness matrix
  f_ele = zeros(n_en*2, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; 
    y_l = 0.0;
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
%笔记10-76的矩阵和pp
      B_a = [Na_x, 0; 
             0, Na_y; 
             Na_y, Na_x]; 

      pp = Dof*(aa-1);
%组装载荷
      f_ele(pp+1) = f_ele(pp+1) + weight(ll) * detJ * f1(x_l, y_l) * Na;
      f_ele(pp+2) = f_ele(pp+2) + weight(ll) * detJ * f2(x_l, y_l) * Na;

      for bb = 1 : n_en
          Nb = Quad(bb, xi(ll), eta(ll));
         [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
%矩阵
        B_b = [Nb_x, 0;
               0, Nb_y; 
               Nb_y, Nb_x];

        for i = 1 : Dof
         e_i = (i==1)*[1,0] + (i==2)*[0,1];
           for j = 1 : Dof
%课上写的qq
           qq = Dof*(bb-1)+j;

           e_j = (j==1)*[1,0]+(j==2)*[0,1];

           BDB = B_a' * D1 * B_b;
%组装刚度阵，D1为平面应力下的矩阵
           k_ele(pp+i, qq) = k_ele(pp+i, qq) + weight(ll) * detJ * e_i * BDB * e_j';

           end % end of j loop
        end % end of i loop
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop


  for aa = 1 : n_en
     for i = 1 : Dof  
    PP = ID(IEN(ee,aa) , i);
    if PP > 0
      F(PP) = F(PP) + f_ele(pp + i);
      
      for bb = 1 : n_en
         for j = 1 : Dof
        QQ = ID(IEN(ee,bb) , j);
         if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(Dof*(aa-1)+i, Dof*(bb-1)+j);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
            % 狄利克雷边界条件dirichlet boundary condition，先把边角节点定住，
            %加网格后再找规律把各个边线的边界条件完善
          end
        end
      end 
    end
    end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 2);

 % modify disp with the g data.
for ii = 1 : n_np
    for dof=1:Dof
  index = ID(ii, dof);
   if index > 0
    disp(ii, dof) = dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
    end
   end
end


for ee = 1 : n_el
    x_ele = x_coor(IEN(ee, 1:n_en));
    y_ele = y_coor(IEN(ee, 1:n_en));
    d_ele = [ disp(IEN(ee, 1:n_en),1), disp(IEN(ee, 1:n_en),2)]';
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
        uh_l   = zeros(2,1); 
        uh_x_l = zeros(2,1); 
        uh_y_l = zeros(2,1);
        for aa = 1 : n_en %this loop is to calculate the value of uh_l
            %计算数值解的位移
            uh_l(i) = uh_l(i) + d_ele(i,aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            %计算数值解位移在x和y方向上的导数用于u_2中
            uh_x_l(i) = uh_x_l(i) + d_ele(i,aa) * (Na_xi * dy_deta - Na_eta * dy_dxi) * reveJ;%后面这坨其实是Na-x
            uh_y_l(i) = uh_y_l(i) + d_ele(i,aa) * (Na_eta * dx_dxi - Na_xi * dx_deta) * reveJ;
        end
        %计算精确解u
        u_l    = [exact(x_l, y_l), exact(x_l, y_l)]';
        u_x_l  = [exact_x(x_l, y_l), exact_x(x_l, y_l)]'; 
        u_y_l  = [exact_y(x_l, y_l), exact_y(x_l, y_l)]';
        u_xx_l = [exact_xx(x_l, y_l), exact_xx(x_l, y_l)]'; 
        u_xy_l = [exact_xy(x_l, y_l), exact_xy(x_l, y_l)]'; 
        u_yy_l = [exact_yy(x_l, y_l), exact_yy(x_l, y_l)]';
        %算误差的平方，表达式在bb提交的作业中有写明
        m1 = uh_l - u_l;
        m2 = uh_x_l - u_x_l;
        m3 = uh_y_l - u_y_l;
        %e_0 = e_0 + weight(ll) * (m1).^2;
        %e_1 = e_1 + weight(ll) * ((m1).^2 + (m2).^2 + (m3).^2);
        %u_2 = u_2 + weight(ll) * (u_l.^2 + u_x_l.^2 + u_xx_l.^2 + 2*u_xy_l.^2 + u_y_l.^2 + u_yy_l.^2);
        ch2 = ch2 + weight(ll) * detJ * m1.^2;
        ch1 = ch1 + weight(ll) * detJ * (m2.^2 + m3.^2);
    end
end

ch2 = sqrt(ch2);
ch1 = sqrt(ch1);

ch22 = [ch22, ch2];
ch11 = [ch11, ch1];
length = [length, hx];%hx=hy
% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "disp", "n_el_x", "n_el_y");
end
%不想用数组拼接了，老是向量不匹配，都是5嘛怎么不匹配了
%后处理，画图。
figure;
loglog(length, ch22(1,:), 'b--o', 'DisplayName', 'L2');
%在双对数坐标系下绘制length与ch22(1,:)的数据点，线型为蓝色虚线带圆圈标记
hold on;
loglog(length, ch11(1,:), 'r:*', 'DisplayName', 'H1');
%线型为红色点线带星号标记
%不用数组拼接，分别计算对数
logL2 = log(ch22(1,:));
logH1 = log(ch11(1,:));
loghx = log(length);
%使用多项式拟合函数polyfit，对loghx和logL2进行一次多项式拟合
coeffL2 = polyfit(loghx, logL2, 1);
coeffH1 = polyfit(loghx, logH1, 1);
%提取coeffL2中的第一个元素，即拟合直线的斜率
slopeL2 = coeffL2(1);
slopeH1 = coeffH1(1);
%使用多项式求值函数polyval，根据拟合系数coeffL2和loghx计算拟合值，然后取指数得到拟合曲线的y值
fitL2 = exp(polyval(coeffL2, loghx));
fitH1 = exp(polyval(coeffH1, loghx));

%怎么老是说我错误使用loglog
loglog(length, fitL2, '--',  'DisplayName', sprintf('L2 Fit (Slope: %.2f)', slopeL2));
loglog(length, fitH1, '-.',  'DisplayName', sprintf('H1 Fit (Slope: %.2f)', slopeH1));

xlabel('diffierent hx');
ylabel('ch');
legend('show');
