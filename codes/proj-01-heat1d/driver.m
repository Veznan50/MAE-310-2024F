% Problem definition
f = @(x) -20*x.^3; % f(x) is the source精确解
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
pp   = 1;              % polynomial degree二次多项式，还要2和3
n_en = pp + 1;         % number of element or local nodes局部map的节点数
n_el = 2;              % number of elements单元数为2，还要改4,6,8,10,12,14,16
n_np = n_el * pp + 1;  % number of nodal points总的节点数
n_eq = n_np - 1;       % number of equations方程数
n_int = 10;            %积分点数量

hh = 1.0 / (n_np - 1); % space between two adjacent nodes单元长度，均匀分布
x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes全局节点的坐标

IEN = zeros(n_el, n_en);%初始化了一个大小为 n_el 行 n_en 列的零矩阵，用于存储每个元素连接的节点编号。

for ee = 1 : n_el      %遍历所有的元素（ee代表当前的元素编号） 
  for aa = 1 : n_en    %遍历该元素的所有节点（aa 代表当前节点的编号）。
    IEN(ee, aa) = (ee - 1) * pp + aa;%计算当前元素的当前节点的全局编号，并将其存储在 IEN 矩阵中。
                                     % 节点的全局编号是通过将元素编号减1（因为元素编号从1开始，而不是0），
                                     % 乘以每个元素的节点数（pp），
                                     % 然后加上当前节点的编号（aa）来计算的。
  end
end

% Setup the ID array for the problem
ID = 1 : n_np;
ID(end) = 0;

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);%高斯积分的节点 xi 和权重 weight。
             % 这里，n_int 是积分点的数量，-1 和 1 分别是积分区间的下限和上限。

% allocate the stiffness matrix
K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);%非零元素的数量估计为 (2*pp+1)*n_eq。
F = zeros(n_eq, 1);%初始化了一个大小为 n_eq 乘以1的零向量 F，用于存储载荷向量。

% Assembly of the stiffness matrix and load vector
for ee = 1 : n_el            %遍历所有的元素
  k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix用于存储元素的局部刚度矩阵。
  f_ele = zeros(n_en, 1);    % allocate a zero element load vector用于存储元素的局部载荷向量。

  x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)从全局坐标数组 x_coor 中提取当前元素的节点坐标

  % quadrature loop
  for qua = 1 : n_int    %遍历所有的积分点（qua 代表当前的积分点编号）
    dx_dxi = 0.0;        %初始化变量 dx_dxi（雅可比行列式）和 x_l（参考坐标 xi 转换到物理坐标 x）。
    x_l = 0.0;           %
    for aa = 1 : n_en    %开始了一个循环，遍历当前元素的所有节点。
      x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
      %计算当前积分点的物理坐标 x_l，通过求和每个节点坐标与相应的形状函数值的乘积。
      %          调用function val = PolyShape(degree, a, xi, der)
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
      %计算雅可比行列式 dx_dxi，通过求和每个节点坐标与相应的形状函数导数的乘积。
    end
    dxi_dx = 1.0 / dx_dxi;%计算雅可比行列式的倒数 dxi_dx。

    for aa = 1 : n_en%开始了另一个循环，遍历当前元素的所有节点。1,2,3差不多
      f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
      %计算元素载荷向量 f_ele 的第 aa 个分量，通过求和每个积分点处形状函数值、
      % 源项函数值 f(x_l)、权重 weight(qua) 和雅可比行列式 dx_dxi 的乘积。
      for bb = 1 : n_en%开始了一个循环，遍历当前元素的所有节点。
        k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
        %计算元素刚度矩阵 k_ele 的第 aa,bb 个分量，
        %通过求和每个积分点处两个形状函数导数的乘积、权重 weight(qua) 和雅可比行列式倒数 dxi_dx 的乘积。
      end
    end
  end
 %形状函数及其导数由 PolyShape 函数计算，该函数接受多项式度数 pp、节点编号 aa 或 bb、参考坐标 xi(qua) 和导数标志
 %（0表示形状函数，1表示导数）。
  % Assembly of the matrix and vector based on the ID or LM data
  for aa = 1 : n_en%遍历当前元素内部的所有节点
    P = ID(IEN(ee,aa));%从 ID 数组中获取当前元素的第 aa 个节点的全局自由度编号。
    if(P > 0)%检查全局自由度编号 P 是否大于0。如果是，这意味着该自由度不是边界条件
      F(P) = F(P) + f_ele(aa);%将元素载荷向量的第 aa 个分量加到全局载荷向量的第 P 个分量。
      for bb = 1 : n_en%%遍历当前元素的所有节点
        Q = ID(IEN(ee,bb));%从 ID 数组中获取当前元素的第 bb 个节点的全局自由度编号。
        if(Q > 0)%检查全局自由度编号 Q 是否大于0。如果是，这意味着该自由度不是边界条件
          K(P, Q) = K(P, Q) + k_ele(aa, bb);%将元素刚度矩阵的第 aa,bb 个分量加到全局刚度矩阵的第 P,Q 个分量。
        else
          F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
          %将元素刚度矩阵的第 aa,bb 个分量乘以边界条件值 g 并从全局载荷向量的第 P 个分量中减去。这处理了狄利克雷边界条件。
        end
      end
    end
  end
end

% ee = 1 F = NA(0)xh
F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;
%将第一个元素的第一个节点的载荷向量分量加上边界条件值 h。这里，ID(IEN(1,1)) 是第一个元素的第一个节点的全局自由度编号。

% Solve Kd = F equation
d_temp = K \ F;
%其中 K 是全局刚度矩阵，d 是未知自由度向量，F 是全局载荷向量。
%使用 MATLAB 的左除运算符 \ 解方程组 Kd = F，得到未知自由度向量 d_temp。
disp = [d_temp; g];
%将解向量 d_temp 和边界条件值 g 组合成一个新向量 disp。
% Postprocessing: visualization
plot(x_coor, disp, '--r','LineWidth',3);
%绘制节点坐标 x_coor 对应的位移 disp，使用红色虚线，线宽为3。
x_sam = 0 : 0.01 : 1;%将会创建一个从0到1的采样点数组，步长为0.01。
y_sam = x_sam.^5;%将会计算采样点的精确解（假设精确解是 x^5）。
hold on;
plot(x_sam, y_sam, '-k', 'LineWidth', 3);
%绘制精确解曲线，默认使用黑色实线，线宽为3。


n_sam = 20;%定义了一个变量 n_sam，表示在每个元素中采样点的数量，注意，是数量
xi_sam = -1 : (2/n_sam) : 1;%创建了一个从-1到1的采样点数组，用于参考坐标系中的采样。注意这里是数组

x_sam = zeros(n_el * n_sam + 1, 1);%初始化了一个零向量 x_sam，用于存储采样点的物理坐标。
y_sam = x_sam; % store the exact solution value at sampling points
% 初始化了一个零向量 y_sam，用于存储采样点处的精确解值。
u_sam = x_sam; % store the numerical solution value at sampling pts
%初始化了一个零向量 u_sam，用于存储采样点处的数值解值。

for ee = 1 : n_el%每个元素都有自己的编号。变量 ee 用于遍历这些元素，从1到n_el（元素总数）。
  x_ele = x_coor( IEN(ee, :) );%从全局坐标数组 x_coor 中提取当前元素的节点坐标。
  u_ele = disp( IEN(ee, :) );%从位移数组 disp 中提取当前元素的节点位移。

  if ee == n_el% 检查当前元素是否是最后一个元素。
    n_sam_end = n_sam+1;%如果是最后一个元素，采样点的数量增加1，以包括元素的端点。
  else                 %如果不是最后一个元素。
    n_sam_end = n_sam;%采样点的数量保持为 n_sam。
  end

  for ll = 1 : n_sam_end%遍历当前元素的所有采样点。
    x_l = 0.0;%物理坐标
    u_l = 0.0;%位移
    for aa = 1 : n_en%
      x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
      %通过求和每个节点坐标与相应的形状函数值的乘积，计算当前采样点的物理坐标 x_l。
      u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
      %通过求和每个节点位移与相应的形状函数值的乘积，计算当前采样点的位移 u_l。
    end

    x_sam( (ee-1)*n_sam + ll ) = x_l;%将计算出的物理坐标 x_l 存储在采样点数组 x_sam 中。
    u_sam( (ee-1)*n_sam + ll ) = u_l;%将计算出的位移 u_l 存储在采样点数组 u_sam 中。
    y_sam( (ee-1)*n_sam + ll ) = x_l^5;%将计算出的精确解 x_l^5 存储在采样点数组 y_sam 中
  end
end


plot(x_sam, u_sam, '-r','LineWidth',3);
hold on;
plot(x_sam, y_sam, '-k','LineWidth',3);
