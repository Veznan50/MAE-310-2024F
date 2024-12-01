% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0
exact = @(x) -20*x.^3;
exact_x = @(x)-60*x.^2;

% Setup the mesh
pp   = 1;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
n_int = 10;

e_L2 = zeros(8, 1);%储存误差
e_H1 = zeros(8, 1);%储存误差

for ii = 1:8
    n_el = 2*ii;              % number of elements
    n_np = n_el * pp + 1;  % number of nodal points
    n_eq = n_np - 1;       % number of equations
    

    hh = 1.0 / (n_np - 1); % space between two adjacent nodes
    x_coor = 0 : hh : 1;   % nodal coordinates for equaiiy spaced nodes

%表格
IEN = zeros(n_el, n_en);

 for ee = 1 : n_el
   for aa = 1 : n_en
    IEN(ee, aa) = (ee - 1)*pp + aa;
   end
 end

% Setup the ID array for the problem
ID = 1 : n_np;
ID(end) = 0;

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

% alocate the stiffness matrix
K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
 F = zeros(n_eq, 1);
 
% Assembly of the stiffness matrix and load vector
for ee = 1 : n_el
  k_ele = zeros(n_en, n_en); 
  f_ele = zeros(n_en, 1);    

  x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(ee, aa)

  % quadrature loop
  for qua = 1 : n_int    
    dx_dxi = 0.0;
    x_l = 0.0;
    for aa = 1 : n_en
      x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
    end
    dxi_dx = 1.0 / dx_dxi;

    for aa = 1 : n_en
      f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
      for bb = 1 : n_en
        k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
      end
    end
  end
 
  % Assembly of the matrix and vector based on the ID or LM data
  
  for aa = 1 : n_en
    P = ID(IEN(ee,aa));
    if P > 0
      F(P) = F(P) + f_ele(aa);
      for bb = 1 : n_en
        Q = ID(IEN(ee,bb));
        if Q > 0
          K(P, Q) = K(P, Q) + k_ele(aa, bb);
        else
          F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
        end
      end
    end
  end
end

% ee = 1 F = NA(0)xh
F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

% Solve Kd = F equation
d_temp = K \ F;

disp = [d_temp; g];

 % e_L2 第一个误差
    numerators = 0.0;     %分子
    denominators = 0.0;   %分母
    for ee = 1 : n_el
        for qua = 1 : n_int   
            dx_dxi = 0.0;
            numerator_1 = 0.0;%数值
            numerator_2 = 0.0;%精确
            denominator = 0.0;
            for aa = 1 : n_en
                dx_dxi = dx_dxi + x_coor(IEN(ee, aa)) * PolyShape(pp, aa, xi(qua), 1);
                numerator_1 = numerator_1 + disp(IEN(ee, aa)) * PolyShape(pp, aa, xi(qua), 0);
                numerator_2 = numerator_2 + x_coor(IEN(ee, aa)) * PolyShape(pp, aa, xi(qua), 0);
                denominator = denominator + x_coor(IEN(ee, aa)) * PolyShape(pp, aa, xi(qua), 0);
            end
            numerator = (numerator_1 - exact(numerator_2))^2 * dx_dxi * weight(qua);
            numerators = numerators + numerator;
            denominator = exact(denominator)^2 * dx_dxi * weight(qua);
            denominators =denominators + denominator;
        end
    end
    e_L2(ii) = sqrt(numerators/denominators);
   

    % e_H1 第二个误差
    numerators = 0.0;
    denominators = 0.0;
    for ee = 1 : n_el
        for qua = 1 : n_int
            dx_dxi = 0.0;
            numerator_1 = 0.0;
            numerator_2 = 0.0;
            denominator = 0.0;
            for aa = 1 : n_en
                dx_dxi = dx_dxi + x_coor(IEN(ee, aa)) * PolyShape(pp, aa, xi(qua), 1);
                numerator_1 = numerator_1 + disp(IEN(ee, aa)) * PolyShape(pp, aa, xi(qua), 1);
                numerator_2 = numerator_2 + x_coor(IEN(ee, aa)) * PolyShape(pp, aa, xi(qua), 0);
                denominator = denominator + x_coor(IEN(ee, aa)) * PolyShape(pp, aa, xi(qua), 0);
            end
            dxi_dx = 1.0/ dx_dxi;
            numerator = (numerator_1 * dxi_dx - exact_x(numerator_2))^2 * dx_dxi * weight(qua);
            numerators = numerators + numerator;
            denominator = (exact_x(denominator))^2 * dx_dxi * weight(qua);
            denominators =denominators + denominator;
        end
    end
    e_H1(ii) = sqrt(numerators/denominators);
end


x_hh = 1./ (2: 2: 16);
plot(log(x_hh), log(e_L2), '-r','LineWidth',3)

hold on

plot(log(x_hh), log(e_H1), '-k','LineWidth',3)

slope_e_L2 = (log(e_L2(end))-log(e_L2(1)))/(log(x_hh(end))-log(x_hh(1)));

slope_e_H1 = (log(e_H1(end))-log(e_H1(1)))/(log(x_hh(end))-log(x_hh(1)));


