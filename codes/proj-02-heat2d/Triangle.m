function val = Triangle(aa, xi, eta)

if aa == 1
    val = 1 - xi - eta;
    elseif aa == 2
    val = xi;
    elseif aa == 3
    val = eta;
    else
    error('Error: value of a should be 1, 2,or 3.');
end
% 好像一个三角形就ok了
%elseif aa == 4
   % val = 1 - xi - eta; % 节点4的形函数，对应于三角形2
%elseif aa == 5
    %val = eta; % 节点5的形函数，对应于三角形2
%elseif aa == 6
    %val = xi; % 节点6的形函数，对应于三角形2
%else
    %error('Error: value of a should be 1, 2, 3, 4, 5, or 6.');