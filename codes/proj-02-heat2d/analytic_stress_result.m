function [sigma_xx, sigma_yy, sigma_xy] = analytic_stress_result(Tx, R, x, y)
    % 极坐标转成xy坐标，计算 r 和 theta
    r = sqrt(x.^2 + y.^2);
    theta = atan2(y, x);
    
    % 两个三角函数 cos(2*theta) 和 sin(2*theta)专门写一下
    cos_2theta = (x.^2 - y.^2) ./ (r.^2);
    sin_2theta = (2 * x .* y) ./ (r.^2);
    
    % 计算解析解的应力分量
    sigma_rr = (Tx/2) * (1 - (R^2 ./ r.^2)) + ...
               (Tx/2) * (1 - 4*(R^2 ./ r.^2) + 3*(R^4 ./ r.^4)) .* cos_2theta;
    sigma_theta_theta = (Tx/2) * (1 + (R^2 ./ r.^2)) - ...
                        (Tx/2) * (1 + 3*(R^4 ./ r.^4)) .* cos_2theta;
    sigma_rtheta = -(Tx/2) * (1 + 2*(R^2 ./ r.^2) - 3*(R^4 ./ r.^4)) .* sin_2theta;
    
    % 转换为笛卡尔坐标下的应力分量
    sigma_xx = sigma_rr .* (x ./ r).^2 + sigma_theta_theta .* (y ./ r).^2 - ...
               sigma_rtheta .* 2 .* (x ./ r) .* (y ./ r);
    sigma_yy = sigma_rr .* (y ./ r).^2 + sigma_theta_theta .* (x ./ r).^2 + ...
               sigma_rtheta .* 2 .* (x ./ r) .* (y ./ r);
    sigma_xy = (sigma_theta_theta - sigma_rr) .* (x ./ r) .* (y ./ r) + ...
               sigma_rtheta .* ((x ./ r).^2 - (y ./ r).^2);
end