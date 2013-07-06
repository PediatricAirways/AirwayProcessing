function  [a_long, b_short, pts_center, pts_t] = getEllipse2DTo3D( score, coeff, xmean, ymean, zmean )  

a = fit_ellipse(score(:,1), score(:,2));
fh = @(x,y) a(1)*x^2 + a(2)*x*y + a(3)*y^2 + a(4)*x + a(5)*y + a(6);
%     figure(2), ezplot(fh, [min(score(:,1)), max(score(:,1)), min(score(:,2)), max(score(:,2))]);
%     axis equal;

% center
xc = (a(2)*a(5) - 2*a(3)*a(4)) / (4*a(1)*a(3) - a(2)^2);
yc = (a(2)*a(4) - 2*a(1)*a(5)) / (4*a(1)*a(3) - a(2)^2);
%     figure(2), plot(xc, yc, 'ro');

% angle
celta = 0.5 * atan(a(2)/(a(1)-a(3)));
x_line = linspace(min(score(:,1)), max(score(:,1)), 10);
y_line = yc + tan(celta)*(x_line-xc);
%     figure(2), plot(x_line, y_line, 'b-*');
y_line = linspace(min(score(:,2)), max(score(:,2)), 10);
x_line = xc + (y_line-yc)/tan(celta+3.1415926/2);
%     figure(2), plot(x_line, y_line, 'm-+');

% long and short axises 
a_long = sqrt(2*(a(1)*xc^2 + a(3)*yc^2 + a(2)*xc*yc - a(6))/(a(1) + a(3) + sqrt((a(1)-a(3))^2 + a(2)^2)));
b_short = sqrt(2*(a(1)*xc^2 + a(3)*yc^2 + a(2)*xc*yc - a(6))/(a(1) + a(3) - sqrt((a(1)-a(3))^2 + a(2)^2)));
if a_long < b_short
    tmp = a_long;
    a_long = b_short;   
    b_short = tmp;
end

% two points on long axis
xa1 = xc + cos(celta)*a_long;
ya1 = yc + sin(celta)*a_long;
%     figure(2), plot(xa1, ya1, 'bo');

xa2 = xc - cos(celta)*a_long;
ya2 = yc - sin(celta)*a_long;
%     figure(2), plot(xa2, ya2, 'bo');

% two points on short axis
xb1 = xc + cos(celta+3.1415926/2)*b_short;
yb1 = yc + sin(celta+3.1415926/2)*b_short;
%     figure(2), plot(xb1, yb1, 'ko');

xb2 = xc - cos(celta+3.1415926/2)*b_short;
yb2 = yc - sin(celta+3.1415926/2)*b_short;
%     figure(2), plot(xb2, yb2, 'ko');

% all points on ellipse
t = 0:0.01:2*3.1415926;
xt = xc + a_long*cos(t)*cos(celta) - b_short*sin(t)*sin(celta);
yt = yc + a_long*cos(t)*sin(celta) + b_short*sin(t)*cos(celta);
%     figure(2), plot(xt, yt, 'r-+');

pts_center = [xc yc 0];
pts_center = pts_center * inv(coeff) + [xmean ymean zmean];
%     figure(1), plot3(pts_center(1), pts_center(2), pts_center(3), 'ro');

pts_a1 = [xa1 ya1 0];
pts_a1 = pts_a1 * inv(coeff) + [xmean ymean zmean];
%     figure(1), plot3(pts_a1(1), pts_a1(2), pts_a1(3), 'go');

pts_a2 = [xa2 ya2 0];
pts_a2 = pts_a2 * inv(coeff) + [xmean ymean zmean];
%     figure(1), plot3(pts_a2(1), pts_a2(2), pts_a2(3), 'go');

pts_b1 = [xb1 yb1 0];
pts_b1 = pts_b1 * inv(coeff) + [xmean ymean zmean];
%     figure(1), plot3(pts_b1(1), pts_b1(2), pts_b1(3), 'ko');

pts_b2 = [xb2 yb2 0];
pts_b2 = pts_b2 * inv(coeff) + [xmean ymean zmean];
%     figure(1), plot3(pts_b2(1), pts_b2(2), pts_b2(3), 'ko');

pts_t(:,1) = xt';
pts_t(:,2) = yt';
pts_t(:,3) = 0;
pts_t = pts_t * inv(coeff);
pts_t(:,1) = pts_t(:,1) + xmean;
pts_t(:,2) = pts_t(:,2) + ymean;
pts_t(:,3) = pts_t(:,3) + zmean;
%     figure(1), plot3(pts_t(:,1), pts_t(:,2), pts_t(:,3), 'r-.');
   