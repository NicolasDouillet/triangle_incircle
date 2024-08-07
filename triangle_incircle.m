function [R, I, r] = triangle_incircle(A, B, C, nb_samples, option_display)
%% triangle_incircle : function to compute and display the incircle of one given triangle
%
% Author : nicolas.douillet (at) free.fr, 2022-2024.
%
% Syntax
%
% triangle_incircle(A, B, C);
% triangle_incircle(A, B, C, nb_samples);
% triangle_incircle(A, B, C, nb_samples, option_display);
% [R, I, r] = triangle_incircle(A, B, C, nb_samples, option_display);
%
% Description
%
% triangle_incircle(A, B, C) computes and displays the incircle of ABC triangle.
% triangle_incircle(A, B, C, nb_samples) uses nb_samples to draw the circle.
% triangle_incircle(A, B, C, nb_samples, option_display) displays the circle when option_display is set either to
% logical true or real numeric 1, and doesn't when it is set to logical false or real numeric 0.
% [R, I, r] = triangle_incircle(A, B, C, nb_samples, option_display) stores the results in [R, I, r] vector.
%
% See also INCENTER CIRCUMCENTER
%
% Input arguments
%
%       [Ax]
% - A = [Ay] : real column vector double. 2 <= numel(A) <= 3. One of the three ABC vertices. 
%       [Az]
%
%       [Bx]
% - B = [By] real column vector double. 2 <= numel(A) <= 3. One of the three ABC vertices.
%       [Bz]
%
%       [Cx]
% - C = [Cy] real column vector double. 2 <= numel(A) <= 3. One of the three ABC vertices.
%       [Cz]
%
% - nb_samples : integer scalar double. The number of samples to draw the
%                incircle. nb_samples >= 3.
%
% - option_display : logical *true(1) / false(0), to enable/disable the display mode.
%
%
% Output arguments
%
%       [- Rx -]
% - R = [- Ry -]: real matrix double. The incircle coordinates. size(R) = [size(A,1), nb_samples].
%       [- Rz -]
%
%       [Ix]
% - I = [Iy] : real column vector double. 2 <= numel(I) <= 3. The incircle centre.
%       [Iz]
%
% - r : real scalar double. the incircle radius.
%
%
% Example
% From a triangle of the 3D space
% A = 2*(rand(3,1)-0.5);
% B = 2*(rand(3,1)-0.5);
% C = 2*(rand(3,1)-0.5);
% nb_samples = 30;
% option_display = true;
% triangle_incircle(A,B,C,nb_samples,option_display);
%
% Example #2
% From a triangle of the 2D space
% A = cat(1,2*(rand(2,1)-0.5),0);
% B = cat(1,2*(rand(2,1)-0.5),0);
% C = cat(1,2*(rand(2,1)-0.5),0);
% [R,I,r] = triangle_incircle(A,B,C);


% Input parsing
assert(nargin > 2, 'Not enought input arguments. Three points required to define one triangle.');
assert(nargin < 6, 'Too many input arguments.');

if nargin < 5
    
    option_display = true;
    
    if nargin < 4
       
        nb_samples = 60;
        
    end
    
end

assert(isequal(size(A),size(B),size(C)),'All inputs points must have the same size.');
assert(isequal(ndims(A),ndims(B),ndims(C),2),'All inputs points must have the same number of dimensions (2).');
assert(isreal(A) && isreal(B) && isreal(C),'All inputs points must contain real numbers only.');
assert(numel(A) > 1 && numel(A) < 4,'Input points must have 2 or 3 elements.');


%% Body

% ABC triangle director and normal vectors computation 
AB = (B-A)/norm(B-A);
AC = (C-A)/norm(C-A);
BC = (C-B)/norm(C-B);
n = cross(AB,AC);

% BAC, ABC, and ACB angles bissectrices computation
A_midangle_vect = 0.5*(AB+AC);
B_midangle_vect = 0.5*(BC-AB);

% Circle centre I computation
I = lines_intersection(A,A_midangle_vect,B,B_midangle_vect)';

% Circle radius r computation
r = point_to_line_distance(I,AB',A');

% Compute the circle point coordinates
angle_step = 2*pi/nb_samples;
theta = linspace(0,2*pi-angle_step,nb_samples);

Cx = r*cos(theta);
Cy = r*sin(theta);
Cz = zeros(1,nb_samples);


% Vector u to rotate around
k = [0 0 1]';


if norm(cross(k,n)) > 1e3*eps
    
    u = cross(k,n)/norm(cross(k,n));
    
    % Angle between k and u
    alpha = atan2(norm(cross(k,n)),dot(k,n));
    
    % 3D rotation matrix around u vector
    Rm = @(delta)[u(1,1)^2+cos(delta)*(1-u(1,1)^2) (1-cos(delta))*u(1,1)*u(2,1)-u(3,1)*sin(delta) (1-cos(delta))*u(1,1)*u(3,1)+u(2,1)*sin(delta);
                  (1-cos(delta))*u(1,1)*u(2,1)+u(3,1)*sin(delta) u(2,1)^2+cos(delta)*(1-u(2,1)^2) (1-cos(delta))*u(2,1)*u(3,1)-u(1,1)*sin(delta);
                  (1-cos(delta))*u(1,1)*u(3,1)-u(2,1)*sin(delta) (1-cos(delta))*u(2,1)*u(3,1)+u(1,1)*sin(delta) u(3,1)^2+cos(delta)*(1-u(3,1)^2)];
    
    R = (Rm(alpha) * cat(1,Cx,Cy,Cz))' + I;
    
else
    
    R = cat(1,Cx,Cy,Cz)' + I;
    
end


%% Display
if option_display
    
    figure
    
    line([A(1,1) B(1,1) C(1,1) A(1,1)],[A(2,1) B(2,1) C(2,1) A(2,1)],[A(3,1) B(3,1) C(3,1) A(3,1)],'Color',[1 0 0],'Linewidth',2), hold on;
    line([R(:,1); R(1,1)],[R(:,2); R(1,2)],[R(:,3); R(1,3)],'Color',[0 0 1],'Linewidth',2), hold on;
    view(3);
    
    axis equal, axis tight;
    
end


end % triangle_incircle


%% lines_intersection subfunction
function [I, rc] = lines_intersection(M1, u1, M2, u2)
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2023.


precision = 1e3*eps;
v = cross(u1,u2);
diff_pts = M2-M1;
catdim = find(size(M1)==1);
n = numel(M1);

% Segment cases
if norm(v) < precision
    
    if norm(cross(u1,diff_pts)) < precision && norm(cross(u2,diff_pts)) < precision
        
        I = M1;
        rc = 2;        
        
    else
        
        I = [];
        rc = 0;        
        
    end
    
else
    
    d = [-v(1) v(2) -v(3)];
    f = find(abs(d) > precision,1);
    
    d_pts = diff_pts(setdiff(1:n,f));    
    dt = det(cat(catdim,d_pts,-u2(setdiff(1:n,f))));
    du = det(cat(catdim,u1(setdiff(1:n,f)),d_pts));
    
    t = dt/d(f);
    u = du/d(f);
    
    if abs(M1+u1*t-M2-u2*u) < precision
        
        I = M1 + u1*t;        
        rc = 1;
        
    else
        
        I = [];
        rc = 0;
        
    end
    
end


end % lines_intersection


%% point_to_line_distance subfunction
function [d2H, H] = point_to_line_distance(P, u, I0)
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2022.


% Body
nb_pts = size(P,1);


t_H = (u(1)*(P(:,1)-repmat(I0(1),[nb_pts,1])) + ...
       u(2)*(P(:,2)-repmat(I0(2),[nb_pts,1])) + ...
       u(3)*(P(:,3)-repmat(I0(3),[nb_pts,1])) ) / ...
       sum(u.^2); 

x_H = I0(1) + t_H*u(1);
y_H = I0(2) + t_H*u(2);
z_H = I0(3) + t_H*u(3);


% Orthogonal projected point
H = zeros(size(P));
H(:,1) = x_H;
H(:,2) = y_H;
H(:,3) = z_H;


% Distance
d2H = vecnorm((P-H)',2)';
H = H(:,1:size(P,2));


end % point_to_line_distance       