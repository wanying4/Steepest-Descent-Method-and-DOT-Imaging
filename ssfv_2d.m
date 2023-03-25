% Solve forward and inverse problems of 2D FV diffusion eqn

clear; close all;clear global;

% common parameters
global L x y dx dy Reff qin itermax tol Q

% system parameters and problem setup
probType = 'inv';       % problem type set to either 'fwd' or 'inv'
filename = 'meas.txt';  % detector readings as input to inverse code
L = 6;                  % cm, thickness of slab
dx = 0.2;               % C.V. size in x-dir
dy = dx;                % C.V. size in y-dir
x = 0.0:dx:L;           % x-coord
y = 0.0:dy:L;           % y-coord
N = size(x,2);          % number of grid points

qin = 1;                % laser power, W/cm^2, as light source on west 
nindex = 1.414;         % refractive index of medium
Reff = -1.4399/nindex^2 + 0.7099/nindex + 0.6681 + 0.0636*nindex; 

tol = 1.e-4;            % tol for inverse problem
itermax = 100000;       % max. iter num for inverse problem

% detector locations and constants for partial fluence
Q = zeros(16,N*N);
Q(1,sub2ind([N,N],1/dx+1,1)) = 0.5*((1-Reff)/(1+Reff));
Q(2,sub2ind([N,N],2/dx+1,1)) = 0.5*((1-Reff)/(1+Reff));
Q(3,sub2ind([N,N],4/dx+1,1)) = 0.5*((1-Reff)/(1+Reff));
Q(4,sub2ind([N,N],5/dx+1,1)) = 0.5*((1-Reff)/(1+Reff));
Q(5,sub2ind([N,N],1,1/dx+1)) = 0.5*((1-Reff)/(1+Reff));
Q(6,sub2ind([N,N],N,1/dx+1)) = 0.5*((1-Reff)/(1+Reff));
Q(7,sub2ind([N,N],1,2/dx+1)) = 0.5*((1-Reff)/(1+Reff));
Q(8,sub2ind([N,N],N,2/dx+1)) = 0.5*((1-Reff)/(1+Reff));
Q(9,sub2ind([N,N],1,4/dx+1)) = 0.5*((1-Reff)/(1+Reff));
Q(10,sub2ind([N,N],N,4/dx+1)) = 0.5*((1-Reff)/(1+Reff));
Q(11,sub2ind([N,N],1,5/dx+1)) = 0.5*((1-Reff)/(1+Reff));
Q(12,sub2ind([N,N],N,5/dx+1)) = 0.5*((1-Reff)/(1+Reff));
Q(13,sub2ind([N,N],1/dx+1,N)) = 0.5*((1-Reff)/(1+Reff));
Q(14,sub2ind([N,N],2/dx+1,N)) = 0.5*((1-Reff)/(1+Reff));
Q(15,sub2ind([N,N],4/dx+1,N)) = 0.5*((1-Reff)/(1+Reff));
Q(16,sub2ind([N,N],5/dx+1,N)) = 0.5*((1-Reff)/(1+Reff));

% switch 
switch probType
    case 'fwd'        
        % small square obejct center at 4cm*4cm of saize 1cm*1cm
        sqr_obj_center_x = find(x==4);
        sqr_obj_center_y = find(y==4);
        sqr_ojb_x_ind = sqr_obj_center_x-floor(.5/dx):sqr_obj_center_x+floor(.5/dx);
        sqr_ojb_y_ind = sqr_obj_center_y-floor(.5/dx):sqr_obj_center_y+floor(.5/dx);

        % optical properties - matrix of mua
        mua = 0.05*ones(N,N);                        % background absorption = 0.05 cm-1
        mua(sqr_ojb_x_ind,sqr_ojb_y_ind) = 0.2;    % square object absorption = 0.2 cm-1
        mua = flipud(mua);
        
        % write detector readings
        pred = solve2DForward_all(N,mua(:));
        fid = fopen(filename,'w');
        fprintf(fid,'%f\n',pred);
        fclose(fid);
        
    case 'inv' 
        % read measurements
        fid = fopen(filename,'r');
        meas = fscanf(fid,'%f');
        fclose(fid);
        
        % initial guess from evolution_strategy.m
        mua(1:N*N) = 0.054916567129257;
        
        % solve inverse problem
        [objFunc_vec,mua_r] = solve2DInverse(N,mua,meas);
        
    otherwise
        warning('unexpected problem type');
end

