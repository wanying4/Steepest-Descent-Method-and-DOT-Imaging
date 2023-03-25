function [objFunc_vec,mua_r] = solve2DInverse(N,mua,meas)
% solve inverse problem with steepest descent method
%SOLVE2DINVERSE solves the 2D inverse model of the diffusion eqn using 
%the finite volume method. The reconstruction is done with the steepest
%descent method.
%   Inputs:
%       N: number of grid points
%       mua: N^2-by-1 vector of the mua profile
%       meas: 64-by-1 measurement data
%   Outputs:
%       objFunc_vec: the objective function values of all iterations
%       u is a fluence vector with dimension N^2-by-1, it needs to be reshaped
%       to get the orignal fluence profile: reshape(u,N,N).
%       mua_r: the solution of the mua profile that yields the lowest
%       objective function

% common parameters
global L x y dx dy Reff qin itermax tol Q

% function evaluation with initial guess
iter = 0;

% store the objFunc from every iteration
objFunc_vec = [];

% save progess data to file
filename='output.mat';
m = matfile(filename, 'Writable', true);

% start finding the mua profile
while iter < itermax
    iter = iter + 1;
    
    % solve adjoint problem
    source_locs = ['W' 'E' 'N' 'S'];
    grad = zeros(N*N,1);
    for i=1:length(source_locs)
        source_loc = source_locs(i);
        [M,u] = solve2DForward(N,mua,source_loc);    % solve forward to get fluence
        P = Q*u;                                     % partial currents
        
        % define volume size for each fluence position
        dAw = [(dx/2)*(dy/2);ones(N-2,1)*(dx/2)*dy;(dx/2)*(dy/2)];
        dAm = [dx*(dy/2);ones(N-2,1)*dx*dy;dx*(dy/2)];
        dAe = dAw;
        dA = [dAw;repmat(dAm,N-2,1);dAe];
        
        % calculate gradient
        eta = inv(M')*Q'*(P./meas((1+(i-1)*16):16*i)-1);
        grad = grad - (u.*eta).*dA;
    end
    
    % initialize parameters for line search
    pred = solve2DForward_all(N,mua);
    objFunc = 0.5*sum((pred./meas-1).^2);
    t = 20;
    mua_temp = mua-(t*grad)';
    pred = solve2DForward_all(N,mua_temp);
    objFunc_new = 0.5*sum((pred./meas-1).^2);
    
    % perform line search
    while objFunc_new > objFunc-(1.0e-4)*t*grad*grad'
        t = 0.5*t;
        mua_temp = mua-(t*grad)';
        pred = solve2DForward_all(N,mua_temp);
        objFunc_new = 0.5*sum((pred./meas-1).^2);
    end
    
    % update
    mua = mua-(t*grad)';
    objFunc_vec(iter) = objFunc_new;
    
    % save progress data to file
    m.muaprofile(iter,1:N*N) = mua;
    m.gradient(iter,1:N*N) = grad';
    m.objFunc(iter,1) = objFunc_new;

    % display progress
    fprintf('iter= %d: f(x) = %e and step size t = %f\n',iter,objFunc,t);
    if ~mod(iter,10)
        FigHandle = figure('Position', [100, 100, 1300, 500]);
        subplot(121);imagesc([dx dy*N],[dy dx*N],reshape(mua,N,N));
        title(sprintf('\\mu_a profile at the %dth iteration',iter));
        xlabel('x direction (cm)');ylabel('y direction (cm)');
        colorbar;colormap(jet);drawnow;
        
        subplot(122);imagesc([dx dy*N],[dy dx*N],reshape(grad,N,N));
        title(sprintf('Gradient at the %dth iteration',iter));
        xlabel('x direction (cm)');ylabel('y direction (cm)');
        colorbar;colormap(jet);drawnow;
    end

end
mua_r = mua;

end