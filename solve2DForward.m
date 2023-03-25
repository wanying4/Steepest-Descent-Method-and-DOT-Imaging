function [M,u] = solve2DForward(N,mua,source_loc)
%SOLVE2DFORWARD solves the 2D forward model of the diffusion eqn using the 
%finite volume method. 
%   Inputs:
%       N: number of grid points
%       mua: N^2-by-1 vector of the mua profile
%       source_loc: string that tells where the source is located
%         'W' = west
%         'E' = east
%         'N' = south
%         'S' = north
%   Outputs:
%       M: the matrix in Mu = b
%       u: a fluence vector with dimension N^2-by-1, it needs to be 
%       reshapedto get the orignal fluence profile: reshape(u,N,N).

global L x dx dy Reff qin itermax tol Q

% construct control volume (C.V.) and calculate other necessary parameters
D = 1/15;               % assume D is constant over the medium
Nt = N*N;               % total number of grid points
delx = dx;              % distance between two adjacent nodes in x-dir
dely = dy;              % distance between two adjacent nodes in y-dir
A = (1+Reff)/(1-Reff);  % reflection correction factor A in Robin B.C.

% initialize variables
M = zeros(Nt,Nt);       % matrix M in Mu=b
b = zeros(Nt,1);        % right-hand side vector b in Mu=b
u = zeros(Nt,1);        % solution vector u
sp = zeros(Nt,1);       % internal source s
qin = zeros(Nt,1);      % boundary source qin [W/cm^2] (e.g., laser illumination)

% place light sources to nodes accordingly
switch source_loc
    case 'W'
        j=1;i=(N-1)/2;          % light source located in the middle of west side
        %sp(1) = 1;
    case 'E'
        j=N;i=(N-1)/2;
    case 'N'
        j=(N-1)/2;i=1;
    case 'S'
        j=(N-1)/2;i=N;
end
qin((j-1)*N+i) = 1.0;   % qin = 1 W/cm^2 directed into the medium

% construct matrix M and vector b in Mu=b according to fomulations 
% in lecture note: ap*Up = ae*Ue+aw*Uw+as*Us+an*Un+b.
for i=1:N
    for j=1:N
        row = (j-1)*N+i; % row index of matrix M       
        % interior node
        if (i>1 && j>1 && i<N && j<N)            
            ae = D*dy/delx;
            aw = D*dy/delx;
            as = D*dx/dely;
            an = D*dx/dely;
            ap = ae + aw + as + an + mua(sub2ind([N,N],i,j))*dx*dy;
            
            col = row;
            M(row,col) = ap; 
            M(row,col+1) = -ae;
            M(row,col-1) = -aw;
            M(row,col+N) = -as;
            M(row,col-N) = -an;            
            b(row) = sp(row)*dx*dy;        
        % east boundary:
        elseif (i==N && j>1 && j<N)
            aw =  D*dy/delx;           
            as = D*(dx/2)/dely;
            an = D*(dx/2)/dely;
            ap = dy/(2*A)+ aw + as + an + mua(sub2ind([N,N],i,j))*(dx/2)*dy;
                        
            col = row;
            M(row,col) = ap; 
            M(row,col-1) = -aw;
            M(row,col+N) = -as;
            M(row,col-N) = -an;            
            b(row) = sp(row)*(dx/2)*dy+4*qin(row)*dy/(2*A);
        % west boundary node
        elseif (i==1 && j>1 && j<N)
            ae = D*dy/delx;
            as = D*(dx/2)/dely;
            an = D*(dx/2)/dely;
            ap = ae + dy/(2*A)+ as + an + mua(sub2ind([N,N],i,j))*(dx/2)*dy;
         
            col = row;
            M(row,col) = ap; 
            M(row,col+1) = -ae;
            M(row,col+N) = -as;
            M(row,col-N) = -an;            
            b(row) = sp(row)*(dx/2)*dy+4*qin(row)*dy/(2*A);                  
        % south boundary node
        elseif (j==N && i>1 && i<N)
            ae = D*(dy/2)/delx;
            aw =  D*(dy/2)/delx;           
            an = D*dx/dely;
            ap = ae + aw + dx/(2*A) + an + mua(sub2ind([N,N],i,j))*dx*dy;
            
            col = row;
            M(row,col) = ap; 
            M(row,col+1) = -ae;
            M(row,col-1) = -aw;
            M(row,col-N) = -an;            
            b(row) = sp(row)*dx*(dy/2)+4*qin(row)*dy/(2*A);        
        % north boundary node
        elseif (j==1 && i>1 && i<N)
            ae = D*(dy/2)/delx;
            aw =  D*(dy/2)/delx;           
            as = D*dx/dely;
            ap = ae + aw + as + dy/(2*A) + mua(sub2ind([N,N],i,j))*dx*(dy/2);
            
            col = row;
            M(row,col) = ap; 
            M(row,col+1) = -ae;
            M(row,col-1) = -aw;
            M(row,col+N) = -as;            
            b(row) = sp(row)*dx*(dy/2)+4*qin(row)*dx/(2*A);            
        % corner nodes
        elseif i==1 && j==1
            ae = D*(dy/2)/delx;
            as = D*(dx/2)/dely;
            ap = ae + (dy/2)/(2*A)+ as + (dx/2)/(2*A) + mua(sub2ind([N,N],i,j))*(dx/2)*(dy/2);
         
            col = row;
            M(row,col) = ap;
            M(row,col+1) = -ae;
            M(row,col+N) = -as;
            b(row) = sp(row)*(dx/2)*(dy/2)+4*qin(row)*(dy/2)/(2*A)...
                                          +4*qin(row)*(dx/2)/(2*A);
        elseif j==1 && i==N
            aw = D*(dy/2)/delx;
            as = D*(dx/2)/dely;
            ap = (dy/2)/(2*A)+ aw + as + (dx/2)/(2*A) + mua(sub2ind([N,N],i,j))*(dx/2)*(dy/2);
         
            col = row;
            M(row,col) = ap;
            M(row,col-1) = -aw;
            M(row,col+N) = -as;
            b(row) = sp(row)*(dx/2)*(dy/2)+4*qin(row)*(dy/2)/(2*A)...
                                          +4*qin(row)*(dx/2)/(2*A);
        elseif j==N && i==1
            ae = D*(dy/2)/delx;
            an = D*(dx/2)/dely;
            ap = ae + (dy/2)/(2*A)+ (dx/2)/(2*A) + an + mua(sub2ind([N,N],i,j))*(dx/2)*(dy/2);
         
            col = row;
            M(row,col) = ap;
            M(row,col+1) = -ae;
            M(row,col-N) = -an;
            b(row) = sp(row)*(dx/2)*(dy/2)+4*qin(row)*(dy/2)/(2*A)...
                                          +4*qin(row)*(dx/2)/(2*A);
        elseif j==N && i==N
            aw = D*(dy/2)/delx;
            an = D*(dx/2)/dely;
            ap = (dy/2)/(2*A) + aw + (dx/2)/(2*A) + an + mua(sub2ind([N,N],i,j))*(dx/2)*(dy/2);
         
            col = row;
            M(row,col) = ap;
            M(row,col-1) = -aw;
            M(row,col-N) = -an;
            b(row) = sp(row)*(dx/2)*(dy/2)+4*qin(row)*(dy/2)/(2*A)...
                                          +4*qin(row)*(dx/2)/(2*A);
        end
    end
end

% solve Ax = b to get fluence
u = M\b;

end

