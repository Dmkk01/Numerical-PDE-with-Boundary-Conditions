clear;

wave2D(100,100,10,10,10,0.5,1,2)

% nx_ - The number of data points in the X-direction
% ny_ - The number of data points in the Y-direction
% Lx_ - Length of X
% Ly_ - Length of Y
% T_ - Duration of the iterations
% CFL_ - The Courant constanst
% c_ - Constant in the wave equation
% boundary - if 1 -> reflective, else -> absorbing
function wave2D(nx_, ny_, Lx_, Ly_, T_, CFL_, c_, boundary)

    nx = nx_; % Number of points in X-direction
    ny = ny_; % Number of points in Y-direction
    Lx = Lx_; % Length of X
    Ly = Ly_; % Length of Y
    x = linspace(0, Lx, nx); % Row vector of nx number points between 0 and Lx
    y = linspace(0, Ly, ny); % Row vector of ny number points between 0 and Ly
    dx = x(2) - x(1); % Delta-x
    dy = y(2) - y(1); % Delta-y
    U = zeros(nx, ny); % Matrix of all the points
    Um1 = U; % U at time n-1
    Up1 = U; % U at time n+1

    % Time and CFL
    T = T_; 
    CFL = CFL_; 
    c = c_;
    dt = CFL * dx/c;
    t = 0;

    while (t < T)
        if (boundary == 1)
            % Reflecting Boundary condition
            U(:,[1 end]) = 0;
            U([1 end], :) = 0;
        else
            % Absorbing Boundary condition
            Up1(1,:)=U(2,:) + ((CFL-1)/(CFL+1))*(Up1(2,:)-U(1,:));
            Up1(end, :) = U(end-1,:) + ((CFL-1)/(CFL+1))*(Up1(end-1, :)-U(end,:));
            Up1(:,1)=U(:,2) + ((CFL-1)/(CFL+1))*(Up1(:,2)-U(:,1));
            Up1(:, end) = U(:,end-1) + ((CFL-1)/(CFL+1))*(Up1(:, end-1)-U(:,end));
        end

        t = t+dt;

        % Saving current and previous arrays
        Um1=U;
        U = Up1;

        %Source
        U(nx/2, ny/2) = dt^2*20*sin(30*pi*t/20);

        for i = 2:nx-1
            for j = 2:ny-1 
                H = (((U(i,j+1)) - 2*U(i,j) + (U(i,j-1)))/(dy^2));
                V = (((U(i+1,j)) - 2*U(i,j) + (U(i-1,j)))/(dx^2));
                Up1(i,j) =(c^2)*(dt^2)*(H+V) + 2*U(i,j) - Um1(i,j);
            end
        end

        clf;
        subplot(2,1,1);
        imagesc(x,y,U'); colorbar; caxis([-0.02 0.02])
        title(sprintf('t = %.2f', t));
        subplot(2,1,2);
        mesh(x,y,U'); colorbar; caxis([-0.02 0.02])
        axis([0 10 0 10 -0.03 0.03]);
        shg; pause(0.01);
    end
end
