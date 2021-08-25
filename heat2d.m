clear all; 
close all;
clc;
heat2D(100,100, 1, 1, 100000, 0.4 , 1.4, 1000, 800, 200, 500)

% nx_ - The number of data points in the X-direction
% ny_ - The number of data points in the Y-direction
% Lx_ - Length of X
% Ly_ - Length of Y
% nt_ - Number of time steps
% CFL_ - The Courant constanst
% alpha_ - Constant in the heat equation
% west_ - West boundary,number or "gradZero" for zero gradient boundary
% south_ - West boundary,number or "gradZero" for zero gradient boundary
% east_ - West boundary,number or "gradZero" for zero gradient boundary
% north_ - West boundary,number or "gradZero" for zero gradient boundary
function heat2D(nx_, ny_, Lx_, Ly_, nt_, CFL_, alpha_, west_, south_, east_, north_)
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size

    nx = nx_; % Number of points in X-direction
    ny = ny_; % Number of points in Y-direction
    Lx = Lx_; % Length of X
    Ly = Ly_; % Length of Y
    x = linspace(0, Lx, nx); % Row vector of nx number points between 0 and Lx
    y = linspace(0, Ly, ny); % Row vector of ny number points between 0 and Ly
    dx = x(2) - x(1); % Delta-x
    dy = y(2) - y(1); % Delta-y
    T = zeros(nx, ny); % Matrix of all the points

    % Defining boundary conditions
    % West wall
    if (string(west_) == "gradZero")
        T(:,1) = T(:,2);
    else
        T(:,1) = west_;
    end
    % South wall
    if (string(south_) == "gradZero")
        T(1,:) = T(2,:);
    else
        T(1,:) = south_; 
    end
    % East wall
    if (string(east_) == "gradZero")
        T(:, nx) = T(:, nx-1);
    else
        T(:, nx) = east_;
    end
    % North wall
    if (string(north_) == "gradZero")
        T(nx,:) = T(nx-1,:);
    else
        T(nx,:) = north_;
    end

    % Time steps and conditions
    alpha = alpha_;
    t = 1e-4;  %%%%
    cfl = CFL_;
    dt = (cfl*(dx^2))/(2 * alpha);
    nt = nt_; %%%%

    Told = T;
    error = 1;

    for k = 1:nt
        for j = 2:nx-1
            for i = 2:ny-1
                H = (((Told(i,j+1)) - 2*Told(i,j) + (Told(i,j-1)))/(dy^2));
                V = (((Told(i+1,j)) - 2*Told(i,j) + (Told(i-1,j)))/(dx^2));
                T(i,j) = Told(i,j) +  alpha*dt*(H+V);
            end
        end
        if error < 1e-6
            break
        end
        error = max(max(abs(Told-T))); %%%%
        Told = T; %%%%

        % Plot the graph every 25 time steps
        if mod(k, 50) == 0
            contourf(x,y,T,'ShowText', 'on')
            colorbar;
            xlabel('X')
            ylabel('Y')
            title(['Time step: ', num2str(k)]);
            drawnow
            frame = getframe(h);
            im = frame2im(frame);
            [~,~] = rgb2ind(im,256);
        end
    end
end
         