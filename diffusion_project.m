clear all; close all; clc
%user input
N=10;
t=0;
ERR = 10^(-3);
%domain
dx = 2*pi/(N-1);    %dist between points
dy = 2*pi/(N-1);    %dist between points
dt = dx^2/4;        
lambda = .5*dt/dx^2;%eigenvalue

ax = 0; bx = 2*pi;      %x bounds
ay = 0; by = 2*pi;      %y bounds

x = ax:dx:bx;       %x grid
y = ay:dy:by;       %y grid

%dirichlet bc
fb = (by-y).^2.*cos(pi*y/by); %given boundary equation 1 
gb = (by-y).^2.*y;            %given boundary equation 2

%The initial solution at t = 0
U_init = zeros(N,N);
U_init(N,:) = fb;       %given bound LEFT
U_init(1,:) = gb;       %given bound RIGHT
U_init(:,N) = fb(N) + (x-ax)/(bx-ax)*(gb(N)-fb(N)); %given bound BOTTOM

tic     %start timer
%Explicit Method
U_solut = U_init;
iteration_count_Explicit=0; % start counting iterations
if exist('checkpoint.mat','file') %checks for checkpoint existence
    load checkpoint.mat; %loads last checkpoint
end
error = 10;
for t = 0:dt:10
    U_init = U_solut;
    for i = 2:N-1
        for j = 1:N-1
        %Neumannn BC with ghost node
        if i >= 2 && j == 1 && i <= N-1
             U_solut(i,j) = lambda*(2* U_init(i, j+1) + U_init(i+1,j) + U_init(i-1,j)) + (1-4*lambda)*(U_init(i,j)); 
        %interior points
        else
             U_solut(i,j) = lambda*(U_init(i,j-1)+U_init(i, j+1) + U_init(i+1,j) + U_init(i-1,j)) + (1-4*lambda)*(U_init(i,j)); 
        end
        end

    end
    iteration_count_Explicit=iteration_count_Explicit+1; %counting iterations for convergence study
    error = abs(mean(mean(U_solut)) - mean(mean(U_init)))/abs(mean(mean(U_solut))); %calculate relative error
    
    if mod(iteration_count_Explicit,1000)==0   % Checkpoint: every 1000 iterations, save
        save checkpoint.mat N t dx U_solut x y %Save variables
    end  
    if error < ERR       %compare error to user desired accuracy
    break
    end
   %Plotting
   surf(x,y,U_solut)
   xlabel('x'); ylabel('y'); zlabel('U');
   title('Solution for U using the Explicit Discretization','fontweight','normal'); 
   rotate3d
   box on
   axis tight
   k =  colorbar;
   k.Label.String = 'U';
   colormap winter;
   drawnow
end
%Statistics
toc
disp('\nThe number of iterations it took to find u using Explicit Discretization is: ')
iteration_count_Explicit
delete checkpoint.mat %removes checkpoints if applicable