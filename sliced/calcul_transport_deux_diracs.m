% Cout de Transport approche SW dans le cas de deux nuages de points
% de deux ements qui forment un parallelogramme
%
% Julien Rabin et Julie Delon, 2013

close all, clear all, clc


%%
choice = {'random uniform','regular','other'};
flag_theta = choice{2};

compute_grad = 0;
plot_contour = 0;


% ensemble d'orientations
N = 100; %nombre de directions utilisees
if strcmpi(flag_theta,choice{1})
  theta = rand(1,N)*pi;
elseif  strcmpi(flag_theta,choice{2})
  theta = [0:N-1]/N*pi; % + rand*2*pi;
else
  theta = ([0 90]+0*45)/180*pi;
  N = size(theta,2);
end

utheta = [cos(theta);sin(theta)];

figure, polar(theta,ones(1,N),'o'), sdf

h=figure, bar([0:.01:pi],hist(theta,[0:.01:pi])),
sdf, xlabel('\theta'), set(gca,'ytick',[]);
set(gca,'XTick',[0 pi])
set(gca,'XTickLabel',['0|\pi'])
axis([0 pi 0 1])
set(h,'OuterPosition',[0 0 600 200])


%distribution Y, fixe
Y = [0,1;0,-1];
Ybis = cat(1,Y(2,:),Y(1,:));

%distribution X, dispose symetriquement selon l'origine
% selon disque
%A=10; R=10; alpha = [0:A-1]/A*2*pi; rho = 2 * 2/pi * [0:R]/R;
% selon rectangle
h=1e2; H = 2*linspace(-1,1,h); 
l=h/2; L = 2*2/pi*linspace(0,1,l); % facteur 2/pi car il y a un point selle en y=0 et x=2/pi
[pos_x,pos_y] = meshgrid(L,H);

tstart = tic;
G = zeros(l,h,2);
T_SW2 = zeros(h,l);
T_W2 = T_SW2; T_SW2_exact = T_SW2;
ell = 2; % constante
for j=1:l
  for i=1:h
    %X = rho(j)*[cos(alpha(i)),sin(alpha);-cos(alpha(i)),sin(alpha)];

    % X = [1,0;-1,0]; cas particulier du carrÃ©
    % X = Y; cas particulier de l'assignement
    % X = 2/pi*[1,0;-1,0]; %cas particulier oÃ¹ le gradient s'annule

    %X = [L(i), H(j); -L(i), -H(j) ];
    X = [pos_x(i,j), pos_y(i,j); -pos_x(i,j), -pos_y(i,j) ];

    
    
    cout = ((X(1,:) - Y(1,:))*utheta).^2 + ((X(2,:) - Y(2,:))*utheta).^2 - ((X(1,:) ...
            - Y(2,:))*utheta).^2 - ((X(2,:) - Y(1,:))*utheta).^2;
    
    if compute_grad    
        gradient1(1) = sum( ( (cout<0).*((X(1,:)-Y(1,:))*utheta) ...
                      +(cout>=0).*((X(1,:)-Y(2,:))*utheta) )*utheta(1,:)');
        gradient1(2) = sum( ( (cout<0).*((X(1,:)-Y(1,:))*utheta) ...
                      +(cout>=0).*((X(1,:)-Y(2,:))*utheta) )*utheta(2,:)');
        gradient2(1) = sum( ( (cout<0).*((X(2,:)-Y(2,:))*utheta) ...
                      +(cout>=0).*((X(2,:)-Y(1,:))*utheta) )*utheta(1,:)');
        gradient2(2) = sum( ( (cout<0).*((X(2,:)-Y(2,:))*utheta) ...
                      +(cout>=0).*((X(2,:)-Y(1,:))*utheta) )*utheta(2,:)');

        G(i,j,:) = 1/N*gradient1; % gradient selon première coordonnée X1
    end        
                
    transport_SW = sum(   ((cout<0).*((X(1,:)-Y(1,:))*utheta)).^2 +  ...
                ((cout>=0).*((X(1,:)-Y(2,:))*utheta)).^2 + ((cout<0).*((X(2,:)-Y(2,:))*utheta)).^2 + ...
                ((cout>=0).*((X(2,:)-Y(1,:))*utheta)).^2 ) ...
                ; %  * 1/N;
    
    T_SW2(i,j) = transport_SW/N;     % SW_2^2
    
    
    T_W2(i,j) = min( norm(vector(X-Y)) , norm(vector(X-Ybis)) ).^2; % W_2^2
    
    % old calculus v1
    % x = pos_y(i,j);
    % y = pos_x(i,j); ay = abs(y);
%     T_SW2_exact(i,j) = 2*pi * (x^2+y^2) - pi*ell*2*x ...
%                  -4*ell*y^2*ay / (x^2+y^2) ...
%                  +4*ell * x ( atan2(ay,x) - ay*x / (x^2+y^2) ));
    % old calculus v2
    % x = pos_y(i,j);
    % y = pos_x(i,j); ay = abs(y);
%     T_SW2_exact(i,j) =  2*pi * (x^2+y^2-2*x) ...
%                      -  8*ay ...
%                      +  8 * x * atan2(ay,x);
    % new calculus
    x = (pos_x(i,j));
    y = (pos_y(i,j));
    r2 = x^2+y^2;
    % T_SW2_exact(i,j) =  r2 +1 - 2*y - 4*x/pi + 4*y/pi * atan2(x,y);
    T_SW2_exact(i,j) =  r2 + 1 - 4/pi * ( x + y * atan2(y,x) );
    
    
  end
end

toc(tstart)


% display
for it=1:3
    switch it
        case 1
            T = T_SW2;
            fig_name = 'approx. SW_2^2'
        case 2 
            T = T_W2;
            fig_name = 'W_2^2'
        case 3 
            T = T_SW2_exact;
            fig_name = 'exact SW_2^2'
    end
    
    %normalization pour affichage
    %T = T ./ max(T(:)) *2;

    fig_num=figure, surf([-flipdim(pos_x,2), pos_x],[pos_y, pos_y],[flipdim(T,2) T])
    shading interp
    xlabel('x'), ylabel('y'), zlabel(fig_name), % zlabel(fig_name,'Interpreter','Latex')
    %material metal
    if exist('sdf'), sdf(), end
    axis normal; % axis equal
    el = 40; az = -60;
    view(az, el);
     alpha('direct')
     alpha(.8)
    set(gca,'ztick',[]);
    axis()
    set(gca,'XTick',[-1 -2/pi 0 2/pi 1])
    set(gca,'YTick',[-2 -1 0 1 2])
    set(gca,'XTickLabel',['-1| |0| |1'])
    %set(gca,'XTickLabel','-1|-2/ \pi |0|2/ \pi |1') %\frac{2}{\pi}

    mytexstr = '-$\frac{2}{\pi}$';
    text('string',mytexstr,'interpreter','latex','fontsize',30,...
                 'units','norm', 'pos',[.82 .07]);
    mytexstr = '$\frac{2}{\pi}$';
    text('string',mytexstr, 'interpreter','latex',...
                 'fontsize',30,'units','norm','pos',[.94 .17]);

    set(fig_num,'OuterPosition',[500*(it-1) 50*it 500 500])
end
%% display other graphs

if plot_contour
    
    for it=1:3
        switch it
            case 1
                T = T_SW2;
                fig_name = 'approx. SW_2^2'
            case 2 
                T = T_W2;
                fig_name = 'W_2^2'
            case 3 
                T = T_SW2_exact;
                fig_name = 'exact SW_2^2'
        end
        %figure, imagesc([-flipdim(L,2) L],H,[flipdim(T,2) T])

        h=figure, contourf([-flipdim(pos_x,2), pos_x],[pos_y, pos_y],[flipdim(T,2) T],30)

        xlabel('x'), ylabel('y'), %zlabel(fig_name), 
        set(gca,'XTick',[-1 -2/pi 0 2/pi 1])
        set(gca,'YTick',[-2 -1 0 1 2])
        set(gca,'XTickLabel',['-1| |0| |1'])
        set(h,'OuterPosition',[500*(it-1) 100*(it+1) 500 500])
        if exist('sdf'), sdf(), end
        mytexstr = '-$\frac{2}{\pi}$';
        text('string',mytexstr,'interpreter','latex',...
                  'fontsize',30,'units','norm','pos',[.23 -.04]);
        mytexstr = '$\frac{2}{\pi}$';
        text('string',mytexstr,'interpreter','latex',...
                  'fontsize',30, 'units','norm', 'pos',[.73 -.04]);

        % figure, contour3([-flipdim(pos_x,2), pos_x],[pos_y, pos_y],[flipdim(T_SW,2) T_SW],20)
        % h = findobj('Type','patch');
        % set(h,'LineWidth',2)
        % title('Twenty Contours of the SW Function')
        
    end
end

return

%% W2 cost function
fig_num2=figure, surf([-flipdim(pos_x,2), pos_x],[pos_y, pos_y],[flipdim(T_W2,2) T_W2])
shading interp
xlabel('x'), ylabel('y'), zlabel('W_2^2'), % zlabel('$W_2^2$','Interpreter','Latex')
%material metal
sdf()
axis equal
el = 40; az = -60;
view(az, el);
 alpha('direct')
 alpha(.8)
set(gca,'ztick',[]);
axis()
set(gca,'XTick',[-1 -.5 0 .5 1])
set(gca,'YTick',[-2 -1 0 1 2])
set(fig_num2,'OuterPosition',[200 0 800 800])

           
%figure, imagesc([-flipdim(L,2) L],H,[flipdim(T,2) T])

h=figure, contourf([-flipdim(L,2) L],H,[flipdim(T_W2,2) T_W2],40)

xlabel('x'), ylabel('y'),
set(gca,'XTick',[-1 -.5 0 .5 1])
set(gca,'YTick',[-2 -1 0 1 2])
set(h,'OuterPosition',[200 0 800 800])
sdf()
set(fig_num,'OuterPosition',[100 100 900 1100])

% display gradient norm
if compute_grad
    norm_G = sum(G.^2,3);
    %figure, imagesc(L,H,norm_G)
end
% display gradient flow
if compute_grad
    Grad_flow = 4*G; %[-flipdim(G,2) G];
    options.subsampling = 1;
    options.linestyle='m';
    figure, 
    %subplot(1,2,1);
    options.display_streamlines = 0;
    plot_vf( Grad_flow, T_SW2, options );
    title('Gradient flow display with arrows'); axis tight;
    figure, 
    %subplot(1,2,2);
    options.display_streamlines = 1;
    plot_vf( Grad_flow, T_SW2, options );
    axis xy; axis equal;
    title('Gradient flow display with streamlines'); axis tight;
end

