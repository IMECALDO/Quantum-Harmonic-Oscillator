clc,clear,close all
format compact

% initial conditions

h = 1;
m = 1;
w = 1;
s  = sqrt(2*h/(m*w));

Ti = 0;    % start of T
Tf = 2*pi; % end of T
Tn = 100;  % divisions of the T array
T_array = linspace(Ti,Tf,Tn+1);

X_lim = 6; % limits of the X array
Xn = 100; % divisions of the X array
X_array = linspace(-X_lim,X_lim,Xn);

iter = 15;  % max amount of Cn

% Toggle what to show
Show_Phi_0         = false;
Show_animation     = true;
Show_mean          = false;
Show_bidimensional = false;
Show_angle         = false;

% Creation of Phi(n)
C_c = (2/(pi*s^2))^(1/4);
Phi_n = @(x,n) C_c * exp(-x.^2/(s^2)).* hermiteH(n,x.*sqrt(2)/s) /sqrt(2^n * factorial(n));

% Creation of Phi(x,0)

    % Case 1

%     C = [1,0.7,-1.2,2];
%     Phi_x0 = @(x,A) A.*(C(1)*Phi_n(x,0) + C(2)*Phi_n(x,1) +C(3)*Phi_n(x,2) + C(4)*Phi_n(x,3));
 
    % Case 2

    % Case 2.1
%     b = 0.7*s;
%     c = 1.5*s;
  
    % Case 2.2
%    b = s;
%    c = 1.5*s;

    % Case 2.3
%    b = 1.5*s;
%    c = 1.5*s;

    % Case #

%     b = 1.5;
%     c = 2;
    
%     Phi_x0 = @(x,A) A.*exp(-(x-c).^2.*(1/b^2));
    
    % Case 3
    
    b = .7*s;
    c = 1.5*s;

%     v = 1
    v = 2
%     v = 3
    Phi_x0 = @(x,A) A.*exp(-(x-c).^2.*(1/b^2)).*exp(-1i.*m.*v.*x.*(1/h));

    % Case 4

%     Phi_x0 = @(x,A) A.*exp(-( ((x-(s/4)).^2).*(1/(s^2)) )).*cos(2.*pi.*x./s);

% Normalization and graph
A = Normalize(Phi_x0,X_array,Show_Phi_0);

%Calculation of Cn
C_array = Get_Cn(Phi_x0,Phi_n,A,iter);
%%
% Calculation of Phi(x,t)
Phi_xt_m = Get_Phi_xt_m(h,w,iter,C_array,Phi_n);

% Cration of Mean Func
mean_x_func = @(x,t,Phi_xt_m,iter,Xn) Phi_xt_Xn(x,t,Phi_xt_m,iter)'.'.*x.*Phi_xt_Xn(x,t,Phi_xt_m,iter);

% Make animation // bidimensional graph // angles

if Show_animation == true
    video = start_animation();
end
if Show_mean == true
    mean_x = zeros(Tn,1);
end

V = zeros(Tn+1,Xn);
for t = 1:Tn+1

    if Show_bidimensional == true || Show_animation == true || Show_angle == true
        V(t,:) = Phi_xt(X_array,T_array(t),Phi_xt_m,iter,Xn);
    end

    if Show_mean == true || Show_animation == true
        mean_x(t) = integral(@(x) mean_x_func(x,T_array(t),Phi_xt_m,iter,Xn),-inf,inf);
    end
    
    if Show_animation == true
        animate(X_array,V(t,:),T_array(t),Tf,X_lim,mean_x(t),video);
    end
    
    disp(T_array(t)*100/Tf + "%")
end

if Show_animation == true
    close_animation(video)
end

if Show_bidimensional == true
    graph_bidimensional(X_lim,Ti,Tf,V)
end

if Show_mean == true
    graph_mean(T_array,mean_x,Ti,Tf,X_lim)
end

if Show_angle == true
    graph_angle(V,X_lim,Ti,Tf)
end

%Check that the normalization is mantained at test_t
% test_t = 0;
% Phi_xt_p2 = @(x) Phi_xt_Xn(x,test_t,Phi_xt_m,iter)'.'.*Phi_xt_Xn(x,test_t,Phi_xt_m,iter);
% disp("Phi(" + test_t + ")^2 = " + integral(Phi_xt_p2,-inf,inf))

% Normalization and graph func
function A = Normalize(Phi_x0,X_array,show)
A_initial = 1;
Phi_x0_p2 = @(x,A) abs(Phi_x0(x,A)'.'.*Phi_x0(x,A));

A = sqrt(1/(integral(@(x) Phi_x0_p2(x,A_initial),-Inf,Inf)));
if show == true
% Plot
figure()
plot(X_array,abs(Phi_x0_p2(X_array,A_initial)),"Color","r")
title("Phi al cuadrado")
hold on
plot(X_array,abs(Phi_x0_p2(X_array,A)),"Color","b")
hold off
legend('Original','Normalized')
xlabel("x")
end

disp("A = " + A)
end 

% Get Cn
function C_array = Get_Cn(Phi_x0,Phi_n,A,iter)
    bra_phi_n_ket_phi_x0 = @(x,n,A) Phi_n(x,n)'.'.*Phi_x0(x,A);
    
    Cn = @(n) integral(@(x) bra_phi_n_ket_phi_x0(x,n,A),-Inf,Inf);
        
    C_array = zeros(1,iter);
    for i = 1:iter
        C_array(i) = Cn(i-1);
    end

    disp("Sum of all |Cn|^2 used = " + sum(C_array'.'.*C_array) )

end

% Get Phi(x,t)
function Phi_xt_m = Get_Phi_xt_m(h,w,iter,C_array,Phi_n)
    E_n = @(n) h*w*(n+1/2);
    Phi_xt_m = cell(iter,1);
    
    for i = 1:iter
        Phi_xt_m{i} = @(x,t) C_array(i)*Phi_n(x,i-1).*exp(-1i*E_n(i-1)*t/h);
    end
    
end

% Phi(x,t) func
function V = Phi_xt(X_array,t,Phi_xt_m,iter,Xn)
   Vt = zeros(iter,Xn);
   for i = 1:iter
        Vt(i,:) = Phi_xt_m{i}(X_array,t);
   end
   V = sum(Vt);
end

% Phi(x,t) func for Integral
function V = Phi_xt_Xn(X_array,t,Phi_xt_m,iter)
   Xn = length(X_array); % Is needed for the integral
   Vt = zeros(iter,Xn);
   for i = 1:iter
        Vt(i,:) = Phi_xt_m{i}(X_array,t);
   end
   V = abs(sum(Vt));
end

% Creates file to start animation
function video = start_animation()
    figure()
    video = VideoWriter('Phi(xt)','MPEG-4'); 
    video.FrameRate = 30;
    open(video);
end

% Graphs and animates
function video = animate(X_array,V,t,Tf,X_lim,mean_x,video)
    V = abs(V'.'.*V);
    plot(X_array,V,"Color","k");
    title("t = " + t*100/Tf + "%")
    ylim([0,1]);
    xlim([-X_lim,X_lim])
    hold on 
    plot([mean_x,mean_x],[1,0],"Color","r")
    hold off

    F = getframe(gcf);
    writeVideo(video, F);
end

% Close file    
function close_animation(video)
    close(video);
end

% Graphs bidimensional Graph
function graph_bidimensional(X_lim,Ti,Tf,V)
    V = abs(V'.'.*V);
    figure()
    x = [-X_lim X_lim];
    y = [Ti Tf];
    Vlims = [0 1];
    imagesc(x,y,V,Vlims)
    colorbar
    title("|Phi(x,t)|^2")
    xlabel("x")
    ylabel("t")
end

% Graphs mean
function graph_mean(T_array,mean_x,Ti,Tf,X_lim)
    figure()
    plot(mean_x,T_array)
    xlim([-X_lim X_lim])
    ylim([Ti Tf])
    xlabel("x")
    ylabel("t")
    title("mean of x")
end

% Graphs angle
function graph_angle(V,X_lim,Ti,Tf)
    figure()
    arg = angle(V);
    x = [-X_lim X_lim];
    y = [Ti Tf];
    imagesc(x,y,arg)
    title("fase de función de estado")
    xlabel("x")
    ylabel("t")
    colorbar
end


