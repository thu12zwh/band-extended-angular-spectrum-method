clc
clear all
close all
%% object parameters and propagation parameters
n0 = 1024;                                    % number of the original signal
r = n0/4;                                     % number of the rectangle aperture
t = zeros(n0,1);
t (n0/2-r+1:n0/2+r,1) = 1;
t = padarray(t,n0/2);                         %  zero-padded signal
lam = 500e-6;                                 % wavelength
k = 2*pi/lam;                                 % wave number
n = size(t,1);                                % number of the zero-padded signal
pitch = 0.001;                                % sampling pitch in spatial domain
l = n*pitch;                  
x = linspace(-l/2,l/2-pitch,n)';              % cordinate in spatial domain
fx = linspace(-1/2/pitch,1/2/pitch-1/l,n)';   % cordinate in spatial frequency domain
figure,plot(x,t);title('object')
z = 300 ;                                     % propagation distance

%% analytical integral 
X = linspace(-l/4,l/4-pitch,n/2)';
uu = zeros(n/2,1);
disp('analytical intergral:')
tic
for j = 1:n/2
      fun = @(xn) 1/2/pi*z./sqrt((X(j)-xn).^2+z^2).*(1./sqrt((X(j)-xn).^2+z^2)...
               -1i*k/pi).*exp(1i*k*sqrt((X(j)-xn).^2+z^2))./sqrt((X(j)-xn).^2+z^2);
      uu(j,1) = integral(fun,-(r-1)*pitch,r*pitch);
end
toc
uu = uu/max(abs(uu)); 
phase_rsi = (angle(uu));
amplitude_rsi = abs(uu);

figure,plot(X,amplitude_rsi);title('Analytical integral amplitude')
% figure,plot(X,(phase_rsi));title('Analytical integral phase ')

%% impulse response method
R = sqrt(x.^2+z^2); 
kernel = 1/2/pi*z./R.*(1./R-1i*k).*exp(1i*k*R)./R;
disp('impulse response menthod:')
tic
t_FT = fftshift(fft(fftshift(t)));
kernel_FT = fftshift(fft(fftshift(kernel)));
t_1 = ifftshift(ifft(ifftshift(kernel_FT.*t_FT)));
toc
t_1  = t_1/max(abs(t_1));
t_1 = t_1(n/2-n/4+1:n/2+n/4,1);
figure,plot(x(n/2-n/4+1:n/2+n/4),abs(t_1));title('impulse response menthod amplitude')
% figure,plot(x(n/2-n/4+1:n/2+n/4),angle(t_pro2));title('impulse response menthod phase')

%% band-limted ASM
fc = n*pitch/lam/z/2;                % f_limit
disp('band-limited ASM:')
tic
H = exp(1i*k*z*sqrt(1-(fx*lam).^2));
H(abs(fx)>fc) = 0;
t_FT = fftshift(fft(fftshift(t)));
t_2 = ifftshift(ifft(ifftshift(H.*t_FT)));
toc
t_2  = t_2/max(abs(t_2 ));
t_2 = t_2(n/2-n/4+1:n/2+n/4,1);
phase_asm_li = (angle(t_2));
amplitude_asm_li = abs(t_2);

figure,plot(x(n/2-n/4+1:n/2+n/4),amplitude_asm_li);title('band-limited ASM amplitude ')
% figure,plot(x(n/2-n/4+1:n/2+n/4),(phase_asm_li));title('band-limited ASM phase ')

%% entended band ASM
iflag = -1;
eps = 10^(-12);           % accuracy of NUFFT
K = n/2/max(fx);
fcn = 1/2*sqrt(n/lam/z);  % f_extend
ss = fcn/max(abs(fx));
zc = n*pitch^2/lam;
if z <zc
    fxn = fx;
else
    fxn = fx*(ss-0.0);
end
Hn = exp(1i*k*(z*sqrt(1-(fxn*lam).^2)));

disp('band-extended ASM:')
tic
t_asmNUFT = nufft1d3(n,x/max(abs(x))*pi,t,iflag,eps,n,(fxn)*K);
t_3 = nufft1d3(n,x/(max(abs(x)))*pi,Hn.*t_asmNUFT,-iflag,eps,n,(fxn)*K);
toc

t_3 = t_3/max(abs(t_3));
t_3 = t_3(n/2-n/4+1:n/2+n/4,1);
phase_asm_ex = (angle(t_3));
amplitude_asm_ex = abs(t_3);

figure,plot(x(n/2-n/4+1:n/2+n/4),(amplitude_asm_ex));title('band-extended ASM amplitude ')
% figure,plot(x(n/2-n/4+1:n/2+n/4),(phase_asm_ex));title('band-extended ASM phase ')
