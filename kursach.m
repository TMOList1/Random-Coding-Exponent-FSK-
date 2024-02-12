close all;
clc;
clear;

SNRdB = [0 5 10 15 20];
gamma_ = [0 0 0 0 0];

for i=1:5
    gamma_(i) = 10 .^(SNRdB(i)/10);
end

%ro = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
%R = [0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
R = 0:0.01:1;
Er = zeros(5, length(R));

R0_1 = log(2) - log(1 + 4*(gamma_(1) + 1)/((gamma_(1) + 2)^2));
R0_2 = log(2) - log(1 + 4*(gamma_(2) + 1)/((gamma_(2) + 2)^2));
R0_3 = log(2) - log(1 + 4*(gamma_(3) + 1)/((gamma_(3) + 2)^2));
R0_4 = log(2) - log(1 + 4*(gamma_(4) + 1)/((gamma_(4) + 2)^2));
R0_5 = log(2) - log(1 + 4*(gamma_(5) + 1)/((gamma_(5) + 2)^2));


% fun1 = @(x) (1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(1))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(1))).^2./2))./2 .* log2((1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(1))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(1))).^2./2))./2);
% C1 = integral(fun1,0,20);
% C1 = -1/2 * log2(2 * pi * exp(1)) - 2 * C1;
% 
% fun2 = @(x) (1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(2))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(2))).^2./2))./2 .* log2((1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(2))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(2))).^2./2))./2);
% C2 = integral(fun2,0,22);
% C2 = -1/2 * log2(2 * pi * exp(1)) - 2 * C2;
% 
% fun3 = @(x) (1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(3))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(3))).^2./2))./2 .* log2((1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(3))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(3))).^2./2))./2);
% C3 = integral(fun3,0,24);
% C3 = -1/2 * log2(2 * pi * exp(1)) - 2 * C3;
% 
% fun4 = @(x) (1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(4))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(4))).^2./2))./2 .* log2((1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(4))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(4))).^2./2))./2);
% C4 = integral(fun4,0,26);
% C4 = -1/2 * log2(2 * pi * exp(1)) - 2 * C4;
% 
% fun5 = @(x) (1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(5))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(5))).^2./2))./2 .* log2((1./sqrt(2.*pi) .* exp(-(x-sqrt(2.*gamma_(5))).^2./2) + 1./sqrt(2.*pi) .* exp(-(x+sqrt(2.*gamma_(5))).^2./2))./2);
% C5 = integral(fun5,0,28);
% C5 = -1/2 * log2(2 * pi * exp(1)) - 2 * C5;


for k=1:5
    for j=1:length(R)
        E_max = 0;
        for ro=0:0.05:1
            tmp = 1;
            l = 0;
%             while (tmp > 0.000001 || tmp < -0.000001)
                fst_s = (1 + ro) * log(2);
                ss = 0;
                for li = 0:fix(1+ro)
                    ss = ss + nchoosek((1 + ro), li) * (1 / (1 + ro + li * gamma_(k))); 
                end
                pro = 2 * (gamma_(k) + 1)*(1 + ro) / (gamma_(k) + 2);
                E0 = fst_s - log(pro * ss);

            E = E0 - ro * R(j);
            
            if (i == 1)
                E_max = E;
            else
                if (E > E_max)
                    E_max = E;
                end
            end
        end

        Er(k,j) = E_max;
        
    end
end

plot(R, Er(5,:), R, Er(4,:), R, Er(3,:), R, Er(2,:), R, Er(1,:) );
legend(['SNR = 20 dB, R_0 = ' num2str(R0_5)], ['SNR = 15 dB, R_0 = ' num2str(R0_4)], ['SNR = 10 dB, R_0 = ' num2str(R0_3)], ['SNR = 5 dB, R_0 = ' num2str(R0_2)], ['SNR = 0 dB, R_0 = ' num2str(R0_1)]);
title({'E_r(R) for Rayleigh channel, binary input, uninterrupted output'})
xlabel('R, бит/символ')
ylabel('E_r(R)')



