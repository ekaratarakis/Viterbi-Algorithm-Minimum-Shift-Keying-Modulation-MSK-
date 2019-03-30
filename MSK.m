clear all
close all
clc

%-------------------------------------------------------------------------%
%                     Name : Karatarakis Evangelos                        %
%                          A.M. : 2013030003                              %
%                   Course : TEL513 - Satellite Links                     %
%            Project 1 : Viterbi algorithm of MSK modulation              %
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%                                PART 1.1                                 %
%-------------------------------------------------------------------------%

N = 1e4; % 10^5
A = 1;
T = 0.01;

% Generate a sequence of symbols {-1,1}
b = (sign(randn(1,2*N))+1)/2;
for i=1:1:length(b) 
    if b(i)==0
        b(i)=-1;
    end
end

% Apply the recursive equations to transform the symbols in QPSK
% symbols
xQ_1 = -1;
x_1 = 1;
xI_n = zeros(1,N);
xQ_n = zeros(1,N);
for n=1:1:N
    if n==1
        xI_n(n) = -xQ_1*x_1;
        xQ_n(n) = -xI_n(n)*b(n);
    elseif n~=1 
        xI_n(n) = -xQ_n(n-1)*b(2*(n-1));
        xQ_n(n) = -xI_n(n)*b(2*n-1);
    end
end

xI_n; % In - Phase component
xQ_n; % Quadrature component

z_n = xI_n + 1i*xQ_n; % Add in - phase and quadrature component to create z_n

% Calculate BER
iter = 1e2; 
SNR_dB = 5;
SNR_lin = 10.^(SNR_dB/10); % 10log_10(SNR_lin)
BER_appr = 0;

for p=1:1:iter
    n = randn(1,N) + 1i*randn(1,N);
    y_n = A*T*z_n + sqrt(A^2*T^2/SNR_lin)*n;
    y_Re = sign(real(y_n)); % returns -1 when a negative number occur and 1 otherwise
    y_Im = sign(imag(y_n));
    BER_appr = BER_appr + sum(y_Im ~= imag(z_n)) + sum(y_Re ~= real(z_n));
end
BER = BER_appr/(N*iter);

D = ['BER = ', num2str(BER)];
disp('------------------------------------------------------------------------------------------------------------------------------------');
disp('An approximation for BER - Bit Error Rate with SNR = 5 dB is:');
disp(D);
disp('------------------------------------------------------------------------------------------------------------------------------------');

SNR_dB_B = 5:12; % Values of SNR in dB
SNR_lin_B = 10.^([5:12]/10); % Values of SNR in decimal

BER_B = zeros(1,length(SNR_lin_B));
BER_appr_B = zeros(1,length(SNR_lin_B));

for k=1:1:length(SNR_lin_B)
    for p=1:1:iter
        n = randn(1,N) + 1i*randn(1,N);
        y_n = A*T*z_n + sqrt(A^2*T^2/SNR_lin_B(k))*n;
        y_Re = sign(real(y_n)); % returns -1 when a negative number occur and 1 otherwise
        y_Im = sign(imag(y_n));
        BER_appr_B(k) = BER_appr_B(k) + sum(y_Im ~= imag(z_n)) + sum(y_Re ~= real(z_n));
    end
    BER_B(k) = BER_appr_B(k)/(2*N*iter);
    Theory_BER(k) = 0.5*erfc(sqrt(SNR_lin_B(k))); % theoretical BER
end

disp(['SNR in dB : ' num2str(SNR_dB_B)]);
disp(['BER Measured : ' num2str(BER_B)]);
disp(['BER Theoretical : ' num2str(Theory_BER)]);
disp('------------------------------------------------------------------------------------------------------------------------------------');

figure(1)
semilogy(SNR_dB_B,BER_B,'b-*') 
xlabel('SNR (dB)')
ylabel('BER - Bit Error Rate Measured')
title('BER - Bit Error Rate vs SNR_{dB}')
grid on

figure(2)
semilogy(SNR_dB_B,BER_B,'b-*') 
hold on
semilogy(SNR_dB_B,Theory_BER,'m-o')
xlabel('SNR (dB)')
ylabel('BER - Bit Error Rate Measured & Theoretical')
title('BER - Bit Error Rate vs SNR_{dB}')
legend('BER - Measured', 'BER - Theoretical')
grid on

%-------------------------------------------------------------------------%
%                                PART 1.2                                 %
%-------------------------------------------------------------------------%

SNR_dB_C = 5;
SNR_lin_C = 10.^(SNR_dB_C/10);

% Creating the symbols x(n) in set {-1,1}
b = (sign(randn(1,N))+1)/2;
for i=1:1:length(b) 
    if b(i)==0
        b(i)=-1;
    end
end

% Constructing the vectors s_1 & s1
% s_1 is for s^(-1)
% s1 is for s^1
s_1 = [-(2*A*sqrt(T)*1i)/pi; (A*sqrt(T)*sqrt(pi^2 - 4))/pi];
s1  = [A*sqrt(T); 0];

beta = (A^2*T)/SNR_lin_C; % From equation (16)
var = 2*beta;
BER_Viterbi = 0; %initialising BER

for p=1:1:iter
    n1_n = sqrt(var)*(randn(1,N) + 1j*randn(1,N)); % Generating the values of the noise vector n1_n
    n2_n = sqrt(var)*(randn(1,N) + 1j*randn(1,N));
    phase(1) = 0;
    for n=1:1:N
        phase(n+1) = phase(n) + b(n)*pi/2;
        if b(n) == -1 
            r_n(:,n) = s_1.*exp(1i*phase(n)) + [n1_n(n); n2_n(n)]; % if x(n) == -1
        elseif b(n) == 1
            r_n(:,n) = s1.*exp(1i*phase(n)) + [n1_n(n); n2_n(n)];  % if x(n) == 1
        end
    end
    x_opt = Viterbi(N,s1,s_1,r_n);
    BER_Viterbi = BER_Viterbi + sum(b~=x_opt); 
end
BER_Viterbi = BER_Viterbi/(2*N*iter);

disp('An approximation for BER - Bit Error Rate with Viterbi aglgorithm and SNR = 5 dB is:');
disp(['BER = ', num2str(BER_Viterbi)])
% disp(['BER Theoretical : ' num2str(Theory_BER)]);
disp('------------------------------------------------------------------------------------------------------------------------------------');

% Creating the BER diagram for SNR = 5,6,7,8,9,10,11,12
SNR_dB_E = 5:12;
SNR_lin_E = 10.^(SNR_dB_E/10);
BER_Viterbi_E = zeros(1,length(SNR_lin_E));
BER_Viterbi_appr_E = zeros(1,length(SNR_lin_E));

% Constructing the vectors s_1 & s1
s_1 = [-(2*A*sqrt(T)*1i)/pi; (A*sqrt(T)*sqrt(pi^2 - 4))/pi];
s1  = [A*sqrt(T); 0];

for k=1:1:length(SNR_lin_E)
    for p=1:1:iter
        n1_n = sqrt(A^2*T/SNR_lin_E(k))*(randn(1,N) + 1j*randn(1,N)); % Generating the values of the noise vector n1_n
        n2_n = sqrt(A^2*T/SNR_lin_E(k))*(randn(1,N) + 1j*randn(1,N));
        phase(1) = 0;
        for n=1:1:N
            phase(n+1) = phase(n) + b(n)*pi/2;
            if b(n) == -1 
                r_n(:,n) = s_1.*exp(1i*phase(n)) + [n1_n(n); n2_n(n)]; % if x(n) == -1
            elseif b(n) == 1
                r_n(:,n) = s1.*exp(1i*phase(n)) + [n1_n(n); n2_n(n)];  % if x(n) == 1
            end
        end
        x_opt = Viterbi(N,s1,s_1,r_n);
        BER_Viterbi_appr_E(k) = BER_Viterbi_appr_E(k) + sum(b~=x_opt);
        Theory_BER_Viterbi(k) = 0.5*erfc(sqrt(SNR_lin_E(k))); % theoretical BER
    end
    BER_Viterbi_E(k) = BER_Viterbi_appr_E(k)/(2*N*iter);
end

disp(['SNR in dB : ' num2str(SNR_dB_E)]);
disp(['BER Measured : ' num2str(BER_Viterbi_E)]);
disp(['BER Theoretical : ' num2str(Theory_BER_Viterbi)]);
disp('------------------------------------------------------------------------------------------------------------------------------------');


figure(3)
semilogy(SNR_dB_E,BER_Viterbi_E,'b-*') 
xlabel('SNR (dB)')
ylabel('BER - Bit Error Rate Measured ')
title('BER - Bit Error Rate vs SNR_{dB} - Viterbi Algorithm')
grid on

figure(4)
semilogy(SNR_dB_E,BER_Viterbi_E,'b-*'),
hold on
semilogy(SNR_dB_B,BER_B,'g-*'),
hold on
semilogy(SNR_dB_B,Theory_BER,'m-o')
xlabel('SNR (dB)')
ylabel('BER - Bit Error Rate Measured & Theoretical')
title('BER - Bit Error Rate vs SNR_{dB}')
legend('BER - Viterbi','BER - Measured', 'BER - Theoretical','Location','southwest')
grid on



