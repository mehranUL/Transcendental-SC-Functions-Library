
clear
format long
N = 256;

sobol = net(sobolset(256), N);
vd(:,1) = vdcorput(N-1,2);
vd(:,2) = vdcorput(N-1,4);
vd(:,3) = vdcorput(N-1,8);
vd(:,4) = vdcorput(N-1,16);
vd(:,5) = vdcorput(N-1,32);
vd(:,6) = vdcorput(N-1,64);
vd(:,7) = vdcorput(N-1,128);
vd(:,8) = vdcorput(N-1,256);
vd(:,9) = vdcorput(N-1,512);
vd(:,10) = vdcorput(N-1,1024);

%z1 = sqrt(-2.*(log(vd(1:end,1)))).*sin(2*pi.*vd(1:end,8));
% z1 = sqrt(-2.*(log(vd(1:end,1)))).*tanh(vd(1:end,1));
% z1 = z1 - floor(z1);

seed_mat = generateBinaryValues(N);
%for jj = 1:1000

%LFSR3([true false true true true false false false false false],N/2,N);

if N == 1024
    %[~,lfval] = LFSR3([true false true true true false false true false false],N/2,N); %N=1024
    %[~,lfval2] = LFSR3_2([false true false true false false false true false true],N/2,N); %N=1024
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
    lfval = lfval/N;
    lfval2 = lfval2/N;
elseif N == 512
    %[~,lfval] = LFSR3([true false true false true false false false false],N/2,N); %N=512
    %[~,lfval2] = LFSR3_2([false false true true true true false false true],N/2,N); %N=512
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
    lfval = lfval/N;
    lfval2 = lfval2/N;
elseif N == 256
    %[~,lfval] = LFSR3([true false true false true false false false],N/2,N); %N=256
    %[~,lfval2] = LFSR3_2([false false false false true true false false],N/2,N); %N=256
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
    lfval = lfval/N;
    lfval2 = lfval2/N;

%LF3(LF3 == -1) = 0;
%lfval = lfval/N;
%lfval2 = lfval2/N;
else
    lfval = rand(1,N);
    lfval2 = rand(1,N);
end

tanh_sobol = zeros(1,N);
abs_tanh_sob = zeros(1,N);
tanh_vdc = zeros(1,N);
abs_tanh_vd = zeros(1,N);
tanh_lfsr = zeros(1,N);
abs_tanh_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X3_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X3_stream_vdc = zeros(N, N);
X4_stream_vdc = zeros(N, N);
X5_stream_vdc = zeros(N, N);
X2_stream_lfsr = zeros(N, N);
X2_stream_lfsr_coeff = zeros(N, N);

for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,1)
        if i/N > vd(k,4) %4
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > sobol(k,220)
        %if i/N > vd(k,10)
            X3_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,5) %1 --- 5
        %if i/N > z1(k)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,4) %4---4
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,1) %4---1
        %if i/N > z1(k)
            X4_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,1) %6
        %if i/N > z1(k)
            X5_stream_vdc(i,k) = 1;
        end
        if i/N > lfval(k)
            X2_stream_lfsr(i,k) = 1;
        end
        if i/N > lfval2(k)
            X2_stream_lfsr_coeff(i,k) = 1;
        end
    end
end

%N=1024, 62/153->415 17/42->414 2/5->409 1/3->341
%N=512, 1/42->12 1/20->26 1/6->85
%N=256, 1/42->6 1/20->13 1/6->42

% n1_s = and(X2_stream_sobol(102,:), circshift(X2_stream_sobol(102,:),3));
% n2_s = not(and(n1_s, X2_stream_sobol(12,:)));
% n3_s = not(and3(X2_stream_sobol(26,:), n2_s, circshift(n1_s,1)));
% n4_s = not(and3(X2_stream_sobol(85,:), n3_s, circshift(n1_s,2)));
% y_s = and(n4_s, circshift(X2_stream_sobol(102,:),6));
% tanh_sobol = sum(y_s)/N
% 
% n1_v = and(X2_stream_vdc(102,:), circshift(X2_stream_vdc(102,:),3));
% n2_v = not(and(n1_v, X2_stream_vdc(12,:)));
% n3_v = not(and3(X2_stream_vdc(26,:), n2_v, circshift(n1_v,1)));
% n4_v = not(and3(X2_stream_vdc(85,:), n3_v, circshift(n1_v,2)));
% y_v = and(n4_v, circshift(X2_stream_vdc(102,:),6));
% tanh_vdc = sum(y_v)/N
% 
% n1_lf = and(X2_stream_lfsr(102,:), circshift(X2_stream_lfsr(102,:),3));
% n2_lf = not(and(n1_lf, X2_stream_lfsr(12,:)));
% n3_lf = not(and3(X2_stream_lfsr(26,:), n2_lf, circshift(n1_lf,1)));
% n4_lf = not(and3(X2_stream_lfsr(85,:), n3_lf, circshift(n1_lf,2)));
% y_lf = and(n4_lf, circshift(X2_stream_lfsr(102,:),6));
% tanh_lfsr = sum(y_lf)/N


for i = 1:N
    %EE(i) = sin((i-1)/N);
    %EE(i) = (i/N) - ((i/N)^3)/3 + 2*((i/N)^5)/15 - 17*((i/N)^7)/315 + 62*((i/N)^9)/2835;
    EE(i) = (i/N) - ((i/N)^3)/3 + 2*((i/N)^5)/15 - 17*((i/N)^7)/315;

%----------------------------------9th order------------------------
%     n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),4));
%     n2_s(i,:) = not(and(n1_s(i,:), circshift(X2_stream_sobol(ceil((62/153)*N),:),0)));
%     n3_s(i,:) = not(and3(circshift(X3_stream_sobol(ceil((17/42)*N),:),0), n2_s(i,:), circshift(n1_s(i,:),1)));
%     n4_s(i,:) = not(and3(circshift(X2_stream_sobol(ceil((2/5)*N),:),0), n3_s(i,:), circshift(n1_s(i,:),2)));
%     n5_s(i,:) = not(and3(circshift(X3_stream_sobol(ceil((1/3)*N),:),0), n4_s(i,:), circshift(n1_s(i,:),0)));
%     y_s(i,:) = and(n5_s(i,:), circshift(X2_stream_sobol(i,:),0));
%     tanh_sobol(i) = sum(y_s(i,:))/N;
%     %abs_tanh_sob(i) = abs(tanh_sobol(i) - EE(i));
%     abs_tanh_sob(i) = (tanh_sobol(i) - EE(i))^2;

%----------------------------------7th order------------------------
    n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),3));
    %n2_s(i,:) = not(and(n1_s(i,:), circshift(X2_stream_sobol(ceil((62/153)*N),:),0)));
    n2_s(i,:) = not(and(circshift(X3_stream_sobol(ceil((17/42)*N),:),0), n1_s(i,:)));
    n3_s(i,:) = not(and3(circshift(X2_stream_sobol(ceil((2/5)*N),:),0), n2_s(i,:), circshift(n1_s(i,:),0)));
    n4_s(i,:) = not(and3(circshift(X3_stream_sobol(ceil((1/3)*N),:),0), n3_s(i,:), circshift(n1_s(i,:),0)));
    y_s(i,:) = and(n4_s(i,:), circshift(X2_stream_sobol(i,:),0));
    tanh_sobol(i) = sum(y_s(i,:))/N;
    %abs_tanh_sob(i) = abs(tanh_sobol(i) - EE(i)); % Absolute error
    abs_tanh_sob(i) = (tanh_sobol(i) - EE(i))^2; % Squarred error


%     input = X2_stream_sobol(i,:);
%     %input = X2_stream_vdc(i,:);
%     %input = X2_stream_lfsr(i,:);
%     n1_v(i,:) = and(input, circshift(input,4));
%     n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(ceil((62/153)*N),:)));
%     n3_v(i,:) = not(and3(X3_stream_vdc(ceil((17/42)*N),:), n2_v(i,:), circshift(n1_v(i,:),0)));
%     n4_v(i,:) = not(and3(X4_stream_vdc(ceil((2/5)*N),:), n3_v(i,:), circshift(n1_v(i,:),0)));
%     n5_v(i,:) = not(and3(X5_stream_vdc(ceil((1/3)*N),:), n4_v(i,:), circshift(n1_v(i,:),0)));
%     y_v(i,:) = and(n5_v(i,:), circshift(input,0));
%     tanh_vdc(i) = sum(y_v(i,:))/N;
%     %abs_tanh_vd(i) = abs(tanh_vdc(i) - EE(i));
%     abs_tanh_vd(i) = (tanh_vdc(i) - EE(i))^2;

    input = X2_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr(i,:);
    n1_v(i,:) = and(input, circshift(input,3));
    %n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(ceil((62/153)*N),:)));
    n2_v(i,:) = not(and(X2_stream_vdc(ceil((17/42)*N),:), circshift(n1_v(i,:),0)));
    n3_v(i,:) = not(and3(X3_stream_vdc(ceil((2/5)*N),:), n2_v(i,:), circshift(n1_v(i,:),0)));
    n4_v(i,:) = not(and3(X4_stream_vdc(ceil((1/3)*N),:), n3_v(i,:), circshift(n1_v(i,:),0)));
    y_v(i,:) = and(n4_v(i,:), circshift(input,0));
    tanh_vdc(i) = sum(y_v(i,:))/N;
    %abs_tanh_vd(i) = abs(tanh_vdc(i) - EE(i));
    abs_tanh_vd(i) = (tanh_vdc(i) - EE(i))^2;


%     n1_lf(i,:) = and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),4));
%     n2_lf(i,:) = not(and(n1_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((62/153)*N),:),1)));
%     n3_lf(i,:) = not(and3(circshift(X2_stream_lfsr_coeff(ceil((17/42)*N),:),1), n2_lf(i,:), circshift(n1_lf(i,:),1)));
%     n4_lf(i,:) = not(and3(circshift(X2_stream_lfsr_coeff(ceil((2/5)*N),:),1), n3_lf(i,:), circshift(n1_lf(i,:),2)));
%     n5_lf(i,:) = not(and3(circshift(X2_stream_lfsr_coeff(ceil((1/3)*N),:),1), n4_lf(i,:), circshift(n1_lf(i,:),3)));
%     y_lf(i,:) = and(n5_lf(i,:), circshift(X2_stream_lfsr(i,:),8));
%     tanh_lfsr(i) = sum(y_lf(i,:))/N;
%     %abs_tanh_lf(i) = abs(tanh_lfsr(i) - EE(i));
%     abs_tanh_lf(i) = (tanh_lfsr(i) - EE(i))^2;

    n1_lf(i,:) = and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),3));
    %n2_lf(i,:) = not(and(n1_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((62/153)*N),:),1)));
    n2_lf(i,:) = not(and(circshift(X2_stream_lfsr_coeff(ceil((17/42)*N),:),1), n1_lf(i,:)));
    n3_lf(i,:) = not(and3(circshift(X2_stream_lfsr_coeff(ceil((2/5)*N),:),1), n2_lf(i,:), circshift(n1_lf(i,:),1)));
    n4_lf(i,:) = not(and3(circshift(X2_stream_lfsr_coeff(ceil((1/3)*N),:),1), n3_lf(i,:), circshift(n1_lf(i,:),2)));
    y_lf(i,:) = and(n4_lf(i,:), circshift(X2_stream_lfsr(i,:),6));
    tanh_lfsr(i) = sum(y_lf(i,:))/N;
    %abs_tanh_lf(i) = abs(tanh_lfsr(i) - EE(i));
    abs_tanh_lf(i) = (tanh_lfsr(i) - EE(i))^2;
end
MSE_sobol = mean(abs_tanh_sob)
MSE_vdc = mean(abs_tanh_vd)
MSE_lfsr = mean(abs_tanh_lf)
%end
%lfsr_lfsr = mean(MSE_lfsr)