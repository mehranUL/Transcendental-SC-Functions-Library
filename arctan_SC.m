
clear
format long
N = 512;

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


seed_mat = generateBinaryValues(N);
%for jj = 1:1024

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

arctan_sob = zeros(1,N);
abs_arctan_sob = zeros(1,N);
arctan_vdc = zeros(1,N);
abs_arctan_vd = zeros(1,N);
arctan_lfsr = zeros(1,N);
abs_arctan_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X3_stream_vdc = zeros(N, N);
X4_stream_vdc = zeros(N, N);
X5_stream_vdc = zeros(N, N);
X2_stream_lfsr = zeros(N, N);
X2_stream_lfsr_coeff = zeros(N, N);

for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,2)%8 for 512, 2 for 1024,256
        if i/N > vd(k,3)%3
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,9) %9
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,3)%3
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,8)%8
        %if i/N > z1(k)
            X4_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,9)
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

%N=1024, 5/21->244 3/5->614 1/3->341
%N=512, 5/21->122 3/5->307 1/3->171
%N=256, 5/21->61 3/5->154 1/3->85

% n1_s = and(X2_stream_sobol(102,:), circshift(X2_stream_sobol(102,:),3));
% n2_s = not(and(n1_s, X2_stream_sobol(12,:)));
% n3_s = not(and3(X2_stream_sobol(26,:), n2_s, circshift(n1_s,1)));
% n4_s = not(and3(X2_stream_sobol(85,:), n3_s, circshift(n1_s,2)));
% y_s = and(n4_s, circshift(X2_stream_sobol(102,:),6));
% sin_sobol = sum(y_s)/N
% 
% n1_v = and(X2_stream_vdc(102,:), circshift(X2_stream_vdc(102,:),3));
% n2_v = not(and(n1_v, X2_stream_vdc(12,:)));
% n3_v = not(and3(X2_stream_vdc(26,:), n2_v, circshift(n1_v,1)));
% n4_v = not(and3(X2_stream_vdc(85,:), n3_v, circshift(n1_v,2)));
% y_v = and(n4_v, circshift(X2_stream_vdc(102,:),6));
% sin_vdc = sum(y_v)/N
% 
% n1_lf = and(X2_stream_lfsr(102,:), circshift(X2_stream_lfsr(102,:),3));
% n2_lf = not(and(n1_lf, X2_stream_lfsr(12,:)));
% n3_lf = not(and3(X2_stream_lfsr(26,:), n2_lf, circshift(n1_lf,1)));
% n4_lf = not(and3(X2_stream_lfsr(85,:), n3_lf, circshift(n1_lf,2)));
% y_lf = and(n4_lf, circshift(X2_stream_lfsr(102,:),6));
% sin_lfsr = sum(y_lf)/N


for i = 1:N
    %EE(i) = sin((i-1)/N);
    EE(i) = i/N - ((i/N)^3)/3 + ((i/N)^5)/5 - ((i/N)^7)/7;

    n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),3));
    n2_s(i,:) = not(and(n1_s(i,:), X2_stream_sobol(ceil((5/21)*N),:)));
    n3_s(i,:) = not(and3(X2_stream_sobol(ceil((3/5)*N),:), n2_s(i,:), circshift(n1_s(i,:),1)));
    n4_s(i,:) = not(and3(X2_stream_sobol(ceil((1/3)*N),:), n3_s(i,:), circshift(n1_s(i,:),2)));
    y_s(i,:) = and(n4_s(i,:), circshift(X2_stream_sobol(i,:),6));
    arctan_sobol(i) = sum(y_s(i,:))/N;
    %abs_arctan_sob(i) = abs(arctan_sobol(i) - EE(i));
    abs_arctan_sob(i) = (arctan_sobol(i) - EE(i))^2;

    input = X2_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr(i,:);
    n1_v(i,:) = and(input, circshift(input,2));%best with 15 Delays
    %n1_v(i,:) = and(X4_stream_vdc(i,:), circshift(X4_stream_vdc(i,:),3));
    n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(ceil((5/21)*N),:)));
    n3_v(i,:) = not(and3(X3_stream_vdc(ceil((3/5)*N),:), n2_v(i,:), circshift(n1_v(i,:),0)));
    n4_v(i,:) = not(and3(X4_stream_vdc(ceil((1/3)*N),:), n3_v(i,:), circshift(n1_v(i,:),0)));
    y_v(i,:) = and(n4_v(i,:), circshift(input,0));
    arctan_vdc(i) = sum(y_v(i,:))/N;
    %abs_arctan_vd(i) = abs(arctan_vdc(i) - EE(i));
    abs_arctan_vd(i) = (arctan_vdc(i) - EE(i))^2;

    n1_lf(i,:) = and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),3));
    n2_lf(i,:) = not(and(n1_lf(i,:), X2_stream_lfsr_coeff(ceil((5/21)*N),:)));
    n3_lf(i,:) = not(and3(X2_stream_lfsr_coeff(ceil((3/5)*N),:), n2_lf(i,:), circshift(n1_lf(i,:),1)));
    n4_lf(i,:) = not(and3(X2_stream_lfsr_coeff(ceil((1/3)*N),:), n3_lf(i,:), circshift(n1_lf(i,:),2)));
    y_lf(i,:) = and(n4_lf(i,:), circshift(X2_stream_lfsr(i,:),6));
    arctan_lfsr(i) = sum(y_lf(i,:))/N;
    %abs_arctan_lf(i) = abs(arctan_lfsr(i) - EE(i));
    abs_arctan_lf(i) = (arctan_lfsr(i) - EE(i))^2;
end
%MAE_sobol = mean(abs_arctan_sob)
MSE_vdc = mean(abs_arctan_vd)
MSE_lfsr = mean(abs_arctan_lf)
%end
%lfsr_lfsr = mean(MSE_lfsr)