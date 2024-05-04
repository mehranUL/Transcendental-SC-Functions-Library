
clear
%for jj = 1:100
format long
N = 1024;

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
%for jj = 1:1024

if N == 2048
    %[~,lfval] = LFSR3([true false true false true true true true false false true],N/2,N); %N=2048
    %[~,lfval2] = LFSR3_2([false true false true false false false true false true false],N/2,N); %N=2048
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
elseif N == 1024
%     [~,lfval] = LFSR3([true false true false true true true false false false],N/2,N); %N=1024
%     [~,lfval2] = LFSR3_2([false true false true false false false true false true],N/2,N); %N=1024
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
elseif N == 512
    %[~,lfval] = LFSR3([true true true true true true true false false],N/2,N); %N=512
    %[~,lfval2] = LFSR3_2([false false true true true true false false true],N/2,N); %N=512
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
elseif N == 256
    %bb = [false true true true true true true true];
    %bb = bb(randperm(numel(bb)));
    %[~,lfval] = LFSR3([false true true true true true true true],N/2,N); %N=256
    %[~,lfval2] = LFSR3_2([false false false false true true false false],N/2,N); %N=256
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
end
lfval = lfval/N;
lfval2 = lfval2/N;

%lfval = rand(1,N);
%lfval2 = rand(1,N);

cos_sobol = zeros(1,N);
abs_cos_sob = zeros(1,N);
cos_vdc = zeros(1,N);
abs_cos_vd = zeros(1,N);
cos_lfsr = zeros(1,N);
abs_cos_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X3_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X3_stream_vdc = zeros(N, N);
X4_stream_vdc = zeros(N, N);
X5_stream_vdc = zeros(N, N);
X2_stream_lfsr_input = zeros(N, N);
X2_stream_lfsr_coeff = zeros(N, N);

for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,4)%4
        if i/N > vd(k,3)
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > sobol(k,9)
        %if i/N > vd(k,10)
            X3_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,3) %3
        %if i/N > z1(k)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,2)%2
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,4)%4
        %if i/N > z1(k)
            X4_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,8)%8
        %if i/N > z1(k)
            X5_stream_vdc(i,k) = 1;
        end
        if i/N > lfval(k)
            X2_stream_lfsr_input(i,k) = 1;
        end
        if i/N > lfval2(k)
            X2_stream_lfsr_coeff(i,k) = 1;
        end
    end
end

%N=1024, 1/56->18 1/30->34 1/12->85 1/2->512
%N=512, 1/56->9 1/30->17 1/12->43 1/2->256
%N=256, 1/56->4 1/30->8 1/12->21 1/2->128




for i = 1:N
%     EE(i) = cos((i-1)/N);
%     if i == N
%         EE(i) = cos(i/N);
%     end
    EE(i) = 1 - ((i/N)^2)/factorial(2) + ((i/N)^4)/factorial(4) - ((i/N)^6)/factorial(6) + ((i/N)^8)/factorial(8);
    %EE(i) = 1 - (((i-1)/N)^2)/factorial(2) + (((i-1)/N)^4)/factorial(4) - (((i-1)/N)^6)/factorial(6) + (((i-1)/N)^8)/factorial(8);

    n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),4));
    n2_s(i,:) = not(and(n1_s(i,:), circshift(X2_stream_sobol(uint32((1/56)*N),:),0)));
    n3_s(i,:) = not(and3(circshift(X3_stream_sobol(uint32((1/30)*N),:),0), n2_s(i,:), circshift(n1_s(i,:),1)));
    n4_s(i,:) = not(and3(circshift(X2_stream_sobol(uint32((1/12)*N),:),0), n3_s(i,:), circshift(n1_s(i,:),2)));
    y_s(i,:) = not(and3(circshift(X3_stream_sobol(uint32((1/2)*N),:),0), n4_s(i,:), circshift(n1_s(i,:),3)));
    cos_sobol(i) = sum(y_s(i,:))/N;
    %abs_cos_sob(i) = abs(cos_sobol(i) - EE(i));
    abs_cos_sob(i) = (cos_sobol(i) - EE(i))^2;

    input = X2_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr_input(i,:);
    n1_v(i,:) = and(input, circshift(input,2));%best with 15 Delays
    n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(uint32((1/56)*N),:)));
    %n2_v(i,:) = n1_v(i,:);
    n3_v(i,:) = not(and3(X3_stream_vdc(uint32((1/30)*N),:), n2_v(i,:), circshift(n1_v(i,:),0)));
    n4_v(i,:) = not(and3(X4_stream_vdc(uint32((1/12)*N),:), n3_v(i,:), circshift(n1_v(i,:),0)));
    y_v(i,:) = not(and3(X5_stream_vdc(uint32((1/2)*N),:), n4_v(i,:), circshift(n1_v(i,:),0)));
    cos_vdc(i) = sum(y_v(i,:))/N;
    %abs_cos_vd(i) = abs(cos_vdc(i) - EE(i));
    abs_cos_vd(i) = (cos_vdc(i) - EE(i))^2;

    n1_lf(i,:) = and(X2_stream_lfsr_input(i,:), circshift(X2_stream_lfsr_input(i,:),4));
    n2_lf(i,:) = not(and(n1_lf(i,:), circshift(X2_stream_lfsr_coeff(uint32((1/56)*N),:),1)));
    n3_lf(i,:) = not(and3(circshift(X2_stream_lfsr_coeff(uint32((1/30)*N),:),1), n2_lf(i,:), circshift(n1_lf(i,:),1)));
    n4_lf(i,:) = not(and3(circshift(X2_stream_lfsr_coeff(uint32((1/12)*N),:),1), n3_lf(i,:), circshift(n1_lf(i,:),2)));
    y_lf(i,:) = not(and3(circshift(X2_stream_lfsr_coeff(uint32((1/2)*N),:),1), n4_lf(i,:), circshift(n1_lf(i,:),3)));
    cos_lfsr(i) = sum(y_lf(i,:))/N;
    %abs_cos_lf(i) = abs(cos_lfsr(i) - EE(i));
    abs_cos_lf(i) = (cos_lfsr(i) - EE(i))^2;
end
%MSE_sobol = mean(abs_cos_sob)
MSE_vdc = mean(abs_cos_vd)
MSE_lfsr = mean(abs_cos_lf)
%end
%lfsr_vdc = mean(MSE_vdc)
%lfsr_lfsr = mean(MSE_lfsr)