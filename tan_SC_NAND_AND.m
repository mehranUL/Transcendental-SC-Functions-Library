
clear
tic
%for jj = 1:50
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
%for jj = 1:100

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

%tan_sobol = zeros(1,N);
abs_tan_sob = zeros(1,N);
tan_vdc = zeros(1,N);
cot_vdc = zeros(1,N);
abs_tan_vd = zeros(1,N);
tan_lfsr = zeros(1,N);
abs_tan_lf = zeros(1,N);

sin_sob = zeros(1,N);
abs_sin_sob = zeros(1,N);
sin_vdc = zeros(1,N);
abs_sin_vd = zeros(1,N);
sin_lfsr = zeros(1,N);
abs_sin_lf = zeros(1,N);

cos_sobol = zeros(1,N);
abs_cos_sob = zeros(1,N);
cos_vdc = zeros(1,N);
abs_cos_vd = zeros(1,N);
cos_lfsr = zeros(1,N);
abs_cos_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X3_stream_sobol = zeros(N, N);
X2_stream_vdc_cos = zeros(N, N);
X3_stream_vdc_cos = zeros(N, N);
X4_stream_vdc_cos = zeros(N, N);
X5_stream_vdc_cos = zeros(N, N);
X2_stream_vdc_sin = zeros(N, N);
X3_stream_vdc_sin = zeros(N, N);
X4_stream_vdc_sin = zeros(N, N);
X5_stream_vdc_sin = zeros(N, N);
X2_stream_lfsr_input = zeros(N, N);
X2_stream_lfsr_coeff = zeros(N, N);


for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,12)%2
        if i/N > vd(k,3)%3
            X2_stream_sobol(i,k) = 1;
        end
        %if i/N > sobol(k,9)%6,8
        if i/N > vd(k,5)%5
            X3_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,8) %4
        %if i/N > z1(k)
            X2_stream_vdc_cos(i,k) = 1;
        end
        if i/N > vd(k,8)%8
        %if i/N > z1(k)
            X3_stream_vdc_cos(i,k) = 1;
        end
        if i/N > vd(k,8)%4
        %if i/N > z1(k)
            X4_stream_vdc_cos(i,k) = 1;
        end
        if i/N > vd(k,8)%7
        %if i/N > z1(k)
            X5_stream_vdc_cos(i,k) = 1;
        end
        if i/N > vd(k,8)%5
        %if i/N > z1(k)
            X2_stream_vdc_sin(i,k) = 1;
        end
        if i/N > vd(k,8)%8
        %if i/N > z1(k)
            X3_stream_vdc_sin(i,k) = 1;
        end
        if i/N > vd(k,4)%4
        %if i/N > z1(k)
            X4_stream_vdc_sin(i,k) = 1;
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
%N=512, 1/42->12 1/20->26 1/6->85
%N=256, 1/42->6 1/20->13 1/6->42


for i = 1:N
    %EE(i) = cos((i-1)/N);
    EE(i) = (i/N) + ((i/N)^3)/3 + (2/15)*((i/N)^5) + (17/315)*((i/N)^7);



    input = X3_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr_input(i,:);
    n1_v_sin(i,:) = and(input, circshift(input,1));
    n2_v_sin(i,:) = not(and(n1_v_sin(i,:), X2_stream_vdc_sin(ceil((1/5040)*N),:)));
    n3_v_sin(i,:) = and(n2_v_sin(i,:), X3_stream_vdc_sin(ceil((1/120)*N),:));
    n4_v_sin(i,:) = not(and(n3_v_sin(i,:), circshift(n1_v_sin(i,:), 0)));%1
    n5_v_sin(i,:) = and(n4_v_sin(i,:), X4_stream_vdc_sin(ceil((1/6)*N),:));
    n6_v_sin(i,:) = not(and(n5_v_sin(i,:), circshift(n1_v_sin(i,:), 0)));%2
    y_v_sin(i,:) = and(n6_v_sin(i,:), circshift(input,2));%2
    sin_vdc(i) = sum(y_v_sin(i,:))/N;


    input = X2_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr_input(i,:);
    n1_v_cos(i,:) = and(input, circshift(input,1));
    n2_v_cos(i,:) = not(and(n1_v_cos(i,:), X3_stream_vdc_cos(ceil((1/40320)*N),:)));
    n3_v_cos(i,:) = and(n2_v_cos(i,:), X3_stream_vdc_cos(ceil((1/720)*N),:));
    n4_v_cos(i,:) = not(and(n3_v_cos(i,:), circshift(n1_v_cos(i,:),0)));
    n5_v_cos(i,:) = and(n4_v_cos(i,:), X3_stream_vdc_cos(ceil((1/24)*N),:));
    n6_v_cos(i,:) = not(and(n5_v_cos(i,:), circshift(n1_v_cos(i,:),0)));
    n7_v_cos(i,:) = and(n6_v_cos(i,:), X3_stream_vdc_cos(ceil((1/2)*N),:));
    y_v_cos(i,:) = not(and(n7_v_cos(i,:), circshift(n1_v_cos(i,:),0)));
    cos_vdc(i) = sum(y_v_cos(i,:))/N;

%----------------------------------Conventional----------------------------------------------------------------    
    n1_lf_sin(i,:) = and(X2_stream_lfsr_input(i,:), circshift(X2_stream_lfsr_input(i,:),1));
    n2_lf_sin(i,:) = not(and(n1_lf_sin(i,:), X2_stream_lfsr_coeff(ceil((1/5040)*N),:)));
    n3_lf_sin(i,:) = and(n2_lf_sin(i,:), X2_stream_lfsr_coeff(ceil((1/120)*N),:));
    n4_lf_sin(i,:) = not(and(n3_lf_sin(i,:), circshift(n1_lf_sin(i,:), 1)));%1
    n5_lf_sin(i,:) = and(n4_lf_sin(i,:), X2_stream_lfsr_coeff(ceil((1/6)*N),:));
    n6_lf_sin(i,:) = not(and(n5_lf_sin(i,:), circshift(n1_lf_sin(i,:), 2)));%2
    y_lf_sin(i,:) = and(n6_lf_sin(i,:), circshift(X2_stream_lfsr_input(i,:),2));%2
    sin_lfsr(i) = sum(y_lf_sin(i,:))/N;
    
    n1_lf_cos(i,:) = and(X2_stream_lfsr_input(i,:), circshift(X2_stream_lfsr_input(i,:),1));
    n2_lf_cos(i,:) = not(and(n1_lf_cos(i,:), X2_stream_lfsr_coeff(ceil((1/40320)*N),:)));
    n3_lf_cos(i,:) = and(n2_lf_cos(i,:), X2_stream_lfsr_coeff(ceil((1/720)*N),:));
    n4_lf_cos(i,:) = not(and(n3_lf_cos(i,:), circshift(n1_lf_cos(i,:),1)));
    n5_lf_cos(i,:) = and(n4_lf_cos(i,:), X2_stream_lfsr_coeff(ceil((1/24)*N),:));
    n6_lf_cos(i,:) = not(and(n5_lf_cos(i,:), circshift(n1_lf_cos(i,:),2)));
    n7_lf_cos(i,:) = and(n6_lf_cos(i,:), X2_stream_lfsr_coeff(ceil((1/2)*N),:));
    y_lf_cos(i,:) = not(and(n7_lf_cos(i,:), circshift(n1_lf_cos(i,:),3)));
    cos_lfsr(i) = sum(y_lf_cos(i,:))/N;

    if sin_lfsr(i) < cos_lfsr(i)
        [~,y_lf_tan(i,:)] = CORLD_DIV(N*sin_lfsr(i), N*cos_lfsr(i), N);
        tan_lfsr(i) = sum(y_lf_tan(i,:))/N;
        %jj = jj+1;
        %abs_tan_vd(i) = abs(tan_vdc(i) - EE(i));
    else
        [~,y_lf_tan(i,:)] = CORLD_DIV(sum(and(y_lf_sin(floor(i/2),:), y_lf_cos(floor(i/2),:))), N*cos_lfsr(i), N);
        tan_lfsr(i) = 2*sum(y_lf_tan(i,:))/N;
    end
    %abs_tan_lf(i) = abs(tan_lfsr(i) - EE(i));
    abs_tan_lf(i) = (tan_lfsr(i) - EE(i))^2;
%--------------------------------------------------------------------------------------------------------------

    if sin_vdc(i) < cos_vdc(i)
        [~,y_v_tan(i,:)] = CORLD_DIV(N*sin_vdc(i), N*cos_vdc(i), N);
        tan_vdc(i) = sum(y_v_tan(i,:))/N;
        %jj = jj+1;
        %abs_tan_vd(i) = abs(tan_vdc(i) - EE(i));
    else
        [~,y_v_tan(i,:)] = CORLD_DIV(sum(and(y_v_sin(floor(i/2),:), y_v_cos(floor(i/2),:))), N*cos_vdc(i), N);
        tan_vdc(i) = 2*sum(y_v_tan(i,:))/N;


%         tmp = sin_vdc(i) - cos_vdc(i);
%         while tmp > cos_vdc(i)
%             tmp = sin_vdc(i) - cos_vdc(i);
%             ii = ii+1;
%         end
%         jj = jj+1;
%         [~,y_v_tan(jj,:)] = CORLD_DIV(N*tmp, N*cos_vdc(jj), N);
%         tan_vdc(jj) = ii + sum(y_v_tan(i,:))/N;
        
    end
    %abs_tan_vd(i) = abs(tan_vdc(i) - EE(i));
    abs_tan_vd(i) = (tan_vdc(i) - EE(i))^2;
%     for xx = 1:jj
%         tan_vdc(jj) = 
%     end
   

%     n1_lf(i,:) = and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),4));
%     n2_lf(i,:) = not(and(n1_lf(i,:), circshift(X2_stream_lfsr(18,:),1)));
%     n3_lf(i,:) = not(and3(circshift(X2_stream_lfsr(34,:),1), n2_lf(i,:), circshift(n1_lf(i,:),1)));
%     n4_lf(i,:) = not(and3(circshift(X2_stream_lfsr(85,:),1), n3_lf(i,:), circshift(n1_lf(i,:),2)));
%     y_lf(i,:) = not(and3(circshift(X2_stream_lfsr(512,:),1), n4_lf(i,:), circshift(n1_lf(i,:),3)));
%     tan_lfsr(i) = sum(y_lf(i,:))/N;
%     abs_tan_lf(i) = abs(tan_lfsr(i) - EE(i));
end
%MAE_sobol = mean(abs_tan_sob)
MSE_vdc = mean(abs_tan_vd)
MSE_lfsr = mean(abs_tan_lf)
%MAE_lfsr = mean(abs_tan_lf)
%end
%lfsr_vdc = mean(MSE_vdc)
%lfsr_lfsr = mean(MSE_lfsr)
toc