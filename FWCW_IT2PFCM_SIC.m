function [C1,C2,UT_upper,UT_low, E1,E2]=FWCW_IT2PFCM_SIC(k,all_pixels,all_pixels_xi,C,gamma_y,gamma_xi,alpha,beta,q,p,m1,m2, ...
    theta1,theta2,Cp,Cf,ETA,t_max,beta_memory)
%
[N,d] = size(all_pixels);

% Initialize cluster weights and feature weights
Z1=ones(1,k)/k;
Z2 = Z1;
W1=ones(k,d)/d;
W2 = W1;

C1 = C;
C2 = C;

% Set initial value
Iter=1;
O_F_old=inf;

% main circulation
while 1
    % a) Update U_upper ºÍ U_low.
    Distence1 = Distance_Function(all_pixels,all_pixels_xi,C1,gamma_y,gamma_xi,alpha,beta);
    Distence2 = Distance_Function(all_pixels,all_pixels_xi,C2,gamma_y,gamma_xi,alpha,beta);
    for j=1:k
        W1_q = transpose(W1(j,:).^q); W1_q(W1_q==inf)=0;
        W2_q = transpose(W2(j,:).^q); W2_q(W2_q==inf)=0;
        dNK1(:,j) = Z1(1,j).^p * reshape(Distence1(j,:,:),[N,d]) * W1_q;
        dNK2(:,j) = Z2(1,j).^p * reshape(Distence2(j,:,:),[N,d]) * W2_q;
    end

    tmp_U1 = zeros(N,k);
    tmp_U2 = zeros(N,k);
    for j=1:k
        tmp2_U1 = (dNK1./repmat(dNK1(:,j),1,k)).^(1/(m1-1));
        tmp2_U1(tmp2_U1==inf)=0; tmp2_U1(isnan(tmp2_U1))=0;
        tmp_U1=tmp_U1+tmp2_U1;

        tmp2_U2 = (dNK2./repmat(dNK2(:,j),1,k)).^(1/(m2-1));
        tmp2_U2(tmp2_U2==inf)=0; tmp2_U2(isnan(tmp2_U2))=0;
        tmp_U2=tmp_U2+tmp2_U2;
    end
    U1 = transpose(1./tmp_U1);
    U1(isnan(U1))=1; U1(U1==inf)=1;

    U2 = transpose(1./tmp_U2);
    U2(isnan(U2))=1; U2(U2==inf)=1;

    U_upper = max(U1,U2);
    U_low = min(U1,U2);

    % b) Update T_upper and T_low.
    for j=1:k
        W1_q = transpose(W1(j,:).^q); W1_q(W1_q==inf)=0;
        dTNK1 =Z1(1,j).^p .* reshape(Distence1(j,:,:),[N,d]).*repmat(W1_q',[N,1]);
        tmp3_1(j,:) = Cp * (theta1.^2) * (theta1.^2 - 1) * sum(dTNK1,2)./repmat(ETA(j,:),[N,1]);

        W2_q = transpose(W2(j,:).^q); W2_q(W2_q==inf)=0;
        dTNK2 =Z2(1,j).^p .* reshape(Distence2(j,:,:),[N,d]).*repmat(W2_q',[N,1]);
        tmp3_2(j,:) = Cp * (theta2.^2) * (theta2.^2 - 1) * sum(dTNK2,2)./repmat(ETA(j,:),[N,1]);
    end

    T1 = (exp(-1.*lambertw(tmp3_1))).^(theta1/(theta1.^2-1));
    si = find (tmp3_1 == Inf); T1(si) = 1;

    T2 = (exp(-1.*lambertw(tmp3_2))).^(theta2/(theta2.^2-1));
    si = find (tmp3_2 == Inf); T2(si) = 1;

    T_upper = max(T1,T2);
    T_low = min(T1,T2);

    % c) Upper bound and Lower bound Membership
    U_T1=(Cf.*U_upper)+(Cp.*T_upper);
    U_T2=(Cf.*U_upper)+(Cp.*T_low);
    U_T3=(Cf.*U_low)+(Cp.*T_upper);
    U_T4=(Cf.*U_low)+(Cp.*T_low);

    DUM1=max(U_T1,U_T2);
    DUM2=max(DUM1,U_T3);
    UT_upper=max(DUM2,U_T4);

    DUM1=min(U_T1,U_T2);
    DUM2=min(DUM1,U_T3);
    UT_low=min(DUM2,U_T4);

    %
    C1_old=C1;
    C2_old=C2;

    % d)Karni-Mendel Alg
    for km=1:k
        UT_mean=(UT_upper+UT_low)/2;
        maxUT_mean = max(UT_mean);

        index_c = find(UT_mean(km, :) == maxUT_mean);
        Y= all_pixels_xi(index_c,:);

        a=UT_low(km,index_c);
        b=UT_upper(km,index_c);
        a=a(:);
        b=b(:);
        F=cat(2,a,b);

        [XLeft,XRight,~,~]=KM_Alg(F,Y);

        C1(km,:)=XLeft';
        C2(km,:)=XRight';
    end

    % e) Defuzzification method
    % f) check termination condition
    E1(Iter) = norm (C1_old - C1, 1);
    E2(Iter) = norm (C2_old - C2, 1);

    fprintf('Iteration count = %d, Termination measure value => Lower= %f    Upper=%f\n', Iter, E1(Iter),E2(Iter));

    if Iter>=t_max || (E1(Iter) < 1e-3 ) && (E2(Iter) < 1e-3 )

        fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
        fprintf('Converging after %d iterations.\n',Iter);
        fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');

        break;

    end
    %
    W1_old=Z1;
    z1_old=W1;
    W2_old=Z2;
    z2_old=W2;

    % Update W
    for j=1:k
        dWkm1(j,:) = UT_low(j,:) * reshape(Distence1(j,:,:),[N,d]);
        dWkm2(j,:) = UT_upper(j,:) * reshape(Distence2(j,:,:),[N,d]);
    end

    tmp1_W1 = zeros(k,d);
    tmp1_W2 = zeros(k,d);
    for j=1:d
        tmp2_W1 = (dWkm1./repmat(dWkm1(:,j),1,d)).^(1/(q-1));
        tmp2_W1(tmp2_W1==inf)=0; tmp2_W1(isnan(tmp2_W1))=0;
        tmp1_W1=tmp1_W1+tmp2_W1;

        tmp2_W2 = (dWkm2./repmat(dWkm2(:,j),1,d)).^(1/(q-1));
        tmp2_W2(tmp2_W2==inf)=0; tmp2_W2(isnan(tmp2_W2))=0;
        tmp1_W2=tmp1_W2+tmp2_W2;
    end
    W1 = 1./tmp1_W1;
    W1(isnan(W1))=1; W1(W1==inf)=1;

    W2 = 1./tmp1_W2;
    W2(isnan(W2))=1; W2(W2==inf)=1;

    if nnz(dWkm1==0)>0
        for j=1:k
            if nnz(dWkm1(j,:)==0)>0
                W1(j,find(dWkm1(j,:)==0)) = 1/nnz(dWkm1(j,:)==0);
                W1(j,find(dWkm1(j,:)~=0)) = 0;
            end
        end
    end

    if nnz(dWkm2==0)>0
        for j=1:k
            if nnz(dWkm2(j,:)==0)>0
                W2(j,find(dWkm2(j,:)==0)) = 1/nnz(dWkm2(j,:)==0);
                W2(j,find(dWkm2(j,:)~=0)) = 0;
            end
        end
    end

    % Update Z
    for j=1:k
        W1_q = transpose(W1(j,:).^q); W1_q(W1_q==inf)=0;
        Dw1(1,j) = transpose(W1_q) * transpose(reshape(Distence1(j,:,:),[N,d])) *  transpose(UT_low(j,:)) ;

        W2_q = transpose(W2(j,:).^q); W2_q(W2_q==inf)=0;
        Dw(1,j) = transpose(W2_q) * transpose(reshape(Distence2(j,:,:),[N,d])) *  transpose(UT_upper(j,:)) ;
    end

    tmp_Z1 = sum((repmat(Dw1,k,1)./transpose(repmat(Dw1,k,1))).^(1/(p-1)));
    tmp_Z1(tmp_Z1==inf)=0; tmp_Z1(isnan(tmp_Z1))=0;
    Z1 = 1./tmp_Z1;
    Z1(isnan(Z1))=1; Z1(Z1==inf)=1;

    tmp_Z2 = sum((repmat(Dw,k,1)./transpose(repmat(Dw,k,1))).^(1/(p-1)));
    tmp_Z2(tmp_Z2==inf)=0; tmp_Z2(isnan(tmp_Z2))=0;
    Z2 = 1./tmp_Z2;
    Z2(isnan(Z2))=1; Z2(Z2==inf)=1;

    if nnz(Dw1==0)>0
        Z1(find(Dw1==0)) = 1/nnz(Dw1==0);
        Z1(find(Dw1~=0)) = 0;
    end

    if nnz(Dw==0)>0
        Z2(find(Dw==0)) = 1/nnz(Dw==0);
        Z2(find(Dw~=0)) = 0;
    end

    % Memory effect.
    Z1=(1-beta_memory)*Z1+beta_memory*W1_old;
    W1=(1-beta_memory)*W1+beta_memory*z1_old;
    Z2=(1-beta_memory)*Z2+beta_memory*W1_old;
    W2=(1-beta_memory)*W2+beta_memory*z1_old;

    Iter=Iter+1;
end
end
