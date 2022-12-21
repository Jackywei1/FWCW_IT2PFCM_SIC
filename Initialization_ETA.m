function [ETA,C]= Initialization_ETA (all_pixels,all_pixels_xi,gamma_y,gamma_xi,alpha,beta,m,k,K)
[N,d] = size(all_pixels);
[C,U] = fcm(all_pixels_xi,k,[NaN 100 0.0001 0]);
mf = U.^m;
Distence = Distance_Function(all_pixels,all_pixels_xi,C,gamma_y,gamma_xi,alpha,beta);
ETA = zeros(k,d);
for j=1:k
    ETA(j,:) = K*sum(repmat((mf(j,:))',[1 d]).*reshape(Distence(j,:,:),[N,d]))./sum(repmat((mf(j,:))',[1 d]));
end
ETA = sum(ETA,2);
end

