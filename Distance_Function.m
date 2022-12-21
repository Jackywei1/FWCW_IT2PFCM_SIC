function distance = Distance_Function (all_pixels,all_pixels_xi,C,gamma_y,gamma_xi,alpha,beta)
[N,d] = size(all_pixels);
[k,~] = size(C);
distance = zeros(k,N,d);
for j=1:k
    distance(j,:,:) = alpha.*((1-exp((-1.*repmat(gamma_y,N,1)).*((all_pixels-repmat(C(j,:),N,1)).^2)))) + ...
           beta.*((1-exp((-1.*repmat(gamma_xi,N,1)).*((all_pixels_xi-repmat(C(j,:),N,1)).^2))));
end
end
