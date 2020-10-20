function [data, V, D, mu] = pca_pablo(A)
    %Take data in tidy format: rows = samples, columns = observations
    
    mu = mean(A,1);             %save mean to add it back
    A = A-mu;                   %mean subtract data on observations
    A = A(:,std(A,1) > 0);      %keep only features with any variance
    cov_mat = (A'*A)/size(A,1); %find the covariance matrix
    corr_mat = corrcov(cov_mat);%This standardizes all variables so that larger variables dont dominate
    [V,D] = eig(corr_mat);      %identify vectors capturing max variance (eignevectors of covariance matrix) with their associated vairance captured (eigenvalues)
    [D,order] = sort(diag(-D)); %give the order from most variance to least (ordering PCs in order of importance), and store variance captured by new order
    D = -D;                     %made D negative just to find descending order, but values should be positive again
    V = V(:,order);             %reorder the PCs in the eigenbassi so that first dimension is first PC, etc.
    data = A*V;                 %project data onto these Principal Axes
end