function [supphat,xhat] = omp(y,A,K)

AA = [];
M = length(y);
N = length(A(1,:));
y_orig = y;

A_nrmlzd = normalize(A,'norm',2);

itr = 1;
if norm(y_orig) ~= 0
    for itr = 1:K
        

%         tmp1 = abs((A./D)'*y);
        tmp1 =abs((A_nrmlzd)'*y);
                    
    % --------------------------------------
        zz = find(max(tmp1)==tmp1);
        if length(zz)>1
            ind_max(itr) = zz(1);
%             itr2 = itr2 + 1;
        else
            ind_max(itr) = zz;
        end
        if sum(ind_max==ind_max(itr))>1
            ind_max(itr) = [];
            break;
        end
          
        AA = [AA,A(:,ind_max(itr))];

        xhat1 = (pinv(AA))*y_orig;

        r = y_orig - AA*xhat1;
            

        y=r;

    end
    
    xhat = zeros(N,1);
    for itr1 = 1:length(ind_max)
        xhat(ind_max(itr1),1) = xhat1(itr1);
    end
    
else
    
    xhat = zeros(N,1);
end
supphat = sort(reshape(ind_max,[numel(ind_max),1]),'descend');
end
