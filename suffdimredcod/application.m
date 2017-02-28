
datasrd = load('datasrdrealfull.tab');
pini = 0;
n1 = size(datasrd, 1);
p = size(datasrd, 2) - 1;
for itr = 1:100
[datasrd0, ind] = datasample(datasrd, (n1 +1)/2, 'Replace', false); 
y0 = datasrd(ind, 1);
x0 = datasrd(ind, 2:(p+ 1));
datasrd1 =datasrd(setdiff(1:size(datasrd, 1),ind), :);
newy = datasrd1(:, 1);
newx = datasrd1(:, 2:(p + 1));
n = size(y0, 1);

x = x0; y = y0; dim = 1;
ss = cov(x,0); mu = mean(x);
z = (x - ones(n,1)*mu)*(inv(ss))^(1/2);

temp =[1, 1.11803278688525, -0.139959016393443, 0.531762295081967, -0.80983606557377, 0.995491803278688];% median(mbeta(:, :, 3))%

temp = repmat(temp, 1, 10);
temp = reshape(temp, p, 10);
betapre = temp;

for i= 1:4

    if i==1
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i), ...
        div1(:, :, i)]=select_main_SIR(z,y,dim,beta0,1, pini, betapre);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n)
%	mbeta(itr, :, i) = betahat(:, :, i)
    end
    if i==2
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i), ...
        div1(:, :, i)]=select_main_SAVE(z,y,dim,beta0,1, pini, betapre);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n)
 %      mbeta(itr, :, i) = betahat(:, :, i)
    end
    if i==3
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i), ...
        div1(:, :, i)]=select_main_DR(z,y,dim,beta0,1, pini, betapre);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n)
  %    mbeta(itr, :, i) = betahat(:, :, i)
    end
    if i==4
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i), ...
        div1(:, :, i)]=select_main_PHD(z,y,dim,beta0,1, pini, betapre);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n)
   %    mbeta(itr, :, i) = betahat(:, :, i)
    end
    i
end
div
otherresult = struct('LIC', LIC, 'ee', ee, 'eflag', eflag, 'dis', ...
                     dis, 'betahat', betahat, 'div1', div1, 'div', div);       

    
%save(sprintf('app_result_3_mave.mat'), 'otherresult'); 

for i = 1:4
cuvx =  x0 * otherresult.betahat(:, 1, i) ;
cuvy = y0;
cuvnewx = newx * otherresult.betahat(:, 1, i)  ;
%fity = csaps(cuvx, cuvy,'xx',  cuvnewx)%fit([cuvx(:, 1)],  cuvy, 'lowess');
[m, m0] = reg_local_linear_smoothing(cuvx, cuvy, 2, 1, cuvnewx)
%cross = median((m0 - newy).^2); 
cross(itr, i) = mean((m0 - newy).^2);

%cuvx =  x0 * otherresult.betahat(:, 2:3, 3);
%cuvy = y0;
%cuvnewx = newx * otherresult.betahat(:, 1:2, i);
%fity = fit([cuvx(:, 1) cuvx(:, 2)],  cuvy, 'poly23');
%fity = feval(fity, [cuvnewx(:, 1), cuvnewx(:, 2)]);
%cross(itr, i) = median((fity - newy).^2);
end
end
%save(sprintf('app_3_10_cross.mat'), 'cross'); 
allOneString = sprintf('%d,' , otherresult.betahat(:, 1, 2))
