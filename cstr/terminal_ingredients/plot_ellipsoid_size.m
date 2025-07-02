for ii = 1:20
    xx(ii)=prod(svd(Zval{ii}([0.4,350])));
end

for ii = 1:20
    yyy=eig(Pval{ii}([0.5,0.3]));
    yy(ii)=yyy(2);
end

% svd(Pval{1})

subplot(2,1,1)
plot(xx,'LineWidth',1.5)
xlim([1 20])
% axis([1 20 0 0.5])
subplot(2,1,2)
plot(yy,'LineWidth',1.5)
xlabel('F-Z iterations')
xlim([1 12])