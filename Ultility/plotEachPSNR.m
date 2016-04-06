% function plotEachPSNR(result,path,subrate);
sub = subrate;

% testImage = [1, 2, 6, 7];
% for mmm = 	i
%     switch testImage(mmm) 
%         case 1
%             path = 'lenna';
%         case 2
%             path = 'barbara';
%         case 3
%             path = 'boats';
%         case 4
%             path = 'house';
%         case 5
%             path = 'peppers';
%         case 6
%             path ='cameraman';
%         case 7
%             path = 'leaves';
%         case 8
%             path = 'parrots';
%         case 9 
%             path = 'mandrill';           
%     end
% end;
% ---------------------------------------
eachPSNR = result.eachPSNR(:);
figure(1);
subplot(1,2,1);
plot(eachPSNR(:))
title([path ' Sub' num2str(sub) ' PSNR Curve']);
xlabel('iteration');
ylabel('PSNR');
grid on;
% ---------------------------------------
eachErr = result.eachErr(:);
% figure(2);R
subplot(1,2,2);
plot(eachErr(:))
title([path ' Sub' num2str(sub) ' Err Curve']);
xlabel('iteration');
ylabel('PSNR');
axis([0 300 0 0.001])
grid on;
% ---------------------------------------