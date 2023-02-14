function affichage_file(filename, type, meastype)
% funtion to reproduce the differeent transition phase in the Paper.
switch lower(type)
    case 'probability'
        load(filename,'outerror','mlist','n');
        siz=size(outerror);
        D=floor(siz(1));
        col={'b--+','r--+','k--+','c--o','m--o'};
        f=figure;
        for i=1:5
            %i=i+(i==3);
            var=1:D;
            data=outerror(var,i);
            plot(mlist(var)./n,data,col{i},'LineWidth',1.5);
            axis tight;
            hold on;
        end
        legend('Wirtinger Flow','Polyak subgradient $\ell_1$','Polyak subgradient $\ell_2^2$','MD random init','MD with spectral init','Interpreter','Latex','FontSize',10);
        xlabel('$\frac{m}{n}$','Interpreter','Latex','FontSize',17), ylabel('Recovery Rate','Interpreter','Latex','FontSize',15);
        title('Probability of success','Interpreter','Latex','FontSize',15);
        %set(f,'Units','Inches');Number of masks $D=
        %pos = get(f,'Position');
        %set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        %print(f,'phase_CDP','-dpdf','-r0');
        hold off;
        
    case 'monte-carlo'
        load(filename,'outerror','nlist','mlist');
        siz=size(outerror);
        D=floor(siz(2));
        %nlist= 8:10:128;
        %mlist=128:12:2*floor(nlist(end)*log(nlist(end)));
        %mstart=floor(nlist(end)/10);
        %var=1:D
        data=outerror; 
        tit={'Wirtinger flow','Polyak SG','MD random init','MD with spectral init'};
	if strcmp(lower(meastype),'real gaussian')	a=1;legtrans='$n\log(n)$';
	elseif strcmp(lower(meastype),'real cdp')	a=2;legtrans='$n\log(n)^2$';
	else disp('Unknown measurement type');
	end
	
        f=figure(1);
        for i=1:4
            %i=i+(i==3);
            subplot(1,4,i)
            h=pcolor(nlist,mlist,data(:,:,i)');shading interp;box on
	    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            colormap(gray);
            %colormap(flipud(gray));
            hold on;
            plot(nlist,nlist.*(log(nlist).^a),'r');%,nlist,2*nlist.*log(nlist),'b');
	    hold off
            legend(legtrans,'Interpreter','Latex','FontSize',10);
	    legend boxoff
            %colorbar;
            xlabel('$n$','Interpreter','Latex','FontSize',15), ylabel('$m$','Interpreter','Latex','FontSize',15);
            title(tit{i},'Interpreter','Latex','FontSize',15);
            %axis tight;
        end
        
end
end




% load(filename,'outerror','nlist','mlist');
% siz=size(outerror);
% D=floor(siz(2));
% %nlist= 8:10:128;
% %mlist=128:12:2*floor(nlist(end)*log(nlist(end)));
% %mstart=floor(nlist(end)/10);
% var=1:D;
% %var(1)=1;
% data=outerror(:,:,:);
% tit={'Wirtinger flow','Polyak SG','Mirror-Descent','Mirror-Descent'};
% f=figure(1);
% for i=1:4
%     i=i+(i==3);
%     subplot(1,3,i-(i==4))
%     pcolor(nlist,mlist(var),data(:,var,i)');
%     shading interp;
%     colormap(flipud(gray));
%     hold on;
%     plot(nlist,2*nlist.*log(nlist),'r');%,nlist,2*nlist.*log(nlist),'b');
%     legend('','$2n\log(n)$','Interpreter','Latex','FontSize',10);
%     %colormap(gray);
%     colorbar;
%     xlabel('$n$','Interpreter','Latex','FontSize',15), ylabel('$m$','Interpreter','Latex','FontSize',15);
%     title(tit{i},'Interpreter','Latex','FontSize',15);
%     axis;
% end
