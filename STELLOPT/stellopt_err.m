function stellopt_err( stel_data,jac_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


dex=length(stel_data.iter);
dex_stel=dex;
jac = jac_data.jac';
Npar = jac_data.n;
Npnt = jac_data.m;
Nfit = jac_data.n;
y_dat = stel_data.TARGETS(dex,:);
y_fit = stel_data.VALS(dex,:);
y_sig = stel_data.SIGMAS(dex,:);
jac   = jac_data.jac';
%for i = 1:Npar
%    jac(:,i) = jac(:,i).*y_sig';
%end
chisq = ((stel_data.TARGETS-stel_data.VALS)./stel_data.SIGMAS).^2;
chisq_tot = sum(chisq,2);
delta_y=(y_dat-y_fit)./y_sig;
weights_sq = (Npnt-Nfit+1)./((delta_y'*delta_y) * ones(Npnt,1));
weights_sq(delta_y' == 0) = 1.0E-18;
JtWJ       = jac' * ( jac.* ( weights_sq * ones(1,Npar)));
covar      = inv(JtWJ);
sigma_p    = abs(sqrt(diag(covar)));
sigma_y    = zeros(Npnt,1);
for i = 1:Npnt
    sigma_y(i) = jac(i,:)*covar*jac(i,:)';
end
sigma_y    = abs(sqrt(sigma_y).*y_sig');
sigma_yp   = abs(sqrt(weights_sq+sigma_y).*y_sig');
%sigma_p = sqrt(diag(Vp));
%sigma_y = sqrt(abs(diag(jac*Vp*jac')));
%sigma_yp = sqrt(W2+diag(jac*Vp*jac'));

subplot(4,2,[1 3]);
plot_stellopt(jac_data);
subplot(4,2,2);
plot(stel_data.iter,100.*(stel_data.VALS-stel_data.TARGETS)./stel_data.TARGETS);
subplot(4,2,4);
plot(stel_data.iter,chisq_tot,'+');
subplot(4,2,5);
errorbar(y_dat,y_sig,'ok');
%plot(y_dat,'o');
hold on;
%plot(y_fit,'+');
err=1.96.*real(sigma_y)';
errorbar(y_fit,err,'+r');
%plot(y_fit+1.96.*real(sigma_y)','-r');
%plot(y_fit-1.96.*real(sigma_y)','-r');
xlim([-5 Nfit+5]);
subplot(4,2,7);
temp=abs(real(sigma_y));
temp(temp ==0) = min(temp(temp>0));
semilogy(real(temp),'-r');
xlim([-5 Nfit+5]);
subplot(4,2,[6 8]);
hist((y_dat-y_fit)./y_sig)

% Now make individual plots
dex  = 1;
dex2 = 1;
cdex = [1.0 0 0];
while dex > 0
    tname=jac_data.target_name{dex};
    dex1=dex;
    dex2=dex1;
    for i=dex+1:size(jac_data.target_name,2)
        if strcmp(tname,jac_data.target_name{i})
            dex2=i;
        end
    end
    fig=figure;
    if (dex1==dex2)
        errorbar(y_dat(dex1:dex2),y_sig(dex1:dex2),'ok');
        hold on
        errorbar(y_fit(dex1),1.96.*real(sigma_y(dex1)),'xr');
    else
        if strcmp(strtrim(tname),'Electron Density') 
            temp=y_dat(dex1:dex2);
            temp2=y_sig(dex1:dex2);
            temp3=sigma_y(dex1:dex2);
            temp4=y_fit(dex1:dex2);
            %[s, idex]=sort(stel_data.NE_S(dex_stel,:),2,'ascend'); % Sort by flux
            [s, idex]=sort(stel_data.NE_R(dex_stel,:),2,'ascend'); % Sort by r
            temp=temp(idex);
            temp2=temp2(idex);
            temp3=temp3(idex);
            temp4=temp4(idex);
            size(temp);
            y_dat(dex1:dex2)=temp;
            y_sig(dex1:dex2)=temp2;
            sigma_y(dex1:dex2)=temp3;
            y_fit(dex1:dex2)=temp4;
        elseif strcmp(strtrim(tname),'Electron Temperature') 
            temp=y_dat(dex1:dex2);
            temp2=y_sig(dex1:dex2);
            temp3=sigma_y(dex1:dex2);
            temp4=y_fit(dex1:dex2);
            %[s, idex]=sort(stel_data.TE_S(dex_stel,:),2,'ascend'); % Sort by flux
            [s, idex]=sort(stel_data.TE_R(dex_stel,:),2,'ascend'); % Sort by r
            temp=temp(idex);
            temp2=temp2(idex);
            temp3=temp3(idex);
            temp4=temp4(idex);
            size(temp);
            y_dat(dex1:dex2)=temp;
            y_sig(dex1:dex2)=temp2;
            sigma_y(dex1:dex2)=temp3;
            y_fit(dex1:dex2)=temp4;
        elseif strcmp(strtrim(tname),'Ion Temperature') 
            temp=y_dat(dex1:dex2);
            temp2=y_sig(dex1:dex2);
            temp3=sigma_y(dex1:dex2);
            temp4=y_fit(dex1:dex2);
            [s, idex]=sort(stel_data.TI_S(dex_stel,:),2,'ascend');
            temp=temp(idex);
            temp2=temp2(idex);
            temp3=temp3(idex);
            temp4=temp4(idex);
            size(temp);
            y_dat(dex1:dex2)=temp;
            y_sig(dex1:dex2)=temp2;
            sigma_y(dex1:dex2)=temp3;
            y_fit(dex1:dex2)=temp4;
        elseif strcmp(strtrim(tname),'Motional Stark Effect Diagnostic') 
            temp=y_dat(dex1:dex2);
            temp2=y_sig(dex1:dex2);
            temp3=sigma_y(dex1:dex2);
            temp4=y_fit(dex1:dex2);
            [s, idex]=sort(stel_data.MSE_S(dex_stel,:),2,'ascend');
            temp=temp(idex);
            temp2=temp2(idex);
            temp3=temp3(idex);
            temp4=temp4(idex);
            size(temp);
            y_dat(dex1:dex2)=temp;
            y_sig(dex1:dex2)=temp2;
            sigma_y(dex1:dex2)=temp3;
            y_fit(dex1:dex2)=temp4;
        else
            s=1:(dex2-dex1+1);
        end
        errorbar(s,y_dat(dex1:dex2),y_sig(dex1:dex2),'ok');
        hold on
        plot(s,y_fit(dex1:dex2),'+b')
        for j = 1:10
            yup=y_fit(dex1:dex2)+1.96.*real(sigma_y(dex1:dex2))'/j;
            ydn=y_fit(dex1:dex2)-1.96.*real(sigma_y(dex1:dex2))'/j;
            x=s;
            patch([x fliplr(x)],[yup fliplr(ydn)],[x fliplr(x)].*0-1 ,'red','FaceAlpha',0.2,'EdgeColor','none');
            %plot(y_fit(dex1:dex2)+1.96.*real(sigma_y(dex1:dex2))'/j,'Color',[1.0 (100-j)/200 (100-j)/200])
            %plot(y_fit(dex1:dex2)-1.96.*real(sigma_y(dex1:dex2))'/j,'Color',[1.0 (100-j)/200 (100-j)/200])
        end
        hold off
        %plot(y_fit(dex1:dex2)+1.96.*real(sigma_y(dex1:dex2))','-r')
        %plot(y_fit(dex1:dex2)-1.96.*real(sigma_y(dex1:dex2))','-r')
    end
    title(tname);
    if (dex2 == size(jac_data.target_name,2))
        dex = -1;
    else
        dex = dex2+1;
    end
end
end

