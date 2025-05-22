%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P_act_opt_one_initial_cond.m ver.250522%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
clear global
close all

%parameter settings
U_para=2;
cu=0.1;
cd=-1;
mu0=0.38;
mu2=-0.01;
sigma2=0.08;
nu=U_para+2i*cu;
gamma=1+1i*cd;

T=5;
dt_num=200;
step_num_max=T*dt_num;
dt=1/dt_num;
t_axis=linspace(0,T,step_num_max+1);
if rem(T,1)~=0
    fprintf('T need to be an integer.\n')
    pause
end

addpath(genpath('DMSUITE'));
n=220;
b=0.239;

T_tilde=5;
index_tilde=[87; 97; 107];
kappa=[-0.3+0.8i; 1.2-0.3i; 0.6+0.2i];
p_tilde=size(index_tilde,1);
if rem(T_tilde,1)~=0
    fprintf('T_tilde need to be an integer.\n')
    pause
end

r_list=[10 5];
p_list=[3 8];
rp_num_max=size(r_list,2);
index_select_manual_list=zeros(max(p_list),rp_num_max);
index_select_manual_list(:,1)=[75; 121; 143; 0; 0; 0; 0; 0];
index_select_manual_list(:,2)=[68; 81; 89; 92; 106; 123; 131; 152];

sigma_min=10^(-10);

disp_every=4;
z_min=-30;
z_max=30;
x_min=-0.5;
x_max=2;

scrsz=get(0,'ScreenSize');
fig_size=[700 600];
font_name='Times New Roman';
font_size=32;
line_width=2.5;
line_width_dashed=2.5;
line_width_thick=9;
marker_size=15;

color_list=zeros(15,3);
color_list(1,:)=[1.00 0.00 0.00];
color_list(2,:)=[1.00 0.65 0.00];
color_list(3,:)=[0.92 0.87 0.20];
color_list(4,:)=[0.50 0.90 0.00];
color_list(5,:)=[0.00 0.70 0.00];
color_list(6,:)=[0.00 1.00 1.00];
color_list(7,:)=[0.00 0.70 0.80];
color_list(8,:)=[0.00 0.10 1.00];
color_list(9,:)=[0.68 0.00 1.00];
color_list(10,:)=[1.00 0.00 1.00];
color_list(11,:)=[1.00 0.70 1.00];
color_list(12,:)=[0.77 0.63 0.49];
color_list(13,:)=[0.68 0.38 0.24];
color_list(14,:)=[0.75 0.75 0.75];
color_list(15,:)=[0.00 0.00 0.00];

%Discretization of the model
[z,D]=herdif(n,2,b);
D1=D(:,:,1);
D2=D(:,:,2);

M=diag(mu0-cu^2+mu2*(z.^2)/2);
A=-nu*D1+gamma*D2+M;
B=zeros(n,n);
for row_num=1:n
    for col_num=1:n
        B(row_num,col_num)=1/sqrt(2*pi*sigma2)*exp(-(z(row_num,1)-z(col_num,1))^2/(2*sigma2));%ガウス分布に従う
    end
end
C=zeros(n,n);
C(1,1)=sqrt((z(2,1)-z(1,1))/2);
for row_num=2:n-1
    C(row_num,row_num)=sqrt((z(row_num+1,1)-z(row_num-1,1))/2);
end
C(n,n)=sqrt((z(n,1)-z(n-1,1))/2);

%Computation of response without input
alpha_tilde=zeros(n,1);
alpha_tilde(index_tilde)=kappa;
x0_nocontrol=expm(A*T_tilde)*B*alpha_tilde;
x_rec_nocontrol=zeros(n,step_num_max+1);
x_rec_nocontrol(:,1)=x0_nocontrol;
e_Adt=expm(A*dt);
for step_num=1:step_num_max
    x_rec_nocontrol(:,step_num+1)=e_Adt*x_rec_nocontrol(:,step_num);
end
y_rec_nocontrol=C*x_rec_nocontrol;
yT_nocontrol_norm=norm(y_rec_nocontrol(:,step_num_max+1));

for rp_num=1:rp_num_max
    %Singular value decomposition
    r=r_list(1,rp_num);
    p=p_list(1,rp_num);
    index_select_manual=index_select_manual_list(1:p,rp_num);
    filename_parameter=sprintf('r%dp%dT%d_',r,p,T);
    phi=C*expm(A*T)*B;
    [U,S,V]=svds(phi,r);
    W=V*S;

    for type_obj=1:3
        close all
        %Actuator selection
        if type_obj==1
            obj=sum(abs(W).^2,2);
            [obj_select,index_select]=maxk(obj,p);
        elseif type_obj==2
            if p<=r
                [obj_select,index_select]=F_act_dg_under1(W,p);
            else
                [obj_select_under,index_select_under]=F_act_dg_under1(W,r);
                [obj_select,index_select]=F_act_dg_over1(W,p,obj_select_under,index_select_under);
            end
        else
            index_select=index_select_manual;
            obj_select=zeros(p,1);
        end
        z_select=z(index_select,1);
        
        %Computation of response with control
        index_select_sort=sort(index_select);
        phis=phi(:,index_select_sort.');
        [Ut,St,Vt]=svd(phis);
        Stp=zeros(p,n);
        for l=1:p
            if St(l,l)>sigma_min
                Stp(l,l)=1/St(l,l);
                l_used_max=l;
            end
        end
        y0T_nocontrol=C*expm(A*T)*x0_nocontrol;
        alphas=-Vt*Stp*Ut'*y0T_nocontrol;
        Bs=B(:,index_select_sort);
        Bs_alphas=Bs*alphas;
        x0_control=x0_nocontrol+Bs_alphas;
        x_rec_control=zeros(n,step_num_max+1);
        x_rec_control(:,1)=x0_control;
        for step_num=1:step_num_max
            x_rec_control(:,step_num+1)=e_Adt*x_rec_control(:,step_num);
        end
        y_rec_control=C*x_rec_control;
        yT_control_norm=norm(y_rec_control(:,step_num_max+1));
        if type_obj==1
            y_rec_control_tr=y_rec_control;
        elseif type_obj==2
            y_rec_control_det=y_rec_control;
        elseif type_obj==3
            y_rec_control_manual=y_rec_control;
        end
    
        %Writing results to text file
        if type_obj==1
            filename_obj_type='tr_';
        elseif type_obj==2
            filename_obj_type='det_';
        else
            filename_obj_type='manual_';
        end
        fileID=fopen(strcat(filename_parameter,filename_obj_type,'result.txt'),'w');
        fprintf(fileID,'Selected points (in selected order)\n');
        for k=1:p
            fprintf(fileID,sprintf('point%d:z=%.15f (index=%d,obj=%.15f)\n',k,z_select(k,1),index_select(k,1),obj_select(k,1)));
        end
        fprintf(fileID,'\n');
        fprintf(fileID,'Singular values of Phis\n');
        for l=1:p
            fprintf(fileID,sprintf('mode%d:sigma=%.15f\n',l,St(l,l)));
        end
        fprintf(fileID,sprintf('l_used_max=%d\n',l_used_max));
        fprintf(fileID,'\n');
        fprintf(fileID,'Norm of input\n');
        fprintf(fileID,sprintf('||alpha||_2=%.15f\n',norm(alphas)));%alphaとalphasの2ノルムは等しい
        fprintf(fileID,'\n');
        fprintf(fileID,'Norm of terminal state\n');
        fprintf(fileID,sprintf('Without control:||y(T)||_2=%.15f\n',yT_nocontrol_norm));
        fprintf(fileID,sprintf('With control:||y(T)||_2=%.15f\n',yT_control_norm));
        fclose(fileID);
    
        %Showing actuator locations with singular vectors
        for k_UVW=1:3
            if k_UVW==1
                fig_position_hor=50;
                UVW=U;
                ylabel_uvw='$\phi_\ell$';
                filename_uvw='phi(all)_';
            elseif k_UVW==2
                fig_position_hor=100;
                UVW=V;
                ylabel_uvw='$\psi_\ell$';
                filename_uvw='psi(all)_';
            else
                fig_position_hor=150;
                UVW=W;
                ylabel_uvw='$\xi_\ell$';
                filename_uvw='xi(all)_';
            end
            for k_ria=1:3
                if k_ria==1
                    fig_position_ver=0.1*scrsz(4);
                    UVW_ria=real(UVW);
                    ylabel_ria='Real part';
                    filename_ria='real';
                elseif k_ria==2
                    fig_position_ver=0.15*scrsz(4);
                    UVW_ria=imag(UVW);
                    ylabel_ria='Imaginary part';
                    filename_ria='imag';
                else
                    fig_position_ver=0.2*scrsz(4);
                    UVW_ria=abs(UVW);
                    ylabel_ria='Absolute value';
                    filename_ria='abs';
                end
                figure('Position',[fig_position_hor fig_position_ver fig_size]);
                hold on
                if rp_num==1
                    p1=plot(z,UVW_ria(:,1),'Color',color_list(1,:),'LineWidth',line_width);
                    p2=plot(z,UVW_ria(:,2),'Color',color_list(2,:),'LineWidth',line_width);
                    p3=plot(z,UVW_ria(:,3),'Color',color_list(3,:),'LineWidth',line_width);
                    p4=plot(z,UVW_ria(:,4),'Color',color_list(4,:),'LineWidth',line_width);
                    p5=plot(z,UVW_ria(:,5),'Color',color_list(5,:),'LineWidth',line_width);
                    p6=plot(z,UVW_ria(:,6),'Color',color_list(6,:),'LineWidth',line_width);
                    p7=plot(z,UVW_ria(:,7),'Color',color_list(7,:),'LineWidth',line_width);
                    p8=plot(z,UVW_ria(:,8),'Color',color_list(8,:),'LineWidth',line_width);
                    p9=plot(z,UVW_ria(:,9),'Color',color_list(9,:),'LineWidth',line_width);
                    p10=plot(z,UVW_ria(:,10),'Color',color_list(10,:),'LineWidth',line_width);
                elseif rp_num==2
                    p1=plot(z,UVW_ria(:,1),'Color',color_list(1,:),'LineWidth',line_width);
                    p2=plot(z,UVW_ria(:,2),'Color',color_list(2,:),'LineWidth',line_width);
                    p3=plot(z,UVW_ria(:,3),'Color',color_list(3,:),'LineWidth',line_width);
                    p4=plot(z,UVW_ria(:,4),'Color',color_list(4,:),'LineWidth',line_width);
                    p5=plot(z,UVW_ria(:,5),'Color',color_list(5,:),'LineWidth',line_width);
                end
                for k=1:p
                    za=z_select(k,1);
                    plot([za za],[x_min x_max],'--','Color','k','LineWidth',line_width_dashed)%全て黒で示したい場合
                end
                xlabel('$z$','Interpreter','latex')
                ylabel(strcat(ylabel_ria,'~of entries in~',ylabel_uvw),'Interpreter','latex');
                xlim([z_min z_max])
                ylim([x_min x_max])
                ax=gca;
                ax.FontName=font_name;
                ax.FontSize=font_size;
                box on
                grid on
                if rp_num==1
                    legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10],{'~$\ell=1$','~$\ell=2$','~$\ell=3$','~$\ell=4$','~$\ell=5$','~$\ell=6$','~$\ell=7$','~$\ell=8$','~$\ell=9$','~$\ell=10$'},'Location','northeast','Interpreter','latex');
                elseif rp_num==2
                    legend([p1 p2 p3 p4 p5],{'~$\ell=1$','~$\ell=2$','~$\ell=3$','~$\ell=4$','~$\ell=5$'},'Location','northeast','Interpreter','latex');
                end
                filename=strcat(filename_parameter,filename_obj_type,filename_uvw,filename_ria,'.fig');
                saveas(gcf,filename)
                filename=strcat(filename_parameter,filename_obj_type,filename_uvw,filename_ria,'.pdf');
                exportgraphics(gcf,filename)
            end
        end
    end
    
    %Showing the initial and terminal output
    for k_0T=1:2
        if k_0T==1
            y0T_nocontrol=y_rec_nocontrol(:,1);
            y0T_control_tr=y_rec_control_tr(:,1);
            y0T_control_det=y_rec_control_det(:,1);
            y0T_control_manual=y_rec_control_manual(:,1);
            ylabel_0T='$y(0)$';
            filename_0T='_y(0)';
        else
            y0T_nocontrol=y_rec_nocontrol(:,step_num_max+1);
            y0T_control_tr=y_rec_control_tr(:,step_num_max+1);
            y0T_control_det=y_rec_control_det(:,step_num_max+1);
            y0T_control_manual=y_rec_control_manual(:,step_num_max+1);
            ylabel_0T='$y(T)$';
            filename_0T='_y(T)';
        end 
        figure('Position',[600 0.1*scrsz(4) fig_size]);
        hold on
        p1=plot(z,abs(y0T_nocontrol),'Color',[0.7 0.7 0.7],'LineStyle',':','LineWidth',line_width_thick);
        p2=plot(z,abs(y0T_control_manual),'Color','g','LineStyle','-','LineWidth',line_width);
        p3=plot(z,abs(y0T_control_tr),'Color','b','LineStyle','-','LineWidth',line_width);
        p4=plot(z,abs(y0T_control_det),'Color','r','LineStyle','-','LineWidth',line_width);
        xlabel('$z$','Interpreter','latex')
        ylabel(strcat('Absolute value of entries in~',ylabel_0T),'Interpreter','latex');
        xlim([z_min z_max])
        ylim([x_min x_max])
        ax=gca;
        ax.FontName=font_name;
        ax.FontSize=font_size;
        box on
        grid on
        legend([p1 p2 p3 p4],{'~w/o input','~Random','~Trace','~Det.'},'Location','northeast','Interpreter','latex');
        filename=strcat(filename_parameter,'trdetman',filename_0T,'.pdf');
        exportgraphics(gcf,filename)
    end
    
    %Showing the time variation of output
    for step_num=0:disp_every:step_num_max
        y_nocontrol=y_rec_nocontrol(:,step_num+1);
        y_control_tr=y_rec_control_tr(:,step_num+1);
        y_control_det=y_rec_control_det(:,step_num+1);
        y_control_manual=y_rec_control_manual(:,step_num+1);
        t=dt*step_num;
        f1=figure('Position',[650 0.1*scrsz(4) fig_size]);
        hold on
        p1=plot(z,abs(y_nocontrol),'Color',[0.7 0.7 0.7],'LineStyle',':','LineWidth',line_width_thick);
        p2=plot(z,abs(y_control_manual),'Color','g','LineStyle','-','LineWidth',line_width);
        p3=plot(z,abs(y_control_tr),'Color','b','LineStyle','-','LineWidth',line_width);
        p4=plot(z,abs(y_control_det),'Color','r','LineStyle','-','LineWidth',line_width);
        title(sprintf('$t=%.2f$',t),'Interpreter','latex')
        xlabel('$z$','Interpreter','latex')
        ylabel('Absolute value of entries in $y(t)$','Interpreter','latex');
        xlim([z_min z_max])
        ylim([x_min x_max])
        ax=gca;
        ax.FontName=font_name;
        ax.FontSize=font_size;
        box on
        grid on
        legend([p1 p2 p3 p4],{'~w/o input','~Random','~Trace','~Det.'},'Location','northeast','Interpreter','latex');
        filename=strcat(filename_parameter,'trdetman_y(t)_',sprintf('t%.2f',t),'.png');
        saveas(gcf,filename,'png');
        close(f1)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%