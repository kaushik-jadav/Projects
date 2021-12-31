clear all
clc
init = 0.001;
last = 1;
xd_init = 1;
xd_last = 25;
N = 1;
tf = 100;
ts = 0.001;
t_nw=[];
hold on
cont = 0;
for ii = 1:100
    u_x_ic = 0;
    u_phi_ic =0;
    x_dot_ic = 0;
    phi_dot_ic = 0;
    x_ic = init + (last+init)*randn(N,1);
    x_dotdot_ic = init + (last+init)*randn(N,1);
    xd_bar = init + (last+init)*randn(N,1);
    phi_ic = init + (last+init)*randn(N,1);
    phi_d_bar = init + (last+init)*randn(N,1);
    K_axp = 0.1;
    gamma_ax = 0.01;
    gamma_aphi = 0.01;
    c_x=init + (last+init)*rand(N,1);
    f_xd=0.25+init + (last+init)*rand(N,1);
    f_phi_d=0.25+init + (last+init)*rand(N,1);
    c_phi=init + (last+init)*rand(N,1);
%     m_x=0.9 + (1)*rand(N,1);
    m_x = init + (last+init)*rand(N,1);
    m_phi = init + (last+init)*rand(N,1);
    l_x = 1 + (1)*rand(N,1);
    %l_phi = 1 + (1)*rand(N,1);
    l_phi =init + (last+init)*rand(N,1);;
    theta_x = [m_x, m_phi, c_x, m_phi*l_x]';
    theta_phi = [m_phi*l_phi, m_phi*l_phi^2, c_phi]';
    g = 9.8;
    a_hat_x_ic =init + (last+init)*rand(N,1);
    a_hat_phi_ic =init + (last+init)*rand(N,1);
    phi_ic=init + (last+init)*rand(N,1);
    s_x = 10;
    s_phi = 10;
    a_x = init + (last+init)*rand(N,1);
    a_phi = init + (last+init)*rand(N,1);
    theta_hat_x_ic=[init + (last+init)*rand(N,1), init + (last+init)*rand(N,1), init + (last+init)*rand(N,1), init + (last+init)*rand(N,1)]';
    theta_hat_phi_ic=[init + (last+init)*rand(N,1), init + (last+init)*rand(N,1), init + (last+init)*rand(N,1)]';
    alpha_x =15;
    beta_x = 15;
    alpha_phi =15;
    beta_phi = 15;
    gamma_x = [0.01, 0.01, 0.01, 0.01];
    gamma_phi = [0.1, 0.1, 0.1];
%     K_AX = [0.1, 0.1];
%     K_ax = diag(K_AX);
%     gamma_x = init + (last+init)*rand(N,1);
%     gamma_phi = init + (last+init)*rand(N,1);
    gamma_x = diag(gamma_x);
    K_clx = 0.001*gamma_x;
    gamma_phi = diag(gamma_phi);
    K_clp = 0.001*gamma_phi;
    t = [0:ts:tf];
    xd = xd_bar*sin(2*pi*f_xd*t);
    xd_dot = 2*(pi*f_xd)*xd_bar*cos(2*pi*f_xd*t);
    xd_dotdot = -(2*pi*f_xd)^2*xd_bar*sin(2*pi*f_xd*t);
    xd_dotdotdot = -(2*pi*f_xd)^3*xd_bar*cos(2*pi*f_xd*t);
    xdslx = timeseries(xd',t');
    xd_dotslx = timeseries(xd_dot',t');
    xd_dotdotslx = timeseries(xd_dotdot',t');
    xd_dotdotdotslx = timeseries(xd_dotdotdot',t');
    phi_d = phi_d_bar*sin(2*pi*f_phi_d*t);
    phi_d_dot = 2*pi*f_phi_d*phi_d_bar*cos(2*pi*f_phi_d*t);
    phi_d_dotdot = -(2*pi*f_phi_d)^2*phi_d_bar*sin(2*pi*f_phi_d*t); 
    phi_d_dotdotdot = -(2*pi*f_phi_d)^3*phi_d_bar*cos(2*pi*f_phi_d*t); 
    pdslx = timeseries(phi_d',t');
    pd_dotslx = timeseries(phi_d_dot',t');
    pd_dotdotslx = timeseries(phi_d_dotdot',t');
    pd_dotdotdotslx = timeseries(phi_d_dotdotdot',t');

    try
     out = sim('cl2')
     Ex(ii,:) = out.e_x_t.Data';
     time = out.e_x_t.Time';
     Rx(ii,:) = out.r_x_t.Data';
     Ep(ii,:) = out.e_phi_t.Data';
     Rp(ii,:) = out.r_phi_t.Data';
     theta_tx1 = out.theta_x_tilda_t.Data';
     theta_tp1 = out.theta_phi_tilda_t.Data';
     a_tx1 = out.a_x_tilda_t.Data';
     a_tp1 = out.a_phi_tilda_t.Data';
     f_forwardx(ii,:) = out.ff_x.Data';
     f_forwardp(ii,:) = out.ff_p.Data';
     f_backx(ii,:) = out.fb_x.Data';
     f_backp(ii,:) = out.fb_p.Data';
     udx(ii,:) = out.ud_x_t.Data';
     udp(ii,:) = out.ud_phi_t.Data';
     total_input(ii,:) = [out.ud_x_t.Data' out.ud_phi_t.Data'];
     total_feedforward(ii,:) = [out.ff_x.Data' out.ff_p.Data'];
     total_feedback(ii,:) = [out.fb_x.Data' out.fb_p.Data'];
     diff_ff(ii,:) = total_input(ii,:) - total_feedforward(ii,:);
     diff_fb(ii,:) = total_input(ii,:) - total_feedback(ii,:);
     u_t_x(ii,:) = out.u_tilda_x_t.Data';
     u_t_p(ii,:) = out.u_tilda_phi_t.Data';
     t_t_x1(ii,:) = out.theta_x_tilda_t.Data(:,1)';
     t_t_x2(ii,:) = out.theta_x_tilda_t.Data(:,2)';
     t_t_x3(ii,:) = out.theta_x_tilda_t.Data(:,3)';
     t_t_x4(ii,:) = out.theta_x_tilda_t.Data(:,4)';
     t_t_p1(ii,:) = out.theta_phi_tilda_t.Data(:,1)';
     t_t_p2(ii,:) = out.theta_phi_tilda_t.Data(:,2)';
     t_t_p3(ii,:) = out.theta_phi_tilda_t.Data(:,3)';
     a_t_x(ii,:) = out.a_x_tilda_t.Data';
     a_t_p(ii,:) = out.a_phi_tilda_t.Data';
     ttt_x(ii,:)=[ t_t_x1(ii,:),  t_t_x2(ii,:), t_t_x3(ii,:),  t_t_x4(ii,:)];
     ttt_p(ii,:)=[ t_t_p1(ii,:),  t_t_p2(ii,:), t_t_p3(ii,:)];
%      t_d_x(ii,:) = ttt_x(ii,:)'.*gamma_x.*ttt_x(ii,:);
%      t_d_p(ii,:) = ttt_p(ii,:)'*gamma_phi.*ttt_p(ii,:);
     

    figure(1)
    subplot(2,2,1)
    hold on
    plot(out.e_x_t.Time,out.e_x_t.Data)
%     norm(ii) = sum(out.e_x_t.Data.^2)^0.5;
    xlabel('time',fontsize =15);
    ylabel('e_x',fontsize =15);
    grid on
    subplot(2,2,2)
    plot(out.r_x_t.Time, out.e_x_t.Data)
    hold on
    xlabel('time',fontsize =15);
    ylabel('r_x',fontsize =15);
    grid on
    subplot(2,2,3)
    hold on
    plot(out.e_phi_t.Time, out.e_phi_t.Data)
    xlabel('time',fontsize =15);
    ylabel('e_\phi',fontsize =15);
    grid on
    subplot(2,2,4)
    hold on
    plot(out.r_phi_t.Time, out.r_phi_t.Data)
    xlabel('time',fontsize =15);
    ylabel('r_\phi',fontsize =15);
    grid on
    saveas(gcf,'e and r.png')



    figure(2)
    subplot(2,2,1)
    hold on
    plot(out.theta_x_tilda_t.Time, out.theta_x_tilda_t.Data)
    xlabel('time',fontsize =15);
    ylabel('theta tilda x',fontsize =15);
    ylim([-4,4]);
    grid on
    subplot(2,2,2)
    hold on
    plot(out.theta_phi_tilda_t.Time, out.theta_phi_tilda_t.Data)
    xlabel('time',fontsize =15);
    ylabel('theta tilda phi',fontsize =15);
    ylim([-4,4]);
    grid on  
    subplot(2,2,3)
    hold on
    plot(out.a_x_tilda_t.Time, out.a_x_tilda_t.Data)
    xlabel('time',fontsize =15);
    ylabel('a tilda x',fontsize =15);
%     ylim([-1000,1000]);
    grid on
    subplot(2,2,4)
    hold on
    plot(out.a_phi_tilda_t.Time, out.a_phi_tilda_t.Data)
    xlabel('time',fontsize =15);
    ylabel('a tilda phi',fontsize =15);
%     ylim([-1000,1000]);
    grid on
    saveas(gcf,'theta and a.png')

    catch
        norm(ii)=1000;
        ii;
        cont = cont +1;
    end
end
hold off
for jj = 1:size(Ex,2)
    EXn(jj) = sum(Ex(:,jj).^2)^0.5;
    RXn(jj) = sum(Rx(:,jj).^2)^0.5;
    EPn(jj) = sum(Ep(:,jj).^2)^0.5;
    RPn(jj) = sum(Rp(:,jj).^2)^0.5;
    E_tx(jj) = sum(theta_tx1(:,jj).^2)^0.5;
    E_tp(jj) = sum(theta_tp1(:,jj).^2)^0.5;
    E_ax(jj) = sum(a_tx1(:,jj).^2)^0.5;
    E_ap(jj) = sum(a_tp1(:,jj).^2)^0.5;
    f_fx(jj) = sum(f_forwardx(:,jj).^2)^0.5;
    f_fp(jj) = sum(f_forwardp(:,jj).^2)^0.5;
    f_bx(jj) = sum(f_backx(:,jj).^2)^0.5;
    f_bp(jj) = sum(f_backp(:,jj).^2)^0.5;
    u_d_x(jj) = sum(udx(:,jj).^2)^0.5;
    u_d_p(jj) = sum(udp(:,jj).^2)^0.5;
    tff(jj) = sum(diff_ff(:,jj).^2)^0.5;
    tfb(jj) = sum(diff_fb(:,jj).^2)^0.5;
    utx(jj) = sum(u_t_x(:,jj).^2)^0.5;
    utp(jj) = sum(u_t_p(:,jj).^2)^0.5;
    ttx1(jj) = sum(t_t_x1(:,jj).^2)^0.5;
    ttx2(jj) = sum(t_t_x2(:,jj).^2)^0.5;
    ttx3(jj) = sum(t_t_x3(:,jj).^2)^0.5;
    ttx4(jj) = sum(t_t_x4(:,jj).^2)^0.5;
    ttp1(jj) = sum(t_t_p1(:,jj).^2)^0.5;
    ttp2(jj) = sum(t_t_p2(:,jj).^2)^0.5;
    ttp3(jj) = sum(t_t_p3(:,jj).^2)^0.5;
    atx(jj) = sum(a_t_x(:,jj).^2)^0.5;
    atp(jj) = sum(a_t_p(:,jj).^2)^0.5;
%     tdx(jj) = sum(t_d_x(:,jj).^2)^0.5;
%     tdp(jj) = sum(t_d_p(:,jj).^2)^0.5;
    tax(jj,:) = [ttx1(jj) ttx2(jj) ttx3(jj) ttx4(jj)];
    tap(jj,:) = [ttp1(jj) ttp2(jj) ttp3(jj)];
    taxn(jj) = sum(ttt_x(:,jj).^2)^0.5;
    tapn(jj) = sum(ttt_p(:,jj).^2)^0.5;
    LFC(jj) = 0.5*EXn(jj)^2 + 0.5*RXn(jj)^2 + 0.5*utx(jj)^2 + 0.5*taxn(jj)^2 + 0.5*atx(jj)^2;
    LFP(jj) = 0.5*EPn(jj)^2 + 0.5*RPn(jj)^2 + 0.5*utp(jj)^2 + 0.5*tapn(jj)^2 + 0.5*atp(jj)^2;
%     mi_x(jj) = 0.5*min(out.e_x_t.Data)^2 + 0.5*min(out.r_x_t.Data)^2 + 0.5*min(out.u_tilda_x_t)^2 + 0.5*min(tax(jj)).^2 + 0.5*min(out.a_x_tilda_t)^2;
%     ma_x(jj) = 0.5*max(out.e_x_t.Data)^2 + 0.5*max(out.r_x_t.Data)^2 + 0.5*max(out.u_tilda_x_t)^2 + 0.5*max(tax(jj)).^2 + 0.5*max(out.a_x_tilda_t)^2;
%     ma_x(jj) = max(LFC(jj));
%     ma_p(jj) = max(LFC(jj));
    LFCUB(jj) = 0.5*EXn(jj)^2 + 0.5*RXn(jj)^2 + 0.5*utx(jj)^2;
    LFPUB(jj) = 0.5*EPn(jj)^2 + 0.5*RPn(jj)^2 + 0.5*utp(jj)^2;

end
% total_input = [u_d_x u_d_p];
% total_feedforward = [f_fx f_fp];
% total_feedback = [f_bx f_bp];

figure(3)
hold on
plot(time, EXn);
plot(time, RXn);
plot(time, EPn);
plot(time, RPn);
legend('e (cart)','r (cart)', 'e (pendulum)','r (pendulum)',fontsize =13);
title('Norm of Tracking errors', fontsize=13);
saveas(gcf,'norm of e and r.png')


figure(4)
hold on
plot(time, E_tx);
plot(time, E_tp);
plot(time, E_ax);
plot(time, E_ap);
legend('theta-tilde (cart)','theta-tilde (pendulum)', 'a-tilde (cart)','a-tilde (pendulum)',fontsize =13);
title('Norm of Estimation errors', fontsize=13);
saveas(gcf,'norm of theta and a tilde.png')


figure(5)
hold on
plot(time, f_fx);
plot(time, f_bx);
plot(time, u_d_x);
legend('feed forward (cart)', 'feedback (cart)', 'ud (cart)' ,fontsize =15);
title({['Total input, Error feedback portion,'] ['Estimated feedforward portion for Collar']}, fontsize=13);
saveas(gcf,'input and ff and fb collar.png')

figure(6)
hold on
plot(time, f_fp);
plot(time, f_bp);
plot(time, u_d_p);
legend('feed forward (pendulum)','feedback (pendulum)','ud (pendulum)',fontsize =13);
title({['Total input, Error feedback portion,'] ['Estimated feedforward portion for Pendulum']}, fontsize=13);
saveas(gcf,'input and ff and fb pendulum.png')

figure(7)
hold on
plot(time,tff);
title({['norm of difference of'] ['total input and estimated feedforward portions']},fontsize =13)
saveas(gcf,'diff input and ff.png')


figure(8)
hold on
plot(time,tfb);
title({['norm of difference of'] ['total input and error feedback portions']},fontsize =13)
saveas(gcf,'diff input and fb.png')

figure(9)
hold on
subplot(2,1,1)
hold on
plot(time, LFC);
plot(time,LFCUB);
% plot(time, mi_x);
% plot(time, ma_x);
legend('Lyapunov Function','Lower Bound',fontsize =13);
title({['Lyapunov function and Bounds'] ['for Collar']},fontsize =13)
saveas(gcf,'LF-collar.png')
%title('min','max')
subplot(2,1,2)
hold on
% plot(time, mi_p);
% plot(time, ma_p);
plot(time, LFP);
plot(time,LFPUB);
legend('Lyapunov Function','Lower Bound',fontsize =13);
title({['Lyapunov function and Bounds'] ['for Pendulum']},fontsize =13)
saveas(gcf,'LF-pendulum.png')