function [] = pzmap_testCL(system,H,gain,feedin,feedout)
% evaluate and plot the pole-zero map for the closed loop system for
% different values of the gain

    [~, n] = size(gain);
    [m1, n1, ~] = size(H);
    [~,n2] = size(feedin);

    figure
    for i = 1:n
        %     if n1 == n2
        system_CL = feedback(system,gain(i)*H,feedin,feedout);

        [P,Z] = pzmap(system_CL);
        plot(real(P(:)),imag(P(:)),'x',real(Z(:)),imag(Z(:)),'o');hold on
        xlabel('Real axis (s^{-1})');ylabel('Imaginary Axis (s^{-1})');
        %         clear P Z
        %     else
        %         system_CL = feedback(system,gain(i)*H(:,1+(i-1)*m1:m1+(i-1)*m1),feedin,feedout);
        %
        %         [P,Z] = pzmap(system_CL);
        %         plot(real(P(:)),imag(P(:)),'x',real(Z(:)),imag(Z(:)),'o');hold on
        %         xlabel('Real axis (s^{-1})');ylabel('Imaginary Axis (s^{-1})');
        %         clear P Z
        %     end
    end
    str = {strcat('gain = ' , num2str(gain(1)))};  % at the end of first loop, z being loop output
    str = [str , strcat('gain = ' , num2str(gain(1)))]; % after 2nd loop
    for i = 2:n
        str = [str , strcat('gain = ' , num2str(gain(i)))]; % after 2nd loop
        str = [str , strcat('gain = ' , num2str(gain(i)))]; % after 2nd loop
    end
    legend(str{:})
end
