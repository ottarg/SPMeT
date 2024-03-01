function status = myOutputFcn(t, y, flag)
    status = 0; % Initialize status to 0 to continue integration
    threshold = 2; % Define a threshold for the first component of y
    
    switch flag
        case 'init' % Initialize plot at the start of integration
%             hold on;
%             xlabel('Time');
%             ylabel('y(1)');
%             title('Solution Progress');
            
        case '' % During integration
%             plot(t, y(1), 'b.'); % Plot the first component of y
%             drawnow;
%             if y(1) > threshold % Check if the first component exceeds the threshold
%                 status = 1; % Set status to 1 to stop integration
%             end
            
            
        case 'done' % At the end of integration
            hold off;
    end
end