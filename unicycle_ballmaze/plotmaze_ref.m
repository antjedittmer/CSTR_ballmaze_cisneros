function plotmaze_ref(Cx,Cy,radius,refx,state_sim,animate_plot,refy,XX,XXref,save_to_file)
% plotmaze_ref plots the unicycle path in the ballmaze. Run Ballmaze_main
% or load DataForPlot if available.

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for presentation at DLR-FT Seminar: Antje Dittmer


if nargin < 10
    save_to_file = 0;
end

figDir = 'figDir';
if ~isfolder(figDir)
    mkdir(figDir)
end

animationDir = fullfile(figDir,'animationDir');
if ~isfolder(animationDir)
    mkdir(animationDir)
end


%% Plot the obstacles
th = 0:pi/50:2*pi;
for i = 1:length(Cx)
    xunit = radius * cos(th) + Cx(i);
    yunit = radius * sin(th) + Cy(i);
    plot(xunit, yunit,'k-','LineWidth',2);
    hold on
end
% plot ([-0.8 -0.8],[-1.1 1.1],'k-','LineWidth',1.5);
% plot ([-0.8 1.1],[-1.1 -1.1],'k-','LineWidth',1.5);
% axis([-0.2 2.2 -0.2 2.2])

cl = lines; % get line color
plot(refx(1:800),refy(1:800),'k--'); % Reference line
plot(state_sim(:,1),state_sim(:,2),'LineWidth',1.5,'color', cl(4,:)); %,'o','MarkerSize',5)

axis([-1.1 1.1 -1.1 1.1])
axis equal

set_title_labels(cl)

%% animate maze

if animate_plot

    filename = fullfile(animationDir,'ballmaze_velin');
    len = 0.1;
    th = 0:pi/50:2*pi;
    for j = 2:2:size(state_sim,1)-1
        for i = 1:length(Cx)
            xunit = radius * cos(th) + Cx(i);
            yunit = radius * sin(th) + Cy(i);
            plot(xunit, yunit,'k-','LineWidth',2);
            hold on
        end

        plot(refx(1:800),refy(1:800),'k--'); % Reference line
        plot(state_sim(:,1),state_sim(:,2),'LineWidth',1.5,'color', cl(4,:)); %,'o','MarkerSize',5)

        axis([-1.1 1.1 -1.1 1.1])
        axis equal
        set_title_labels(cl)


        fwx = state_sim(j,1) + len*cos(state_sim(j,4));
        bkx = state_sim(j,1) - len*cos(state_sim(j,4));
        fwy = state_sim(j,2) + len*sin(state_sim(j,4));
        bky = state_sim(j,2) - len*sin(state_sim(j,4));


        plot([bkx fwx],[bky fwy],'LineWidth',8)
        plot(XX(10:9:end,j),XX(11:9:end,j),'bo','MarkerFaceColor','blue','Markersize',4)
        %     plot(XX(6:5:end,j),XX(7:5:end,j),'bo','MarkerFaceColor','blue','Markersize',4)
        %     plot(XX(11:10:end,j),XX(12:10:end,j),'*')

        plot(XXref(1:9:end,j),XXref(2:9:end,j),'o','Color',[0.8 0.8 0.8],'Markersize',4);
        %     plot(XXref(1:5:end,j),XXref(2:5:end,j),'o','Color',[0.8 0.8 0.8],'Markersize',4);

        plot(refx(j),refy(j),'rx','Markersize',10);
        %     plot ([-0.8 -0.8],[-1 1],'k-','LineWidth',1.5);
        %     plot ([-0.8 1],[-1 -1],'k-','LineWidth',1.5);

        hold off
        axis([-1.1 1.1 -1.1 1.1])
        axis equal
        if (save_to_file == 1)
            frame = getframe(1);
            im = frame2im(frame);

            imwrite(im,strcat(filename,num2str(j/2,'%04d'),'.png'));

        end
        pause(0.025)
    end
end

if animate_plot && save_to_file
    movieFile = fullfile(figDir, 'output_movie.avi');

    filePattern = fullfile(animationDir, '*.png'); % Pattern for PNG-files
    pngFiles = dir(filePattern); % Alle PNG-Dateien im Verzeichnis finden

    % Create a video object to write the video
    outputVideo = VideoWriter(movieFile); % Use MPEG-4 format
    outputVideo.FrameRate = 10; % Set frame rate

    open(outputVideo); % Open the video object

    % Loop through all PNG files and add them to the video
    for k = 1:length(pngFiles)
        % Read each image
        baseFileName = pngFiles(k).name;
        fullFileName = fullfile(animationDir, baseFileName);
        img = imread(fullFileName);

        % Check if the image is colored
        if size(img, 3) == 1
            disp(['Image ', baseFileName, ' is grayscale.']);
        end

        % Write the image to the video
        writeVideo(outputVideo, img);
    end

    % Close the video object
    close(outputVideo);
end


strFig = 'Fig03_ballmaze';
print(fullfile(figDir,strFig), '-dpng');
%print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-depsc');

function set_title_labels(cl)
xlabel('x (m)'); ylabel('y (m)');

tileSt = {
    'Reference trajectory (--)', ...
    'Vehicle closed loop trajecory (-)', ...
    };

% Define colors (rows correspond to strings)
cl1 = [
    0     0     0;   % black
    cl(4,:) %violet
    ];

% Helper function to convert RGB to LaTeX-style string
rgbStr = @(v) sprintf('%.2f,%.2f,%.2f', v);

% Build LaTeX string
titleCell = cell(numel(tileSt),1);
for i = 1:numel(tileSt)
    aStr = ['{\color[rgb]{', rgbStr(cl1(i,:)), '}', tileSt{i}, '}'];
    titleCell{i} = aStr;
end
title(titleCell);

set(findall(gcf,'-property','FontSize'),'FontSize',11.5);