function DuneShrubMovie(ShrubDomain,Steps,Type)

% Create and save movie frames of island domain over time
% Type: use 1 for Shrub Age, 2 for Percent Cover

% Ian Reeves 19 Sep 18

axi = 10;
BarLabel = 'Shrub Age (yr)';
FrameType = 'AgeFrame';

if Type == 2
    axi = 100;
    BarLabel = 'Shrub Percent Cover';
    FrameType = 'CoverFrame';
end

for t = 1:Steps

    Tframe = squeeze(ShrubDomain(t,:,:));
    Tframe = rot90(Tframe,3);
    imagesc(Tframe); figure(gcf);
    set(gcf,'rend','painters','pos',[150,300,1500,150]); % Position and size of frame
    set(gca,'Ydir','normal')
    ylabel('Island Width');
    xlabel('Island Length');
    bar = colorbar;
    caxis([0 axi]);
    ylabel(bar,BarLabel,'Rotation', 270);
    time = strcat(['t = ',num2str(t),' yr']);
    annotation('textbox',[.13 .05 .05 .15],'String',time,'FitBoxToText','off','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
    pmovie(t) = getframe;
    set(gcf, 'InvertHardCopy', 'off','PaperPositionMode','auto');
    outputfilename = ['C:/Reeves/DuneShrub/MovieFrames/' FrameType num2str(t)]; % Save frames  !!Hardwired to IRBR desktop computer!!
    print('-dpng',outputfilename)

end

% % Play movie
% movie(pmovie,1,3) % Plays 1 time, 2 fps
% save(['C:/Reeves/DuneShrub/MovieFrames/pmovie.mat'],'pmovie')


end

