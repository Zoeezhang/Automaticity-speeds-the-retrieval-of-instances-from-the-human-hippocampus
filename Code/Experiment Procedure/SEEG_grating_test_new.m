function SEEG_grating_test
addpath('io64');
run config_io.m
Screen('Preference', 'SkipSyncTests', 1)
p.SessStart = GetSecs;
Screen('CloseAll')
close all;
clc;
warning('off','MATLAB:dispatcher:InexactMatch');                            % turn off case-mismatch manager (it's annoying)
%key setup
EscapeKey = KbName('q');
continuekey=KbName('space');
funsetting = 0;obj=1;
% Dimensions for Screen
p.sub=inputdlg({'编号','姓名','姓别(1男/2女)','年龄','记忆负荷(1/2)','Group1/2','Practice(1yes/0no)'});
p.subID=p.sub{1};
p.subName=p.sub{2};
p.load=str2num(p.sub{5});
p.group=str2num(p.sub{6});
p.practice=str2num(p.sub{7});
p.vDist = 57;                                                         % 57; % viewing distance (cm)
rect=Screen('Rect', max(Screen('Screens')));                                %resolution，0 0 2560 1440 pixel
[width, height]=Screen('DisplaySize', max(Screen('Screens')));              %553 311 mm
p.px = (width/rect(3))/10;                                                  % 0.0216 pixel size (cm) for BenQ monitors

if p.practice
    fNameall=[pwd,'\Subject Data\subID-',p.subID,'-',p.subName,'-load',num2str(p.load),'-prac.mat'] ;
    fNameall_retrieval=[pwd,'\Subject Data\subID-',p.subID,'-',p.subName,'-load',num2str(p.load),'-retri-prac.mat'] ;
else
    fNameall=[pwd,'\Subject Data\subID-',p.subID,'-',p.subName,'-load',num2str(p.load),'.mat'] ;
    fNameall_retrieval=[pwd,'\Subject Data\subID-',p.subID,'-',p.subName,'-load',num2str(p.load),'-retri.mat'] ;
end

if exist( fNameall_retrieval ,'file')
    error('The data file exists! Please enter a different subject ID.');
    Screen('CloseAll');
end
load (fNameall);
stim1 = stim;
clear stim
stim.tgOriraw = stim1.tgOriraw;
stim.tg2Oriraw = stim1.tg2Oriraw;
%--------------------------------------------------------------------------
% General Experimental Parameters
%--------------------------------------------------------------------------
p.nTrials = 30;%p.nStim; % number of trials per block   
p.nBlocks = 6;%p.nTotal/p.nStim;  % number of blocks
if p.practice; p.nBlocks = 1;p.nTrials=6; end % 1 block if practice1
p.windowed =0; % 1 = small,    
p.port = 1;

% setup addresss
address = hex2dec('CEFC');
ioObj=io64;
if funsetting==1
    obj= serial('COM4');
    fopen(obj);
    fwrite(obj,char(0));
end
if p.port
    status=io64(ioObj);
    io64(ioObj,address,0);
end

%--------------------------------------------------------------------------
% Build the stimulus display
%--------------------------------------------------------------------------
AssertOpenGL;                                                               % Make sure we're on an OpenGL-compatible machine
s = max(Screen('Screens'));                                                 % Grab a handle for the display ID (should return 0 for single-monitor setups)
p.nFrames = 60;%  round(FrameRate(s));                                      % grab and save the frame rate

if p.windowed                                                               % if we're debugging, open up a small window
    [w,p.sRect] = Screen('OpenWindow',s,p.foreCol,[0 0 800 600],[],[],[],[20]);%[120 120 800 600]
else
    [w,p.sRect] = Screen('OpenWindow',s,p.foreCol,[],[],[],[],[20]);                         % otherwise, go full-screen
    HideCursor;
    Priority(MaxPriority(w));                                               % set priority to max to discourage interruptions
end

% Compute and store the center of the screen
p.xCenter = (p.sRect(3)-p.sRect(1))/2;
p.yCenter = (p.sRect(4)-p.sRect(2))/2;

% Foreground rectangle
p.foreRect = p.sRect; % set to the full screen

% Fixation rectangle
p.fixRect = CenterRect([0 0 p.fixSize p.fixSize],p.foreRect);

% Make rect for gratings
p.gratingRect = CenterRect([0 0 p.gratingDiameter p.gratingDiameter],p.foreRect);

% New Color Wheel Prefs
p.colorWheelRect = CenterRect([0 0 p.colorWheelDiameter p.colorWheelDiameter],p.foreRect);
p.colorWheelRectouter = CenterRect([0 0 p.colorWheelDiameter+p.colorWheelThickness p.colorWheelDiameter+p.colorWheelThickness],p.foreRect);
p.colorWheelRectinner = CenterRect([0 0 p.colorWheelDiameter-p.colorWheelThickness p.colorWheelDiameter-p.colorWheelThickness],p.foreRect);

% make a dummy call to GetSecs to load the .dll before we need it
dummy = GetSecs; clear dummy;

clear x TempImage youxunhuan
x=[pwd '\youxunhuan.jpg'];
TempImage=imread(x);
youxunhuan=Screen('MakeTexture',w,TempImage); 

clear x TempImage nexttrial
x=[pwd '\nexttrial.jpg']; 
TempImage=imread(x);
nexttrial=Screen('MakeTexture',w,TempImage); 
        
clear x TempImage huiyi
x=[pwd '\session.jpg']; 
TempImage=imread(x);
session=Screen('MakeTexture',w,TempImage);
        
for b = 1:p.nBlocks
    fName=[pwd,'\Subject Data\subID-',p.subID,'-',p.subName,'-load',num2str(p.load),'-b',num2str(b),'-retri.mat'] ;
    % make grating stimulus
    [tg,bg] = makeSmoothGrating(w,p,p.white,1);
    [tg2,~] = makeSmoothGrating(w,p,p.white,1);
    % probe grating ori
    stim.tgOri(b,:)=Shuffle([zeros(1,p.nTrials/2),90.*ones(1,p.nTrials/2)]);% begin with 0 or 90 deg
    stim.tg2Ori(b,:)=Shuffle([zeros(1,p.nTrials/2),90.*ones(1,p.nTrials/2)]);
        
    stim.cue(b,:)=Shuffle([ones(1,p.nTrials/2),2.*ones(1,p.nTrials/2)]);
    %------------------------------------------------------- 
    % Begin Block Loop
    %-------------------------------------------------------
        Screen('FillRect',w,p.foreCol,p.foreRect);            % Draw the foreground window
        Screen('DrawTexture',w,youxunhuan,[],p.foreRect);
        Screen('FillOval',w,p.fixCol,p.fixRect);           % Draw the fixation point    
        Screen('TextSize', w, [55]);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        WaitSecs(.5);

    % Wait for a spacebar press to continue
       while 1
        [keyIsDown,secs,keyCode]=KbCheck;
            if keyIsDown && keyCode(continuekey)
                break;
            elseif  keyIsDown && keyCode(EscapeKey)
                if funsetting == 1
                    fclose(obj);
                end
                Screen('CloseAll')
                return
            end
       end
        
    %-------------------------------------------------------
    % Begin Trial Loop
    %-------------------------------------------------------
    cnt=1;
    for t = 1:p.nTrials;
        
        if ~p.windowed;HideCursor; end

        % ITI
        %----------------------------------------------
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('DrawingFinished',w);
        Screen('Flip',w);
        WaitSecs(p.ITI(randsample(length(p.ITI),1)));   

        % fixation
        %----------------------------------------------
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fixCol,p.fixRect);
        Screen('DrawingFinished',w);
        Screen('Flip',w);  
        if funsetting==1
            fwrite(obj,char(0));
        end
        if p.port ==1
%            outp(address,b);
           io64(ioObj,address,0);
        end
        WaitSecs(p.fixDur);
        
        % second grating
        %----------------------------------------------
        if p.load==2
            % Draw cue
            %---------------------------------------------- 
            text1 = num2str(stim.cue(b,t));
            tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-20];
            Screen('FillRect',w,p.foreCol,p.foreRect);            % Draw the foreground window
            Screen('DrawText', w, text1, tCenter1(1), tCenter1(2), p.white);
            Screen('DrawingFinished', w);
            Screen('Flip', w);
            WaitSecs(p.cue);
            
            % Draw delay
            %---------------------------------------------- 
            Screen('FillRect',w,p.foreCol,p.foreRect);
    %         Screen('FillOval',w,p.fixCol,p.fixRect);
            Screen('DrawingFinished',w);
            Screen('Flip',w);
            WaitSecs(p.cuedelay);
        end
      
        % Get response
        %----------------------------------------------
        [p,stim] = posProbe(p,stim,t,b,w,ioObj,address,funsetting,obj);
 
        % next trial
        %----------------------------------------------
        Screen('DrawTexture',w,nexttrial,[],p.foreRect); 
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        
        % Wait for a spacebar press to continue
       while 1
        [keyIsDown,secs,keyCode]=KbCheck;
            if keyIsDown && keyCode(continuekey)
                break;
            elseif  keyIsDown && keyCode(EscapeKey) 
                if funsetting == 1
                    fclose(obj);
                end
                Screen('CloseAll')
                return
            end
       end

        % save data file at the end of each trial
        save(fName,'p','stim');
        cnt=cnt+1;
    end  % end of trial loop

        % get 2min rest within block
        %----------------------------------------------
        if p.practice~=1 
            Screen('DrawTexture',w,session,[],p.foreRect); 
            Screen('DrawingFinished', w);
            Screen('Flip', w);
            
            WaitSecs(30);
        end
end
% ListenChar(0);
if funsetting == 1
    fclose(obj);
end
p.SessEnd = GetSecs;
save(fNameall_retrieval,'p','stim'); 
Screen('CloseAll');
ShowCursor('Arrow');

function [p, stim] = posProbe(p, stim, t, b, w, ioObj, address, funsetting, obj)
    [tg, ~] = makeSmoothGrating(w, p, p.white, 1); 
    
    stim.rtStart(b, t) = GetSecs;
    stim.noresp(b, t) = 0;
    HideCursor;
    
    % ========== 新增变量 ========== %
    stim.movementStartTime(b, t) = NaN;  % 鼠标开始移动时间
    stim.movementDuration(b, t) = NaN;   % 移动持续时间
    stim.movementEndTime(b, t) = NaN;   % 移动结束时间
    movementThreshold = 5;  % 移动检测阈值(像素)
    % ============================= %
    
    if p.load == 2
        if stim.cue(b, t) == 1
            test_tgOri = stim.tgOri(b, t);
            raw_tgOri = stim.tgOriraw;
        else
            test_tgOri = stim.tg2Ori(b, t);
            raw_tgOri = stim.tg2Oriraw;
        end
    else
        test_tgOri = stim.tgOri(b, t);
        raw_tgOri = stim.tgOriraw;
    end

    % 设置初始位置
    [pointX, pointY] = polar2xy(test_tgOri, p.colorWheelRadius, p.xCenter, p.yCenter);
    SetMouse(pointX, pointY, w);
    prevPos = [pointX, pointY];  % 记录初始位置

    % 呈现probe
    Screen('FillRect', w, p.foreCol, p.foreRect);
    Screen('DrawTexture', w, tg, [], p.gratingRect, test_tgOri, 1, [], []);
    Screen('DrawingFinished', w); 
    Screen('Flip', w);
    
    if funsetting == 1
        fwrite(obj, char(4));
    end
    if p.port == 1
        io64(ioObj, address, 4);
    end
    
    % 主循环
    buttons = [];
    oldAngle = -1;
    while ~any(buttons)
        Screen('DrawTexture', w, tg, [], p.gratingRect, test_tgOri, 1, [], []);

        [curX, curY, buttons] = GetMouse(w);
        curPos = [curX, curY];
        
        % ===== 新增移动检测逻辑 ===== %
        if isnan(stim.movementStartTime(b, t)) && ...
           sqrt(sum((curPos - prevPos).^2)) > movementThreshold
            stim.movementStartTime(b, t) = GetSecs;
        end
        prevPos = curPos;
        % ========================== %
        
        curAngle = xy2polar(curX, curY, p.xCenter, p.yCenter);
       
        Screen('DrawTexture', w, tg, [], p.gratingRect, curAngle, 1, [], []);
        Screen('Flip', w);
        
        tempSecs=GetSecs;
%         if tempSecs-stim.rtStart(b,t)>=3 % no response within 3s
%             stim.noresp(b,t)=1;
%             break;
%         end
    end
    
    % 记录点击时间
    if funsetting == 1
        fwrite(obj, char(0));
    end    
    
    tgOriDiff = curAngle - raw_tgOri;
    if tgOriDiff > 180
        tgOriDiff = tgOriDiff - 360;
    elseif tgOriDiff < -180
        tgOriDiff = tgOriDiff + 360;
    end
    
    if tgOriDiff < -90
        tgOriDiff = tgOriDiff + 180;
    elseif tgOriDiff > 90
        tgOriDiff = tgOriDiff - 180;
    end
    
    stim.tgOriDiff(b, t) = tgOriDiff;
    stim.rtEnd(b, t) = GetSecs;
    stim.movementEndTime(b, t) = stim.rtEnd(b, t);
    stim.movementDuration(b, t) = stim.movementEndTime(b, t) - stim.movementStartTime(b, t);
    stim.rt(b, t) = stim.rtEnd(b, t) - stim.rtStart(b, t);
    
    Screen('FillRect', w, p.foreCol, p.foreRect);
    Screen('DrawingFinished', w); 
    Screen('Flip', w);
    WaitSecs(floor(3-stim.rt(b, t)));
    
    text1 = num2str(stim.tgOriDiff(b, t));
    tCenter1 = [p.xCenter - RectWidth(Screen('TextBounds', w, text1))/2, p.yCenter - 120];
    Screen('FillRect', w, p.foreCol, p.foreRect);
    Screen('FillOval', w, p.fixCol, p.fixRect);
    Screen('DrawText', w, text1, tCenter1(1), tCenter1(2), p.fixCol);
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    WaitSecs(1);

    
%--------------------------------------------------------------------------
function pix = deg2pix(deg,vDist,pixSize)
% Convert degrees of visual angle to pixels for easy specification of
% stimulus size for PsychToolbox. Returns the size of a stimulus in
% pixels:
%
% INPUTS:
% deg: desired stim size in degrees of visual angle.
% vDist: viewing distance.
% pxSize: pixel size (in same units as viewing distance).
%
rad = deg2rad(deg); % convert visual angle from degrees to radians
sz = vDist*tan(rad); % size of stimulus (in same units as p.triVand pxSize)
pix = round(sz/pixSize); % convert to pixels

function [angle, radius] = xy2polar(h,v,centerH,centerV)

  % get polar coordinates
  hdist   = h-centerH;
  vdist   = v-centerV;
  radius     = sqrt(hdist.*hdist + vdist.*vdist)+eps;
  
  % determine angle using cosine (hyp will never be zero)
  angle = acos(vdist./radius)./pi*180;
 
  
  % correct angle depending on quadrant
  angle(hdist == 0 & vdist > 0) = 180;
  angle(hdist == 0 & vdist < 0) = 0;
  angle(vdist == 0 & hdist > 0) = 90;
  angle(vdist == 0 & hdist < 0) = 270;
 
  angle(vdist < 0 & hdist > 0)=180-angle(vdist < 0 & hdist > 0);%一
  angle(vdist > 0 & hdist > 0)=180-angle(vdist > 0 & hdist > 0);%二
    angle(vdist > 0 & hdist < 0)=180+angle(vdist > 0 & hdist < 0);%三
     angle(vdist < 0 & hdist < 0)=180+angle(vdist < 0 & hdist < 0);%四
  
function [x, y] = polar2xy(angle,radius,centerH,centerV) 
  y = round(centerV - radius.*cosd(angle));
  x = round(centerH + radius.*sind(angle));

