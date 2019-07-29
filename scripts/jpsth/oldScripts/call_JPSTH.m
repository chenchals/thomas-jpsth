function [OUTSTRUCT]=call_JPSTH(CellXIn,CellYIn,Align_On,Pre_Time,Post_Time,TLIST,BINSIZE,varargin)
%only 2 JPSTHs
%function [OUTSTRUCT]=call_JPSTH(CellX,CellY,Align_On,Pre_Time,Post_Time,TLIST,BINSIZE [,Saccade_])
%call_JPSTH(DSP03b,DSP04b,Target_,0,300,RF5TrialsE,25,Saccade_);
%*Empty matrix or 0 must be passed for TLIST : for all trials
% NBINS no. of integral bins from Pre_Time to Post_Time
%Check if # trials in cellX and cellY are same
DRAWRAST=0;

%Make sure that the selected range of times is evenly divisible by BinSize
BINSIZE=ceil(BINSIZE);
SelRange=range([Pre_Time Post_Time]);
Remain=rem(SelRange,BINSIZE);
if Remain % binsize is not evenly divisible into the selected range.
   disp('Increasing the selected Pre_time and Post_Time for selected binsize')
   MaxRange=(ceil(SelRange/BINSIZE))*BINSIZE;
   PadRange=(MaxRange-SelRange)/2;
   Pre_Time=Pre_Time-PadRange;
   Post_Time=Post_Time+PadRange;
end
Range=range([Pre_Time Post_Time]);
NBINS=Range; %effectively binwidth =1;
%*****************Check no of arg to see if End event is passed****************
if nargin==8
   EndEventName=inputname(8);
   EndEventName=EndEventName(1:findstr(EndEventName,'_')-1);
   EndEvent=varargin{1};
   EndEventFlag=1;
else
   EndEventName=[];
   EndEventFlag=0;
end
if EndEventFlag
   [EndEvent]=GetEventTime(EndEvent,size(CellXIn,1));
   if ischar(EndEvent)
      error([EndEvent,' ',EndEventName])  
   end
end
%*****************Check no of arg to see if End event is passed****************
%****************Set flag for dropping spikes AFTER End Event*****************
if EndEventFlag
   DropSpikesAfterEndEvent=1;
else
   DropSpikesAfterEndEvent=1;
end
%****************Set flag for dropping spikes AFTER End Event*****************
%****************Set flag for dropping TRIALS AFTER End Event*****************
if DropSpikesAfterEndEvent
   DropTrials=1;
else
   DropTrials=1;
end   
%****************Set flag for dropping TRIALS AFTER End Event*****************
%****************Set flag for smoothing Coincidence histogram*****************
smoothing=1;%Do smmothing of the Histogram
SMOOTHPLOT=1;%Plot smoothed Histogram
%****************Set flag for smoothing Coincidence histogram*****************
%*************************For Manual Scaling***********************************
ManualScale=0;% for manuallly scaling Xcoins;
ScaleXCOINS=[-1.0 1.0];%must span equaaly on both neg and pos side
%ScaleHIST=[0 2.60];%for 20 ms bin width corresponds to 130 Hz =(2.6*1000/20)
ScaleHIST=[0 1.60];%for 20 ms bin width corresponds to 80 Hz =(1.6*1000/20)
allmean=-0.0091;%0.0074(117vs121);%got by mean from all images to be compared
allsd=0.1394;%0.1169;
nSD=5;
ScaleJPST=[allmean-(nSD*allsd) allmean+(nSD*allsd)];
%*************************For Manual Scaling***********************************
%******************Global variables from External Call*************************
global XCellInfo YCellInfo RFPOS XCorrFlag % params from call_eggermont1
global LowerCoincidenceTime UpperCoincidenceTime COLOR SAVEFIG figname
global sizescale
%******************Global variables from External Call*************************
%*****************Test for coincidence histogram (for DEBUGGUNG)***************
TEST=0;
if TEST==1
   UpperCoincidenceTime=((Post_Time-Pre_Time)/NBINS)*10;
   LowerCoincidenceTime=(-((Post_Time-Pre_Time)/NBINS))*10;
end
%*****************Test for coincidence histogram (for DEBUGGUNG)***************
%***************************Set coincidence limits*****************************
if LowerCoincidenceTime, XFiresBeforeY=LowerCoincidenceTime; else XFiresBeforeY=0;end %Above main Diagonal (cellX fires before cellY fires)
if UpperCoincidenceTime, XFiresAfterY=UpperCoincidenceTime; else XFiresAfterY=0;end %Below Main Diagonal (cellX fires after cellY fires)
%***************************Set coincidence limits*****************************
%*************SAVEFIGURE / DO CROSS CORR **************************************
if (~isempty(SAVEFIG) & SAVEFIG==1), FigFile=figname; else FigFile=''; end
if ~isempty(XCorrFlag) & XCorrFlag==1, CallXCorr=1;else, CallXCorr=0;end
%*************SAVEFIGURE / DO CROSS CORR **************************************
%*********************INLINE functions*****************************************
global n2s n2s2
n2s=inline(['num2str(v,','''%4.3f'')'],'v');
n2s2=inline(['num2str(v,','''%4.2f'')'],'v');
%*********************INLINE functions*****************************************
%***************************Variables for figure**********************
CellX=CellXIn;CellY=CellYIn;
sizescale=1.0;
fontscale=sizescale;
linescale=0.05;
delete(allchild(0));


CellXID=inputname(1);CellYID=inputname(2);
%***************************Variables for figure**********************
%******************Check Align variable*******************************
% If input is a named variable 
%if ~isempty(inputname(3))
if ischar(inputname(3))
   Align=inputname(3);
   LL=findstr(Align,'_');
   if ~isempty(LL)
      Align=Align(1:findstr(Align,'_')-1);
   else
      Align=Align;
   end
else  %The 3rd argument could be a scalar or a vector
   if prod(size(Align_On))==0 % empty
      Align='Scalar 0 (ms)';
   elseif prod(size(Align_On))==1 %Scalar
      Align=Align_On;
   elseif prod(size(Align_On))>1
      if prod(size(Align_On))==size(CellX,1) %vector length = no of trials
         Align='Vector (ms)'
      else
         disp('No of Align Times must equal No of Trials')
         disp('Changing Align_On variables to ZERO')
         Align_On=0;
         Align='Scalar 0 (ms)';
      end
   else
   end
end
[Align_On]=GetEventTime(Align_On,size(CellX,1));
if ischar(Align_On)
   error([Align_On])  
end
%***********Check Trial list and no of trials*********************
NTrials=size(CellX,1);%rows
if NTrials~=size(CellY,1)
   error('ROW ERROR: No of Trials (rows) for simultaneously recorded cells must be SAME!')
end
TLIST=nonzeros(TLIST);
if(isempty(TLIST))
   TLIST=(1:size(CellX,1))';
end
if max(TLIST(:))>NTrials
   error('ROW ERROR: Max. Trial no. in TLIST exceeds no. of trials (rows) in the Spike matrices')
end
%***********Check Trial list and no of trials*********************
%***Extract Align Time and EndEventTime for all trials in the Trial list******
NTrials=length(TLIST);
Align_On=Align_On(TLIST);
if EndEventFlag
   EndEvent=EndEvent(TLIST);
else
   %if no end event specified then set End event to be the Post time
   EndEvent=repmat(Post_Time,size(Align_On,1),1);
end

XCellNo=[];
if ~isstruct(XCellInfo),XCellNo=CellXID;else,XCellNo=['#',num2str(XCellInfo.CellNo)];end
YCellNo=[];
if ~isstruct(YCellInfo),YCellNo=CellYID;else,YCellNo=['#',num2str(YCellInfo.CellNo)];end
%***********Check Trial list and extract spike matrix of Trial List*******
COLOR=1;

%Call figure template
invgr=[];gr=[];
NoOFJPSTHs=2;
jpst_fig4(NoOFJPSTHs);
setcolormap(COLOR);
global H_Figure H_JPST H_JPSTX H_JPSTY H_JPSTX2 H_COIN H_CORREL H_INFO
COINS_MaxY=0;


for DIRECTION=1:2
   PP=1;   
   JJ=DIRECTION;
   %***Extract Spike Matrix
   CellX=CellXIn(TLIST,:);
   CellY=CellYIn(TLIST,:);
   
   disp(' ')
   disp(['Doing DIRECTION ',n2s(DIRECTION)])
   disp('')
   %SAWAP ALIGN, END, PRE, POST for backward alignment
   if DIRECTION==2;
      TEMP=[];
      TEMP=Align_On;Align_On=EndEvent;EndEvent=TEMP;clear TEMP;
      TEMP=Pre_Time; Pre_Time=-(Post_Time);Post_Time=-(TEMP); clear TEMP;
      TEMP=Align;Align=EndEventName;EndEventName=TEMP;clear TEMP
   end
   
   %Aligning and sorting
   CX=[];CY=[];
   
   for JJJ=1:NTrials
      Xtemp=[];Ytemp=[];Atime=[];Etime=[];
      Xtemp=CellX(JJJ,:);Ytemp=CellY(JJJ,:);
      Atime=Align_On(JJJ);Etime=(EndEvent(JJJ)-Atime);
      Xtemp=Xtemp-Atime;Ytemp=Ytemp-Atime;
      xtemp=[];ytemp=[];
      %get all spikes between pre and post time
      xtemp=Xtemp(find(Xtemp>=Pre_Time & Xtemp<=Post_Time));
      ytemp=Ytemp(find(Ytemp>=Pre_Time & Ytemp<=Post_Time));
      if isempty(xtemp), xtemp=0; end
      if isempty(ytemp), ytemp=0; end
      CX(JJJ,1:length(xtemp),PP)=xtemp;
      CY(JJJ,1:length(ytemp),PP)=ytemp;
   end
   clear Xtemp Ytemp Atime Etime
   %sort all Variables based on diff times
   SI=[];ST=[];
   ST=abs(EndEvent-Align_On);
   [ST,SI]=sort(ST);
   EndEvent=EndEvent(SI);Align_On=Align_On(SI);
   CX=CX(SI,:,:);CY=CY(SI,:,:);TLIST=TLIST(SI);
   CX(find(CX==0))=nan;CY(find(CY==0))=nan;
   
   %******************Check Align variable*******************************
   %*******************Time range of JPSTH,NBINs and BINWIDTH**************
   maxRange=abs(max(Post_Time-Pre_Time));
   if(NBINS==0 | NBINS >= maxRange)
      NBINS=ceil(maxRange);
   elseif NBINS < 10
      NBINS=ceil(0.10*maxRange);
   else
   end
   BinWidth1ms=(Post_Time-Pre_Time)/NBINS;%MUST be 1 before calling getJPSTH
   disp(['Initial Bin width for spike histogram = ',n2s(BinWidth1ms)]);
   %*******************Time range of JPSTH,NBINs and BINWIDTH**************
   %********************DIAGONALS USED FOR COINCIDENCE*********************
   AboveMain=round(XFiresBeforeY/BINSIZE);
   BelowMain=round(XFiresAfterY/BINSIZE);
   %********************DIAGONALS USED FOR COINCIDENCE***********************
   %Make bins for binning spike times between Pre and Post time 
   %BinC=(Pre_Time+(BinWidth/2):BinWidth:Post_Time-(BinWidth/2));
   BinC1ms=(Pre_Time:BinWidth1ms:Post_Time);
   x_axis=[Pre_Time:BINSIZE:Post_Time]';;%for PSTH will be in user specified binsize
   HISTAxis=x_axis; 
   %make a matrix of Ntrials when spikes are dropped off after Saccade event
   %temp=[]; NTrialsMatrix=[];tempH=[];
   %temp=EndEvent-Align_On;
   %*****************************************************************************
   %for rasters
   MaxRasters=40;
   NoOfRasters=ceil(NTrials/MaxRasters)*MaxRasters;
   %*************************************************************************
   CellX=[];CellY=[];NTRIALS=[];
   CellX=CX(:,:,1);
   CellY=CY(:,:,1);
   NTRIALS=NTrials;
   CellXHist=[];CellYHist=[];
   CellXRast1ms=[];CellYRast1ms=[];
   CellXRast1ms=(hist(CellX',BinC1ms))';
   CellYRast1ms=(hist(CellY',BinC1ms))';
   if BinWidth1ms==1 % will enter this loop 
      CellXRast1ms(find(CellXRast1ms>1))=1;
      CellYRast1ms(find(CellYRast1ms>1))=1;
   end
   %count spikes
   XSpkCount=0;YSpkCount=0;
   XSpkCount=sum(CellXRast1ms(:));
   YSpkCount=sum(CellYRast1ms(:));
   CellXHist=(hist(CellX(:),HISTAxis))./NTRIALS;
   CellYHist=(hist(CellY(:),HISTAxis))./NTRIALS;
	% different formats of calling JPSTH   
   %[RAW,N_JPST]=GetJPSTH(CellXRast1ms,CellYRast1ms,BINSIZE);
   %[RAW,N_JPST,P_JPST,E_JPST,I_JPST]=GetJPSTH(CellXRast1ms,CellYRast1ms,BINSIZE);
   %[RAW,N_JPST,P_JPST,E_JPST,I_JPST]=GetJPSTH(CellXRast1ms,CellYRast1ms,BINSIZE);
   [RAW,N_JPST]=GetJPSTH(CellXRast1ms,CellYRast1ms,BINSIZE);
   %NOW CHANGE THE BINWIDTH from 1ms to the user specified
   BinWidth=BINSIZE;
   %find the XCoincidence **************************************************************
   XCOIN=[];rXCOIN=[];XCORREL=[];rXCORREL=[];
   [XCOIN,XCORREL]=GetCoin(N_JPST,XFiresBeforeY,XFiresAfterY,BinWidth);
   [rXCOIN,rXCORREL]=GetCoin(RAW,XFiresBeforeY,XFiresAfterY,BinWidth);
   XCOINS=[];rXCOINS=[];XCORRELS=[];rXCORRELS=[];temp=[];
   XCOINS=XCOIN;rXCOINS=rXCOIN;
   XCORRELS=XCORREL;rXCORRELS=rXCORREL;
   if(smoothing)
      temp=gsmooth(XCOINS(:,2),'gaussian');
      XCOINS=[[XCOINS(:,1)][temp(:)]];
      temp=[];
      temp=gsmooth(rXCOINS(:,2),'gaussian');
      rXCOINS=[[rXCOINS(:,1)][temp(:)]];
      temp=[];
      temp=gsmooth(XCORRELS(:,2),'5point');
      XCORRELS=[[XCORRELS(:,1)][temp(:)]];
      temp=[];
      temp=gsmooth(rXCORRELS(:,2),'5point');
      rXCORRELS=[[rXCORRELS(:,1)][temp(:)]];
   end
   
   
   %OUTPUTS
   OUTSTRUCT.RAW{DIRECTION}=RAW;
   OUTSTRUCT.N_JPST{DIRECTION}=N_JPST;
   if exist('E_JPST'), OUTSTRUCT.E_JPST{DIRECTION}=E_JPST; end
   if exist('I_JPST'), OUTSTRUCT.I_JPST{DIRECTION}=I_JPST; end
   OUTSTRUCT.XCOINS{DIRECTION}=XCOINS;
   OUTSTRUCT.CellXRast{DIRECTION}=CellXRast1ms;
   OUTSTRUCT.CellYRast{DIRECTION}=CellYRast1ms;
   OUTSTRUCT.CellXHist{DIRECTION}=CellXHist;
   OUTSTRUCT.CellYHist{DIRECTION}=CellYHist;
   OUTSTRUCT.HISTAxis{DIRECTION}=HISTAxis;
   OUTSTRUCT.CellX{DIRECTION}=CellX;
   OUTSTRUCT.CellY{DIRECTION}=CellY;
   
   %%%DRAWING ROUTINE**********************************************************************
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Abscissa
   %x_lim=[min(x_axis)-(BinWidth/2) max(x_axis)+(BinWidth/2)];
   x_axis=HISTAxis;
   x_lim=[min(HISTAxis) max([max(HISTAxis) 1])];
   % always make the x-axis into 4 parts
   x_tick=[];step=[];
   step=(range(x_lim))/4; x_tick=([0:4]'.*step)+(min(x_lim));
   maxy=max(max(CellXHist),max(CellYHist));
   maxy=maxy*1.10;
   miny=min(min(CellXHist),min(CellYHist));
   if miny==0, miny=-0.01*maxy; else, miny=miny-(0.05*miny);end
   if maxy, maxy=maxy*1.05; else, maxy=1; end
   y_lim=[miny maxy];
   if ManualScale
      miny=min(ScaleHIST);
      maxy=max(ScaleHIST);
   end
   
   axes(H_JPSTX(JJ))
   set(bar(x_axis,CellXHist,1),'facecolor',[0.0 0.75 0.25],'edgecolor','none')
   %SWAPAXES=0;
   hold on
   %**************************************************************************************
   stepsize=(0.95*maxy)/NoOfRasters;
   SWAPAXES=0;
   DrawRasters(CellX,Align_On,EndEvent,stepsize,SWAPAXES,DRAWRAST,DIRECTION);
   %**************************************************************************************
   set(gca,'xdir','reverse','cameraupvector',[0 -1  0],'xtick',x_tick,'xticklabel',[],...
      'xaxislocation','top','yticklabel',[]);
   %NEW FOR LOOP*********************************************************************************
   for LL=1:length(x_tick)
      set(line([x_tick(LL) x_tick(LL)],y_lim),'linestyle','--','color',[0 0 0])
   end   
   %set(line([0 0],y_lim),'linestyle','-','color',[0 1 1])
   if(0>=min(x_lim) & 0 <=max(x_lim))
      text(0,max(y_lim)*1.08,'*','fontweight','bold','fontsize',8,'horizontalalignment','center')
   end
   %END NEW FOR LOOP******************************************************************************
   
   set(gca,'xlim',x_lim,'ylim',y_lim,'layer','top');
   text(min(x_lim)-(0.05*range(x_lim)), max(y_lim)*0.95, ['Cell ',XCellNo],...
      'fontname','Arial','fontweight','bold','fontangle','italic','fontsize',fontscale*6,...
      'Horizontalalignment','right','rotation',0);
   %text min time
   text(min(x_lim)-(0.025*range(x_lim)),max(y_lim)*0.5,[num2str(Pre_Time),' ms'],...
      'fontname','Arial','fontsize',fontscale*5,'fontweight','bold',...
      'Horizontalalignment','center','verticalalignment','middle','rotation',90);
   %text max time
   text(max(x_lim)+(0.025*range(x_lim)),max(y_lim)*0.5,[num2str(Post_Time),' ms'],...
      'fontname','Arial','fontsize',fontscale*5,'fontweight','bold',...
      'Horizontalalignment','center','verticalalignment','middle','rotation',90);
   hold off
   %Ordinate
   axes(H_JPSTY(JJ))
   set(bar(x_axis,CellYHist,1),'facecolor',[1.0 0.5 0.5],'edgecolor','none')
   SWAPAXES=1;
   hold on
   DrawRasters(CellY,Align_On,EndEvent,stepsize,SWAPAXES,DRAWRAST,DIRECTION);
   set(gca,'xdir','normal','cameraupvector',[1 0 0],'xticklabel',[],'yticklabel',[],'xtick',x_tick);
   %NEW FOR LOOP*********************************************************************************
   for LL=1:length(x_tick)
      set(line([x_tick(LL) x_tick(LL)],y_lim),'linestyle','--','color',[0 0 0])
   end   
   %set(line([0 0],y_lim),'linestyle','-','color',[0 1 1])
   if(0>=min(x_lim) & 0 <=max(x_lim))
      text(0,max(y_lim)*1.08,'*','fontweight','bold','fontsize',8,'horizontalalignment','center')   
   end
   %END NEW FOR LOOP******************************************************************************
   set(gca,'xlim',x_lim,'ylim',y_lim,'layer','top');
   text(min(x_lim)-(0.09*range(x_lim)), max(y_lim)*1.0, ['Cell ',YCellNo],... 
      'fontname','Arial','fontweight','bold','fontangle','italic','fontsize',fontscale*6,...
      'Horizontalalignment','left','rotation',0);
   %text min time
   text(min(x_lim)-(0.04*range(x_lim)),max(y_lim)*0.5,[num2str(Pre_Time),' ms'],...
      'fontname','Arial','fontsize',fontscale*5,'fontweight','bold',...
      'Horizontalalignment','center','verticalalignment','middle','rotation',0);
   %text max time
   text(max(x_lim)+(0.04*range(x_lim)),max(y_lim)*0.5,[num2str(Post_Time),' ms'],...
      'fontname','Arial','fontsize',fontscale*5,'fontweight','bold',...
      'Horizontalalignment','center','verticalalignment','middle','rotation',0);
   %Text for ALIGNMENT EVENT
   text(min(x_lim)-(0.23*range(x_lim)), max(y_lim)*0.65, [{'Time from'; upper(Align)}],...
      'fontname','Arial','fontweight','bold','fontsize',fontscale*8,...
      'Horizontalalignment','center','rotation',0);
   hold off
   %other plots
   %Cross coincidence
   axes(H_COIN(JJ))
   temp=0;
   if(~SMOOTHPLOT),temp=rXCOIN;else,temp=XCOINS;,end
   x_axis_CCH=temp(:,1);
   temp=temp(:,2);
   %MANUALLY SCALE the Y-AXIS for XCOINS
   if ManualScale
      y_limc=ScaleXCOINS;
   else
      y_limc=[-max(abs(temp))*1.0 max(abs(temp))*1.0+0.0001];
   end
   
   set(bar(x_axis_CCH,temp,1),'facecolor',[1.0 0.0 1.0],'edgecolor','none')
   x_lim_CCH=[min(x_axis_CCH) max(x_axis_CCH)];
   set(gca,'xlim',x_lim_CCH)
   set(gca,'ylim',y_limc)
   
   set(gca,'xticklabel',[],'yticklabel',[])
   set(line(x_lim_CCH,[0 0]),'color',[0 0 0],'linewidth',linescale*1.2);
   set(line([max(x_lim_CCH) max(x_lim_CCH)],[0 max(y_limc)]),'color',[0 0 0],'linewidth',linescale*1.0);
   set(line([max(x_lim_CCH) max(x_lim_CCH)-(4*BinWidth)],[max(y_limc) max(y_limc)]),'color',[0 0 0],'linewidth',linescale*1.0);
   
   text(max(x_lim_CCH)*0.96,max(y_limc)*0.96,num2str(max(y_limc),'%2.3f '),...
      'horizontalalignment','right','verticalalignment','bottom',...
      'fontname','Arial','fontsize',fontscale*6,'fontweight','normal','rotation',45);
   pbaspect([3 1 1])
   camzoom(sqrt(2));
   camorbit(-45,0)
   set(gca,'layer','bottom')
   %Cell X Hist second time
   axes(H_JPSTX2(JJ))
   set(bar(x_axis,CellXHist,1),'facecolor',[0.0 0.75 0.25],'edgecolor','none')
   SWAPAXES=0;
   hold on
   DrawRasters(CellX,Align_On,EndEvent,stepsize,SWAPAXES,DRAWRAST,DIRECTION);
   set(gca,'xdir','reverse','cameraupvector',[0 -1  0],'xtick',x_tick,'xticklabel',[],...
      'xaxislocation','top','yticklabel',[]);
   %NEW FOR LOOP*********************************************************************************
   for LL=1:length(x_tick)
      set(line([x_tick(LL) x_tick(LL)],y_lim),'linestyle','--','color',[0 0 0])
   end   
   %set(line([0 0],y_lim),'linestyle','-','color',[0 1 1])
   if(0>=min(x_lim) & 0 <=max(x_lim))
      text(0,max(y_lim)*1.08,'*','fontweight','bold','fontsize',8,'horizontalalignment','center')   
   end
   %END NEW FOR LOOP******************************************************************************
   set(gca,'xdir','reverse','cameraupvector',[0 -1  0]);
   set(gca,'xaxislocation','top','xlim',x_lim,'ylim',y_lim,'yticklabel',[],'xticklabel',[],'layer','top')
   hold off   
   
   %DRAW  CROSS CORRELATION
   axes(H_CORREL(JJ))
   %x_axisCorrel=[-(Post_Time-Pre_Time)/2:BinWidth:(Post_Time-Pre_Time)/2];
   x_axisCorrel=XCORRELS(:,1);
   
   %always have a total of 9 ticks for x axis (including 0)
   xc_tick=[];step=[];
   step=(range(x_axisCorrel))/8;xc_tick=([-4:1:4])'.*step;
   temp=0;tmp=0;
   if(~SMOOTHPLOT),tmp=rXCORREL;else,tmp=XCORRELS;,end
   temp=tmp(:,2);
   minmax=[min(temp) max(temp)];
   y_axis=[];
   y_axis=minmax;
   if max(temp)<0 ,y_axis(2)=0;end
   if min(temp)>0 ,y_axis(1)=0;end
   yscale=0;
   if abs(min(temp))> abs(max(temp))
      yscale=min(temp);
   else
      yscale=max(temp);
   end
   
   %set(bar(x_axisCorrel+BinWidth,temp,1),'facecolor',[1 0.25 1],'edgecolor','none')
   set(bar(tmp(:,1),tmp(:,2),1),'facecolor',[1 0.25 1],'edgecolor','none')
   
   set(gca,'yticklabel',[],'xticklabel',[],'xtick',xc_tick)
   set(gca,'xlim',[-range(x_axisCorrel) range(x_axisCorrel)],...
      'ylim',y_axis);
   set(line([-range(x_axisCorrel) range(x_axisCorrel)],[0 0]),'color',[0 0 0],'linewidth',linescale*1.0);
   set(line([-range(x_axisCorrel) -range(x_axisCorrel)],[0 yscale]),'color',[0 0 0],'linewidth',linescale*1.0);
   set(line([-range(x_axisCorrel) -range(x_axisCorrel)+4*BinWidth],[yscale yscale]),'color',[0 0 0],'linewidth',linescale*1.0);
   text(-range(x_axis)*0.9,yscale,num2str(yscale,'%2.3f '),...
      'horizontalalignment','left','verticalalignment','middle',...
      'fontname','Arial','fontsize',fontscale*6,'fontweight','normal','rotation',-45);
   set(line([(AboveMain+0.5)*BinWidth (BelowMain+0.5)*BinWidth], [yscale*1.05 yscale*1.05] ),'color',[0 0 0],'linewidth',linescale*1.0)
   set(line([(AboveMain+0.5)*BinWidth (AboveMain+0.5)*BinWidth], [yscale*1.05 yscale*0.99] ),'color',[0 0 0],'linewidth',linescale*1.0)
   set(line([(BelowMain+0.5)*BinWidth (BelowMain+0.5)*BinWidth], [yscale*1.05 yscale*0.99] ),'color',[0 0 0],'linewidth',linescale*1.0)
   %set(line([0 0], [max(temp)*1.07 max(temp)*1.12] ),'color',[0 0 0],'linewidth',linescale*1.0)
   set(line([0 0], [0 yscale*1.15]),'color',[0 0 0],'linewidth',linescale*1.0)
   pbaspect([5 1 1]);
   camzoom(sqrt(2))
   camorbit(45,0)
   
   %DRAW JPST IMAGES
   axes(H_JPST(JJ))
   temp=[];infotxt=[];
   temp=N_JPST;
   infotxt='ALL Spikes';
   if ~ManualScale
      ScaleJPST=[-max(temp(:)) max(temp(:))];
   end
   imagesc(temp,ScaleJPST);   
   clear temp
   %draw the intensity bar
   cp=get(gca,'position');
   cp(2)=cp(4)+cp(2)+0.10;
   cp(4)=0.2*cp(3);
   %scalefactor
   scalefactor=(cp(4)/cp(3))/0.05;
   set(gca,'xticklabel',[],'yticklabel',[],'Xtick',[],'Ytick',[]);
   h=axes('position',cp);
   h1=colorbar(h);
   set(h1,'xticklabel',[],'yticklabel',[])
   %text(0,max(get(h1,'ylim')),n2s(min(temp(:))),...
   text(0,max(get(h1,'ylim')),n2s(min(ScaleJPST)),...
      'fontname','Arial','fontsize',fontscale*6,'rotation',0,...
      'horizontalalignment','left','verticalalignment','bottom');
   text(max(get(h1,'xlim')),max(get(h1,'ylim')),n2s(max(ScaleJPST)),...
      'fontname','Arial','fontsize',fontscale*6,'rotation',0,...
      'horizontalalignment','right','verticalalignment','bottom');
   x=max(get(h1,'xlim'))/2;
   y=max(get(h1,'ylim'))*(10/scalefactor);
   text(x,y*1.5,infotxt,'horizontalalignment','center',...
      'fontname','Arial','fontweight','bold','fontsize',fontscale*11);
   minx=min(get(gca,'xlim'));
   off=(range(get(gca,'xlim')))/4;
   xpos=[minx+(0.6*off) minx+(2*off) minx+(3.2*off)];
   %ypos=[-24 -18 -12 -6];
   ypos=[-24 -20 -16 -12 -8];
   ypos=ypos./scalefactor;
   
   text(xpos(2),ypos(5),['Cell ',XCellNo],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   text(xpos(3),ypos(5),['Cell ',YCellNo],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   
   text(xpos(1),ypos(4),['Total Trials'],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   text(xpos(2),ypos(4),[num2str(NTrials)],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   text(xpos(3),ypos(4),[num2str(NTrials)],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   
   text(xpos(1),ypos(3),['Total Spikes'],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   text(xpos(2),ypos(3),[num2str(XSpkCount)],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   text(xpos(3),ypos(3),[num2str(YSpkCount)],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   
   text(xpos(1),ypos(2),['Bin Width (ms)'],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   text(xpos(2),ypos(2),[n2s2(BinWidth)],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   text(xpos(3),ypos(2),[n2s2(BinWidth)],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   
   text(xpos(1),ypos(1),['Coincidence (ms)'],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   text(xpos(3),ypos(1),[n2s2(AboveMain*BinWidth),' to ',n2s2(BelowMain*BinWidth)],'fontsize',6,'fontweight','normal',...
      'horizontalalignment','center')
   % draw the grid for all 3 sub JPSTHs
   gridlines=1;
   if gridlines
      axes(H_JPSTX2(JJ));
      for jj=1:4
         text(x_tick(jj+1),0,repmat('-',1,round(15*jj)),'color',[0.5 0.5 0.5],'rotation',90,'fontname','Arial');
      end
      axes(H_JPSTY(JJ));
      m=max(get(gca,'ylim'));
      for jj=1:4
         text(x_tick(jj+1),-(m*2.85),repmat('-',1,round(15*jj)+1),'color',[0.5 0.5 0.5],'rotation',0,'fontname','Arial');
      end   
   end
   %FIND max scale for scaling the Y axis of 
   COINS_MaxY=max(COINS_MaxY,max(abs(XCOINS(:,2))));
end %for DIRECTION=1:2

%*********************************************************************
axes(H_INFO)
XCellI=[];YCellI=[];FX=[];FY=[];F=[];

if isstruct(XCellInfo)
   FX=fieldnames(XCellInfo);
   for jj=1:size(FX,1)
      tmp=num2str(eval(['XCellInfo.'FX{jj}]));
      XCellI(jj,1:length(tmp))=tmp; tmp=[];
   end
else
   XCellI=XCellNo;
end

if isstruct(YCellInfo)
   FY=fieldnames(YCellInfo);
   for jj=1:size(FY,1)
      tmp=num2str(eval(['YCellInfo.'FY{jj}]));
      YCellI(jj,1:length(tmp))=tmp; tmp=[];
   end
else
   YCellI=YCellNo;
end

if ~isempty(FX)
   F=FX;
elseif ~isempty(FY)
   F=FY;
else
   F='Input Cells';
end
xpos=[0.6 0.75 0.9];
ypos=0.92;%was 0.95
text(0.75,0.95,[num2str(XCellNo),' vs ',num2str(YCellNo)],'fontsize',12,'fontweight','bold','horizontalalignment','center','verticalalignment','middle')
text(xpos(1),ypos,F,'fontsize',8,'fontweight','normal','horizontalalignment','center','verticalalignment','top');
text(xpos(2),ypos,char(XCellI),'fontsize',8,'fontweight','normal','horizontalalignment','center','verticalalignment','top');
text(xpos(3),ypos,char(YCellI),'fontsize',8,'fontweight','normal','horizontalalignment','center','verticalalignment','top');
text(0.6,0.12,['RF Used'],'fontsize',8,'fontweight','bold','horizontalalignment','center','verticalalignment','top');
text((0.75+0.9)/2,0.12,['['num2str(RFPOS(:)'),']'],'fontsize',8,'fontweight','bold','horizontalalignment','center','verticalalignment','top');
text(0.6,0.08,['Ref. Cell'],'fontsize',8,'fontweight','bold','horizontalalignment','center','verticalalignment','top');
text((0.75+0.9)/2,0.08,['[',YCellNo,']'],'fontsize',8,'fontweight','bold','horizontalalignment','center','verticalalignment','top');
%FLAG for callin X-CORR
if CallXCorr==1
   %*************Call X-CORRELATION ******************************************
   %function [ZCORR,CCH,SCCH,SIGNIF]=spike_corr(RefCell,TargetCell,alignT,FromT,ToT,CG_time,CG_bin,varargin)
   %when calling from eggermont1:
   %RefCell=CellY, TargetCell=CellX
   %alignT=0, FromT=Pre_Time, ToT=Post_Time, CG_time=[Post_Time-Pre_Time]/2,  CG_bin=BinWidth or 1(better)
   %varargin is to pass list of trials; else TrlList will be [1:no of rows in RefCell]
   CrossCorr=spike_corr(CellY,CellX,0,Pre_Time,Post_Time,[Post_Time-Pre_Time]/2,BinWidth/BinWidth);
   H_CrossCorr=axes('parent',H_Figure,'position',[0.57 0.07 0.16 0.36],'box','on');
   axes(H_CrossCorr)
   %fill(CrossCorr(:,1),CrossCorr(:,2),'r')
   temp=gsmooth(CrossCorr(:,2),'5point');
   CrossCorr(1:length(temp),2)=temp;temp=[];
   temp=gsmooth(CrossCorr(:,3),'5point');
   CrossCorr(1:length(temp),3)=temp;temp=[];
   plot(CrossCorr(:,1),CrossCorr(:,2),'r')
   ylabel('Z-Score','fontsize',8,'fontweight','bold')
   xlabel('Time delay (ms)','fontsize',8,'fontweight','bold')
   hold on
   plot(CrossCorr(:,1),CrossCorr(:,3),'b')
   h(1)=line([min(CrossCorr(:,1)) max(CrossCorr(:,1))],[3 3]);
   h(2)=line([min(CrossCorr(:,1)) max(CrossCorr(:,1))],[-3 -3]);
   set(h,'linestyle','--');
   title('Crosscorrelation','fontsize',10,'fontweight','bold')
   hold off
end

axes(H_INFO)
text(0.4,0.02,FigFile,'interpreter','none','fontsize',4,'fontweight','normal')
text(0.85,0.02,[datestr(datenum(now),1),', ',datestr(datenum(now),13)],'fontsize',4,'fontweight','normal');
if CallXCorr==1, axes(H_CrossCorr), end



%*************SUB-FUNCTIONS******************************************
%GetCoin -------------> moved to a function getcoin
%*************SUB-FUNCTIONS******************************************
%***********************SMOOTHING FUNCTION****************************
function out=gsmooth(in,KERNEL)
%smooth psthistogram with a gaussian of sd SIGMA bins
%generate kernel for smoothing
switch KERNEL
case 'gaussian'
   Sigma=1.5;
   disp(['Convolving with GAUSSIAN sigma(bins) = ',num2str(Sigma)])
   Kernel=[-3*Sigma:3*Sigma];
   GBinSize=length(Kernel);
   Half_BW=(GBinSize-1)/2;
   Kernel=[-GBinSize/2:GBinSize/2];
   Factor=1.00/sqrt(Sigma*2*pi);
   Kernel=Factor*(exp(-(0.5*((Kernel./Sigma).^2))));
   Kernel=Kernel(:);
case '5point'
   %disp(['5 point smoothing'])
   Kernel=[0.1;0.2;0.5;0.2;0.1];
   in=in(:);
otherwise
end
%Now smooth the histograms
in=in(:);
out=convn(in,Kernel,'same');
%***********************GET TOTAL SPIKE COUNT****************************
%***************************GET EVENT TIME************************
function [RETURN_]=GetEventTime(InEvent,NTrls)
global n2s2
[r,c,p]=size(InEvent);
if p>1, RETURN_=['Matrix cannot be 3-dimensional '], end
if (r==1 & c==1) %Scalar Align on time
   RETURN_=repmat(InEvent,NTrls,1);
elseif (r==1 & c>1) | (r>1 & c==1) %Row vector or a column vector
   RETURN_=InEvent(:);
elseif (r>1 & c==2) % Event variable from file
   RETURN_=InEvent(:,1);
else
   RETURN_=['CHECK matrix .....'];
end
clear r c p
%********************************************************
function DrawRasters(SPKMatrix,AlignTimes,EndTimes,stepsize,SWAP,DR,DIRECT)
SM=SPKMatrix;AT=AlignTimes;ET=EndTimes;
NT=size(SM,1);
%set times LT x_axis to nan and GT x-axis to nan
XLIM=get(gca,'xlim');
AM=0;AM=AT-AT;EM=ET-AT;
if DIRECT==1
	EM(find(EM>max(XLIM)))=nan;   
else
	EM(find(EM<min(XLIM)))=nan;   
end
SO=stepsize*(NT*0.05);

M='I';
FS='fontsize';RFS=2;MFS=2;
FW='fontweight';RFW='Normal';MFW='bold';
if SWAP
   ROT='Rotation';RROT=90;MROT=0;
else
   ROT='Rotation';RROT=0;MROT=90;
end
FC='color';RFC=[0.6 0.6 0.6];AMFC='R';EMFC='b';
VA='VerticalAlignment';va='Middle';
HA='HorizontalAlignment';ha='Center';
if DR==1
   for KK=1:NT
      R=0;YP=0;
      R=SM(KK,:);YP=repmat(stepsize*KK,1,length(R));
      text(R,YP+SO,M,FS,RFS,FW,RFW,FC,RFC,ROT,RROT,VA,va,HA,ha);
   end
end
YPOS=[1:NT]';
YPOS=(YPOS.*stepsize)+SO;
text(AM,YPOS,M,FS,MFS,FW,MFW,FC,AMFC,ROT,MROT,VA,va,HA,ha);
text(EM,YPOS,M,FS,MFS,FW,MFW,ROT,MROT,FC,EMFC,VA,va,HA,ha);

%*************************************************
function setcolormap(COLOR)
if isempty(COLOR)
   COLOR=0;
end
if(COLOR==0)
   colormap('gray');
   gr=colormap;
   for LL=1:64;invgr(LL,:)=gr(65-LL,:);end
   colormap(invgr);
else
   CC=getcolors;
   colormap(CC./400)
end
%***************************************
function CC=getcolors
%Color Table
CC=[0	0	225
   0	0	250
   0	0	275
   0	0	300
   0	0	325
   0	0	350
   0	0	375
   0	0	400
   0	25	400
   0	50	400
   0	75	400
   0	100	400
   0	125	400
   0	150	400
   0	175	400
   0	200	400
   0	225	400
   0	250	400
   0	275	400
   0	300	400
   0	325	400
   0	350	400
   0	375	400
   0	400	400
   25	400	400
   50	400	375
   75	400	350
   100	400	325
   125	400	300
   150	400	275
   175	400	250
   200	400	225
   225	400	200
   250	400	175
   275	400	150
   300	400	125
   325	400	100
   350	400	75
   375	400	50
   400	400	25
   400	400	0
   400	375	0
   400	350	0
   400	325	0
   400	300	0
   400	275	0
   400	250	0
   400	225	0
   400	200	0
   400	175	0
   400	150	0
   400	125	0
   400	100	0
   400	75	0
   400	50	0
   400	25	0
   400	0	0
   375	0	0
   350	0	0
   325	0	0
   300	0	0
   275	0	0
   250	0	0
   225	0	0
];