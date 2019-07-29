%create figre template for JPSTH
delete(allchild(0))
global H_Figure H_JPST H_JPSTX H_JPSTY H_JPSTX2 H_COIN H_CORREL H_INFO
global sizescale
if isempty(sizescale)
   sizescale=1.0;
end

set(0,'units','pixels')
set(0,'defaulttextfontsize',6*sizescale,...
   'defaulttextfontname','Arial',...
   'defaultaxesfontsize',6*sizescale,...
   'defaultaxeslinewidth',0.05)
screen_size=get(0,'ScreenSize');
w=screen_size(3);h=screen_size(4);margin=40;
scale=sizescale;
w=w*scale;h=h*scale;

FigPos=[margin margin w-(2*margin) h-(2*margin)];
Name=(['   JSTH   ']);
%Main figure window
H_Figure=figure('Position',FigPos,...
   'color',[1 1 1],'numbertitle','off','renderer','painters',...
   'renderermode','manual','menubar','none','name',Name,...
   'Tag','H_Figure');
orient landscape
set(H_Figure,'units','normalized')

OutlinePos=[0.002 0.002 0.996 0.996];
H_Outline=axes('parent',H_Figure,'position',OutlinePos,'box','on',...
   'xtick',[],'ytick',[],'xcolor',[0 0 0],'ycolor',[0 0 0],'Tag','H_Outline');

%Scaling
yaspect=w/h;offset=0.01;
xoffset=0.025;
for ii=1:3
   switch ii
   case 1
      FigPos1=[0.1-xoffset 0.6 0.15 0.15*yaspect];
   case 2
      FigPos1=[0.6-xoffset 0.6 0.15 0.15*yaspect];
   case 3
      FigPos1=[0.1-xoffset 0.1 0.15 0.15*yaspect];
   end
   
   %JPSH figure
   H_JPST(ii)=axes('parent',H_Figure,'position',FigPos1,'box','on',...
      'layer','top','Tag','H_JPSTH');
   set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
end

%Cell On Abscissa (XCell)
%FigPos=[0.2 0.04 0.25 0.115*yaspect];
for ii=1:3
   switch ii
   case 1
      FigPos=[0.1-xoffset 0.52 0.15 0.055*yaspect];
   case 2
      FigPos=[0.6-xoffset 0.52 0.15 0.055*yaspect];
   case 3
      FigPos=[0.1-xoffset 0.02 0.15 0.055*yaspect];
   end
   H_JPSTX(ii)=axes('parent',H_Figure,'position',FigPos,'box','on',...
      'layer','top','Tag','H_JPSTX');
   set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
end

%Cell On Ordinate
%FigPos=[0.08 0.2 0.115 0.25*yaspect];
for ii=1:3
   switch ii
   case 1
      FigPos=[0.1-((2/300)+0.055)-xoffset 0.6 0.055 0.15*yaspect];
   case 2
      FigPos=[0.6-((2/300)+0.055)-xoffset 0.6 0.055 0.15*yaspect];
   case 3
      FigPos=[0.1-((2/300)+0.055)-xoffset 0.1 0.055 0.15*yaspect];
   end
   H_JPSTY(ii)=axes('parent',H_Figure,'position',FigPos,'box','on',...
      'layer','top','Tag','H_JPSTY');
   set(gca,'cameraupvector',[1 0 0])
   set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
end

%Cell On Abscissa(repeat)
for ii=1:3
   switch ii
   case 1
      FigPos=[0.1+0.15+0.01-xoffset 0.52 0.15 0.055*yaspect];
   case 2
      FigPos=[0.6+0.15+0.01-xoffset 0.52 0.15 0.055*yaspect];
   case 3
      FigPos=[0.1+0.15+0.01-xoffset 0.02 0.15 0.055*yaspect];
   end
   H_JPSTX2(ii)=axes('parent',H_Figure,'position',FigPos,'box','on',...
      'layer','top','Tag','H_JPSTX');
   set(gca,'xtick',[],'ytick',[],'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5]);
   
end



%DIAG figure
%FigPoscoin=[0.5 0.2 0.25 0.25*yaspect];
for ii=1:3
   switch ii
   case 1
      FigPoscoin=[0.1+0.15+0.01-xoffset 0.6 0.15 0.15*yaspect];
   case 2
      FigPoscoin=[0.6+0.15+0.01-xoffset 0.6 0.15 0.15*yaspect];
   case 3
      FigPoscoin=[0.1+0.15+0.01-xoffset 0.1 0.15 0.15*yaspect];
   end
   H_COIN(ii)=axes('parent',H_Figure,'position',FigPoscoin,'box','off',...
      'xtick',[],'ytick',[],'yaxislocation','right','xcolor',[1 1 1],...
      'ycolor',[1 1 1],'layer','bottom','Tag','H_COIN');
   set(gca,'nextplot','add');
   %ylim([-1 +1]);
end
%CORRELOGRAM fig
%FigPoscorr=[0.5+(0.25*0.67) 0.2+(0.25*yaspect*0.67) 0.25 0.25*yaspect];
for ii=1:3
   switch ii
   case 1
      FigPoscorr=[0.315+(0.15*0.275)-xoffset 0.6+(0.15*yaspect*0.64) 0.15 0.15*yaspect];
   case 2
      FigPoscorr=[0.815+(0.15*0.275)-xoffset 0.6+(0.15*yaspect*0.64) 0.15 0.15*yaspect];
   case 3
      FigPoscorr=[0.315+(0.15*0.275)-xoffset 0.1+(0.15*yaspect*0.64) 0.15 0.15*yaspect];
   end
   H_CORREL(ii)=axes('parent',H_Figure,'position',FigPoscorr,'box','off',...
      'xtick',[],'ytick',[],'yaxislocation','left','xcolor',[1 1 1],...
      'ycolor',[1 1 1],'layer','bottom','Tag','H_CORREL');
   set(gca,'nextplot','add');
%   ylim([-1 +1]);
end
InfoPos=[0.52 0.02 0.46 0.46];
H_INFO=axes('parent',H_Figure,'position',InfoPos,'box','on',...
   'xtick',[],'ytick',[]);
set(gca,'linewidth',0.01);
