
function mpr(data,OV,varargin)
% MPR(volume,overlay_volume,'tag1',value1,...) - Creates MPR view of the given data set
%	volume  = 3D matrix
%	overlay
%
% Navigation:
%   Left mouse button - navigation
%   Scroll button - slice scrolling
%   Right button - window and level
%   Shift-Right (axial-only) - Create an oblique view
%   Spacebar - Turn slice indicators on/off
%
% Tags/Values			Default			Function
% ===============================================================================
% 'clim',[num num]		[min max]		Set colorbar limits
% 'slices',[x y z]		[nx/2 ny/2 nz/2]	Set slices to display
% 'oblique',[bx,by,ex,ey]	[]			Draw Oblique cut in LR
% 'thickness',num		1			Thickness for slices/oblique

% Web Stayman, August 2010 - April 2012

% USER OVERRIDES
for j=1:length(varargin)/2,
	if iscell(varargin{j*2}),
		eval(['clear ',varargin{j*2-1},';']);
		for k=1:length(varargin{j*2}),
% 			disp([varargin{j*2-1} '{',num2str(k),'} = ''' varargin{j*2}{k} ''';']);
			eval([varargin{j*2-1} '{',num2str(k),'} = ''' varargin{j*2}{k} ''';']);
		end
	else

% 		disp([varargin{j*2-1} ' = ' mat2str(varargin{j*2}) ';']);
		eval([varargin{j*2-1} ' = ' mat2str(varargin{j*2}) ';']);
	end
end

[nx,ny,nz] = size(data); m = max([nx ny nz]);
cx = ceil(nx/2); cy = ceil(ny/2); cz = ceil(nz/2); cm = ceil(m/2);

if ~exist('OV','var'), OV={}; end
if ~exist('clim','var'), clim = [min(data(:)), max(data(:))]; end
if ~exist('slices','var'), slices(1)=cx; slices(2)=cy; slices(3)=cz; end
if ~exist('oblique','var'), oblique=[]; end
if ~exist('thickness','var'), thickness=1; end

zoom off;
set(gcf,'UserData',[]);

[nx,ny,nz] = size(data); m = max([nx ny nz]);
cx = ceil(nx/2); cy = ceil(ny/2); cz = ceil(nz/2); cm = ceil(m/2);

cx = slices(1); cy = slices(2); cz = slices(3);

MPRdata.nx = nx; MPRdata.ny = ny; MPRdata.nz = nz;
MPRdata.cx = cx; MPRdata.cy = cy; MPRdata.cz = cz;
MPRdata.data = data;
MPRdata.OV = OV;
MPRdata.handlebuttondown = 0;
MPRdata.thickness = thickness;

cols = [1 0 0; 0 0 1; 0 1 0; 1 1 0; 1 0 1; 0 1 1;];
for j=1:length(MPRdata.OV), 
	ind = mod(j-1,length(MPRdata.OV))+1;
	MPRdata.OV1{j} = cat(3,cols(ind,1)*ones(nx,ny),cols(ind,2)*ones(nx,ny),cols(ind,3)*ones(nx,ny));
	MPRdata.OV2{j} = cat(3,cols(ind,1)*ones(nz,nx),cols(ind,2)*ones(nz,nx),cols(ind,3)*ones(nz,nx));
	MPRdata.OV3{j} = cat(3,cols(ind,1)*ones(nz,ny),cols(ind,2)*ones(nz,ny),cols(ind,3)*ones(nz,ny));
end

msubplot(2,2,1,0,0,0,0,0,0); set(gca,'Clim',clim); mpr_slice(MPRdata,1); 
hold on; MPRdata.h1=plot([1 1]*MPRdata.cy,[1 MPRdata.nx],'g'); hold off;
hold on; MPRdata.h2=plot([1 MPRdata.ny],[1 1]*MPRdata.cx,'g'); hold off;
mtext(0.01*MPRdata.ny,0.03*MPRdata.nx,sprintf('%d/%d',MPRdata.cz,MPRdata.nz),10,'g');

msubplot(2,2,2,0,0,0,0,0,0); set(gca,'Clim',clim); mpr_slice(MPRdata,2); 
hold on; MPRdata.h3=plot([1 MPRdata.nx],[1 1]*MPRdata.cz,'g'); hold off;
hold on; MPRdata.h4=plot([1 1]*MPRdata.cx,[1 MPRdata.nz],'g'); hold off;
mtext(0.01*MPRdata.nx,0.03*MPRdata.nz,sprintf('%d/%d',MPRdata.cy,MPRdata.ny),10,'g');

msubplot(2,2,3,0,0,0,0,0,0); set(gca,'Clim',clim); mpr_slice(MPRdata,3);
hold on; MPRdata.h5=plot([1 1]*MPRdata.cy,[1 MPRdata.nz],'g'); hold off;
hold on; MPRdata.h6=plot([1 MPRdata.ny],[1 1]*MPRdata.cz,'g'); hold off;
mtext(0.01*MPRdata.ny,0.03*MPRdata.nz,sprintf('%d/%d',MPRdata.cx,MPRdata.nx),10,'g');

msubplot(2,2,4,0,0,0,0,0,0); imagesc(0);
%vol3d('CData',data(1:2:end,1:2:end,1:2:end)); axis image;
%alphamap([zeros(1,4),linspace(0,.14,34),0.9*ones(1,26)]); set(gca,'ZDir','reverse');
h=colorbar('East'); set(h,'XColor','g','YColor','g');
axis image; set(gca,'Clim',clim); axis off;
colormap gray(256);

MPRdata.h7=[]; 
if ~isempty(oblique),
	mpr_oblique(oblique(1),oblique(2),oblique(3),oblique(4),MPRdata,clim);
	MPRdata = get(gcf,'UserData');
else
	MPRdata.oblique = [];
end

set(gcf,'WindowScrollWheelFcn',@mpr_scroll);
set(gcf,'WindowButtonDownFcn',@mpr_buttondown);
set(gcf,'WindowButtonUpFcn',@mpr_buttonup);
set(gcf,'KeyReleaseFcn',@mpr_keyrelease);
set(gcf,'DoubleBuffer','on','UserData',MPRdata,'Color',[0 0 0]);
%pos = get(gcf,'Position'); set(gcf,'Position',[pos(1) pos(2) max(pos(3:4)) max(pos(3:4))]);

if exist('l_off','var'),
	if l_off, 
		evnt.Character=' ';
		mpr_keyrelease(0,evnt);
	end
end

function mpr_scroll(src,evnt)
	MPRdata = get(gcf,'UserData');
	xy = get(0,'PointerLocation');
	pos = get(gcf,'Position');
	clim = get(gca,'Clim');

	if strcmp(get(MPRdata.h1,'Visible'),'on'), viz=1; else viz=0; end

	if xy(1)>pos(1)+pos(3)/2,
		if xy(2)>pos(2)+pos(4)/2,	% UR
			MPRdata.cy = max(min(MPRdata.cy+evnt.VerticalScrollCount,MPRdata.ny),1);
			msubplot(2,2,2,0,0,0,0,0,0); mpr_slice(MPRdata,2); 
			mtext(0.01*MPRdata.ny,0.03*MPRdata.nz,sprintf('%d/%d',MPRdata.cy,MPRdata.ny),10,'g');
			hold on; MPRdata.h3=plot([1 MPRdata.nx],[1 1]*MPRdata.cz,'g'); hold off;
			hold on; MPRdata.h4=plot([1 1]*MPRdata.cx,[1 MPRdata.nz],'g'); hold off;
			msubplot(2,2,1,0,0,0,0,0,0); delete(MPRdata.h1);
			hold on; MPRdata.h1=plot([1 1]*MPRdata.cy,[1 MPRdata.nx],'g'); hold off;
			msubplot(2,2,3,0,0,0,0,0,0); delete(MPRdata.h5);
			hold on; MPRdata.h5=plot([1 1]*MPRdata.cy,[1 MPRdata.nz],'g'); hold off;
		else						% LR
		end
	else
		if xy(2)>pos(2)+pos(4)/2,	% UL
			MPRdata.cz = max(min(MPRdata.cz+evnt.VerticalScrollCount,MPRdata.nz),1);
			msubplot(2,2,1,0,0,0,0,0,0); mpr_slice(MPRdata,1); 
			mtext(0.01*MPRdata.ny,0.03*MPRdata.nx,sprintf('%d/%d',MPRdata.cz,MPRdata.nz),10,'g');
			hold on; MPRdata.h1=plot([1 1]*MPRdata.cy,[1 MPRdata.nx],'g'); hold off;
			hold on; MPRdata.h2=plot([1 MPRdata.ny],[1 1]*MPRdata.cx,'g'); hold off;
			if ~isempty(MPRdata.oblique),
				bx=MPRdata.oblique(1); by=MPRdata.oblique(2); ex=MPRdata.oblique(3); ey=MPRdata.oblique(4);
				hold on; MPRdata.h7=plot([by ey],[bx ex],'x-'); set(MPRdata.h7,'Color',[1 1 0]); hold off;
			end
			msubplot(2,2,2,0,0,0,0,0,0); delete(MPRdata.h3);
			hold on; MPRdata.h3=plot([1 MPRdata.nx],[1 1]*MPRdata.cz,'g'); hold off;
			msubplot(2,2,3,0,0,0,0,0,0); delete(MPRdata.h6);
			hold on; MPRdata.h6=plot([1 MPRdata.ny],[1 1]*MPRdata.cz,'g'); hold off;
		else						% LL
			MPRdata.cx = max(min(MPRdata.cx+evnt.VerticalScrollCount,MPRdata.nx),1);
			msubplot(2,2,3,0,0,0,0,0,0); mpr_slice(MPRdata,3);
			mtext(0.01*MPRdata.ny,0.03*MPRdata.nz,sprintf('%d/%d',MPRdata.cx,MPRdata.nx),10,'g');
			hold on; MPRdata.h5=plot([1 1]*MPRdata.cy,[1 MPRdata.nz],'g'); hold off;
			hold on; MPRdata.h6=plot([1 MPRdata.ny],[1 1]*MPRdata.cz,'g'); hold off;
			msubplot(2,2,1,0,0,0,0,0,0); delete(MPRdata.h2);
			hold on; MPRdata.h2=plot([1 MPRdata.ny],[1 1]*MPRdata.cx,'g'); hold off;
			msubplot(2,2,2,0,0,0,0,0,0); delete(MPRdata.h4);
			hold on; MPRdata.h4=plot([1 1]*MPRdata.cx,[1 MPRdata.nz],'g'); hold off;
		end
	end
	if viz==0, 
		set([MPRdata.h1,MPRdata.h2,MPRdata.h3,MPRdata.h4,MPRdata.h5,MPRdata.h6],'Visible','off'); 
		if ~isempty(MPRdata.h7), set(MPRdata.h7,'Visible','off'); end
	end
	

	set(gcf,'UserData',MPRdata);
end

function mpr_buttondown(src,evnt)
	MPRdata = get(gcf,'UserData');
	xy = get(0,'PointerLocation');
	pos = get(gcf,'Position');
	clim = get(gca,'Clim');

	switch get(gcf,'SelectionType'),
	case 'normal'
		pt = [(xy(1)-pos(1)) (xy(2)-pos(2))];
		set(gcf,'CurrentPoint',pt);
		pt = get(gca,'CurrentPoint');
		u = round(pt(1))+1; v = round(pt(3))-1;
		if xy(1)>pos(1)+pos(3)/2,
			if xy(2)>pos(2)+pos(4)/2,	% UR
				MPRdata.cz = max(min(v,MPRdata.nz),1);
				MPRdata.cx = max(min(u,MPRdata.nx),1);
			else						% LR
			end
		else
			if xy(2)>pos(2)+pos(4)/2,	% UL
				MPRdata.cx = max(min(v,MPRdata.nx),1);
				MPRdata.cy = max(min(u,MPRdata.ny),1);
			else						% LL
				MPRdata.cz = max(min(v,MPRdata.nz),1);
				MPRdata.cy = max(min(u,MPRdata.ny),1);
			end
		end

		if strcmp(get(MPRdata.h1,'Visible'),'on'), viz=1; else viz=0; end

		msubplot(2,2,1,0,0,0,0,0,0); mpr_slice(MPRdata,1); 
		hold on; MPRdata.h1=plot([1 1]*MPRdata.cy,[1 MPRdata.nx],'g'); hold off;
		hold on; MPRdata.h2=plot([1 MPRdata.ny],[1 1]*MPRdata.cx,'g'); hold off;
		if ~isempty(MPRdata.oblique),
			bx=MPRdata.oblique(1); by=MPRdata.oblique(2); ex=MPRdata.oblique(3); ey=MPRdata.oblique(4);
			hold on; MPRdata.h7=plot([by ey],[bx ex],'x-'); set(MPRdata.h7,'Color',[1 1 0]); hold off;
		end
		mtext(0.01*MPRdata.ny,0.03*MPRdata.nx,sprintf('%d/%d',MPRdata.cz,MPRdata.nz),10,'g');

		msubplot(2,2,2,0,0,0,0,0,0); mpr_slice(MPRdata,2); 
		hold on; MPRdata.h3=plot([1 MPRdata.nx],[1 1]*MPRdata.cz,'g'); hold off;
		hold on; MPRdata.h4=plot([1 1]*MPRdata.cx,[1 MPRdata.nz],'g'); hold off;
		mtext(0.01*MPRdata.nx,0.03*MPRdata.nz,sprintf('%d/%d',MPRdata.cy,MPRdata.ny),10,'g');

		msubplot(2,2,3,0,0,0,0,0,0); mpr_slice(MPRdata,3); 
		hold on; MPRdata.h5=plot([1 1]*MPRdata.cy,[1 MPRdata.nz],'g'); hold off;
		hold on; MPRdata.h6=plot([1 MPRdata.ny],[1 1]*MPRdata.cz,'g'); hold off;
		mtext(0.01*MPRdata.ny,0.03*MPRdata.nz,sprintf('%d/%d',MPRdata.cx,MPRdata.nx),10,'g');

		if viz==0, 
			set([MPRdata.h1,MPRdata.h2,MPRdata.h3,MPRdata.h4,MPRdata.h5,MPRdata.h6],'Visible','off'); 
			if ~isempty(MPRdata.h7), set(MPRdata.h7,'Visible','off'); end
		end

		set(gcf,'UserData',MPRdata);
	case 'alt'
		MPRdata.handlebuttondown = 1; set(gcf,'UserData',MPRdata);

		cl = get(gca,'Clim');
		width  = cl(2)-cl(1);
		center = mean(cl);
		x0 = xy(1); y0=xy(2); j=0;

		while 1,
			[xy]=get(0,'PointerLocation');
			x = (x0-xy(1))/1280; y=(xy(2)-y0)/1024;

			c0 = center+width*y;
			w0 = width*max(x+1,1e-8);

			msubplot(2,2,1,0,0,0,0,0,0); 
			set(gca,'Clim',[c0-w0/2 c0+w0/2]);
			msubplot(2,2,2,0,0,0,0,0,0); 
			set(gca,'Clim',[c0-w0/2 c0+w0/2]);
			msubplot(2,2,3,0,0,0,0,0,0); 
			set(gca,'Clim',[c0-w0/2 c0+w0/2]);
			msubplot(2,2,4,0,0,0,0,0,0);
			set(gca,'Clim',[c0-w0/2 c0+w0/2]);

			pause(0.02);
			j=j+1;
			if j>2000/2, disp('timeout'); break; end
			MPRdata = get(gcf,'UserData');
			if ~MPRdata.handlebuttondown, break; end
		end
	case 'extend'
		MPRdata.handlebuttondown = 1; set(gcf,'UserData',MPRdata);

		pt = [(xy(1)-pos(1)) (xy(2)-pos(2))];
		set(gcf,'CurrentPoint',pt);
		pt = get(gca,'CurrentPoint');
		u = round(pt(1))+1; v = round(pt(3))-1;
		if xy(1)>pos(1)+pos(3)/2,
			if xy(2)>pos(2)+pos(4)/2,	% UR
				return
			else						% LR
				return
			end
		else
			if xy(2)>pos(2)+pos(4)/2,	% UL
				bx = max(min(v,MPRdata.nx),1);
				by = max(min(u,MPRdata.ny),1);
			else						% LL
				return
			end
		end

		j=0; hold on; h=plot([by by],[bx bx],'x-'); hold off;
		while 1,
			xy = get(0,'PointerLocation');
			pt = [(xy(1)-pos(1)) (xy(2)-pos(2))];
			set(gcf,'CurrentPoint',pt);
			pt = get(gca,'CurrentPoint');
			u = round(pt(1))+1; v = round(pt(3))-1;

			ex =  max(min(v,MPRdata.nx),1);
			ey =  max(min(u,MPRdata.ny),1);

			delete(h);
			hold on; h=plot([by ey],[bx ex],'x-'); set(h,'Color',[1 1 0]); hold off; drawnow;

			pause(0.2);

			j=j+1;
			if j>2000/2, disp('timeout'); break; end

			MPRdata = get(gcf,'UserData');
			if ~MPRdata.handlebuttondown, break; end
		end

		if ~isempty(MPRdata.h7), 
			delete(MPRdata.h7); MPRdata.h7=[]; drawnow;
		end

		delete(h);
		mpr_oblique(bx,by,ex,ey,MPRdata,clim);
	end
end

function mpr_slice(MPRdata,dim)
	clim = get(gca,'Clim');
	if MPRdata.thickness==1,
		switch dim,
		case 1, imagesc(MPRdata.data(:,:,MPRdata.cz)); hold on;
				for j=1:length(MPRdata.OV), set(imagesc(MPRdata.OV1{j}),'AlphaData',MPRdata.OV{j}(:,:,MPRdata.cz)); end
		case 2, imagesc(squeeze(MPRdata.data(:,MPRdata.cy,:)).'); hold on;
				for j=1:length(MPRdata.OV), set(imagesc(MPRdata.OV2{j}),'AlphaData',squeeze(MPRdata.OV{j}(:,MPRdata.cy,:)).'); end
		case 3, imagesc(squeeze(MPRdata.data(MPRdata.cx,:,:)).'); hold on;
				for j=1:length(MPRdata.OV), set(imagesc(MPRdata.OV3{j}),'AlphaData',squeeze(MPRdata.OV{j}(MPRdata.cx,:,:)).'); end
		end
	else
		rr = ceil((-MPRdata.thickness/2):(MPRdata.thickness/2));
		switch dim,
		case 1, imagesc(mean(MPRdata.data(:,:,min(max(MPRdata.cz+rr,1),MPRdata.nz)),3)); hold on;
				for j=1:length(MPRdata.OV), set(imagesc(MPRdata.OV1{j}),'AlphaData',min(sum(MPRdata.OV{j}(:,:,min(max(MPRdata.cz+rr,1),MPRdata.nz)),3),1)); end
		case 2, imagesc(squeeze(mean(MPRdata.data(:,min(max(MPRdata.cy+rr,1),MPRdata.ny),:),2)).'); hold on;
				for j=1:length(MPRdata.OV), set(imagesc(MPRdata.OV2{j}),'AlphaData',min(squeeze(sum(MPRdata.OV{j}(:,min(max(MPRdata.cy+rr,1),MPRdata.ny),:),2)).',1)); end
		case 3, imagesc(squeeze(mean(MPRdata.data(min(max(MPRdata.cx+rr,1),MPRdata.nx),:,:),1)).'); hold on;
				for j=1:length(MPRdata.OV), set(imagesc(MPRdata.OV3{j}),'AlphaData',min(squeeze(sum(MPRdata.OV{j}(min(max(MPRdata.cx+rr,1),MPRdata.nx),:,:),1)).',1)); end
		end
	end	
	axis image; set(gca,'Clim',clim,'XTick',[],'YTick',[],'XColor','k','Ycolor','k');  axis off;
end

function mpr_oblique(bx,by,ex,ey,MPRdata,clim)
	MPRdata.oblique=[bx,by,ex,ey];
	msubplot(2,2,1,0,0,0,0,0,0); 
	hold on; MPRdata.h7=plot([by ey],[bx ex],'x-'); set(MPRdata.h7,'Color',[1 1 0]); hold off;
	set(gcf,'UserData',MPRdata);

	d = round(sqrt((bx-ex)^2+(by-ey)^2));
	if d<10, delete(h); return; end

	thk = MPRdata.thickness;
	dx = (ey-by)/d; dy = (bx-ex)/d;
	msubplot(2,2,1,0,0,0,0,0,0); 
	if thk>1,
		hold on; plot([by ey]-dy*thk/2,[bx ex]-dx*thk/2,'yo--'); hold off;
		hold on; plot([by ey]+dy*thk/2,[bx ex]+dx*thk/2,'yo--'); hold off;
	end

	rx = linspace(bx,ex,d);
	ry = linspace(by,ey,d);
	[Ri,Rz]=ndgrid(1:d,1:MPRdata.nz);
	for j=1:thk,
		Rx(:,:,j) = rx(Ri)+(j-(thk+1)/2)*dx;
		Ry(:,:,j) = ry(Ri)+(j-(thk+1)/2)*dy;
		Rz(:,:,j) = Rz(:,:,1);
	end

	msubplot(2,2,4,0,0,0,0,0,0); 
	tmp = mean(interpn(MPRdata.data,Rx,Ry,Rz),3)';
	imagesc(tmp);
	cols = [1 0 0; 0 0 1; 0 1 0; 1 1 0; 1 0 1; 0 1 1;];
	for j=1:length(MPRdata.OV),
		hold on;
		ind = mod(j-1,length(MPRdata.OV))+1;
		overlay = cat(3, cols(ind,1)*ones(size(tmp)), cols(ind,2)*ones(size(tmp)), cols(ind,3)*ones(size(tmp)));
		set(imagesc(overlay),'AlphaData',min(sum(interpn(MPRdata.OV{j},Rx,Ry,Rz),3).',1));
	end

	h=colorbar('East'); set(h,'XColor','g','YColor','g');
	axis image; set(gca,'Clim',clim); axis off;
	colormap gray; hold off; drawnow;
	clear Rx Ry Rz tmp;
end

function mpr_buttonup(src,evnt)
	MPRdata = get(gcf,'UserData');
	MPRdata.handlebuttondown = 0; set(gcf,'UserData',MPRdata);
end

function mpr_keyrelease(src,evnt)
	MPRdata = get(gcf,'UserData');
	switch evnt.Character,
	case ' ',
		if strcmp(get(MPRdata.h1,'Visible'),'on'),
			set([MPRdata.h1,MPRdata.h2,MPRdata.h3,MPRdata.h4,MPRdata.h5,MPRdata.h6],'Visible','off');
			if ~isempty(MPRdata.h7), set(MPRdata.h7,'Visible','off'); end
		else
			set([MPRdata.h1,MPRdata.h2,MPRdata.h3,MPRdata.h4,MPRdata.h5,MPRdata.h6],'Visible','on');
			if ~isempty(MPRdata.h7), set(MPRdata.h7,'Visible','on'); end
		end
	case {'=','+'}
		MPRdata.thickness = MPRdata.thickness+1;
		disp(sprintf('Thickness is now = %d',MPRdata.thickness))
	case '-'
		MPRdata.thickness = max(MPRdata.thickness-1,1);
		disp(sprintf('Thickness is now = %d',MPRdata.thickness))
	end
	set(gcf,'UserData',MPRdata);
end

function cam=mpr_getcamera()
	cam.h1=get(gca,'CameraPosition');
	cam.h2=get(gca,'CameraTarget');
	cam.h3=get(gca,'CameraUpVector');
	cam.h4=get(gca,'CameraViewAngle');
end
function mpr_setcamera(cam)
	set(gca,'CameraPosition',cam.h1);
	set(gca,'CameraTarget',cam.h2);
	set(gca,'CameraUpVector',cam.h3);
	set(gca,'CameraViewAngle',cam.h4);
end

end


function h=msubplot(r,c,n,hs,vs,tb,bb,lb,rb)
% MSUBPLOT(r,c,n,hs,vs,tb,bb,lb,rb)
% Modified subplot utility
%      r - rows
%      c - columns
%      n - plot number
%      hs - horizontal spacing
%      vs - vertical spacing
%      tb - top border
%      bb - bottom border
%      lb - left border
%      rb - right border

% Web Stayman, January 1998

yi = ceil(n/c);
xi = mod(n-1,c)+1;

if ~exist('tb','var'), tb = 0.05; end
if ~exist('bb','var'), bb = 0.05; end
if ~exist('lb','var'), lb = 0.05; end
if ~exist('rb','var'), rb = 0.05; end
if ~exist('hs','var'), hs = 0.03; end
if ~exist('vs','var'), vs = 0.03; end

xd = (1-lb-rb-(c-1)*hs)/c;
yd = (1-tb-bb-(r-1)*hs)/r;

x = lb+(xi-1)*(xd+hs);
y = bb+(r-yi)*(yd+vs);

ih=subplot('Position',[x,y,xd,yd]);

if (nargout > 0)
	h=ih;
end
end


function [hh] = mtext(x,y,str,sz,col)
% [h] - Text(x,y,str,ptsz,col) - Better text function.
%   Inputs:
%      x,y - text position
%      str - string for label
%      sz  - point size text (default=10)
%      col - color of text (default=[0 0 0])
%   Outputs:
%      h   - Handle to label for current axes.

% Web Stayman, January 1998
% Needs to be extended for all text options.

if exist('str','var'),
	h = text(x,y,str);

	if ~exist('sz','var'),  sz=10; end
	if ~exist('col','var'), col=[0 0 0]; end
	
	set(h, 'FontAngle',  get(h, 'FontAngle'), ...
	    'FontName',   get(h, 'FontName'), ...
	    'FontSize',   get(h, 'FontSize'), ...
	    'FontWeight', get(h, 'FontWeight'), ...
	    'Rotation',   0, ...
	    'FontSize',   sz, ...
	    'Color',      col, ...
	    'string',     str);
end

if nargout>0,
	hh=h;
end
end
