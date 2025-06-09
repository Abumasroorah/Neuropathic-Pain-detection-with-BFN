function plot_topo(layout, data, label)
% Version 0.0
% Created and Modified by Yubo Wang
%%%%%%%%%%%%%%%%%%%%%%%%
% At least provied 2 funcitons
% label = 0, only generate a head model wiht sensor locaiton
% label = 1, amplitude mapping
% 
% The data is automatically classified into 'Connectivity' and 'AmMap' mode
% by its dimension.


% Data type Checking
  temp = size(data);
  if temp(1,1) >1
     datatype = 'Connectivity';
  else datatype = 'AmMap';
  end
  
% In order to generate multiplot at same figure
% following code should be modified
hpos = 0;
vpos = 0;
width         = 1 ;
height        = 1;
gridscale     = 67; % 67 in original
shading       = 'flat';
% interplim     = ft_getopt(varargin, 'interplim',      'electrodes');
interpmethod  = 'v4'; % for funciton of gridata.
mask = layout.mask;
outline = layout.outline;
isolines = 6; % used for contour plot, fixed to 6 Layer.
% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

%%%%
chanX =layout.pos(:,1);
chanY =layout.pos(:,2);

chanX = chanX(:) * width  + hpos;
chanY = chanY(:) * height + vpos;

% Get the limit for plot
  hlim = [inf -inf];
  vlim = [inf -inf];
  for i=1:length(mask)
    hlim = [min([hlim(1); mask{i}(:,1)*width+hpos]) max([hlim(2); mask{i}(:,1)*width+hpos])];
    vlim = [min([vlim(1); mask{i}(:,2)*width+vpos]) max([vlim(2); mask{i}(:,2)*width+vpos])];
  end

  % check if all mask point are inside the limits otherwise redefine mask
newpoints = [];
if length(mask)==1
  % which channels are outside
  outside = false(length(chanX),1);
  inside  = inside_contour1([chanX chanY], mask{1});
  outside = ~inside;
  newpoints = [chanX(outside) chanY(outside)];
end
  

% convert the mask into a binary image
  maskimage = false(gridscale);
  %hlim      = [min(chanX) max(chanX)];
  %vlim      = [min(chanY) max(chanY)];
  xi        = linspace(hlim(1), hlim(2), gridscale);   % x-axis for interpolation (row vector)
  yi        = linspace(vlim(1), vlim(2), gridscale);   % y-axis for interpolation (row vector)
  [Xi,Yi]   = meshgrid(xi', yi);
  
  for i=1:length(mask)
    mask{i}(:,1) = mask{i}(:,1)*width+hpos;
    mask{i}(:,2) = mask{i}(:,2)*height+vpos;
    mask{i}(end+1,:) = mask{i}(1,:);                   % force them to be closed
    maskimage(inside_contour1([Xi(:) Yi(:)], mask{i})) = true;
  end
  
  
  % plot the outline of the head, ears and nose

  for i=1:length(outline)
     xval = outline{i}(:,1) * width  + hpos;
     yval = outline{i}(:,2) * height + vpos;
     % This funciton is from Fieldtrip Toolbox
     ft_plot_vector(xval, yval, 'Color','g', 'LineWidth',0.6)
  end
  
  
 if strcmp(datatype,'AmMap') 
  % for amplitude plot
   % the data used to interpolate based on the grid
      xi = linspace(hlim(1), hlim(2), gridscale);       % x-axis for interpolation (row vector)
      yi = linspace(vlim(1), vlim(2), gridscale);       % y-axis for interpolation (row vector)
     [Xi,Yi,Zi] = griddata(chanX', chanY, data, xi', yi, interpmethod); % interpolate the topographic data
  if ~isempty(maskimage)
  % apply anatomical mask to the data, i.e. that determines that the interpolated data outside the circle is not displayed
     Zi(~maskimage) = NaN;
  end

  
   % Creat contour 
   [cont, h] = contour(Xi,Yi,Zi,isolines,'k');
   % Plot Color Surface
     deltax = xi(2)-xi(1); % length of grid entry
     deltay = yi(2)-yi(1); % length of grid entry
     h = surf(Xi-deltax/2,Yi-deltay/2,zeros(size(Zi)), Zi, 'EdgeColor', 'interp', 'FaceColor', shading); % replaced none --> interp : allows to change color of edge
 end % Amplitude Map
 
 
    % For Sensor Position plot
   
       % convert between true/false/yes/no etc. statements
       % point   = istrue(point);
       % box     = istrue(box);
%         label   = istrue(label);
%         mask    = istrue(mask);
%         outline = istrue(outline);
%         verbose = istrue(verbose);
        % Location Data Per
        X      = layout.pos(:,1) + hpos;
        Y      = layout.pos(:,2) + vpos;
        Width  = layout.width;
        Height = layout.height;
        Lbl    = layout.label;
        
        % Plot sensor Location
        plot(X, Y, 'marker','.','color','k','markersize',8,'linestyle','none'); % marker, marker size and color changed
 
 
 
 % Connectivity Pattern
  if strcmp(datatype,'Connectivity') 
     %  For connectivity pattern, the data matrix contain the paried
     %  locations. Add a line to connect each location
     % Version 0.0 for unweighted map only
     
     % Select the channle to be linkd
     [pair_s,pair_e]=find(data==1);
     
     xbeg = layout.pos(pair_s,1);
     ybeg = layout.pos(pair_s,2);
     xend = layout.pos(pair_e,1);
     yend = layout.pos(pair_e,2);
     
      x = [xbeg xend]';
      y = [ybeg yend]';
      h = line(x, y);
      C = [0 1 0]; % color
      h = patch(x, y, C);
    color =  set(h,'EdgeColor','g');
%       h = patch(x,y,'g','Tag','PatchBorder');
%       set(findobj('Tag','PatchBorder'),'FaceColor','w');
     
  end
 
 
 

        

    
   % For Label Names

   if label
      labeloffset = 0;
      Lbl = Lbl(1:length(Lbl)-2);
%       text(X+labeloffset, Y+(labeloffset*1.5), Lbl,'fontsize',0.5);   %%%%% labels font size changed 
   end
    % Refine the plot 
     for i=1:length(layout.outline)
        if ~isempty(layout.outline{i})
      X = layout.outline{i}(:,1) + hpos;
      Y = layout.outline{i}(:,2) + vpos;
      h = line(X, Y);
      c = [0 0 0];
      set(h, 'color', c);                         %%%%%%%% Outer look color modified
      set(h, 'linewidth', 1);
        end
     end
     
     
%       X = layout.mask{1}(:,1) + hpos;
%       Y = layout.mask{1}(:,2) + vpos;
%       % the polygon representing the mask should be closed
%       X(end+1) = X(1);
%       Y(end+1) = Y(1);
%       h = line(X, Y);
%       set(h, 'color', 'k');
%       set(h, 'linewidth', 1.5);
%       set(h, 'linestyle', ':');

     
axis auto
%axis fill
axis equal
axis([-0.55 0.555 -0.55 0.6]); % added 
axis off
end % end of the function


