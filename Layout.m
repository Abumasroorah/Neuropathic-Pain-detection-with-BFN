function lay = Layout(config)
% Head Model
 % Read Layout File % Currently only support EEG10/20 system
  fprintf('reading layout from file %s\n', config.layout);
  lay = readlay(config.layout);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % check whether outline and mask are available
  % if not, add default "circle with triangle" to resemble the head
  % in case of "circle with triangle", the electrode positions should also be
  % scaled
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    rmax  = 0.5;
    l     = 0:2*pi/100:2*pi;
    HeadX = cos(l).*rmax;
    HeadY = sin(l).*rmax;
    NoseX = [0.18*rmax 0 -0.18*rmax];
    NoseY = [rmax-.004 rmax*1.15 rmax-.004];
    EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
    EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
    % Scale the electrode positions to fit within a unit circle, i.e. electrode radius = 0.45
    ind_scale = strmatch('SCALE', lay.label);
    ind_comnt = strmatch('COMNT', lay.label);
    sel = setdiff(1:length(lay.label), [ind_scale ind_comnt]); % these are excluded for scaling
    x = lay.pos(sel,1);
    y = lay.pos(sel,2);
    % Only leaves me the location that being used in plot
    lay = rmfield(lay,'pos'); 
    lay.pos(:,1)= x;
    lay.pos(:,2) = y;
    
    xrange = range(x);
    yrange = range(y);
    % First scale the width and height of the box for multiplotting
    lay.width  = lay.width./xrange;
    lay.height = lay.height./yrange;
    % Then shift and scale the electrode positions
    lay.pos(:,1) = 0.9*((lay.pos(:,1)-min(x))/xrange-0.5);
    lay.pos(:,2) = 0.9*((lay.pos(:,2)-min(y))/yrange-0.5);
    % Define the outline of the head, ears and nose
    lay.outline{1} = [HeadX(:) HeadY(:)];
    lay.outline{2} = [NoseX(:) NoseY(:)];
    lay.outline{3} = [ EarX(:)  EarY(:)];
    lay.outline{4} = [-EarX(:)  EarY(:)];
    % Define the anatomical mask based on a circular head
    lay.mask{1} = [HeadX(:) HeadY(:)];

     
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% read the layout information from the ascii file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lay = readlay(filename)
if ~exist(filename, 'file')
  error(sprintf('could not open layout file: %s', filename));
end
[chNum,X,Y,Width,Height,Lbl,Rem] = textread(filename,'%f %f %f %f %f %q %q');

if length(Rem)<length(Lbl)
  Rem{length(Lbl)} = [];
end

for i=1:length(Lbl)
  if ~isempty(Rem{i})
    % this ensures that channel names with a space in them are also supported (i.e. Neuromag)
    Lbl{i} = [Lbl{i} ' ' Rem{i}];
  end
end
lay.pos    = [X Y];
lay.width  = Width;
lay.height = Height;
lay.label  = Lbl;
return % function readlay
