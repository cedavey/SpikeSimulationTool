% SETTOOLTIPS Assigns each item in the application its tooltip (description
% on mouse hover).
%
% Syntax: 
%     h = setTooltips(h); % Sets all the default tooltips to all the
%                               default elements
%     h = setTooltips(h, {'object1' 'object2'}, {'String 1', 'String 2'});
%
% Inputs:
%     h        - GUI handles
%     object   - (Optional) updates object's tooltip with input string
%                 'str'. Can be a single object's name or a vector of names
%     str      - (Optional) Input string to update object's tooltip. It can
%                 be a single string or a cell vector.
%     
% Outputs:
%     h        - GUI handles
%
% Artemio - 11/July/2019
function h = setSstTooltips(h, varargin)
   if nargin == 1
      object = {'nAxonsTextBox', 'totalTimeTextBox', 'samplingRateTextBox',...
      'maxSpikeRateTextBox', 'overlapCheckBox', 'repeatTemplateCheckBox',...
      'snrTextBox', 'samplingRateTextBox', 'hasNoiseCheckBox',...
      'hasDriftCheckBox', 'preNoiseCheckBox', 'doFilterCheckBox',...
      'lowBandTextBox', 'highBandTextBox', 'plotCheckBox', 'runButton'};
      
      str = getSstTooltips(object);

   elseif nargin ~= 3
      str = sprintf('\tInvalid number of parameters ''getTooltips()''\n');
      cprintf('Errors', str);
      return      
      
   else
      object = varargin{1};
      str = varargin{2};
   end
   
   for i = 1:min(length(object), length(str))
      set(h.(object{i}), 'TooltipString', str{i});
   end   
end