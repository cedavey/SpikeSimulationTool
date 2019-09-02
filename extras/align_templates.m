% ALIGN_TEMPLATES Aligns to the peak value
%
% Syntax:
%     d = align_templates(templates);
%
% Inputs:
%     templates - matrix with different spike shapes to align to the
%                 maximum value
%     normalize - (Optional) Binary. If true, it normalizes the amplitude
%                 of the templates
%
% Outputs:
%     d         - Aligned <and normalized> templates
%
% Artemio Soto-Breceda | 15 - August - 2019
function d = align_templates(templates, varargin)
   if nargin > 1
      normalize = varargin{1};
   else 
      normalize = true;
   end
   
   if normalize
      templates = normalize_templates(templates);
   end
   
   [~, idx] = max(templates);
   mi = max(idx);
   for i = 1:size(templates, 2)
      d_ = templates(:,i);
      [~, idx_] = max(d_);
      if idx_ < mi
         d_ = [zeros(mi - idx_,1); d_(1:end - (mi - idx_))];
      end
      d(:,i) = d_;
   end
end

function d = normalize_templates(templates)
   for i = 1:size(templates, 2)
      dd = templates(:,i);
      dd = dd / abs(max(dd));
      d(:,i) = dd;
   end
end

% Removes extra samples in a template (leaving the spikes with info only
% within zero cross)
function d = clean_templates(templates)
   d = templates; % Assign the name of the variable to 'd'
   new_templates = zeros(size(d));
   maxLength = 0;
   for i = 1:size(d,2)
      dd = d(:,i)>0;
      ddd = diff(dd);
      dddd = find(ddd == 1);
      try
         d_ = (d(dddd(1)+1 : dddd(2), i));
      catch E
         d_ = d(:,i);
      end

      new_templates(1:numel(d_), i) = d_;

      if numel(d_) > maxLength, maxLength = numel(d_); end
   end

   % Replace the original variable
   d = new_templates(1:maxLength, :);
end

function saveFile(vsim)
   %% Save file
   [file,path] = uiputfile(['simulations' filesep 'sim_voltage.mat'],'Save file name');
   if file
      file_name = [path filesep file];
      save(file_name, 'vsim');
   else
      fprintf('\tUser didn''t chose a file location. The simulation wasn''t saved.\n');
   end
end