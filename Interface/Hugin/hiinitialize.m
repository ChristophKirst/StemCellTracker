function hiinitialize()
%
% hiinitialize()
%
% description:
%      checks existence of the needed command line tools
%
% todo: windows version, include tools into the software

global hitools

if isunix()

   [~, hitools] = hipath();

   tools = hitools.keys;
   cmds = hitools.values;
   for i = 1:hitools.Count
      fprintf('hiinitialize: tool %s -> %s\n', tools{i}, cmds{i});
   end
   
   fprintf('hiinitialize: hugin tools initialized!\n');
   
else
   
   error('hiinitialize: Hugin interface for Windows not implemented!')
   
end
   
end