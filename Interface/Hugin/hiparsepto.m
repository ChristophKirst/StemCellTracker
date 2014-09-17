function pto = hiparsepto(ptofile)
%
% pto = hiparsepto(ptofile)
%
%
% description:
%    parse hugin pto file for relevant parameter
%
% input:
%   ptofile   the pto filename
%
% output:
%   pto       strut with parameter
%
% note: we only parse the shift parameter here

if ~isfile(ptofile)
   error('hiparsepto: no such pto file %s', ptofile);
end

% get image spec lines starting with i
fid = fopen(ptofile);
imgspecs = {}; i = 1;

ln = fgetl(fid);
while ischar(ln)
   if ~isempty(ln) && ln(1) == 'i'
      imgspecs{i} = ln; %#ok<AGROW>
      i = i + 1;
   end
   ln = fgetl(fid);
end
fclose(fid);

ni = length(imgspecs);
pto = struct();
for i = ni:-1:1
   ln = imgspecs{i};
   pw = regexp(ln, 'w');
   ph = regexp(ln, 'h');
   phe= regexp(ln, 'f');

   pv = regexp(ln, 'v');
   pve= regexp(ln, 'Ra');
   px = regexp(ln, 'TrX');
   py = regexp(ln, 'TrY');
   pz = regexp(ln, 'TrZ');
   pe = regexp(ln ,'j');
   
   pto(i).w   = str2double(ln(pw+1 : ph-1));
   pto(i).h   = str2double(ln(ph+1 : phe-1)); 
   pto(i).v   = str2double(ln(pv+1 : pve-1));
   pto(i).TrX = str2double(ln(px+3 : py-1));
   pto(i).TrY = str2double(ln(py+3 : pz-1));
   pto(i).TrZ = str2double(ln(pz+3 : pe-1));
   
   
   vs = {'Va', 'Vb', 'Vc', 'Vd', 'Vx', 'Vy', 'Vm'};
   for v = 1:length(vs)-1
      pv1 = regexp(ln, vs{v});
      pv2 = regexp(ln, vs{v+1});
      pto(i).(vs{v}) = str2double(ln(pv1+2 : pv2-1));
   end
   
   e1 = regexp(ln,'Eev');
   e2 = regexp(ln, 'Er');
   
   pto(i).Eev = str2double(ln(e1+3:e2-1));
  
   %pn = regexp(ln ,'n');
   %pto(i).filename = ln(pn+2:end-1);
end

end