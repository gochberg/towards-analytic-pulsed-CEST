function wave = read_phased_rf_shape(file_name);

fid = fopen(file_name,'r');
if fid == -1
    str = sprintf('Can not open file %s',file_name);
    error(str);
end

d = fread(fid);
dd = char(d');

star_points = find(dd == '*');
last = max([star_points 0]);

ddd = dd(last+1:length(dd));
dddd = str2num(ddd);


wave = exp(i*dddd(:,1)*pi/180).*dddd(:,2);

% if abs(abs(real(wave))-abs(wave)) <= .0001*abs(real(wave))
%    wave = real(exp(i*dddd(:,1)*pi/180).*dddd(:,2));
% end

%wave = real(exp(i*dddd(:,1)*pi/180).*dddd(:,2));
% NOTE: real above is not correct if phase changing pulse, but works for
% real pulses that change sign. e.g. sinc.  And it avoids having a very
% small imag part.

fclose(fid);
