function M = readbin(fid,Nd)
%M = readbin(fid,Nd)
%
%  For use when postprocessing binary data in Matlab
%	AUTHOR:  Trevor Wood, Boston University
%	DATE: c. 2000
%
%   Reads array M (of Nd dimensions) from binary file, fid.
%   If Nd is not present, then one recond containing one integer is read.
%
   dum = fread(fid, 1,'int');
   if nargin < 2
      M   = fread(fid, 1,'int');
   else
      if(isstr(Nd))
         M   = fread(fid, 1, Nd);
      else
         N   = fread(fid,Nd,'int');
         M   = fread(fid,prod(N),'double');
         if(length(N) > 1)
            M = reshape(M,N');
         end
      end
   end
   dum = fread(fid, 1,'int');

