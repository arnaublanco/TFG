function itd = itmark(itd, itdm);
% Usage:  
%  itd = itmark(itd, itdm);
%            itdm: total iterations expected
%                  0 for not known
%
%  %Inisialization: 
%    itd = itmark(0, 10);
%    itd = itmark(itd, 10);
%  %End: 
%    itd = itmark(-1, 10); 
%
%  Alex Perera, Feb 2003
  
  
  
  mv = {'\\','|','/','-'};
  mt = {'.','-','o','*','O'};
  mvl = length(mv);
  i = rem(itd,mvl)+1;

  if itd==-1,
    if itdm <1,
      fprintf( '\b' );
    else
      fprintf( '\b\b');
    end
  else
    if itdm <1,
      if itd ==0,
        fprintf(mv{i} );
      else
        fprintf( ['\b'  mv{i} ]);
      end
    else 
      if itd ==0,
        fprintf([ mt{1} mv{i} ]);
      else
        t = ceil(length(mt)*itd/itdm);
        fprintf( ['\b\b'  mt{min(t, length(mt))}  mv{i} ]);
      end
    end
    itd=itd+1;
  end
    
    
  
