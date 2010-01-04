function K=Kxy(s,index)
  global Kx_array;
  global Ky_array;
  global s_array;

  if(index == 1) 
    K_array = Kx_array;
  elseif (index == 2)
    K_array = Ky_array;
  else
    printf("Kxy error: you stink! index is %d\n",index);
    K = 0;
    return;
  endif

  for i=1:length(s)
    if (s(i) < s_array(1)) || (s(i) > s_array(length(s_array)))
      K(i) = 0.0;
    else
      K(i) = interp1(s_array,K_array,s(i),"linear");
    endif
  endfor

endfunction
