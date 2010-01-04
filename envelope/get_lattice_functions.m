function lattice=get_lattice_functions(madfile_name,lat_name,kin_energy_GeV,...
				       num_markers,particle,latticecalc_path)
  
  if (nargin <5 )
    particle ="proton";
  endif
  if (nargin <6 )
    latticecalc_path="../latticecalc";
  endif

  out_file = sprintf("/tmp/latticecalc_%s",getenv("USER"));

  command = sprintf("%s %s %s %g %d 16 %s> %s",latticecalc_path,madfile_name,...
		    lat_name, kin_energy_GeV,num_markers,particle,out_file);
  printf("command = %s\n",command);
  system(command);

  try
    lattice = load(out_file);
  catch
    lattice = 0;
    printf("%s failed. Output was:\n",latticecalc_path);
    printf("---------------------------------------------------------\n");
    try
      system(sprintf("cat %s",out_file));
    catch
      printf("<absolutely nothing. latticecalc produced no output.>\n");
    end_try_catch
    printf("---------------------------------------------------------\n");
  end_try_catch


endfunction
