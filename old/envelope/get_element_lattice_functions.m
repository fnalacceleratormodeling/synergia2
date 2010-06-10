function lattice=get_element_lattice_functions(madfile_name,lat_name,...
					       kin_energy_GeV,...
					       elementcalc_path)
  
  if (nargin < 4)
    elementcalc_path="../elementcalc"
  endif

  out_file = sprintf("/tmp/elementcalc_%s",getenv("USER"));

  command = sprintf("%s %s %s %g | grep -v error > %s",elementcalc_path,madfile_name,...
		    lat_name, kin_energy_GeV, out_file);

  system(command);

  lattice = load(out_file);

endfunction
