#/bin/sh

mods="beam_parameters computational_grid density_plot field function_cache gourmet job_manager loadfile matching memory options physics_constants processor_grid syn2_diagnostics_impact syn2_diagnostics tracker diagnostics bunch"

cp $1 translate.tmp1
for mod in $mods
do
    grep -v "import $mod" translate.tmp1 |     sed "s/$mod\\./synergia./g" > translate.tmp2
    mv translate.tmp2 translate.tmp1
done
mv translate.tmp1 $1
