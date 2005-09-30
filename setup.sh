top_srcdir=/home3/amundson/work/fnal/branches/jfa1
for dir in \
    $top_srcdir/physics_toolkit/src/.libs \
    $top_srcdir/beamline/src/.libs \
    $top_srcdir/bmlfactory/src/.libs \
    $top_srcdir/basic_toolkit/src/.libs \
    $top_srcdir/mxyzptlk/src/.libs
do
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$dir
done
