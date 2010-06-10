import octapy
import numarray

def test_get(oct,name):
    print "testing getting",name,"from octave:"
    value = oct.get_value(name)
    print name,"=",value

def test_set(oct,name,value):
    print "testing setting",name,"in octave:"
    print "Python:",name,"=",value
    oct.set_value(name,value)
    test_get(oct,name)

o = octapy.Octave()

print "********** Execute some code in Octave:"
cmd = """a = 1.1;
b = 2.2;
c = sqrt(-4.0);"""
print cmd
o.execute(cmd)

print "\n********** See which variables are currently defined in Octave:"
tlv = o.top_level_vars()
print "top level variables =", tlv

print "\n********** Get the values of the variables back from Octave:"
test_get(o,"a")
test_get(o,"b")
test_get(o,"c")

print "\n********** Set those variables directly from Python..."
test_set(o,"a",17)
test_set(o,"b","strings work, too")
print "\n********** and a new one..."
test_set(o,"d",complex(3.0,4.0))

print "\n********** Do something a little more complicated in Octave:"
cmd = """my_x = [1.0,2.0,7.0,9.0];
my_xx = [11.0,12.0,13.0;21.0,22.0,23.0];
function [avg,dev] = get_stats(x)
    avg = mean(x);
    dev = std(x);
endfunction
[my_avg,my_dev] = get_stats(my_x);
"""
print cmd
o.execute(cmd)
test_get(o,"my_avg")
test_get(o,"my_dev")

print "\n********** Vectors and matrices are handled using numarray:"
test_get(o,"my_x")
test_get(o,"my_xx")
my_xx = o.get_value("my_xx")
print "my_xx.is_c_array() =", my_xx.is_c_array()
print "my_xx.is_f_array() =", my_xx.is_f_array()
test_set(o,"x1",numarray.array([1,2,3]))
test_set(o,"x2",numarray.array([[1,2],[3,4]]))


