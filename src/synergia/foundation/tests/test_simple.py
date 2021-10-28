import sys
sys.stdout = open("test_simple.out", "w")
print(f"sys.path is: {sys.path}")
import synergia
print("importing synergia succeeded")

def test_simple():
    assert True 


