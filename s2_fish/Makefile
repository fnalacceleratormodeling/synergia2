%.o : %.cc
	g++ -c -g $<

space_charge: scalar_field.o main.o Double_tensor.o
	g++ -o $@ $^

clean:
	rm -rf space_charge *.so *.o 
