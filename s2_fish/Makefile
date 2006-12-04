%.o : %.cc
	g++ -c -g $<

space_charge: scalar_field.o main.o Double_tensor.o
	g++ -o $@ $^

scalar_field.o: scalar_field.h Double_tensor.h

clean:
	rm -rf space_charge *.so *.o 
