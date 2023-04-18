all : lib Rfuncs

Rfuncs: src/C/def.h src/C/dmatrix.h src/C/gen_lasso.h src/C/filtering.h src/C/util.h
	R CMD SHLIB src/C/Rfuncs.cpp src/C/gen_lasso.cpp src/C/filtering.cpp src/C/util.cpp src/C/def.c src/C/dmatrix.c -o lib/Rfuncs.so

clean: 
	rm -f src/C/*.o src/C/*.so lib/*.so

lib: 
	mkdir lib
