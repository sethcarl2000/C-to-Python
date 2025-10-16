# We build two versions of the shared library
# Nothing special is required to load the C++ libraries into python using ROOT
# If we use ctypes and load the functions by hand in python, they must be
# linked using the extern "C" setting.  Otherwise python will not be able to
# recognize the names of the functions (c.f. name mangling)


all: mylib.so mylib.soc mycpplib.so mycpplib.soc

mylib.so: mylib.cpp mylib.h # build shared library
	g++ -O3 -Wall -fPIC -shared -omylib.so mylib.cpp 

mylib.soc: mylib.cpp mylib.h  # build "so" w/ C style linkage
	g++ -O3 -Wall -fPIC -DEXTERNC -shared -omylib.soc mylib.cpp 

mycpplib.so: mycpplib.cpp mycpplib.h
	g++ -O3 -Wall -fPIC -shared -o mycpplib.so mycpplib.cpp

mycpplib.soc: mycpplib.cpp mycpplib.h
	g++ -O3 -Wall -fPIC -DEXTERNC -shared -o mycpplib.soc mycpplib.cpp

clean:
	rm -f mylib.so mylib.soc mycpplib.so mycpplib.soc *~

cleanall: clean
	rm -rf __pycache__ .ipynb_checkpoints
