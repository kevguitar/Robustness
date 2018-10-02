# Define what target we normally want to make
default: robustness_main
	
# You may need to adjust the compiler, compiler flags or library paths. -lgsl
# and -lgslcblas refer to the Gnu Scientific Library, which must be installed,
# as well as the openSSL library, which is linked via -lssl and -lcrypto.

SHELL = /bin/sh

S = ./source/

# Settings for C++ compiler.
CPP = g++
CPPFLAGS = -g -O0 
CPPINCLUDES = -I/usr/local/include

SOURCE = $(S)robustness_main.cc \
         $(S)utils.cc \
         $(S)print_to_file.cc \
         $(S)lin_gauss.cc \
         $(S)rob_analyt.cc \
         $(S)parse_from_init_file.cc \
         $(S)md5hash.cc

OBJ = $(S)robustness_main.o \
      $(S)utils.o \
      $(S)print_to_file.o \
      $(S)lin_gauss.o \
      $(S)rob_analyt.o \
      $(S)parse_from_init_file.o \
      $(S)md5hash.o

# Linker settings
LINKER_CMD = g++
AR_CMD = ar
LDFLAGS = $(CPPFLAGS) -L/usr/local/lib 
LDLIBS = -lgsl -lgslcblas -lm -lssl -lcrypto

%.o: %.cc
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) -c $*.cc -o $*.o

# Make the main program
robustness_main : $(SOURCE)
	$(LINKER_CMD) $(LDFLAGS) $(SOURCE) -o robustness_main $(LDLIBS)

clean:
	rm -f robustness_main $(S)*.o 
