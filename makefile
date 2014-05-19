# specify compilers
FC = gfortran # /usr/bin/gfortran
CC = gcc # /usr/bin/gcc

# build the program

PROG = flow
PROG_OBJS = fluid.o

PROG_DIR = .

PGPLOT_DIR = $(MESASDK_ROOT)/pgplot
X11_LIB_DIR = /opt/X11/lib 
# /opt/X11/lib for hermes
# /usr/lib64/X11 for fedora

FCopenmp = -fopenmp

$(PROG) : $(PROG_OBJS)
	$(FC) $(FCbasic) $(FCopenmp) -o $(PROG_DIR)/$(PROG) $(PROG_OBJS) \
	$(LOAD_OTHER) -lz -L $(PGPLOT_DIR) -lpgplot -L $(X11_LIB_DIR) -lX11

#################################################################

FC_free_preprocess = -x f95-cpp-input
FCfree = -ffree-form $(FC_free_preprocess)
MY_FC_FLAGS = $(FCfree)
SRC_DIR = .

%.o: $(SRC_DIR)/%.f
	 $(FC) -g $(FCbasic) $(MY_FC_FLAGS) -c $<

clean:
	-@rm -f *.o *.mod flow
