CXX=g++

## make dir to output obj files
OBJDIR := obj

DEBUG ?= 0
ifeq ($(DEBUG), 1)
    CXXFLAGS=-std=c++11 -g -O0
else ifeq ($(DEBUG), 2)
		CXXFLAGS=-std=c++11 -g -O2
else
    CXXFLAGS=-std=c++11 -Wall -W -O2
endif

## specify the eigen and fwdpp directory
CPPFLAGS=-I/users/pyjiang/expr_evol/cell_type/eigen-3.3.1 -I/users/pyjiang/expr_evol/cell_type/fwdpp/fwdpp-0.5.1

## create rule for generating .o files from .cc files (because the source files are in a different dirrectory, cannot use
## implicit rule)
$(OBJDIR)/%.o : src/%.cc
	$(CXX) -c  $(CXXFLAGS)  $(CPPFLAGS) $< -o $@


# ## list all the obj files
IO_OBJ := $(OBJDIR)/ct_IO_high_level.o
BASIC_SIMU_OBJ := $(OBJDIR)/ct_haploid_save_state_new.o
BASIC_RESUME_OBJ := $(OBJDIR)/ct_haploid_resume_new.o
RAND_B_OBJ := $(OBJDIR)/ct_haploid_save_rand_known_geno_target.o
FIX_ALPHA_INIT_OBJ 	:= $(OBJDIR)/ct_haploid_fixed_rob_alpha_init.o
FIX_ALPHA_OBJ := $(OBJDIR)/ct_haploid_fixed_diff_rob_alpha.o
EVOL_ALPHA_OBJ := $(OBJDIR)/ct_haploid_diff_init_rob_geno_evol_two_loc_no_h.o

$(IO_OBJ)   $(BASIC_SIMU_OBJ)  $(BASIC_RESUME_OBJ)  : | $(OBJDIR)

## create obj dir
$(OBJDIR):
	mkdir $(OBJDIR)

## add new which make all the targets that are needed
all: haploid_save_new   haploid_resume_new   haploid_fixed_rob_alpha_init haploid_diff_fixed_rob_alpha   haploid_diff_init_rob_evol

## basic haploid simualtion code with fixed robustness determined by gamma
haploid_save_new:  $(BASIC_SIMU_OBJ) $(IO_OBJ)
	$(CXX) -o ct_haploid_save_state_new  $(BASIC_SIMU_OBJ) $(IO_OBJ)  -lgsl -lgslcblas -lz

## resuming running basic haploid
haploid_resume_new:  $(BASIC_RESUME_OBJ)  $(IO_OBJ)
	$(CXX) -o ct_haploid_resume_new $(BASIC_RESUME_OBJ)  $(IO_OBJ) -lgsl -lgslcblas -lz

## set up simulation with initial b vector generated from a random genotype vector
haploid_save_rand_b_geno: $(RAND_B_OBJ)   $(IO_OBJ)
	$(CXX) -o haploid_save_rand_b_geno  $(RAND_B_OBJ)   $(IO_OBJ) -lgsl -lgslcblas -lz

## this function is used to simulate initial genotype, geno-pheno map which will be used to apply different alpha
## for setting up different robustness by alpha
haploid_fixed_rob_alpha_init:  $(FIX_ALPHA_INIT_OBJ) $(IO_OBJ)
	$(CXX) -o ct_haploid_fixed_rob_alpha_init   $(FIX_ALPHA_INIT_OBJ) $(IO_OBJ)   -lgsl -lgslcblas -lz

## set up simulations for fixed different alpha (after init)
haploid_diff_fixed_rob_alpha:  $(FIX_ALPHA_OBJ)  $(IO_OBJ)
	$(CXX) -o ct_haploid_fixed_diff_rob_alpha  $(FIX_ALPHA_OBJ) $(IO_OBJ)   -lgsl -lgslcblas -lz

## initial population with different degrees of robustness and evolve robustness
haploid_diff_init_rob_evol : $(EVOL_ALPHA_OBJ)  $(IO_OBJ)
	$(CXX) -o ct_haploid_diff_init_rob_evol  $(EVOL_ALPHA_OBJ)  $(IO_OBJ)   -lgsl -lgslcblas -lz


clean: | $(OBJDIR)
		rm -f $(OBJDIR)/*.o
