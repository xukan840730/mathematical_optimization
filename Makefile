#################################################################
#################################################################

CC := g++

EIGEN_INC := eigen

CFLAGS := -std=c++11 -I$(EIGEN_INC) -Wno-deprecated-declarations -g 

#################################################################
## output directory for everything.
#################################################################

MKDIR_P := "mkdir" -p
OUTPUT_DIR := output

#################################################################
## COMMON
#################################################################
COMMON_LIB_NAME := common_lib
COMMON_LIB_DIR := common
COMMON_LIB_OUTPUT_DIR := $(OUTPUT_DIR)/$(COMMON_LIB_DIR)
COMMON_LIB_LIB := $(COMMON_LIB_OUTPUT_DIR)/lib$(COMMON_LIB_NAME).a

COMMON_LIB_SRC := $(wildcard $(COMMON_LIB_DIR)/*.cpp)
COMMON_LIB_OBJ := $(patsubst %,$(COMMON_LIB_OUTPUT_DIR)/%, $(notdir $(COMMON_LIB_SRC:.cpp=.o)))

#################################################################
## LINEAR_ALGEBRA
#################################################################
LINEAR_ALGEBRA_NAME := linear_algebra
LINEAR_ALGEBRA_DIR := linear_algebra
LINEAR_ALGEBRA_OUTPUT_DIR := $(OUTPUT_DIR)/$(LINEAR_ALGEBRA_DIR)
LINEAR_ALGEBRA_LIB := $(LINEAR_ALGEBRA_OUTPUT_DIR)/lib$(LINEAR_ALGEBRA_NAME).a

LINEAR_ALGEBRA_SRC := $(wildcard $(LINEAR_ALGEBRA_DIR)/*.cpp)
LINEAR_ALGEBRA_OBJ := $(patsubst %,$(LINEAR_ALGEBRA_OUTPUT_DIR)/%, $(notdir $(LINEAR_ALGEBRA_SRC:.cpp=.o)))

#################################################################
## LINE SEARCH
#################################################################
LINE_SEARCH_NAME := line_search
LINE_SEARCH_DIR := line_search
LINE_SEARCH_OUTPUT_DIR := $(OUTPUT_DIR)/$(LINE_SEARCH_NAME)
LINE_SEARCH_LIB := $(LINE_SEARCH_OUTPUT_DIR)/lib$(LINE_SEARCH_NAME).a

LINE_SEARCH_SRC := $(wildcard $(LINE_SEARCH_DIR)/*.cpp)
LINE_SEARCH_OBJ := $(patsubst %,$(LINE_SEARCH_OUTPUT_DIR)/%, $(notdir $(LINE_SEARCH_SRC:.cpp=.o)))

#################################################################
## UNCONSTRAINED OPTIMIZATION
#################################################################
UNCONSTRAINED_OPTIMIZATION_NAME := unconstr_opt
UNCONSTRAINED_OPTIMIZATION_DIR := unconstr_opt
UNCONSTRAINED_OPTIMIZATION_OUTPUT_DIR := $(OUTPUT_DIR)/$(UNCONSTRAINED_OPTIMIZATION_NAME)
UNCONSTRAINED_OPTIMIZATION_LIB := $(UNCONSTRAINED_OPTIMIZATION_OUTPUT_DIR)/lib$(UNCONSTRAINED_OPTIMIZATION_NAME).a

UNCONSTRAINED_OPTIMIZATION_SRC := $(wildcard $(UNCONSTRAINED_OPTIMIZATION_DIR)/*.cpp)
UNCONSTRAINED_OPTIMIZATION_OBJ := $(patsubst %,$(UNCONSTRAINED_OPTIMIZATION_OUTPUT_DIR)/%, $(notdir $(UNCONSTRAINED_OPTIMIZATION_SRC:.cpp=.o)))

#################################################################
## CONSTRAINED OPTIMIZATION
#################################################################
CONSTRAINED_OPTIMIZATION_NAME := constr_opt
CONSTRAINED_OPTIMIZATION_DIR := constr_opt
CONSTRAINED_OPTIMIZATION_OUTPUT_DIR := $(OUTPUT_DIR)/$(CONSTRAINED_OPTIMIZATION_NAME)
CONSTRAINED_OPTIMIZATION_LIB := $(CONSTRAINED_OPTIMIZATION_OUTPUT_DIR)/lib$(CONSTRAINED_OPTIMIZATION_NAME).a

CONSTRAINED_OPTIMIZATION_SRC := $(wildcard $(CONSTRAINED_OPTIMIZATION_DIR)/*.cpp)
CONSTRAINED_OPTIMIZATION_OBJ := $(patsubst %,$(CONSTRAINED_OPTIMIZATION_OUTPUT_DIR)/%, $(notdir $(CONSTRAINED_OPTIMIZATION_SRC:.cpp=.o)))


#################################################################
## TEST CASES
#################################################################
TEST_CASE_NAME := test_case
TEST_CASE_DIR := test_case
TEST_CASE_OUTPUT_DIR := $(OUTPUT_DIR)/$(TEST_CASE_DIR)

TEST_CASE_SRC := $(wildcard $(TEST_CASE_DIR)/*.cpp)
TEST_CASE_OBJ := $(patsubst %,$(TEST_CASE_OUTPUT_DIR)/%, $(notdir $(TEST_CASE_SRC:.cpp=.o)))

PROGRAM := $(OUTPUT_DIR)/main

LFLAGS := -l$(COMMON_LIB_NAME) -L$(COMMON_LIB_OUTPUT_DIR) -l$(LINEAR_ALGEBRA_NAME) -L$(LINEAR_ALGEBRA_OUTPUT_DIR) -l$(LINE_SEARCH_NAME) -L$(LINE_SEARCH_OUTPUT_DIR) -l$(UNCONSTRAINED_OPTIMIZATION_NAME) -L$(UNCONSTRAINED_OPTIMIZATION_OUTPUT_DIR)

all : directories $(PROGRAM)

$(PROGRAM) : $(TEST_CASE_OBJ) $(COMMON_LIB_LIB) $(LINEAR_ALGEBRA_LIB) $(LINE_SEARCH_LIB) $(UNCONSTRAINED_OPTIMIZATION_LIB) $(CONSTRAINED_OPTIMIZATION_LIB)
	$(CC) -o $@ $^ $(LFLAGS)
	
directories : $(OUTPUT_DIR) $(COMMON_LIB_OUTPUT_DIR)  $(LINEAR_ALGEBRA_OUTPUT_DIR) $(LINE_SEARCH_OUTPUT_DIR) $(UNCONSTRAINED_OPTIMIZATION_OUTPUT_DIR) $(CONSTRAINED_OPTIMIZATION_OUTPUT_DIR) $(TEST_CASE_OUTPUT_DIR)

$(TEST_CASE_OUTPUT_DIR)/%.o : $(TEST_CASE_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
	
#################################################################
## create directories
#################################################################
$(OUTPUT_DIR) : 
	$(MKDIR_P) $(OUTPUT_DIR)

$(COMMON_LIB_OUTPUT_DIR) :
	$(MKDIR_P) $(COMMON_LIB_OUTPUT_DIR)
	
$(LINEAR_ALGEBRA_OUTPUT_DIR) :
	$(MKDIR_P) $(LINEAR_ALGEBRA_OUTPUT_DIR)
	
$(LINE_SEARCH_OUTPUT_DIR) :
	$(MKDIR_P) $(LINE_SEARCH_OUTPUT_DIR)
	
$(UNCONSTRAINED_OPTIMIZATION_OUTPUT_DIR) : 
	$(MKDIR_P) $(UNCONSTRAINED_OPTIMIZATION_OUTPUT_DIR)
	
$(CONSTRAINED_OPTIMIZATION_OUTPUT_DIR) : 
	$(MKDIR_P) $(CONSTRAINED_OPTIMIZATION_OUTPUT_DIR)
	
$(TEST_CASE_OUTPUT_DIR) :
	$(MKDIR_P) $(TEST_CASE_OUTPUT_DIR)
	
#################################################################
## Program Libs
#################################################################

$(COMMON_LIB_OUTPUT_DIR)/%.o : $(COMMON_LIB_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(COMMON_LIB_LIB) : $(COMMON_LIB_OBJ)
	$(AR) cr $@ $^

$(LINEAR_ALGEBRA_LIB) : $(LINEAR_ALGEBRA_OBJ)
	$(AR) cr $@ $^
	
$(LINEAR_ALGEBRA_OUTPUT_DIR)/%.o : $(LINEAR_ALGEBRA_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(LINE_SEARCH_LIB) : $(LINE_SEARCH_OBJ)
	$(AR) cr $@ $^
	
$(LINE_SEARCH_OUTPUT_DIR)/%.o : $(LINE_SEARCH_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(UNCONSTRAINED_OPTIMIZATION_LIB) : $(UNCONSTRAINED_OPTIMIZATION_OBJ)
	$(AR) cr $@ $^
	
$(UNCONSTRAINED_OPTIMIZATION_OUTPUT_DIR)/%.o : $(UNCONSTRAINED_OPTIMIZATION_DIR)/%.cpp	
	$(CC) $(CFLAGS) -c $< -o $@

$(CONSTRAINED_OPTIMIZATION_LIB) : $(CONSTRAINED_OPTIMIZATION_OBJ)
	$(AR) cr $@ $^
	
$(CONSTRAINED_OPTIMIZATION_OUTPUT_DIR)/%.o : $(CONSTRAINED_OPTIMIZATION_DIR)/%.cpp	
	$(CC) $(CFLAGS) -c $< -o $@
	
#################################################################
## utility
#################################################################

.PHONY : help clean

clean :
	$(RM) $(PROGRAM) $(COMMON_LIB_LIB) $(LINE_SEARCH_LIB) $(LINE_SEARCH_OBJ) $(LINEAR_ALGEBRA_LIB) $(LINEAR_ALGEBRA_OBJ) $(UNCONSTRAINED_OPTIMIZATION_LIB) $(UNCONSTRAINED_OPTIMIZATION_OBJ) $(CONSTRAINED_OPTIMIZATION_LIB) $(CONSTRAINED_OPTIMIZATION_OBJ) $(TEST_CASE_OBJ)

help :
	#@echo "LINE_SEARCH_SRC:=$(LINE_SEARCH_SRC)"
	#@echo "LINE_SEARCH_OBJ_DIR:=$(LINE_SEARCH_OUTPUT_DIR)"
	#@echo "LINE_SEARCH_OBJ:=$(LINE_SEARCH_OBJ)"
	#@echo "LINE_SEARCH_LIB:=$(LINE_SEARCH_LIB)"
	@echo "CONSTR_OPT_SRC:=$(CONSTRAINED_OPTIMIZATION_SRC)"
	@echo "CONSTR_OPT_OBJ:=$(CONSTRAINED_OPTIMIZATION_OBJ)"
	
