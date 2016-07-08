#################################################################
#################################################################

CC := g++

EIGEN_INC := eigen

CFLAGS := -std=c++11 -I$(EIGEN_INC) 

#################################################################
## output directory for everything.
#################################################################

MKDIR_P := "mkdir" -p
OUTPUT_DIR := output

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
## NEWTONS METHOD
#################################################################
NEWTONS_METHOD_NAME := newtons_method
NEWTONS_METHOD_DIR := newtons_method
NEWTONS_METHOD_OUTPUT_DIR := $(OUTPUT_DIR)/$(NEWTONS_METHOD_NAME)
NEWTONS_METHOD_LIB := $(NEWTONS_METHOD_OUTPUT_DIR)/lib$(NEWTONS_METHOD_NAME).a

NEWTONS_METHOD_SRC := $(wildcard $(NEWTONS_METHOD_DIR)/*.cpp)
NEWTONS_METHOD_OBJ := $(patsubst %,$(NEWTONS_METHOD_OUTPUT_DIR)/%, $(notdir $(NEWTONS_METHOD_SRC:.cpp=.o)))

#################################################################
## TEST CASES
#################################################################
TEST_CASE_NAME := test_case
TEST_CASE_DIR := test_case
TEST_CASE_OUTPUT_DIR := $(OUTPUT_DIR)/$(TEST_CASE_DIR)

TEST_CASE_SRC := $(wildcard $(TEST_CASE_DIR)/*.cpp)
TEST_CASE_OBJ := $(patsubst %,$(TEST_CASE_OUTPUT_DIR)/%, $(notdir $(TEST_CASE_SRC:.cpp=.o)))

PROGRAM := $(OUTPUT_DIR)/main

LFLAGS := -l$(LINEAR_ALGEBRA_NAME) -L$(LINEAR_ALGEBRA_OUTPUT_DIR) -l$(LINE_SEARCH_NAME) -L$(LINE_SEARCH_OUTPUT_DIR) -l$(NEWTONS_METHOD_NAME) -L$(NEWTONS_METHOD_OUTPUT_DIR)

all : directories $(PROGRAM)

$(PROGRAM) : $(TEST_CASE_OBJ) $(LINEAR_ALGEBRA_LIB) $(LINE_SEARCH_LIB) $(NEWTONS_METHOD_LIB)
	$(CC) -o $@ $^ $(LFLAGS)
	
directories : $(OUTPUT_DIR) $(LINEAR_ALGEBRA_OUTPUT_DIR) $(LINE_SEARCH_OUTPUT_DIR) $(NEWTONS_METHOD_OUTPUT_DIR) $(TEST_CASE_OUTPUT_DIR)

$(TEST_CASE_OUTPUT_DIR)/%.o : $(TEST_CASE_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
	
#################################################################
## create directories
#################################################################
$(OUTPUT_DIR) : 
	$(MKDIR_P) $(OUTPUT_DIR)
	
$(LINEAR_ALGEBRA_OUTPUT_DIR) :
	$(MKDIR_P) $(LINEAR_ALGEBRA_OUTPUT_DIR)
	
$(LINE_SEARCH_OUTPUT_DIR) :
	$(MKDIR_P) $(LINE_SEARCH_OUTPUT_DIR)
	
$(NEWTONS_METHOD_OUTPUT_DIR) : 
	$(MKDIR_P) $(NEWTONS_METHOD_OUTPUT_DIR)
	
$(TEST_CASE_OUTPUT_DIR) :
	$(MKDIR_P) $(TEST_CASE_OUTPUT_DIR)
	
#################################################################
##
#################################################################

$(LINEAR_ALGEBRA_LIB) : $(LINEAR_ALGEBRA_OBJ)
	$(AR) cr $@ $^
	
$(LINEAR_ALGEBRA_OUTPUT_DIR)/%.o : $(LINEAR_ALGEBRA_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(LINE_SEARCH_LIB) : $(LINE_SEARCH_OBJ)
	$(AR) cr $@ $^
	
$(LINE_SEARCH_OUTPUT_DIR)/%.o : $(LINE_SEARCH_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(NEWTONS_METHOD_LIB) : $(NEWTONS_METHOD_OBJ)
	$(AR) cr $@ $^
	
$(NEWTONS_METHOD_OUTPUT_DIR)/%.o : $(NEWTONS_METHOD_DIR)/%.cpp	
	$(CC) $(CFLAGS) -c $< -o $@

	
#################################################################
## utility
#################################################################

.PHONY : help clean

clean :
	$(RM) $(TARGET) $(LINE_SEARCH_LIB) $(LINE_SEARCH_OBJ) $(LINEAR_ALGEBRA_LIB) $(LINEAR_ALGEBRA_OBJ) $(NEWTONS_METHOD_LIB) $(NEWTONS_METHOD_OBJ) $(TEST_CASE_OBJ)

help :
	@echo "LINE_SEARCH_SRC:=$(LINE_SEARCH_SRC)"
	@echo "LINE_SEARCH_OBJ_DIR:=$(LINE_SEARCH_OUTPUT_DIR)"
	@echo "LINE_SEARCH_OBJ:=$(LINE_SEARCH_OBJ)"
	@echo "LINE_SEARCH_LIB:=$(LINE_SEARCH_LIB)"
	