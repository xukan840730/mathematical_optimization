
# TO BE CONTINUED

TARGET = main
SOURCE = $(wildcard *.cpp)
OBJECT = $(SOURCE:.cpp=.o)

help :
	@echo  "target: $(TARGET)"
