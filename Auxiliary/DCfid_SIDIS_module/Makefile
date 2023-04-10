CXX      := -c++
CXXFLAGS := -pedantic-errors -Wall -Wextra -Werror -std=c++17
LDFLAGS  := -L/usr/lib -lstdc++ -lm -L/opt/homebrew/Cellar/root/6.26.06_2/lib/root
ROOTCFLAGS:= `root-config --cflags`
ROOTLFLAGS:= `root-config --glibs`
CFLAGS   := -c -g -Wall `root-config --cflags`
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps
TARGET   := program
INCLUDE  := -Iinclude/  -I/opt/homebrew/Cellar/root/6.26.06_2/include/root
SRC      :=                      \
   $(wildcard src/*.cpp)         \

OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.so)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.so: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(ROOTCFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) $(ROOTLFLAGS) -shared -Wl -o $(APP_DIR)/$(TARGET) $(OBJECTS)

.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*



