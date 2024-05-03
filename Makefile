CXX      := g++
CXXFLAGS := -pthread -O3
# REPR_FLAGS:= -D __FLOAT_MP__ -D __CSTATE_VCHAR__
# REPR_FLAGS:= -D __FLOAT_MP__ 
# REPR_FLAGS:= -D __CSTATE_VCHAR__
LDFLAGS  := -Llocal/lib -lgmp -lgmpxx
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps
SHELL	 := /bin/bash
TARGET   := Feynman
INCLUDE  := -IFeynman_MCSimulator/simulators/ -IFeynman_MCSimulator/ -IFeynman_MCSimulator/pcg_random/ -Ilocal/include
SRC      :=                      \
   $(wildcard Feynman_MCSimulator/simulators/*.cpp) \
   $(wildcard Feynman_MCSimulator/*.cpp)         \

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES \
         := $(OBJECTS:.o=.d)

all:	build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	module load gcc-11.3 ; \
	$(CXX) $(CXXFLAGS) $(REPR_FLAGS) $(INCLUDE) -c $< -MMD -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	@cp Feynman_MCSimulator/*.sh $(APP_DIR)
	module load gcc-11.3 ; \
	$(CXX) $(CXXFLAGS) $(REPR_FLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug release info

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"

