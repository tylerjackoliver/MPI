# Helpful macros
TARGET_EXEC ?= edge_to_image
OBJ_DIR ?= ./objects
SRC_DIRS ?= ./source

# Macros for CC
FC = mpif90
CFLAGS = -pedantic -Wall

# Directories
SRCS := $(shell find $(SRC_DIRS) -name *.f90)
OBJS := $(SRCS:%=$(OBJ_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

# Include Macros
INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

# Compile the source
$(TARGET_EXEC): $(OBJS)
	$(FC) $(OBJS) $(CFLAGS) $(LIBS) -o $@

# Create object files
$(OBJ_DIR)/%.f90.o: %.f90
	$(MKDIR_P) $(dir $@)
	$(FC) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) -r $(OBJ_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p
