HIPC=hipcc
radius=4
vec=2

# Directories
SRC_DIR=src
INC_DIR=include
BUILD_DIR=build

# GPU Architecture: can be specified via command line (e.g., make ARCH=gfx1100)
# If not specified, auto-detect from available GPU using rocminfo
# Falls back to gfx90a (MI250X) if detection fails
ifndef ARCH
  DETECTED_ARCH := $(shell rocminfo 2>/dev/null | grep -o -m1 'gfx[0-9a-z]*')
  ifeq ($(DETECTED_ARCH),)
    ARCH := gfx90a
    $(info Auto-detection failed, using default: $(ARCH))
  else
    ARCH := $(DETECTED_ARCH)
    $(info Auto-detected GPU architecture: $(ARCH))
  endif
else
  $(info Using specified GPU architecture: $(ARCH))
endif

# Build configuration suffix
CONFIG=R$(radius)_vec$(vec)_$(ARCH)

# Build subdirectory for this configuration
BUILD_SUBDIR=$(BUILD_DIR)/$(CONFIG)

# Optional: save temp files (assembly, IR) for debugging
# Usage: make SAVE_TEMPS=1
ifdef SAVE_TEMPS
  TEMPS_FLAG=--save-temps
else
  TEMPS_FLAG=
endif

# Compiler flags
HIPFLAGS=-O3 -g --offload-arch=$(ARCH) -std=c++17 -I$(CURDIR)/$(INC_DIR) -DRADIUS=$(radius) -DVEC_EXP=$(vec) $(TEMPS_FLAG)

# Output binary (includes architecture in name)
TARGET=baseline_$(CONFIG).x

# Object files
OBJS=$(BUILD_SUBDIR)/main.o $(BUILD_SUBDIR)/baseline_kernel.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(HIPC) $(HIPFLAGS) $(OBJS) -o $@

# Compile from within the build subdirectory so temp files end up there
$(BUILD_SUBDIR)/main.o: $(SRC_DIR)/main.cpp | $(BUILD_SUBDIR)
	cd $(BUILD_SUBDIR) && $(HIPC) $(HIPFLAGS) -c $(CURDIR)/$< -o main.o

$(BUILD_SUBDIR)/baseline_kernel.o: $(SRC_DIR)/baseline_kernel.hip | $(BUILD_SUBDIR)
	cd $(BUILD_SUBDIR) && $(HIPC) $(HIPFLAGS) -c $(CURDIR)/$< -o baseline_kernel.o

$(BUILD_SUBDIR):
	mkdir -p $(BUILD_SUBDIR)

clean:
	rm -rf $(BUILD_DIR) *.x

clean-config:
	rm -rf $(BUILD_SUBDIR) $(TARGET)

.PHONY: all clean clean-config
