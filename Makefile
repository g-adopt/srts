# Main Makefile for TOMOFILT
# Modernized build system for tomographic filtering tools
# Author: Auto-generated modernization
# Date: September 2025

# Compiler settings
FF = gfortran
CC = gcc
FFLAGS = -ffixed-line-length-none -fno-automatic -O3
CFLAGS = -std=c90 -O3

# Directories
ROOT_DIR = $(shell pwd)
LIB_DIR = $(ROOT_DIR)/lib
SRC_DIR = $(ROOT_DIR)/src
BIN_DIR = $(ROOT_DIR)/bin

# Library directories
LIBGEO_DIR = $(LIB_DIR)/libgeo
LIBHJVH_DIR = $(LIB_DIR)/libHJVH
LIBS20_DIR = $(LIB_DIR)/libS20
LIBSPHFUNC_DIR = $(LIB_DIR)/libSphFunc

# Library files
LIBGEO = $(LIB_DIR)/libgeo.a
LIBHJVH = $(LIB_DIR)/libHJVH.a
LIBS20 = $(LIB_DIR)/libS20.a
LIBSPHFUNC = $(LIB_DIR)/libSphFunc.a

# All libraries (in dependency order)
LIBRARIES = $(LIBGEO) $(LIBHJVH) $(LIBS20) $(LIBSPHFUNC)

# Executables that will be built
EXECUTABLES = mkexpmatxy invexpandxy mk3d_res mk3d_res_ESEP sphadd spt2sph \
              sphexp linegrep depmap raw2xyz-eqdist sph2v_input pwrspecsph correlatorsph

# Default target
.PHONY: all
all: check-env create-bin libraries executables
	@echo "Build completed successfully!"
	@echo "Executables are available in: $(BIN_DIR)"

# Check environment
.PHONY: check-env
check-env:
	@echo "Checking build environment..."
	@if [ -z "$(TOMOFILT)" ]; then \
		echo "WARNING: TOMOFILT environment variable not set."; \
		echo "Consider setting: export TOMOFILT=$(ROOT_DIR)"; \
	else \
		echo "TOMOFILT environment variable set to: $(TOMOFILT)"; \
	fi
	@which $(FF) > /dev/null || (echo "ERROR: $(FF) not found" && exit 1)
	@which $(CC) > /dev/null || (echo "ERROR: $(CC) not found" && exit 1)
	@echo "Environment check passed."

# Create bin directory if it doesn't exist
.PHONY: create-bin
create-bin:
	@mkdir -p $(BIN_DIR)

# Build all libraries
.PHONY: libraries
libraries: $(LIBRARIES)
	@echo "All libraries built successfully."

# Individual library targets
$(LIBGEO): $(LIBGEO_DIR)/Makefile $(wildcard $(LIBGEO_DIR)/*.f $(LIBGEO_DIR)/*.F)
	@echo "Building libgeo..."
	@if [ ! -d "$(LIBGEO_DIR)" ]; then \
		echo "ERROR: libgeo directory $(LIBGEO_DIR) not found"; \
		exit 1; \
	fi
	@cd $(LIBGEO_DIR) && $(MAKE) -f Makefile ../libgeo.a
	@echo "libgeo built successfully."

$(LIBHJVH): $(LIBHJVH_DIR)/Makefile $(wildcard $(LIBHJVH_DIR)/*.f $(LIBHJVH_DIR)/*.F)
	@echo "Building libHJVH..."
	@if [ ! -d "$(LIBHJVH_DIR)" ]; then \
		echo "ERROR: libHJVH directory $(LIBHJVH_DIR) not found"; \
		exit 1; \
	fi
	@cd $(LIBHJVH_DIR) && $(MAKE) -f Makefile ../libHJVH.a
	@echo "libHJVH built successfully."

$(LIBS20): $(LIBS20_DIR)/Makefile $(wildcard $(LIBS20_DIR)/*.f $(LIBS20_DIR)/*.F)
	@echo "Building libS20..."
	@if [ ! -d "$(LIBS20_DIR)" ]; then \
		echo "ERROR: libS20 directory $(LIBS20_DIR) not found"; \
		exit 1; \
	fi
	@cd $(LIBS20_DIR) && $(MAKE) -f Makefile ../libS20.a
	@echo "libS20 built successfully."

$(LIBSPHFUNC): $(LIBSPHFUNC_DIR)/Makefile $(wildcard $(LIBSPHFUNC_DIR)/*.f $(LIBSPHFUNC_DIR)/*.F)
	@echo "Building libSphFunc..."
	@if [ ! -d "$(LIBSPHFUNC_DIR)" ]; then \
		echo "ERROR: libSphFunc directory $(LIBSPHFUNC_DIR) not found"; \
		exit 1; \
	fi
	@cd $(LIBSPHFUNC_DIR) && $(MAKE) -f Makefile ../libSphFunc.a
	@echo "libSphFunc built successfully."

# Build executables (depends on libraries)
.PHONY: executables
executables: $(LIBRARIES) create-bin
	@echo "Building executables..."
	@$(MAKE) -C $(SRC_DIR) all
	@echo "All executables built successfully."

# Individual executable targets (for selective building)
.PHONY: $(EXECUTABLES)
$(EXECUTABLES): $(LIBRARIES) create-bin
	@$(MAKE) -C $(SRC_DIR) $@

# Individual library build targets (for selective building)
.PHONY: libgeo libHJVH libS20 libSphFunc
libgeo: $(LIBGEO)
libHJVH: $(LIBHJVH)
libS20: $(LIBS20)
libSphFunc: $(LIBSPHFUNC)

# Clean targets
.PHONY: clean
clean: clean-libs clean-src clean-bin
	@echo "Cleaned all build artifacts."

.PHONY: clean-libs
clean-libs:
	@echo "Cleaning libraries..."
	@rm -f $(LIB_DIR)/*.a
	@for dir in $(LIBGEO_DIR) $(LIBHJVH_DIR) $(LIBS20_DIR) $(LIBSPHFUNC_DIR); do \
		if [ -d "$$dir" ]; then \
			echo "Cleaning $$dir..."; \
			find "$$dir" -name "*.o" -delete 2>/dev/null || true; \
			if [ -f "$$dir/Makefile" ]; then \
				cd "$$dir" && $(MAKE) clean 2>/dev/null || true; \
			fi; \
		fi; \
	done

.PHONY: clean-src
clean-src:
	@echo "Cleaning source directory..."
	@find $(SRC_DIR) -name "*.o" -delete 2>/dev/null || true

.PHONY: clean-bin
clean-bin:
	@echo "Cleaning executables..."
	@rm -rf $(BIN_DIR)

# Deep clean (including backup files)
.PHONY: distclean
distclean: clean
	@echo "Performing deep clean..."
	@find . -name "*~" -delete 2>/dev/null || true
	@find . -name "*.f~" -delete 2>/dev/null || true
	@find . -name "*.F~" -delete 2>/dev/null || true
	@find . -name "*.h~" -delete 2>/dev/null || true

# Rebuild everything from scratch
.PHONY: rebuild
rebuild: clean all

# Install target (copies to system location if desired)
.PHONY: install
install: all
	@if [ -n "$(PREFIX)" ]; then \
		echo "Installing to $(PREFIX)/bin..."; \
		mkdir -p $(PREFIX)/bin; \
		cp $(BIN_DIR)/* $(PREFIX)/bin/; \
		echo "Installation completed."; \
	else \
		echo "To install, run: make install PREFIX=/path/to/install"; \
	fi

# Test build (just compile, don't link)
.PHONY: test-compile
test-compile:
	@echo "Testing compilation without linking..."
	@for dir in $(LIBGEO_DIR) $(LIBHJVH_DIR) $(LIBS20_DIR) $(LIBSPHFUNC_DIR); do \
		if [ -d "$$dir" ]; then \
			echo "Testing compilation in $$dir..."; \
			$(MAKE) -C "$$dir" -n > /dev/null; \
		fi; \
	done
	@echo "Compilation test passed."

# Parallel build support
.PHONY: parallel
parallel:
	@echo "Building with parallel jobs..."
	@$(MAKE) -j$(shell nproc 2>/dev/null || echo 4) all

# Help target
.PHONY: help
help:
	@echo "TOMOFILT Build System"
	@echo "===================="
	@echo ""
	@echo "Main targets:"
	@echo "  all          - Build everything (default)"
	@echo "  libraries    - Build only the libraries"
	@echo "  executables  - Build only the executables (requires libraries)"
	@echo "  clean        - Remove all build artifacts"
	@echo "  rebuild      - Clean and build everything"
	@echo "  distclean    - Deep clean including backup files"
	@echo "  parallel     - Build with parallel jobs"
	@echo ""
	@echo "Installation:"
	@echo "  install PREFIX=/path - Install executables to specified path"
	@echo ""
	@echo "Individual libraries:"
	@echo "  libgeo       - Build $(LIBGEO)"
	@echo "  libHJVH      - Build $(LIBHJVH)"
	@echo "  libS20       - Build $(LIBS20)"
	@echo "  libSphFunc   - Build $(LIBSPHFUNC)"
	@echo ""
	@echo "Individual executables:"
	@echo "  $(EXECUTABLES)"
	@echo ""
	@echo "Environment:"
	@echo "  Set TOMOFILT environment variable for proper operation"
	@echo "  Current TOMOFILT: $(TOMOFILT)"

# Check for common issues
.PHONY: doctor
doctor: check-env
	@echo "Running build system diagnostics..."
	@echo "Checking library directories..."
	@for dir in $(LIBGEO_DIR) $(LIBHJVH_DIR) $(LIBS20_DIR) $(LIBSPHFUNC_DIR); do \
		if [ ! -d "$$dir" ]; then \
			echo "ERROR: Directory $$dir not found"; \
			exit 1; \
		fi; \
		if [ ! -f "$$dir/Makefile" ]; then \
			echo "ERROR: Makefile not found in $$dir"; \
			exit 1; \
		fi; \
	done
	@echo "Checking source directory..."
	@if [ ! -d "$(SRC_DIR)" ]; then \
		echo "ERROR: Source directory $(SRC_DIR) not found"; \
		exit 1; \
	fi
	@if [ ! -f "$(SRC_DIR)/Makefile" ]; then \
		echo "ERROR: Makefile not found in $(SRC_DIR)"; \
		exit 1; \
	fi
	@echo "All diagnostics passed!"

# Show build status
.PHONY: status
status:
	@echo "Build Status Report"
	@echo "=================="
	@echo "Libraries:"
	@for lib in $(LIBRARIES); do \
		if [ -f "$$lib" ]; then \
			echo "  ✓ $$lib (built)"; \
		else \
			echo "  ✗ $$lib (not built)"; \
		fi; \
	done
	@echo ""
	@echo "Executables:"
	@for exe in $(EXECUTABLES); do \
		if [ -f "$(BIN_DIR)/$$exe" ]; then \
			echo "  ✓ $$exe (built)"; \
		else \
			echo "  ✗ $$exe (not built)"; \
		fi; \
	done
