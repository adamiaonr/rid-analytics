# Set the main compiler here. Options e.g.: 'gcc', 'g++'
CC := g++

# Special directories
SRCDIR := src
BUILDDIR := build
TARGET_SANITY := sanity-check
TARGET_SENSITIVITY := sensitivity-analysis

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS_SANITY := $(filter-out $(BUILDDIR)/$(TARGET_SENSITIVITY).o, $(OBJECTS))
OBJECTS_SENSITIVITY := $(filter-out $(BUILDDIR)/$(TARGET_SANITY).o, $(OBJECTS))

# use -ggdb for GNU debugger and -std=c++11 for tree.hh
#CFLAGS := -g -ggdb -Wall -std=c++11 -stdlib=libc++
CFLAGS := -g -ggdb -O3 -Wall -std=c++11 -stdlib=libc++

LIB := -lm
INC := -Iinclude

all: $(TARGET_SANITY) $(TARGET_SENSITIVITY)
	@echo " Doing nothing..."

$(TARGET_SANITY): $(OBJECTS_SANITY)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET_SANITY) $(LIB)"; $(CC) $^ -o $(TARGET_SANITY) $(LIB)

$(TARGET_SENSITIVITY): $(OBJECTS_SENSITIVITY)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET_SENSITIVITY) $(LIB)"; $(CC) $^ -o $(TARGET_SENSITIVITY) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	$(RM) -r $(BUILDDIR) $(TARGET_SANITY) $(TARGET_SENSITIVITY) *~

.PHONY: clean
