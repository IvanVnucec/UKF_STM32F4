.SECONDARY:

PROJECT = renode-example
BUILD_DIR = build
Q ?= @

CC = arm-none-eabi-gcc
LD = arm-none-eabi-ld
OCPY = arm-none-eabi-objcopy
MKDIR = mkdir
GIT=git
ECHO=@echo
CAT=cat
PYTHON ?= python

GIT_SHA := \"$(shell $(GIT) rev-parse --short HEAD)\"


SRCS_APP = \
  src/app.c \
  src/clock.c \
  src/gpio.c \
  src/usart.c \
  src/syscalls.c \
  src/ukfLib/lib/mtxLib.c \
  src/ukfLib/lib/ukfLib.c \
  src/ukfLib/cfg/ukfCfg.c \
  src/ukfLib/cfg/ukfCfg1.c \
  src/ukfLib/cfg/ukfCfg2.c


INCLUDES = \
  src \
  src/ukfLib/cfg \
  src/ukfLib/lib

RENODE_REPO = renode

DEFINES += \
	STM32F4 \
	GIT_SHA=$(GIT_SHA) \

CFLAGS += \
  -mcpu=cortex-m4 \
  -mfloat-abi=hard \
  -mfpu=fpv4-sp-d16 \
  -mthumb \
  -Wall \
  -Werror \
  -std=gnu11 \
  -O0 \
  -g \
  -ffunction-sections \
  -fdata-sections

LDFLAGS += \
  -static \
  -nostartfiles \
  -specs=nosys.specs \
  -Wl,--start-group -lc -lgcc -lnosys -Wl,--end-group \
  -Wl,-Map=$(BUILD_DIR)/$(PROJECT).map \
  -u _printf_float

# remove "-lm" if you want to exclude <math.h> library
MATH_LIB = -lm

LDFLAGS_APP = $(LDFLAGS) -T stm32f429i-discovery.ld

OPENCM3_PATH = ./libopencm3
OPENCM3_INCLUDES = $(OPENCM3_PATH)/include
OPENCM3_LIB = $(OPENCM3_PATH)/lib/libopencm3_stm32f4.a

INCLUDES += $(OPENCM3_INCLUDES)
CFLAGS += $(foreach i,$(INCLUDES),-I$(i))
CFLAGS += $(foreach d,$(DEFINES),-D$(d))
LDSCRIPT = stm32f429i-discovery.ld

.PHONY: all test_docker test_local renode
all: $(BUILD_DIR)/$(PROJECT).elf

$(BUILD_DIR)/$(PROJECT).elf: $(SRCS_APP) $(OPENCM3_LIB)
	$(ECHO) "  LD		$@"
	$(Q)$(MKDIR) -p $(BUILD_DIR)
	$(Q)$(CC) $(CFLAGS) $(LDFLAGS_APP) $^ -o $@ $(MATH_LIB)

$(RENODE_REPO):
	$(ECHO) "renode not found, cloning it..."
	$(Q)$(GIT) clone https://github.com/renode/renode.git 2>1

$(OPENCM3_LIB):
	$(ECHO) "Building libopencm3"
	$(Q)$(MAKE) -s -C $(OPENCM3_PATH) TARGETS=stm32/f4

test_docker:
	./docker-test.sh

test_local: $(RENODE_REPO)
	./run_tests.sh

start_renode:
	./start.sh

.PHONY: clean
clean:
	$(ECHO) "  CLEAN		rm -rf $(BUILD_DIR)"
	$(Q)rm -rf $(BUILD_DIR)
	$(Q)make -C $(OPENCM3_PATH) TARGETS=stm32/f4 clean