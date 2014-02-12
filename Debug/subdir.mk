################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../bntseq.c \
../utils.c 

CC_SRCS += \
../lea_index.cc \
../lea_index_reglen.cc \
../lea_utility.cc \
../main.cc 

OBJS += \
./bntseq.o \
./lea_index.o \
./lea_index_reglen.o \
./lea_utility.o \
./main.o \
./utils.o 

C_DEPS += \
./bntseq.d \
./utils.d 

CC_DEPS += \
./lea_index.d \
./lea_index_reglen.d \
./lea_utility.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


