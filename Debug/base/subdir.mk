################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../base/kernels.c \
../base/utils.c 

OBJS += \
./base/kernels.o \
./base/utils.o 

C_DEPS += \
./base/kernels.d \
./base/utils.d 


# Each subdirectory must supply rules for building sources it contributes
base/%.o: ../base/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -std=gnu99 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


