################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../test/eng_heat_model.c \
../test/reconstruction_test.c \
../test/synthetic_heat_test.c 

OBJS += \
./test/eng_heat_model.o \
./test/reconstruction_test.o \
./test/synthetic_heat_test.o 

C_DEPS += \
./test/eng_heat_model.d \
./test/reconstruction_test.d \
./test/synthetic_heat_test.d 


# Each subdirectory must supply rules for building sources it contributes
test/%.o: ../test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -std=c99 -std=gnu99 -DNDEBUG -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


