################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../neighboring/cell_list.c 

OBJS += \
./neighboring/cell_list.o 

C_DEPS += \
./neighboring/cell_list.d 


# Each subdirectory must supply rules for building sources it contributes
neighboring/%.o: ../neighboring/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -std=gnu99 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


