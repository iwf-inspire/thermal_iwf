################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../structures/particle.c \
../structures/singleton_geometry.c \
../structures/singleton_physics.c 

OBJS += \
./structures/particle.o \
./structures/singleton_geometry.o \
./structures/singleton_physics.o 

C_DEPS += \
./structures/particle.d \
./structures/singleton_geometry.d \
./structures/singleton_physics.d 


# Each subdirectory must supply rules for building sources it contributes
structures/%.o: ../structures/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -std=c99 -std=gnu99 -DNDEBUG -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


