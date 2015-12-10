
.origin 0
.entrypoint START

#define PRU0_ARM_INTERRUPT 19

#define GPIO1 				0x4804c000		// The adress of the GPIO1 
#define GPIO_CLEARDATAOUT 	0x190
#define GPIO_SETDATAOUT 	0x194
#define GPIO_DATAIN         0x138

START:
    LBCO r0, C4, 4, 4					// Load Bytes Constant Offset (?)
    CLR  r0, r0, 4						// Clear bit 4 in reg 0
    SBCO r0, C4, 4, 4					// Store Bytes Constant Offset

    MOV r1, 0xF00000
    MOV r2, 1<<28
BLINK:
    MOV r3, GPIO1 | GPIO_SETDATAOUT
    MOV r0, 8
    SBBO r2, r3, 0, 4
DELAY:
    SUB r0, r0, 1
    QBNE DELAY, r0, 0
//    ADD r0, r0, 1  // [1]
//    ADD r0, r0, 1  // [1]

    MOV r3, GPIO1 | GPIO_DATAIN
    LBBO r2, r3, 0, 4
    MOV r2, 1<<28


    MOV r3, GPIO1 | GPIO_CLEARDATAOUT
    SBBO r2, r3, 0, 4

    MOV r0, 4
DELAY2:
    SUB r0, r0, 1
    QBNE DELAY2, r0, 0
//    ADD r0, r0, 1  // [2]

    SUB r1, r1, 1
    QBNE BLINK, r1, 0

    MOV R31.b0, PRU0_ARM_INTERRUPT+16   // Send notification to Host for program completion
HALT

