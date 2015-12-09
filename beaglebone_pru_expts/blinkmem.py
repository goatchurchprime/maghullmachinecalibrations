from mmap import mmap
import time, struct

# given in p182 table 2-3 of the AM335x Sitara processor reference manual 
GPIO1_offset = 0x4804c000
GPIO1_size = 0x4804cfff-GPIO1_offset
USR1 = 1<<28   # for pin P9_28

# values from p4877 section 25.4.1
GPIO_OUTPUTENABLE = 0x134  
GPIO_SETDATAOUT = 0x194
GPIO_CLEARDATAOUT = 0x190
GPIO_DATAIN = 0x138

f = open("/dev/mem", "r+b" )
mem = mmap(f.fileno(), GPIO1_size, offset=GPIO1_offset)

reg = struct.unpack("<L", mem[GPIO_OUTPUTENABLE:GPIO_OUTPUTENABLE+4])[0]


if False:   # read case
  mem[GPIO_OUTPUTENABLE:GPIO_OUTPUTENABLE+4] = struct.pack("<L", reg | USR1)
  while(True):
    reg = struct.unpack("<L", mem[GPIO_DATAIN:GPIO_DATAIN+4])[0]
    print(reg & USR1)
    time.sleep(0.1)

else:  # blink case
  mem[GPIO_OUTPUTENABLE:GPIO_OUTPUTENABLE+4] = struct.pack("<L", reg & ~USR1)
  while True:
    mem[GPIO_SETDATAOUT:GPIO_SETDATAOUT+4] = struct.pack("<L", USR1)
    #time.sleep(0.2)
    mem[GPIO_CLEARDATAOUT:GPIO_CLEARDATAOUT+4] = struct.pack("<L", USR1)
    #time.sleep(0.2)

