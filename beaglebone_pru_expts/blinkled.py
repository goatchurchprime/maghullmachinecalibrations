import subprocess, os
p = subprocess.Popen("/usr/bin/pasm -b blinkled.p", shell=True)
pid, sts = os.waitpid(p.pid, 0)

# git clone https://bitbucket.org/intelligentagent/pypruss.git 
# cd pypruss
# python setup.py install
# export LD_LIBRARY_PATH=/usr/local/lib
# do this on the command line at start up if the device needs to be enbabled
#    echo BB-BONE-PRU-01 > /sys/devices/bone_capemgr.9/slots

import pypruss

pypruss.modprobe()                                  # This only has to be called once pr boot
pypruss.init()                                      # Init the PRU
pypruss.open(0)                                     # Open PRU event 0 which is PRU0_ARM_INTERRUPT
pypruss.pruintc_init()                              # Init the interrupt controller
pypruss.exec_program(0, "./blinkled.bin")           # Load firmware "blinkled.bin" on PRU 0
pypruss.wait_for_event(0)                           # Wait for event 0 which is connected to PRU0_ARM_INTERRUPT
pypruss.clear_event(0)                              # Clear the event
pypruss.pru_disable(0)                              # Disable PRU 0, this is already done by the firmware
pypruss.exit()                                      # Exit, don't know what this does.

