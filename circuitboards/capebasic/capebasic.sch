EESchema Schematic File Version 2
LIBS:power
LIBS:device
LIBS:transistors
LIBS:conn
LIBS:linear
LIBS:regul
LIBS:74xx
LIBS:cmos4000
LIBS:adc-dac
LIBS:memory
LIBS:xilinx
LIBS:special
LIBS:microcontrollers
LIBS:dsp
LIBS:microchip
LIBS:analog_switches
LIBS:motorola
LIBS:texas
LIBS:intel
LIBS:audio
LIBS:interface
LIBS:digital-audio
LIBS:philips
LIBS:display
LIBS:cypress
LIBS:siliconi
LIBS:opto
LIBS:atmel
LIBS:contrib
LIBS:valves
LIBS:beaglebone
LIBS:capebasic-cache
EELAYER 27 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title ""
Date "4 mar 2016"
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L BEAGLEBONE U1
U 5 1 5649CF8A
P 6700 6500
F 0 "U1" H 6750 6450 60  0000 C CNN
F 1 "BEAGLEBONE" H 6950 6350 60  0000 C CNN
F 2 "" H 6700 6500 60  0000 C CNN
F 3 "" H 6700 6500 60  0000 C CNN
	5    6700 6500
	1    0    0    -1  
$EndComp
$Comp
L SW_PUSH SW1
U 1 1 5649D004
P 4400 3900
F 0 "SW1" H 4550 4010 50  0000 C CNN
F 1 "SW_PUSH" H 4400 3820 50  0000 C CNN
F 2 "~" H 4400 3900 60  0000 C CNN
F 3 "~" H 4400 3900 60  0000 C CNN
	1    4400 3900
	-1   0    0    1   
$EndComp
$Comp
L +3.3V #PWR01
U 1 1 5649D06D
P 4500 4250
F 0 "#PWR01" H 4500 4210 30  0001 C CNN
F 1 "+3.3V" H 4500 4360 30  0000 C CNN
F 2 "" H 4500 4250 60  0000 C CNN
F 3 "" H 4500 4250 60  0000 C CNN
	1    4500 4250
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR02
U 1 1 5649D083
P 4750 4300
F 0 "#PWR02" H 4750 4300 30  0001 C CNN
F 1 "GND" H 4750 4230 30  0001 C CNN
F 2 "" H 4750 4300 60  0000 C CNN
F 3 "" H 4750 4300 60  0000 C CNN
	1    4750 4300
	1    0    0    -1  
$EndComp
$Comp
L +3.3V #PWR03
U 1 1 5649DB80
P 3150 5900
F 0 "#PWR03" H 3150 5860 30  0001 C CNN
F 1 "+3.3V" H 3150 6010 30  0000 C CNN
F 2 "" H 3150 5900 60  0000 C CNN
F 3 "" H 3150 5900 60  0000 C CNN
	1    3150 5900
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR04
U 1 1 5649E7C5
P 3100 6100
F 0 "#PWR04" H 3100 6100 30  0001 C CNN
F 1 "GND" H 3100 6030 30  0001 C CNN
F 2 "" H 3100 6100 60  0000 C CNN
F 3 "" H 3100 6100 60  0000 C CNN
	1    3100 6100
	1    0    0    -1  
$EndComp
$Comp
L C C1
U 1 1 564A03A1
P 3000 5800
F 0 "C1" H 3000 5900 40  0000 L CNN
F 1 "C" H 3006 5715 40  0000 L CNN
F 2 "~" H 3038 5650 30  0000 C CNN
F 3 "~" H 3000 5800 60  0000 C CNN
	1    3000 5800
	1    0    0    -1  
$EndComp
$Comp
L 7404 U3
U 1 1 564A07FE
P 9050 1700
F 0 "U3" H 8550 2400 60  0000 L CNN
F 1 "7404" H 8550 1000 60  0000 L CNN
F 2 "" H 9050 2150 60  0000 C CNN
F 3 "" H 9050 2150 60  0000 C CNN
	1    9050 1700
	1    0    0    -1  
$EndComp
$Comp
L C C2
U 1 1 564A0B24
P 9500 1000
F 0 "C2" H 9500 1100 40  0000 L CNN
F 1 "C" H 9506 915 40  0000 L CNN
F 2 "~" H 9538 850 30  0000 C CNN
F 3 "~" H 9500 1000 60  0000 C CNN
	1    9500 1000
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR05
U 1 1 564A0C94
P 9400 850
F 0 "#PWR05" H 9400 850 30  0001 C CNN
F 1 "GND" H 9400 780 30  0001 C CNN
F 2 "" H 9400 850 60  0000 C CNN
F 3 "" H 9400 850 60  0000 C CNN
	1    9400 850 
	1    0    0    -1  
$EndComp
$Comp
L +5V #PWR06
U 1 1 564A0CA3
P 9700 1050
F 0 "#PWR06" H 9700 1140 20  0001 C CNN
F 1 "+5V" H 9700 1140 30  0000 C CNN
F 2 "" H 9700 1050 60  0000 C CNN
F 3 "" H 9700 1050 60  0000 C CNN
	1    9700 1050
	1    0    0    -1  
$EndComp
$Comp
L RR8 RR1
U 1 1 5649DDD6
P 650 6450
F 0 "RR1" H 700 7000 70  0000 C CNN
F 1 "RR8" V 680 6450 70  0000 C CNN
F 2 "~" H 650 6450 60  0000 C CNN
F 3 "~" H 650 6450 60  0000 C CNN
	1    650  6450
	-1   0    0    1   
$EndComp
$Comp
L CONN_2 P1
U 1 1 564A21A5
P 1300 950
F 0 "P1" V 1250 950 40  0000 C CNN
F 1 "CONN_2" V 1350 950 40  0000 C CNN
F 2 "" H 1300 950 60  0000 C CNN
F 3 "" H 1300 950 60  0000 C CNN
	1    1300 950 
	-1   0    0    1   
$EndComp
$Comp
L +5V #PWR07
U 1 1 564A224B
P 1700 1050
F 0 "#PWR07" H 1700 1140 20  0001 C CNN
F 1 "+5V" H 1700 1140 30  0000 C CNN
F 2 "" H 1700 1050 60  0000 C CNN
F 3 "" H 1700 1050 60  0000 C CNN
	1    1700 1050
	0    1    1    0   
$EndComp
NoConn ~ 4850 4550
NoConn ~ 4850 5850
NoConn ~ 4850 5950
NoConn ~ 4850 6050
NoConn ~ 4850 6150
NoConn ~ 4850 6250
NoConn ~ 4850 6350
NoConn ~ 4850 6450
NoConn ~ 6350 6450
NoConn ~ 6350 6350
NoConn ~ 6350 6150
NoConn ~ 6350 6050
NoConn ~ 6350 5950
NoConn ~ 6350 5850
NoConn ~ 6350 5750
NoConn ~ 6350 5250
NoConn ~ 6350 5150
NoConn ~ 6350 5050
NoConn ~ 6350 4650
NoConn ~ 6350 4550
NoConn ~ 6350 4350
NoConn ~ 6350 4250
Text Label 4200 4850 0    60   ~ 0
HOME BB
Text Label 4200 4950 0    60   ~ 0
PROBE BB
Text Label 4200 5050 0    60   ~ 0
X DRIVER BB
Text Label 4200 5150 0    60   ~ 0
Y DRIVER BB
$Comp
L GND #PWR08
U 1 1 564A17B9
P 4700 3950
F 0 "#PWR08" H 4700 3950 30  0001 C CNN
F 1 "GND" H 4700 3880 30  0001 C CNN
F 2 "" H 4700 3950 60  0000 C CNN
F 3 "" H 4700 3950 60  0000 C CNN
	1    4700 3950
	1    0    0    -1  
$EndComp
Text Label 4200 4750 0    60   ~ 0
ESTOP BB
$Comp
L GND #PWR09
U 1 1 564A496C
P 8250 2150
F 0 "#PWR09" H 8250 2150 30  0001 C CNN
F 1 "GND" H 8250 2080 30  0001 C CNN
F 2 "" H 8250 2150 60  0000 C CNN
F 3 "" H 8250 2150 60  0000 C CNN
	1    8250 2150
	1    0    0    -1  
$EndComp
Text Label 8250 2150 0    60   ~ 0
GND
$Comp
L +5V #PWR010
U 1 1 564A4DD9
P 4300 4250
F 0 "#PWR010" H 4300 4340 20  0001 C CNN
F 1 "+5V" H 4300 4340 30  0000 C CNN
F 2 "" H 4300 4250 60  0000 C CNN
F 3 "" H 4300 4250 60  0000 C CNN
	1    4300 4250
	1    0    0    -1  
$EndComp
$Comp
L BEAGLEBONE U1
U 4 1 5649CF57
P 5200 6500
F 0 "U1" H 5250 6450 60  0000 C CNN
F 1 "BEAGLEBONE" H 5450 6350 60  0000 C CNN
F 2 "" H 5200 6500 60  0000 C CNN
F 3 "" H 5200 6500 60  0000 C CNN
	4    5200 6500
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR011
U 1 1 564A5EFE
P 1750 900
F 0 "#PWR011" H 1750 900 30  0001 C CNN
F 1 "GND" H 1750 830 30  0001 C CNN
F 2 "" H 1750 900 60  0000 C CNN
F 3 "" H 1750 900 60  0000 C CNN
	1    1750 900 
	1    0    0    -1  
$EndComp
$Comp
L +5V #PWR012
U 1 1 564A66F4
P 900 7050
F 0 "#PWR012" H 900 7140 20  0001 C CNN
F 1 "+5V" H 900 7140 30  0000 C CNN
F 2 "" H 900 7050 60  0000 C CNN
F 3 "" H 900 7050 60  0000 C CNN
	1    900  7050
	1    0    0    -1  
$EndComp
Text Label 1200 6100 0    60   ~ 0
ESTOP
Text Label 1200 6200 0    60   ~ 0
HOME
Text Label 1200 6300 0    60   ~ 0
PROBE
Text Label 1200 6400 0    60   ~ 0
X DRIVER
Text Label 1200 6500 0    60   ~ 0
Y DRIVER
Text Label 1200 6600 0    60   ~ 0
Z DRIVER
Text Label 3250 6800 0    60   ~ 0
QEPB BB
Text Label 3250 6900 0    60   ~ 0
QEPA BB
Text Label 7950 1200 0    60   ~ 0
X STEP BB
Text Label 7950 1500 0    60   ~ 0
Y STEP BB
Text Label 7950 1800 0    60   ~ 0
Z STEP BB
Text Label 9550 1350 0    60   ~ 0
X DIR BB
Text Label 9550 1650 0    60   ~ 0
Y DIR BB
Text Label 9550 1950 0    60   ~ 0
Z DIR BB
Text Label 5950 5350 0    60   ~ 0
X STEP BB
Text Label 5950 5650 0    60   ~ 0
X DIR BB
Text Label 5950 5450 0    60   ~ 0
Y STEP BB
Text Label 4400 5650 0    60   ~ 0
Y DIR BB
Text Label 5950 5550 0    60   ~ 0
Z STEP BB
Text Label 4400 5750 0    60   ~ 0
Z DIR BB
Entry Wire Line
	4300 5750 4400 5650
Entry Wire Line
	4300 5850 4400 5750
Entry Wire Line
	5850 5550 5950 5650
Entry Wire Line
	5850 5450 5950 5550
Entry Wire Line
	5850 5350 5950 5450
Entry Wire Line
	5850 5250 5950 5350
Entry Wire Line
	7850 1300 7950 1200
Entry Wire Line
	7850 1600 7950 1500
Entry Wire Line
	7850 1900 7950 1800
Entry Wire Line
	9950 1950 10050 2050
Entry Wire Line
	9950 1650 10050 1750
Entry Wire Line
	9950 1350 10050 1450
Entry Bus Bus
	7850 2550 7950 2450
Entry Wire Line
	7300 1850 7400 1950
Entry Wire Line
	7300 1550 7400 1650
Entry Wire Line
	7300 1250 7400 1350
Entry Wire Line
	10300 2100 10400 2000
Entry Wire Line
	10300 1800 10400 1700
Entry Wire Line
	10300 1500 10400 1400
Entry Bus Bus
	7200 700  7300 800 
Entry Wire Line
	6150 1300 6250 1200
Entry Wire Line
	6150 1400 6250 1300
Entry Wire Line
	6150 2200 6250 2100
Entry Wire Line
	6150 2300 6250 2200
Entry Wire Line
	6150 3050 6250 2950
Entry Wire Line
	6150 3150 6250 3050
Text Label 7950 1350 0    60   ~ 0
X STEP
Text Label 7950 1650 0    60   ~ 0
Y STEP
Text Label 7950 1950 0    60   ~ 0
Z STEP
Text Label 9550 1500 0    60   ~ 0
X DIR
Text Label 9550 1800 0    60   ~ 0
Y DIR
Text Label 9550 2100 0    60   ~ 0
Z DIR
Text Label 5800 1300 0    60   ~ 0
X STEP
Text Label 5800 1400 0    60   ~ 0
X DIR
Text Label 5800 2200 0    60   ~ 0
Y STEP
Text Label 5800 2300 0    60   ~ 0
Y DIR
Text Label 5800 3050 0    60   ~ 0
Z STEP
Text Label 5800 3150 0    60   ~ 0
Z DIR
$Comp
L 74LVC245 U2
U 1 1 5649DB71
P 2300 6450
F 0 "U2" H 1800 7150 60  0000 L CNN
F 1 "74LVC245" H 1800 5750 60  0000 L CNN
F 2 "" H 2300 6350 60  0000 C CNN
F 3 "" H 2300 6350 60  0000 C CNN
	1    2300 6450
	1    0    0    -1  
$EndComp
Entry Wire Line
	3800 4850 3900 4750
Entry Wire Line
	3800 4950 3900 4850
Entry Wire Line
	3800 5050 3900 4950
Entry Wire Line
	3800 5150 3900 5050
Entry Wire Line
	3800 5250 3900 5150
Entry Wire Line
	3800 5350 3900 5250
$Comp
L GND #PWR013
U 1 1 564CA83F
P 1550 7000
F 0 "#PWR013" H 1550 7000 30  0001 C CNN
F 1 "GND" H 1550 6930 30  0001 C CNN
F 2 "" H 1550 7000 60  0000 C CNN
F 3 "" H 1550 7000 60  0000 C CNN
	1    1550 7000
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR014
U 1 1 564CA917
P 2900 5650
F 0 "#PWR014" H 2900 5650 30  0001 C CNN
F 1 "GND" H 2900 5580 30  0001 C CNN
F 2 "" H 2900 5650 60  0000 C CNN
F 3 "" H 2900 5650 60  0000 C CNN
	1    2900 5650
	1    0    0    -1  
$EndComp
Text Label 3250 6200 0    60   ~ 0
ESTOP BB
Text Label 3250 6300 0    60   ~ 0
HOME BB
Text Label 3250 6400 0    60   ~ 0
PROBE BB
Text Label 3250 6500 0    60   ~ 0
X DRIVER BB
Text Label 3250 6600 0    60   ~ 0
Y DRIVER BB
Text Label 3250 6700 0    60   ~ 0
Z DRIVER BB
Text Label 1200 6700 0    60   ~ 0
QEPB
Text Label 1250 2050 0    60   ~ 0
ESTOP
Text Label 1250 2250 0    60   ~ 0
PROBE
Text Label 1250 2350 0    60   ~ 0
X DRIVER
Text Label 1250 2450 0    60   ~ 0
Y DRIVER
Text Label 1250 2550 0    60   ~ 0
Z DRIVER
Text Label 1250 2650 0    60   ~ 0
QEPB
Text Label 1250 2750 0    60   ~ 0
QEPA
Text Label 1250 2150 0    60   ~ 0
HOME
$Comp
L CONN_8 P2
U 1 1 56D985BC
P 850 2400
F 0 "P2" V 800 2400 60  0000 C CNN
F 1 "CONN_8" V 900 2400 60  0000 C CNN
F 2 "" H 850 2400 60  0000 C CNN
F 3 "" H 850 2400 60  0000 C CNN
	1    850  2400
	-1   0    0    -1  
$EndComp
Entry Wire Line
	2700 2050 2800 2150
Entry Wire Line
	2700 2150 2800 2250
Entry Wire Line
	2700 2250 2800 2350
Entry Wire Line
	2700 2350 2800 2450
Entry Wire Line
	2700 2450 2800 2550
Entry Wire Line
	2700 2550 2800 2650
Entry Wire Line
	2700 2650 2800 2750
Entry Wire Line
	2700 2750 2800 2850
$Comp
L +5V #PWR015
U 1 1 56D9882F
P 4950 1200
F 0 "#PWR015" H 4950 1290 20  0001 C CNN
F 1 "+5V" H 4950 1290 30  0000 C CNN
F 2 "" H 4950 1200 60  0000 C CNN
F 3 "" H 4950 1200 60  0000 C CNN
	1    4950 1200
	1    0    0    -1  
$EndComp
$Comp
L +5V #PWR016
U 1 1 56D9883E
P 4950 2100
F 0 "#PWR016" H 4950 2190 20  0001 C CNN
F 1 "+5V" H 4950 2190 30  0000 C CNN
F 2 "" H 4950 2100 60  0000 C CNN
F 3 "" H 4950 2100 60  0000 C CNN
	1    4950 2100
	1    0    0    -1  
$EndComp
$Comp
L +5V #PWR017
U 1 1 56D9884D
P 4950 2900
F 0 "#PWR017" H 4950 2990 20  0001 C CNN
F 1 "+5V" H 4950 2990 30  0000 C CNN
F 2 "" H 4950 2900 60  0000 C CNN
F 3 "" H 4950 2900 60  0000 C CNN
	1    4950 2900
	1    0    0    -1  
$EndComp
NoConn ~ 6350 4950
NoConn ~ 6350 4850
NoConn ~ 6350 4750
$Comp
L CONN_3X2 P3
U 1 1 56D98D90
P 5350 1450
F 0 "P3" H 5350 1700 50  0000 C CNN
F 1 "CONN_3X2" V 5350 1500 40  0000 C CNN
F 2 "" H 5350 1450 60  0000 C CNN
F 3 "" H 5350 1450 60  0000 C CNN
	1    5350 1450
	1    0    0    -1  
$EndComp
$Comp
L CONN_3X2 P4
U 1 1 56D98DA9
P 5350 2350
F 0 "P4" H 5350 2600 50  0000 C CNN
F 1 "CONN_3X2" V 5350 2400 40  0000 C CNN
F 2 "" H 5350 2350 60  0000 C CNN
F 3 "" H 5350 2350 60  0000 C CNN
	1    5350 2350
	1    0    0    -1  
$EndComp
$Comp
L CONN_3X2 P5
U 1 1 56D98DB8
P 5350 3200
F 0 "P5" H 5350 3450 50  0000 C CNN
F 1 "CONN_3X2" V 5350 3250 40  0000 C CNN
F 2 "" H 5350 3200 60  0000 C CNN
F 3 "" H 5350 3200 60  0000 C CNN
	1    5350 3200
	1    0    0    -1  
$EndComp
NoConn ~ 5750 3250
NoConn ~ 4950 3250
NoConn ~ 4950 2400
NoConn ~ 5750 2400
NoConn ~ 4950 1500
NoConn ~ 5750 1500
Text Label 1200 6800 0    60   ~ 0
QEPA
Entry Wire Line
	1100 6700 1200 6800
Entry Wire Line
	1100 6600 1200 6700
Entry Wire Line
	1100 6500 1200 6600
Entry Wire Line
	1100 6400 1200 6500
Entry Wire Line
	1100 6300 1200 6400
Entry Wire Line
	1100 6200 1200 6300
Entry Wire Line
	1100 6100 1200 6200
Entry Wire Line
	1100 6000 1200 6100
Text Label 4200 5250 0    60   ~ 0
Z DRIVER BB
Entry Wire Line
	3800 5650 3900 5550
NoConn ~ 4850 5450
NoConn ~ 4850 5350
Text Label 4200 5550 0    60   ~ 0
QEPB BB
Wire Wire Line
	4850 4350 4500 4350
Wire Wire Line
	4500 4350 4500 4250
Wire Wire Line
	4750 4300 4750 4250
Wire Wire Line
	4750 4250 4850 4250
Wire Wire Line
	3000 6000 3150 6000
Wire Wire Line
	3000 5600 2900 5600
Wire Wire Line
	2900 5600 2900 5650
Wire Wire Line
	3150 6000 3150 5900
Wire Wire Line
	3900 5250 4850 5250
Wire Wire Line
	3900 5150 4850 5150
Wire Wire Line
	3900 5050 4850 5050
Wire Wire Line
	3900 4950 4850 4950
Wire Wire Line
	3900 4850 4850 4850
Wire Wire Line
	3900 4750 4850 4750
Wire Wire Line
	9500 800  9400 800 
Wire Wire Line
	9400 800  9400 850 
Wire Wire Line
	8350 2100 8250 2100
Wire Wire Line
	8250 2100 8250 2150
Wire Wire Line
	4300 4250 4300 4450
Wire Wire Line
	4300 4450 4850 4450
Wire Wire Line
	4100 4650 4850 4650
Wire Wire Line
	1650 850  1750 850 
Wire Wire Line
	1750 850  1750 900 
Wire Wire Line
	1650 1050 1700 1050
Wire Wire Line
	3000 6100 3000 6050
Wire Wire Line
	3000 6050 3100 6050
Wire Wire Line
	3100 6050 3100 6100
Wire Wire Line
	900  7050 1000 7050
Wire Wire Line
	4850 5750 4400 5750
Wire Wire Line
	6350 5650 5950 5650
Wire Wire Line
	6350 5550 5950 5550
Wire Wire Line
	6350 5450 5950 5450
Wire Wire Line
	6350 5350 5950 5350
Wire Wire Line
	4850 5650 4400 5650
Wire Wire Line
	9500 1350 9950 1350
Wire Wire Line
	9500 1650 9950 1650
Wire Wire Line
	9500 1950 9950 1950
Wire Wire Line
	7950 1800 8350 1800
Wire Wire Line
	7950 1500 8350 1500
Wire Wire Line
	7950 1200 8350 1200
Wire Bus Line
	7850 1300 7850 3950
Wire Bus Line
	7850 3950 5850 3950
Wire Bus Line
	7950 2450 9300 2450
Wire Bus Line
	9300 2450 9300 2300
Wire Bus Line
	9300 2300 10050 2300
Wire Bus Line
	10050 2300 10050 1450
Wire Wire Line
	7400 1950 8350 1950
Wire Wire Line
	8350 1650 7400 1650
Wire Wire Line
	8350 1350 7400 1350
Wire Wire Line
	9500 1500 10300 1500
Wire Wire Line
	9500 1800 10300 1800
Wire Wire Line
	9500 2100 10300 2100
Wire Bus Line
	10400 700  10400 2000
Wire Bus Line
	6250 700  10400 700 
Wire Wire Line
	5750 1300 6150 1300
Wire Wire Line
	6150 1400 5750 1400
Wire Wire Line
	6150 2200 5750 2200
Wire Wire Line
	5750 2300 6150 2300
Wire Wire Line
	5750 3050 6150 3050
Wire Wire Line
	5750 3150 6150 3150
Wire Bus Line
	4300 5750 4300 6800
Wire Bus Line
	4300 6800 5850 6800
Wire Bus Line
	5850 6800 5850 3950
Wire Wire Line
	9500 1200 9700 1200
Wire Wire Line
	9700 1200 9700 1050
Wire Wire Line
	4700 3950 4700 3900
Wire Wire Line
	4100 3900 4100 4650
Wire Wire Line
	1000 6100 1600 6100
Wire Wire Line
	1000 6200 1600 6200
Wire Wire Line
	1000 6300 1600 6300
Wire Wire Line
	1000 6400 1600 6400
Wire Wire Line
	1000 6500 1600 6500
Wire Wire Line
	1000 6600 1600 6600
Wire Wire Line
	1000 6700 1600 6700
Wire Wire Line
	1000 6800 1600 6800
Wire Wire Line
	1600 6900 1550 6900
Wire Wire Line
	1550 6900 1550 7000
Wire Wire Line
	1200 2250 2700 2250
Wire Wire Line
	1200 2350 2700 2350
Wire Wire Line
	1200 2450 2700 2450
Wire Wire Line
	1200 2550 2700 2550
Wire Wire Line
	1200 2650 2700 2650
Wire Wire Line
	1200 2750 2700 2750
Wire Wire Line
	1200 2050 2700 2050
Wire Wire Line
	1200 2150 2700 2150
Wire Wire Line
	4950 2900 4950 3150
Connection ~ 4950 3050
Wire Wire Line
	4950 2100 4950 2300
Connection ~ 4950 2200
Wire Wire Line
	4950 1200 4950 1400
Connection ~ 4950 1300
Wire Bus Line
	7300 800  7300 1850
Wire Bus Line
	6250 3050 6250 700 
Wire Wire Line
	1000 7050 1000 6900
Wire Bus Line
	1100 6700 1100 3400
Wire Bus Line
	1100 3400 2800 3400
Wire Bus Line
	2800 3400 2800 2150
Wire Bus Line
	3800 6000 3800 4850
Wire Wire Line
	4850 5550 3900 5550
Wire Wire Line
	3000 6200 3600 6200
Wire Wire Line
	3000 6300 3600 6300
Wire Wire Line
	3000 6400 3600 6400
Wire Wire Line
	3000 6500 3600 6500
Wire Wire Line
	3000 6600 3600 6600
Wire Wire Line
	3000 6700 3600 6700
Wire Wire Line
	3000 6800 3600 6800
Wire Wire Line
	6000 6900 3000 6900
Entry Wire Line
	3600 6200 3700 6100
Entry Wire Line
	3600 6300 3700 6200
Entry Wire Line
	3600 6400 3700 6300
Entry Wire Line
	3600 6500 3700 6400
Entry Wire Line
	3600 6600 3700 6500
Entry Wire Line
	3600 6700 3700 6600
Entry Wire Line
	3600 6800 3700 6700
Wire Bus Line
	3700 6000 3800 6000
Wire Wire Line
	6350 6250 6000 6250
Wire Wire Line
	6000 6250 6000 6900
Wire Bus Line
	3700 6700 3700 6000
Text Label 6000 6250 0    60   ~ 0
QEPA BB
Wire Wire Line
	1600 6000 1600 5850
$Comp
L +5V #PWR018
U 1 1 56D9A390
P 6200 4400
F 0 "#PWR018" H 6200 4490 20  0001 C CNN
F 1 "+5V" H 6200 4490 30  0000 C CNN
F 2 "" H 6200 4400 60  0000 C CNN
F 3 "" H 6200 4400 60  0000 C CNN
	1    6200 4400
	1    0    0    -1  
$EndComp
Wire Wire Line
	6200 4400 6200 4450
Wire Wire Line
	6200 4450 6350 4450
$Comp
L +3.3V #PWR019
U 1 1 56D9AAB5
P 1600 5850
F 0 "#PWR019" H 1600 5810 30  0001 C CNN
F 1 "+3.3V" H 1600 5960 30  0000 C CNN
F 2 "" H 1600 5850 60  0000 C CNN
F 3 "" H 1600 5850 60  0000 C CNN
	1    1600 5850
	1    0    0    -1  
$EndComp
$EndSCHEMATC
