#!/usr/bin/env python3 
import cmath
import pandas as texout

DEBUG=False     #To turn on/off debug messages. 
PRINT_INTERMEDIATE_RESULTS_TO_TXT_FILE=True

def replusZ(Z1, Z2):
    return (Z1*Z2)/(Z1+Z2)

def GMRfromA(A):
    return 0.78*cmath.sqrt(A/1000/1000/cmath.pi)

def self_reactance(GMD,GMR):
    return 0.145*cmath.log10(GMD/GMR)

def mutual_reactance(earth_return_depth,GMD):
    return 0.145*cmath.log10(earth_return_depth/GMD)

def current_divider_ratio(Z1,Z2):
    Ir1=Z2/(Z1+Z2)
    return Ir1

def quadratic_formula_roots(a, b, c):
    D=(b*b-4*a*c)
    x1=(-b+cmath.sqrt(D))/(2*a)
    x2=(-b-cmath.sqrt(D))/(2*a)
    return (x1, x2)

class PowerLine():
    def __init__(self, tower_count):
        self.U = 76.21                      #Phase voltage of the power line.
        self.Izsubst = complex(0, -31.5)    #Substation short circuit current
        self.Izref = self.Izsubst           #Reference short circuit current.
        self.Zgrid = self.U /self.Izsubst   #Impedance of grid 
        self.Uphwline = None
        self.Ustwline = None
        self.Uz = None          #Short circuit voltage of power line
        self.Zz = None          #Short circuit impedance of power line
        self.Iz = self.Izref    #Calculated short circuit current.
        self.line_length = 0    #To calculate overall line lenght.
        self.towers = None      #To store Tower objects
        self.Igmax = 0          #Maximum static wire current
        self.Igmaxtwrid = 0     #Index of span with max. static wire current
        self.create_powerline(tower_count)
        self.calculate_spans_self_and_mutual_impedances()
        for i in range(2):  
            #First pass to calculate real Iz.
            #Second pass to calculate real curents distribution. 
            self.set_short_circuit_current_for_towers()
            self.calculate_towers_thevenin_voltages_and_reactances()
            self.calculate_spans_static_wire_currents()
            self.calculate_overall_short_circuit_voltage()
            self.calculate_short_circuit_current()
            self.calculate_line_length_static_wire_max_current_and_location()
            if DEBUG:
                self.print_towers_intermediate_results_to_file()
    def create_powerline(self, tower_count):
        self.towers = []              
        for i in range(tower_count):
            tower = TowerAndSpan()
            self.towers.append(tower) 
        #Fist tower to substation span can be shorter than others. 
        self.towers[-1].Dtowers = 0.10  
        #To set earth return depth half at span of faulty tower. 
        #self.towers[0].Dgrr = self.towers[0].Dgrr / 2 
    def calculate_spans_self_and_mutual_impedances(self):
        for tower in self.towers:
            tower.calculate_self_and_mutual_impedances()
    def calculate_towers_thevenin_voltages_and_reactances(self):
        next_tower = None
        for tower in reversed(self.towers):
            tower.calculate_thevenin_voltage_and_reactance(next_tower)
            next_tower = tower
    def calculate_spans_static_wire_currents(self):
        previous_tower = None
        for tower in self.towers:
            tower.calculate_static_wire_currents(previous_tower)
            previous_tower = tower
    def calculate_overall_short_circuit_voltage(self):
        self.Uz = self.towers[-1].Istwsum * self.towers[-1].Rsubstgr
        for tower in self.towers:
            tower.calculate_span_voltage()
            self.Uz =self.Uz + tower.Uspan
    def calculate_short_circuit_current(self):
        self.Zz = self.Uz / self.Iz
        self.Iz = self.U / (self.Zz + self.Zgrid )
    def calculate_line_length_static_wire_max_current_and_location(self):
        self.line_length = 0
        self.Igmax = 0
        for idx, tower in enumerate(self.towers):
            self.line_length = self.line_length + int(tower.Dtowers *1000)
            if (abs(self.Igmax) < abs(tower.Istwsum)):
                self.Igmax = abs(tower.Istwsum)
                self.Igmaxtwrid = idx
    def set_short_circuit_current_for_towers(self):
        self.Izref = self.Iz
        for tower in self.towers:
            tower.Izref = self.Izref
    def print_towers_intermediate_results_to_file(self):
        msgtxt = (f'{"#":>3}{"Táv":>4}'
            f'{"Zstw":>13}{"Zphw":>13}'
            f'{"Zmut":>13}{"Ustphw":>13}'
            f'{"Ze":>13}{"Ue":>13}'
            f'{"Istwind":>13}{"Istwcon":>13}{"Ist":>13}'
            f'{"Uphstw":>16}'
            f'{"Uspan":>16}{"Izreftwr":>16}'
            f'\n')
        for idx, tower in enumerate(self.towers):
            msgtxt = msgtxt +(
                    f'{idx:3d}{tower.Dtowers*1000:4.0f}'
                    f'{tower.Zstw:13.2f}{tower.Zphw:13.2f}' 
                    f'{tower.Zphstmut:13.2f}{(tower.Ustphw):13.2f}'
                    f'{(tower.Ze):13.2f}{(tower.Ue):13.2f}'
                    f'{tower.Istwind:13.2f}{tower.Istwcond:13.2f}'
                    f'{(tower.Istwind+tower.Istwcond):13.2f}'
                    f'{(tower.Uphstw):16.2f}'
                    f'{tower.Uspan:16.2f}{tower.Izref:16.2f}' 
                    f'\n')
        with open("towermsg.txt", 'w') as f:
            f.write(msgtxt)
    def print_line_intermediate_results_to_file(self):
        msgtxt = (f'{"Twr c":>6}{"Izref":>13}'
            f'{"Uz":>22}{"Zz":>22}{"Iz":>22}\n')
        msgtxt = msgtxt + (
            f'{len(self.towers):6}{self.Izref:13.2f}'
            f'{self.Uz:16.3f}({abs(self.Uz):4.2f})'
            f'{self.Zz:16.2f}({abs(self.Zz):4.2f})'
            f'{self.Iz:16.2f}({abs(self.Iz):4.2f})\n') 
        with open("linemsg.txt", "w") as f:
            f.write(msgtxt)

class TowerAndSpan():
    def __init__(self):
        self.Izref=None     #reference value of earth fault current in kA 
        self.Dtowers=0.3    #general span lengt in kilometers
        self.Rtowergr=5     #earth resistance of the current pole in ohm
        self.Rsubstgr=0.05  #earth resistance of the substation in ohm
        self.Astw=95+50     #cross section of the static wire in mm2
        self.Aphw=250+40    #cross section of the phase wire in mm2
        self.Rphw=0.1154    #resistivity of the phase wire in ohm/km
        self.Rstw=0.2992    #resistivity of the static wire in ohm/km
        self.Rgrr=0         #resistivity of earth return in ohm/km
        self.Dphstw=10      #distance of ground and static wire in meter
        self.Dgrr=1000      #depth of earth return current in meter
        self.Zstw=None      #static wire self impedance ohm/span
        self.Zphw=None      #phase wire self impedance ohm/span
        self.Zphstmut=None  #phase-static mutual impedance ohm/span 
        self.Ustphw=None    #static wire votgate induced by phase current
        self.Uphstw=None    #phase wire voltega induced by static current
        self.Ue=None        #Thevenin's voltage 
        self.Ze=None        #Thevenin's impedance
        self.Istwind=None   #induced part of static wire current
        self.Istwcond=None  #conductive part of static wire current
        self.Istwsum=None   #oversall static wire current
        self.Ustwspan=None  #overall static wire voltage in span
        self.Uphwspan=None  #overall phase wire voltage in span
        self.Uspan=None     #span Kirchhoff's loop voltage  
    def calculate_self_and_mutual_impedances(self):
        self.Zstw=self.Dtowers*complex(
            self.Rstw,self_reactance(self.Dgrr,GMRfromA(self.Astw)))  
        self.Zphw = self.Dtowers * complex(
            self.Rphw, self_reactance(self.Dgrr, GMRfromA(self.Aphw))) 
        self.Zphstmut = self.Dtowers * complex(
            self.Rgrr, mutual_reactance(self.Dgrr, self.Dphstw))
    def calculate_thevenin_voltage_and_reactance(self, next):
        self.Ustphw=self.Zphstmut*self.Izref
        if next is None:    #This case I'am next to the substation.
            self.Ze=self.Zstw+self.Rsubstgr  
            self.Ue=self.Ustphw
            return
        self.Ze=self.Zstw+replusZ(self.Rtowergr,next.Ze)
        self.Ue=(self.Ustphw+next.Ue*self.Rtowergr 
            /(self.Rtowergr+next.Ze))
    def calculate_static_wire_currents(self,previous):
        current_ratio=current_divider_ratio(self.Ze,self.Rtowergr)
        if DEBUG:
            dbgmsg = f'Current ratio in static wire: {current_ratio:5.2f}'
            print(dbgmsg)
        if previous is None:    #Means I'am at the faulty tower. 
            self.Istwcond=self.Izref*current_ratio
            self.Istwind=self.Ue/(self.Rtowergr+self.Ze)
            self.Istwsum=(self.Istwcond+self.Istwind)
            self.Uphstw=self.Zphstmut*self.Istwsum
            return
        self.Istwcond=previous.Istwcond*current_ratio
        self.Istwind=(current_ratio*previous.Istwind
            +self.Ue/(self.Rtowergr+self.Ze))
        self.Istwsum=(self.Istwcond+self.Istwind)
        self.Uphstw=self.Zphstmut*(self.Istwsum)
    def calculate_span_voltage(self):
        #positive values are in the Iz direction
        self.Ustwspan=self.Istwsum*self.Zstw-self.Izref*self.Zphstmut
        self.Uphwspan=self.Izref*self.Zphw-self.Istwsum*self.Zphstmut
        self.Uspan=self.Uphwspan+self.Ustwspan
        if DEBUG:
            dbgtxt = (f'Ustwspan = Istw * Zstw - Izref * Zphstmut \n'
                f'{self.Ustwspan:.2f}={self.Istwsum:.2f}*{self.Zstw:.2f}'
                f'-{self.Izref:.2f}*{self.Zphstmut:.2f}\n'
                f'Uphwspan = Izref * Zphw - Istw * Zphstmut \n'
                f'{self.Uphwspan:.2f}={self.Izref:.2f}*{self.Zphw:.2f}'
                f'-{self.Istwsum:.2f}*{self.Zphstmut:.2f}\n'
                f'Uspan={self.Uspan:.2} |Uspan|={abs(self.Uspan)}\n'
                f'{abs(self.Uphwspan/self.Uspan):5.2f}'
                f'{abs(self.Ustwspan/self.Uspan):5.2f}'
                f'{abs(self.Ustwspan/self.Izref):5.2f}'
                f'{abs(self.Uphwspan/self.Izref):5.2f}'
                f'{abs((self.Uphwspan+self.Ustwspan)/self.Izref):5.2f}')
            print(dbgtxt)

def export_towers_results_to_latex(towers):
    tower_idx = []
    tower_distance = []
    Zphw = []
    Zstw = []
    Zstphw = []
    Istwind = []
    Istwcond = []
    Istwsum = []
    for idx, tower in enumerate(towers):
        tower_idx.append(idx)
        tower_distance.append(int(tower.Dtowers*1000))
        Zphw.append(abs(tower.Zphw))
        Zstw.append(abs(tower.Zstw))
        Zstphw.append(abs(tower.Zphstmut))
        Istwind.append(abs(tower.Istwind))
        Istwcond.append(abs(tower.Istwcond))
        Istwsum.append(abs(tower.Istwsum))
    data = {
            "Oszlop": tuple(tower_idx),
            "Oszlopköz": tuple(tower_distance),
            "Zphw": tuple(Zphw),
            "Zstw": tuple(Zstw),
            "Zstphw": tuple(Zstphw),
            "Ig kond.": tuple(Istwcond),
            "Ig indukt.": tuple(Istwind),
            "Ig": tuple(Istwsum),
        }
    filename = (f'./out/tavvezetek{len(towers):03d}py.tex')
    write_data_to_latex_file(data, filename)

def write_data_to_latex_file(data, filename):
    formatteddata = texout.DataFrame(data)
    latextxt = formatteddata.to_latex(
        index = False,
        position = "h",
        float_format="%.2f",
        ).replace(".", ",")
    print(latextxt)
    with open(filename, "w") as f:
        f.write(latextxt)
    
def calculate_set_of_power_lines(linesets):
    tower_counts = []
    line_lengths = []
    Uz = []
    Igmax = []
    Igmaxtwrid = []
    Iz = []
    Izrel = []
    Igmaxrelline = []
    Igmaxrelsubst = []
    for tower_count in linesets:
        powerline = PowerLine(tower_count)
        if PRINT_INTERMEDIATE_RESULTS_TO_TXT_FILE:
            powerline.print_towers_intermediate_results_to_file()
            powerline.print_line_intermediate_results_to_file()
        export_towers_results_to_latex(powerline.towers)
        tower_counts.append(len(powerline.towers))
        line_lengths.append(powerline.line_length)
        Igmax.append(powerline.Igmax)
        Igmaxtwrid.append(powerline.Igmaxtwrid)
        Uz.append(abs(powerline.Uz))
        Iz.append(abs(powerline.Iz))
        Izrel.append(abs(powerline.Iz) / abs(powerline.Izsubst))
        Igmaxrelline.append(abs(powerline.Igmax / powerline.Iz))
        Igmaxrelsubst.append(abs(powerline.Igmax / powerline.Izsubst))
    data = {
        "Oszlopszám" : tuple(tower_counts),
        "Vonalhossz": tuple(line_lengths),
        "Igmax": tuple(Igmax),
        "Igmax oszl": tuple(Igmaxtwrid),
        "Uz": tuple(Uz),
        "Iz": tuple(Iz),
        "Izrel": tuple(Izrel),
        "Igmaxrel": tuple(Igmaxrelline),
        "Igmaxrelsubst":tuple(Igmaxrelsubst)
    }
    filename = f'./out/tavvezsorozat.tex'
    write_data_to_latex_file(data, filename)
    
def main():
    linesets = tuple(range(1, 26)) + (50, 100, 200)
    calculate_set_of_power_lines(linesets)
    

if __name__=="__main__":
    main()
