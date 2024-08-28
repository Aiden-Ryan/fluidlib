import CoolProp.CoolProp as CP

class Material:
    def __init__(self,density,hardness=0, TS = 0, YS = 0, E=0, Charpy =0, Izod =0, 
                 resitivity=0, magneticPermeability =0,
                 specificHeatCapacity=0, k = 0, meltingPoint = 0, CTE20C =0, CTE250C = 0, CTE500C = 0):
        self.density = density #kg/m3
        self.hardness = hardness
        self.TS = TS #MPa
        self.YS = YS #MPa
        self.E = E #GPa
        self.Charpy = Charpy #J
        self.Izod = Izod #J
        self.resistivity = resitivity #ohm-cm
        self.magneticPermeability = magneticPermeability
        self.specificHeatCapacity = specificHeatCapacity #J/kg-K
        self.k = k #W/m-K
        self.meltingPoint = meltingPoint #K
        self.CTE20C = CTE20C #micrometers/m-C
        self.CTE250C = CTE250C
        self.CTE500C = CTE500C

'''MATERIAL PROPERTIES SOURCED FROM MATWEB'''

SS316 = Material(800,79, 580, 290, 193, 105, 129, 7.4e-5, 1.008, 500, 16.3, CTE20C = 16, CTE250C = 16.2, CTE500C = 17.5)


