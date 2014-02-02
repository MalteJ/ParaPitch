###########################################
#                ParaPitch                #
###########################################
#                                         #
#   Copyright 2013 by Malte Janduda       #
#   all rights reserved                   #
#                                         #
#   feel free to contact me!              #
#   mail@malte-janduda.de                 #
#                                         #
###########################################

import numpy as np
from matplotlib import *

class LinAlg:
    def rotate2d(vector, angle):
        a = angle
        R = np.array(((np.cos(a), -np.sin(a)),
                      (np.sin(a),  np.cos(a))))
        return np.dot(R, vector)


class Polar:
    polar = []

    def readCsv(self, csvfilename):
        import csv
        with open(csvfilename, 'rt') as csvfile:
            rdr = csv.reader(csvfile, delimiter=',')
            i=1
            for row in rdr:
                if i>1: # first line is header
                    arr = np.array((float(row[0])/180*np.pi, float(row[1]), float(row[2]), float(row[3])))
                    self.polar.append(arr)
                i=i+1

    def getAlpha(self, alpha):
        polar = self.polar
        assert len(polar) > 0
        assert alpha >= polar[0][0]
        assert alpha <= polar[len(polar)-1][0]

        for i in range(0,len(polar)):
            if alpha >= polar[i][0] and alpha <= polar[i+1][0]:
                return polar[i] * (polar[i+1][0] - alpha) / (polar[i+1][0] - polar[i][0]) \
                     + polar[i+1] * (alpha - polar[i][0]) / (polar[i+1][0] - polar[i][0])
            
    def getCL(self, alpha):
        return self.getAlpha(alpha)[1]
            
    def getCD(self, alpha):
        return self.getAlpha(alpha)[2]

    def getCM(self, alpha):
        return self.getAlpha(alpha)[3]


class AeroDevice:
    alpha = 0.0
    u     = 0.0
    rho   = 0.0
    A     = 0.0
    cd    = 0.0
    cl    = 0.0
    cm    = 0.0

    # mass and moment of inertia
    m   = 0.0
    moi = 1.0

    # position in parent system
    x_rel = np.zeros(2)

    def getCD(self):
        return self.cd

    def getCL(self):
        return self.cl

    def getCM(self):
        return self.cm

    def calculateForce(self):
        # calculate and rotate back to parent coordinate system
        return LinAlg.rotate2d(0.5 * self.rho * self.u**2 * np.array((-self.getCD(), self.getCL())) * self.A, \
                               -self.alpha)

    def calculateTorque(self):
        # calculate and rotate back to parent coordinate system
        return 0.5 * self.rho * self.u**2 * self.getCM() * self.A


class Lines(AeroDevice):
    def getCD(self):
        # lines have a full cd when alpha = 0. if alpha = pi/2 (90 deg) the cd should be 0
        # due to no projected area.
        return np.cos(self.alpha) * self.cd


class Paraglider(AeroDevice):
    def setPolar(self, polar):
        self.polar = polar

    def getCL(self):
        return self.polar.getCL(self.alpha)

    def getCD(self):
        #print("alpha = %f and cl = %f cd = %f (tan^-1 cd/cl = %f) cm = %f v = %f" % \
        #      (self.alpha*180/np.pi, self.polar.getCL(self.alpha), self.polar.getCD(self.alpha), \
        #       np.arctan(self.polar.getCD(self.alpha)/self.polar.getCL(self.alpha))*180/np.pi, \
        #       self.polar.getCM(self.alpha), self.u))
        return self.polar.getCD(self.alpha)

    def getCM(self):
        return self.polar.getCM(self.alpha)


class MovingSystem:
    x      = np.zeros(2)    # position
    phi    = 0.0            # rotation
    v      = np.zeros(2)    # velocity
    omega  = 0.0            # rot. velocity
    dv     = np.zeros(2)    # transl. acceleration
    domega = 0.0            # rot. acceleration

    m      = 0.0            # mass
    moi    = 1.0            # moment of inertia

    def getIntegrationVariables(self):
        return ((self.x, self.v), (self.phi, self.omega),(self.v, self.dv), (self.omega, self.domega))

    def setIntegratedVariables(self, arr):
        # must be in the same order as in "getIntegrationVariables"
        (self.x, self.phi, self.v, self.omega) = arr


class ParagliderSystem(MovingSystem):
    g   = 9.80665 # m/s^2
    rho = 1.225 # kg/m^3

    def setParaglider(self, p):
        self.paraglider = p

    def setLines(self, l):
        self.lines = l

    def setPilot(self, p):
        self.pilot = p

    def adjustCenterOfGravity(self):
        cg = ( self.paraglider.m * self.paraglider.x_rel \
             + self.lines.m    * self.lines.x_rel \
             + self.pilot.m    * self.pilot.x_rel \
             ) / (self.paraglider.m + self.lines.m + self.pilot.m)
        self.paraglider.x_rel = self.paraglider.x_rel - cg
        self.lines.x_rel = self.lines.x_rel - cg
        self.pilot.x_rel = self.pilot.x_rel - cg
        cg = ( self.paraglider.m * self.paraglider.x_rel \
             + self.lines.m    * self.lines.x_rel \
             + self.pilot.m    * self.pilot.x_rel \
             ) / (self.paraglider.m + self.lines.m + self.pilot.m)


    def calculateState(self):
        d = 1
        # calculating values for the paraglider
        paragliderV = self.v + self.omega*d * np.array((-self.paraglider.x_rel[1], self.paraglider.x_rel[0]))
        self.paraglider.rho = self.rho
        self.paraglider.u   = np.linalg.norm(paragliderV) 
        self.paraglider.alpha = -np.arctan2(paragliderV[1], paragliderV[0]) + self.phi

        # calculating values for the lines
        linesV = self.v + self.omega*d * np.array((-self.lines.x_rel[1], self.lines.x_rel[0]))
        self.lines.rho      = self.rho
        self.lines.u        = np.linalg.norm(linesV)
        self.lines.alpha    = -np.arctan2(linesV[1], linesV[0]) + self.phi

        # calculating values for the pilot
        pilotV = self.v + self.omega*d * np.array((-self.pilot.x_rel[1], self.pilot.x_rel[0]))
        self.pilot.rho      = self.rho
        self.pilot.u        = np.linalg.norm(pilotV)
        self.pilot.alpha    = -np.arctan2(pilotV[1], pilotV[0]) + self.phi

        # summing up all masses to get the system's mass
        self.m = self.paraglider.m + self.lines.m + self.pilot.m

        # using parallel axis theorem to get the moment of inertia
        self.moi = self.paraglider.moi + self.paraglider.m * np.dot(self.paraglider.x_rel, self.paraglider.x_rel) + \
                   self.lines.moi      + self.lines.m      * np.dot(self.lines.x_rel, self.lines.x_rel) + \
                   self.pilot.moi      + self.pilot.m      * np.dot(self.pilot.x_rel, self.pilot.x_rel)

        # calculating accelerations
        a_paraglider_abs = self.domega * np.array((-self.paraglider.x_rel[1], self.paraglider.x_rel[0])) \
                           - self.domega**2 * self.paraglider.x_rel + self.dv
        a_lines_abs = self.domega * np.array((-self.lines.x_rel[1], self.lines.x_rel[0])) \
                      - self.domega**2 * self.lines.x_rel + self.dv
        a_pilot_abs = self.domega * np.array((-self.pilot.x_rel[1], self.pilot.x_rel[0])) \
                      - self.domega**2 * self.pilot.x_rel + self.dv

        # gathering forces and summing up
        FG_norm = LinAlg.rotate2d(np.array((0,-1))*self.g, -self.phi)
        # gravitational force + aerodynamic force + centripetal force
        Fpar_local = FG_norm * self.paraglider.m + self.paraglider.calculateForce() \
                     - (a_paraglider_abs - self.dv) / self.paraglider.m
        Flin_local = FG_norm * self.lines.m + self.lines.calculateForce() - (a_lines_abs - self.dv) / self.lines.m
        Fpil_local = FG_norm * self.pilot.m + self.pilot.calculateForce() - (a_pilot_abs - self.dv) / self.pilot.m
        Fsum_local = Fpar_local + Flin_local + Fpil_local

        # gathering torques and summing up
        Fpar_abs = LinAlg.rotate2d(Fpar_local, self.phi)
        Flin_abs = LinAlg.rotate2d(Flin_local, self.phi)
        Fpil_abs = LinAlg.rotate2d(Fpil_local, self.phi)
        Fsum_abs = Fpar_abs + Flin_abs + Fpil_abs
        #print("self.v = %s self.omega = %s paraglider.u = %f" % (self.v, self.omega, self.paraglider.u))
        #print("Fpar_abs = %s" % Fpar_abs)
        #print("Flin_abs = %s" % Flin_abs)
        #print("Fpil_abs = %s" % Fpil_abs)
        #print("Fpar_local = %s" % Fpar_local)
        #print("Flin_local = %s" % Flin_local)
        #print("Fpil_local = %s" % Fpil_local)

        Mpar = self.paraglider.calculateTorque() - Fpar_local[0] * self.paraglider.x_rel[1] \
               + Fpar_local[1] * self.paraglider.x_rel[0]
        Mlin = self.paraglider.calculateTorque() - Flin_local[0] * self.lines.x_rel[1] \
               + Flin_local[1] * self.lines.x_rel[0]
        Mpil = self.paraglider.calculateTorque() - Fpil_local[0] * self.pilot.x_rel[1] \
               + Fpil_local[1] * self.pilot.x_rel[0]
        Msum = Mpar + Mlin + Mpil
        #print("Msum = Mpar + Mlin + Mpil = %f + %f + %f = %f" % (Mpar, Mlin, Mpil, Msum))
        #print("moi = %f" % self.moi)

        # calculating accelerations
        self.dv     = Fsum_abs / self.m
        self.domega = Msum / self.moi


    def __str__(self):
        return "x  = %s\nv  = %s\nphi   = %s\nomega = %s\nGZ = %f (GW = %f [°])" % \
                (self.x, self.v, self.phi*180/np.pi, self.omega*180/np.pi, -self.v[0] / self.v[1], \
                 np.arctan(self.v[1]/self.v[0]*(-1)) * 180/np.pi)

class IntegrationMethod:
    time = 0.0
    timestep = 1.0

    def setSystem(self, sys):
        self.system = sys

    def iterate(self):
        pass

    def __str__(self):
        return "#%i\tt = %f\tdt = %f" % (round(self.time/self.timestep), self.time, self.timestep)


class ForwardEuler(IntegrationMethod):
    def iterate(self):
        dt = self.timestep
        self.time = self.time + dt
      
        intvars = self.system.getIntegrationVariables()
        res = []
        for (fn, df) in intvars:
            res.append(fn + dt * df)
        self.system.setIntegratedVariables(res)
        self.system.calculateState()

class ClassicalRungeKutta(IntegrationMethod):
    def iterate(self):
        dt = self.timestep
        self.time = self.time + dt

        # 1st step: half step predictor

        intvarsN = self.system.getIntegrationVariables()
        res = []
        fsn12 = []
        for (fn, df) in intvarsN:
            fsn12.append(fn + dt/2 * df)
        
        self.system.setIntegratedVariables(fsn12)

        # 2nd step: half step corrector

        self.system.calculateState()
        intvarsNS12 = self.system.getIntegrationVariables()
        fssn12 = []
        for i in range(0, len(intvarsNS12)):
            fn = intvarsN[i][0]
            df = intvarsNS12[i][1]
            fssn12.append(fn + dt/2 * df)

        self.system.setIntegratedVariables(fssn12)

        # 3rd step: rectangle method as predictor

        self.system.calculateState()
        intvarsNSS12 = self.system.getIntegrationVariables()
        fsn1 = []
        for i in range(0, len(intvarsNSS12)):
            fn = intvarsN[i][0]
            df = intvarsNSS12[i][1]
            fsn1.append(fn + dt * df)

        self.system.setIntegratedVariables(fsn1)

        # 4th step: simpson rule as final corrector

        self.system.calculateState()
        intvarsNS1 = self.system.getIntegrationVariables()
        fn1 = []
        for i in range(0, len(intvarsNS1)):
            fn      = intvarsN[i][0]
            dfn     = intvarsN[i][1]
            dfns12  = intvarsNS12[i][1]
            dfnss12 = intvarsNSS12[i][1]
            dfns1   = intvarsNS1[i][1]
            fn1.append(fn + dt/6 * (dfn + 2*dfns12 + 2*dfnss12 + dfns1))

        self.system.setIntegratedVariables(fn1)
        

if __name__== "__main__":
    paraglider     = Paraglider()
    paraglider.A   = 20.0
    paraglider.m   = 8.0
    paraglider.moi = 1.0
    paraglider.x_rel = np.array((0.8,7))

    polar = Polar()
    polar.readCsv("gleitschirm.csv")
    paraglider.setPolar(Polar())

    lines     = Lines()
    lines.cl  = 0.0
    lines.cd  = 0.35
    lines.A   = 1.0
    lines.m   = 1.5
    lines.moi = 1.0
    lines.x_rel = np.array((0,3.5))

    pilot     = AeroDevice()
    pilot.cl  = 0.0
    pilot.cd  = 0.5
    pilot.A   = 0.5
    pilot.m   = 100.0
    pilot.moi = 1.0
    pilot.x_rel = np.array((0,0))

    sys = ParagliderSystem()
    sys.rho = 1.2

    sys.setParaglider(paraglider)
    sys.setLines(lines)
    sys.setPilot(pilot)


    #integrator = ForwardEuler()
    integrator = ClassicalRungeKutta()
    integrator.timestep = 0.05
    integrator.setSystem(sys)

    # initial conditions
    integrator.time = -250.0
    sys.v = np.array((11.1,-1))
    sys.phi = 3.1/180*np.pi
    
    # Preparations for Solving
    sys.adjustCenterOfGravity()
    sys.calculateState()

    print("t,u,v,alpha,phi,gamma")
    for i in range(0,10000):
        integrator.iterate()
        print("%f,%f,%f,%f,%f,%f" % (integrator.time, sys.v[0], sys.v[1], paraglider.alpha* 180/np.pi, sys.phi* 180/np.pi, np.arctan(sys.v[1]/sys.v[0]*(-1)) * 180/np.pi))
        if integrator.time >= 0:
            ## Pilot Flaeche verdoppeln
            # pilot.A = 1

            # Gas geben: Pilot innerhalb von 10 Sekunden 50cm nach vorne schieben
            paraglider.x_rel = np.array((0.8-min(0.5,0.5*(integrator.time/10)),7))
            sys.adjustCenterOfGravity()

