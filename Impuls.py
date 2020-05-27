import tkinter as tk
from tkinter import ttk
import numpy as np
import math

from Graphs import *
class Impuls:

    def __init__(self, impulsFrame, num_a, num_b, settings, radSelected, resolution, step_size):
        self.impulsFrame = impulsFrame

        self.a = self.convert(num_a)
        self.b = self.convert(num_b)
        self.settings = self.convert_dict(settings)
        self.radSelected = radSelected
        self.resolution = float(resolution.get())
        self.step = int((1./self.resolution)*float(step_size.get()))
        Graphs.step_init(Graphs, self.step)
        self.input_singal = []
        self.output_signal = []
        self.time = [0]
        self.max_ab = max(self.a + self.b)

        self.calculations()

        self.w, self.mag, self.phase = 0, 0, 0
        self.bode()

    
    def calculations(self) :
        # wyznaczanie czasu
        sum = 0
        for _ in range(int(self.settings['duration']) * (int(1.0 / self.resolution)) - 1) :
            sum += self.resolution
            self.time.append(sum)

        # wyznaczanie sygnału wejściowego
        if self.radSelected == 1 :  # prostokatny
            self.input_singal = [
                self.settings["amplitude"] if t % self.settings["period"] < (self.settings["fulfillment"] / 100.) *
                                              self.settings["period"]
                else -self.settings["amplitude"] for t in self.time]

        elif self.radSelected == 2 :  # skok
            self.input_singal = [self.settings["amplitude"] if i >= self.settings["start"] else 0 for i in self.time]

        elif self.radSelected == 3 :  # sinusoida
            self.input_singal = [self.settings["amplitude"] * math.sin(2 * math.pi * (1 / self.settings["period"]) * i)
                                 for i in self.time]

        elif self.radSelected == 4 :  # trojkatny
            sum = 0
            for t in self.time :
                if t % (self.settings["period"]) <= 0.5 * self.settings["period"] :
                    self.input_singal.append(
                        (t % self.settings["period"]) * (self.settings["amplitude"] / (0.5 * self.settings["period"])))
                else :
                    self.input_singal.append((t % self.settings["period"] - self.settings["period"]) * (
                            -self.settings["amplitude"] / (0.5 * self.settings["period"])))

        # wyznaczanie sygnału wyjściowego
        self.output_signal = [0 for _ in range(math.floor(len(self.input_singal) / self.step))]
        v = [0, 0, 0, 0]
        v[3] = self.input_singal
        for _ in range(int(self.settings["duration"] / self.resolution)) :
            v[2] = self.integrate(v[3])
            v[1] = self.integrate(v[2])
            v[0] = self.integrate(v[1])
            v = self.amplifying(0, v)
            v[3] = self.subtracting(v[2], v[1], v[0])
        i = len(v) - 1
        while i >= 0 :
            if i != 3 : v[i] = self.integrate(v[i + 1])
            i -= 1
        v = self.amplifying(1, v)
        self.output_signal = self.adding(self.output_signal, v[3], v[2], v[1], v[0])


    def integrate(self, data):
        sum = 0
        integral = [data[0]]
        dx = self.time[self.step] - self.time[0]
        for f in data[1:]:
            sum += f * dx /2.
            integral.append(sum)
        return integral
    
    def subtracting(self, v2, v1, v0):
        result = [0]
        i = 0
        for i in range(len(self.input_singal[0::self.step])):
            result.append(self.input_singal[0::self.step][i] - v2[i] - v1[i] - v0[i])
        return result

    def adding(self, out, v3, v2, v1, v0):
        result = []
        for i in range(len(out)) :
            result.append(out[i] + v3[i] + v2[i] + v1[i] + v0[i])
        return result

    def multiplying(self, p, v):
        result = []
        for i in range(len(v)) :
            result.append(p * v[i])
        return result

    def amplifying(self, in_or_out, v):
        if in_or_out == 0 :
            for i in range(len(self.a)) :
                v[i] = self.multiplying(self.a[i], v[i])
        else :
            for i in range(len(self.b)) :
                v[i] = self.multiplying(self.b[i], v[i])
        return v

    # funkcja zwracajaca czy uklad jest stabilny 
    def is_stable(self):
        if self.a[0] <= 0 or self.a[1] <= 0 or self.a[2] <= 0:
             return 0
        elif self.a[2]*self.a[1] - self.a[0] == 0:
            return 1
        elif self.a[2] > 0 and self.a[2]*self.a[1] - self.a[0] > 0:
            return 2

    # funkcja konwertujaca StringVar to int
    def convert(self, num_x):
        a = []
        for x in num_x:
            if self.is_number(x.get()):
                a.append(float(x.get()))
            else:
                a.append(0)
        return a

    # funkcja konwertujaca StringVar to int w słowniku
    def convert_dict(self, set):
        d = {}
        for k,v in set.items():
            if self.is_number(v.get()):
                d[k] = float(v.get())
            else:
                d[k] = 0
        return d

    def is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def calculate_pha(self,a,b,w,dt):
        s = b[1] * w - b[3] * (w ** 3)
        p = b[0] - b[2] * (w ** 2)
        r = a[1] * w - 1 * (w ** 3)
        t = a[0] - a[2] * (w ** 2)
        value1 =0
        try:
            value1 = np.rad2deg(np.arctan(s / p) - np.arctan(r / t))
        except:
            s = b[1] * w - b[3] * ((w - dt) ** 3)
            p = b[0] - b[2] * ((w - dt) ** 2)
            r = a[1] * w - 1 * ((w - dt) ** 3)
            t = a[0] - a[2] * ((w - dt) ** 2)
            value1 = np.rad2deg(np.arctan(s / p) - np.arctan(r / t))

        return value1;


    def calculate(self,a,b,w,dt):
        s = b[1] * w - b[3] * (w ** 3)
        p = b[0] - b[2] * (w ** 2)
        r = a[1] * w - (w ** 3)
        t = a[0] - a[2] * (w ** 2)
        try:
            value1 = math.log10((s ** 2 + p ** 2) ** 10)
        except:
            try:
                s = b[1] * w - b[3] * ((w-dt)** 3)
                p = b[0] - b[2] * ((w-dt) ** 2)
                value1 = math.log10((s ** 2 + p ** 2) ** 10)
            except:
                value1 = 0

        try:
            value2 = math.log10( (r**2 + t**2 )**10 )
        except:
            try:
                r = a[1] * w - ((w - dt) ** 3)
                t = a[0] - a[2] * ((w - dt) ** 2)
                value2 = math.log10((r ** 2 + t ** 2) ** 10)
            except:
                value2 = 0

        return value1 - value2


    def bode(self):
        self.phase = []
        self.mag = []
        dt = 0.01
        self.w = np.arange(dt, 100.0, dt)
        for i in self.w:
            self.phase.append(self.calculate_pha(self.a, self.b, i, dt))
            self.mag.append(self.calculate(self.a, self.b, i, dt))

