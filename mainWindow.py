import tkinter as tk
from tkinter import ttk

from DataEntry import *
from Graphs import *

class Window():
    def __init__(self, master):
        self.master = master
        master.title("")
        master.minsize(1415, 810)
        
        self.dataFrame = tk.LabelFrame(master, text = "Wprowadź dane", relief="groove")
        self.dataFrame.place(x = 5, y = 5, width = 300, height = 800)
        self.dataEntry = DataEntry(self.dataFrame)
        

        # tworzenie ramki do wyświetlania czy układ jest stabilny
        self.stabilityFrame = tk.LabelFrame(self.dataFrame, relief="groove")
        self.stabilityFrame.pack(side="top", fill="x", padx = 2)
        self.stability_label = tk.Label(self.stabilityFrame, text = "Czy jest stabilny?", font = "Arial 16")
        self.stability_label.pack(fill="both", ipady = 5, pady = 2, padx = 2)
              
        # przycisk do wykreslenia wykresow
        self.draw_button = tk.Button(self.dataFrame, text = "Rysuj", font = "Arial 16", relief="groove", command = self.draw_function)
        self.draw_button.pack(side="bottom", fill="x", ipady = 5, pady = 2, padx = 2)

        self.drawFrame = tk.LabelFrame(master, text = "Wykresy", relief="groove")
        self.drawFrame.place(x = 310, y = 5, width = 1100, height = 800)

    def draw_function(self):
        for widget in self.drawFrame.winfo_children():
            widget.destroy()

        self.impulse = self.dataEntry.get_impuls()
        self.impulse.wyswietl_wspolczynniki()

        if self.impulse.is_stable() == 0:
            self.stability_label.config(text = "NIESTABILNY", compound = "center", foreground = "red")
        elif self.impulse.is_stable() == 1:
            self.stability_label.config(text = "NA GRANICY STABILNOŚCI", compound="center", foreground="orange")
        elif self.impulse.is_stable() == 2:
            self.stability_label.config(text = "STABILNY", compound = "center", foreground = "green")

        self.Graphs = Graphs(self.drawFrame, self.impulse)

        


if __name__ == "__main__":
    root = tk.Tk()
    app = Window(root)
    root.mainloop()