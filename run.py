#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tkinter as tk
from tkinter.filedialog import askdirectory

import os
import sys


class ChooseFolder(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        
        self.master = master
        
        self.gui()
        
        
    def gui(self):
        def BrowseFolder(self):
            self.ChosenFolder = askdirectory(initialdir = os.path.dirname(os.path.abspath(__file__)))
            self.Folder.config(text=os.path.basename(self.ChosenFolder))
            
        def Quit(self):
            sys.exit()
            
        def RUN(self):
            folder = self.ChosenFolder
            os.system("python " + folder + "/setup.py build_ext --inplace")
            os.system("python " + folder + "/setup.py build_ext --inplace")
            os.system("python MINOTAUR.py " + folder)
            
            sys.exit()
        
        self.master.title("Folder for Expressions")
        self.master.geometry("400x75")
        self.grid(row=1, column=2)
        
        self.Folder = tk.Button(self, text="Browse folder containing expressions", command = lambda: BrowseFolder(self))

        self.Start = tk.Button(self, text="Launch MINOTAUR", command = lambda: RUN(self))
        self.Quit = tk.Button(self, text="Quit", command = lambda: Quit(self))

        self.Folder.grid(row=0, column=0, columnspan=2)
        self.Start.grid(row=1, column=0)
        self.Quit.grid(row=1, column=1)
        

if __name__ == "__main__":
    if len(sys.argv) > 4:
        print('Too many arguments.')
        sys.exit()
        
    if len(sys.argv) > 1:
        folder = str(sys.argv[1])
        os.system("python " + folder + "/setup.py build_ext --inplace")
        os.system("python " + folder + "/setup.py build_ext --inplace")
        
        if len(sys.argv) == 2:
            os.system(f"python MINOTAUR.py {folder}")
            sys.exit()
            
        if len(sys.argv) == 3:
            os.system(f"python MINOTAUR.py {folder} {sys.argv[2]}")
            sys.exit()
            
        if len(sys.argv) == 4:
            os.system(f"python MINOTAUR.py {folder} {sys.argv[2]} {sys.argv[3]}")
            sys.exit()
            
    root1 = tk.Tk()
    ChooseFolder(root1)
    root1.mainloop()
