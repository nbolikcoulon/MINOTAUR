#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 20:18:19 2019

@author: nbc
"""

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
            self.ChosenFolder = askdirectory()
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
    root1 = tk.Tk()
    ChooseFolder(root1)
    root1.mainloop()