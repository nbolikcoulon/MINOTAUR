load 1D3Z.pdb
hide
show cartoon, all

#Color the whole structure in grey as default color
color grey, 1D3Z

#Define a color in the blue to red range for each residue we have data for
set_color resi_2 = [1.0, 0.0, 0.0]
set_color resi_3 = [0.6134386027082572, 0.0, 0.38656139729174277]
set_color resi_4 = [0.0, 0.0, 1.0]
set_color resi_5 = [0.7474590627594241, 0.0, 0.2525409372405759]

#Coloring each residue according to its color
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
