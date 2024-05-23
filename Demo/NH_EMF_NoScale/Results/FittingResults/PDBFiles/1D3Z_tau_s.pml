load 1D3Z.pdb
hide
show cartoon, all

#Color the whole structure in grey as default color
color grey, 1D3Z

#Define a color in the blue to red range for each residue we have data for
set_color resi_2 = [1.0, 0.0, 0.0]
set_color resi_3 = [0.3126521174017098, 0.0, 0.6873478825982902]
set_color resi_4 = [0.7113339158750932, 0.0, 0.2886660841249068]
set_color resi_5 = [0.0, 0.0, 1.0]

#Coloring each residue according to its color
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
