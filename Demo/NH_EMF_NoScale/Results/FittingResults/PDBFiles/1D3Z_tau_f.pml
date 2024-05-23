load 1D3Z.pdb
hide
show cartoon, all

#Color the whole structure in grey as default color
color grey, 1D3Z

#Define a color in the blue to red range for each residue we have data for
set_color resi_2 = [0.5244943626288757, 0.0, 0.4755056373711243]
set_color resi_3 = [1.0, 0.0, 0.0]
set_color resi_4 = [0.08483505700667271, 0.0, 0.9151649429933273]
set_color resi_5 = [0.0, 0.0, 1.0]

#Coloring each residue according to its color
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
