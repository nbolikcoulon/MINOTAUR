load 1D3Z.pdb
hide
show cartoon, all

#Color the whole structure in grey as default color
color grey, 1D3Z

#Define a color in the blue to red range for each residue we have data for
set_color resi_2 = [1.0, 0.0, 0.0]
set_color resi_3 = [0.8564287837087986, 0.0, 0.1435712162912014]
set_color resi_4 = [0.3615393632802415, 0.0, 0.6384606367197585]
set_color resi_5 = [0.0, 0.0, 1.0]

#Coloring each residue according to its color
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
color resi_{AA}, resi {AA}
