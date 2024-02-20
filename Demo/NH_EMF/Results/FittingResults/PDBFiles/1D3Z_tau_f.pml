load 1D3Z.pdb

hide

show cartoon, all

#Color the whole structure in grey as default color
color grey, 1D3Z

#Define a color in the blue to red range for each residue we have data for
set_color resi2 = [0.7090197398726281, 0.0, 0.2909802601273719]
set_color resi3 = [1.0, 0.0, 0.0]
set_color resi4 = [0.27811350220692777, 0.0, 0.7218864977930722]
set_color resi5 = [0.0, 0.0, 1.0]

#Coloring each residue according to its color
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
