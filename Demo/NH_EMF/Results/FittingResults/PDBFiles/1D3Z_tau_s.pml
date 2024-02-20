load 1D3Z.pdb

hide

show cartoon, all

#Color the whole structure in grey as default color
color grey, 1D3Z

#Define a color in the blue to red range for each residue we have data for
set_color resi2 = [1.0, 0.0, 0.0]
set_color resi3 = [0.0, 0.0, 1.0]
set_color resi4 = [0.5751872021678162, 0.0, 0.42481279783218384]
set_color resi5 = [0.4833857424688618, 0.0, 0.5166142575311382]

#Coloring each residue according to its color
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
