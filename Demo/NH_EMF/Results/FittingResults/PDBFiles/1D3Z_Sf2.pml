load 1D3Z.pdb

hide

show cartoon, all

#Color the whole structure in grey as default color
color grey, 1D3Z

#Define a color in the blue to red range for each residue we have data for
set_color resi2 = [0.5139619957365054, 0.0, 0.48603800426349464]
set_color resi3 = [0.8694331213277879, 0.0, 0.13056687867221206]
set_color resi4 = [1.0, 0.0, 0.0]
set_color resi5 = [0.0, 0.0, 1.0]

#Coloring each residue according to its color
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
