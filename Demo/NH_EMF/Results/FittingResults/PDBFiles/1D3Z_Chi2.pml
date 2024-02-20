load 1D3Z.pdb

hide

show cartoon, all

#Color the whole structure in grey as default color
color grey, 1D3Z

#Define a color in the blue to red range for each residue we have data for
set_color resi2 = [1.0, 0.0, 0.0]
set_color resi3 = [0.022699568778640673, 0.0, 0.9773004312213593]
set_color resi4 = [0.03951850794966897, 0.0, 0.9604814920503311]
set_color resi5 = [0.0, 0.0, 1.0]

#Coloring each residue according to its color
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
