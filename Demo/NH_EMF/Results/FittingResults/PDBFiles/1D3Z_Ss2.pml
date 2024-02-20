load 1D3Z.pdb

hide

show cartoon, all

#Color the whole structure in grey as default color
color grey, 1D3Z

#Define a color in the blue to red range for each residue we have data for
set_color resi2 = [0.0, 0.0, 1.0]
set_color resi3 = [1.0, 0.0, 0.0]
set_color resi4 = [0.10349981358670443, 0.0, 0.8965001864132955]
set_color resi5 = [0.9345184207891791, 0.0, 0.06548157921082087]

#Coloring each residue according to its color
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
color resi{AA}, resi {AA}
