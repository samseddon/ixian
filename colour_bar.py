import matplotlib.pyplot as plt
import matplotlib as mpl

"""
Created on Thu Jan 12 09:32:48 2023

@author: samseddon

NonLinCdict will make a colour bar for whatever set of colours (list) and 
position between 0 and 1 that you give it (with the first and last element of
the list having to be 0 and 1). Mostly we want linear of course, which is what
linear_dist returns when you parse it a list of colours
"""

def NonLinCdict(steps, hexcol_array):
    cdict = {"red": (), "green": (), "blue": ()}
    for s, hexcol in zip(steps, hexcol_array):
        rgb = mpl.colors.hex2color(hexcol)
        cdict["red"]   = cdict["red"] + ((s, rgb[0], rgb[0]),)
        cdict["green"] = cdict["green"] + ((s, rgb[1], rgb[1]),)
        cdict["blue"]  = cdict["blue"] + ((s, rgb[2], rgb[2]),)
    return cdict


def linear_dist(cbar):
    cbdist = []
    for _ in range(len(cbar)):
        cbdist.append(_/(len(cbar)-1))
    return cbdist


# Yellow through Purple through blue, made by samseddon
bladerunner = ["#003f5c",
               "#374c80",
               "#7a5195",
               "#bc5090",
               "#ef5675",
               "#ff764a",
               "#ffa600"]
# Gwyddion colours
code_V = ["#05059A",
          "#0024C0",
          "#005CD5",
          "#01A2CF",
          "#00D0A6",
          "#01D354",
          "#44D600",
          "#9DDD00",
          "#E7E200",
          "#F5A201",
          "#FC5B00",
          "#E90E0A"]
sky = ["#000000", 
       "#262965", 
       "#4B6477", 
       "#CA7A3F", 
       "#FCD356", 
       "#FFFFFF"]          
gwyddion = ["#000000",
            "#A8280F",
            "#F3C25D",
            "#FFFFFF"]
# Adapted colour bar for TUD concept colours from code V:
code_TUD = [#"#00305d",
            "#59358c",
            "#00a1d9",
            "#00aca9",
            "#65b32e",
            "#E7B500",
            "#C94B17",
            "#cd1519"]
#Homemade by Malsch and Milde, for MFM:
mag_sky_replacer = ["#000000",
                    "#10041F",
                    "#230C42",
                    "#361265",
                    "#461785", 
                    "#561AA4", 
                    "#5D31AD", 
                    "#6344B3", 
                    "#6855B8",
                    "#6A66BB", 
                    "#6C75BB", 
                    "#7185AF", 
                    "#74949E", 
                    "#76A483", 
                    "#76B35F", 
                    "#77C03F", 
                    "#8BC74C", 
                    "#9ECD59", 
                    "#B2D365", 
                    "#C5D971",
                    "#D5DF7F",
                    "#DEE695", 
                    "#E7EDAC", 
                    "#F1F4C5", 
                    "#FBFAE4", 
                    "#FFFFFF"]
bladerunner_cbar = mpl.colors.LinearSegmentedColormap("bladerunner",
                        NonLinCdict(linear_dist(bladerunner), bladerunner))
code_V_cbar = mpl.colors.LinearSegmentedColormap("code_V",
                        NonLinCdict(linear_dist(code_V), code_V))
mag_sky_replacer_cbar = mpl.colors.LinearSegmentedColormap("mag_sky_replacer",
                        NonLinCdict(linear_dist(mag_sky_replacer), mag_sky_replacer))
sky_cbar = mpl.colors.LinearSegmentedColormap("sky",
                        NonLinCdict(linear_dist(sky), sky))
code_TUD_cbar = mpl.colors.LinearSegmentedColormap("code_TUD",
                        NonLinCdict(linear_dist(code_TUD), code_TUD))
gwyddion_cbar = mpl.colors.LinearSegmentedColormap("gwyddion",
                        NonLinCdict(linear_dist(gwyddion), gwyddion))

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    # NOTE this will plot will colour bar, put desired colour bar here as cmap = 
    cmap = code_TUD_cbar
 #   cmap = mag_sky_replacer_cbar
    norm = mpl.colors.Normalize(vmin=5, vmax=10)  
    
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                    norm=norm,    
                                    orientation="horizontal")
    cb1.set_label("Some Units")
    plt.show()                              

