import DataFrames
import Gadfly

Phi_0_ = 0.1
phi = 10

# read results
table = DataFrames.readtable("results_good.csv")

# select subset
subset = table[DataFrames.array(table[:Phi_0_] .== Phi_0_, 0) & DataFrames.array(table[:phi] .== phi, 0), :]

# color palete
colors1 = [Gadfly.color(c) for c in ["\#e41a1c", "\#377eb8", "\#4daf4a", "\#984ea3", "\#ff7f00", "\#ffff33", "\#a65628", "\#f781bf"]]

# TODO: Color order
# better
# tableaufriction
colors2 = [Gadfly.color(c) for c in ["\#d62728", # red: 214,39,40
									"\#1f77b4", # blue: 31,119,180
									"\#2ca02c", # green: 44,160,44
									"\#ff7f0e", # orange: 255,127,14
									"\#9467bd", # purple: 148,103,189
									"\#bcbd22", # greenish: 188,189,34
									"\#17becf", # blueish: 23,190,207
									"\#e377c2", # rose: 227,119,194
									"\#8c564b", # brown: 140,86,75
									"\#7f7f7f"]] # grey: 127,127,127

# layers
layers = Gadfly.Layer[]

# split table into rows with same dVoxel, method
for plotdata in DataFrames.groupby(subset, [:dVoxel, :method])
	# add plot layer
	push!(layers, Gadfly.layer(plotdata,
					x ="theta",
					y ="pressuredrop",
					Gadfly.Geom.point, Gadfly.Geom.line,
					# FIXME: Unterstriche in Legende machen Probleme
					color=["$(plotdata[:method][1]) $(plotdata[:dVoxel][1])"])...)
end

# plot
p = Gadfly.plot(layers,
	Gadfly.Scale.color_discrete_manual(colors2...),
	#Gadfly.Scale.x_continuous(minvalue=90, maxvalue=30),
	Gadfly.Theme(grid_color=Gadfly.color("\#888888")),
	#Gadfly.Guide.xticks(ticks=[1:10]),
	Gadfly.Guide.xlabel("𝜃"),
	Gadfly.Guide.ylabel("pressuredrop \\Delta p"),
	Gadfly.Guide.title("Rotated Filter: Permeability Scaling for mixed Voxels\n \\phi = $phi"))

# save image to file
# TODO: in PGF: deaktiviere Font-Wahl, text momentan in $\text{}$
# \definecolor{mycolorD0D0E0}{rgb}{0.5,0.5,0.5}
# reduziere y Komponente in Titel um ihn nach oben zu verschieben
image = Gadfly.PGF("Phi_0_$(Phi_0_)_phi_$(phi).tex", 12Gadfly.cm, 7.5Gadfly.cm)
Gadfly.draw(image, p)

file = open("Phi_0_$(Phi_0_)_phi_$(phi).tex","r")
tex = readall(file)
close(file)
file = open("Phi_0_$(Phi_0_)_phi_$(phi).tex","w")
# ersetze "\selectfont" durch ""
tex = replace(tex, "\\selectfont", "")
# ersetze "\fontspec{PT Sans Caption}" durch ""
tex = replace(tex, "\\fontspec{PT Sans Caption}", "")
# ersetze "\fontspec{PT Sans}" durch ""
tex = replace(tex, "\\fontspec{PT Sans}", "")
# "text=mycolor6C606B," durch ""
tex = replace(tex, "text=mycolor6C606B,", "")
# "text=mycolor4C404B," durch ""
tex = replace(tex, "text=mycolor4C404B,", "")
# "text=mycolor564A55," durch ""
tex = replace(tex, "text=mycolor564A55,", "")
# "text=mycolor362A35," durch ""
tex = replace(tex, "text=mycolor362A35,", "")
write(file, tex)
close(file)
