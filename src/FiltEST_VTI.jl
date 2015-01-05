module FiltEST_VTI

export Mt, Material, DataArray, FiltEST_VTIFile, add_data, write_file


# Materials
# baremodule Mt
# 	const Solid = 0
# 	const Fluid = 5
# 	const Porous = 10
# 	const Inflow = 60
# 	const Outflow = 80
# 	const Wall = 101
# end

typealias Material Int64

Mt = Dict{ASCIIString, Material}()
Mt = {"Solid" => 0,
	"Fluid" => 5,
	"Porous" => 10,
	"Inflow" => 60,
	"Outflow" => 80,
	"Wall" => 101}

MtName = Dict{Material, ASCIIString}()
MtName = { 0 => "Solid",
	5 => "Fluid",
	10 => "Porous",
	60 => "Inflow",
	80 => "Outflow",
	101 => "Wall"}


# DataArray
type DataArray
	datatype :: DataType
	components :: Int64
	data
end

DataArray() = DataArray(Any, 0, [])


# FiltEST-VTIFile
type FiltEST_VTIFile
	traversedirection
	units
	dimension
	origin
	spacing
	voxeldata :: Dict{ASCIIString, DataArray}
end

FiltEST_VTIFile() = FiltEST_VTIFile("ZYX", "mm", 3, [0.0,0.0,0.0], 0.0, {"none" => DataArray()})


# TODO: name::UTF8 möglich?
function add_data(vti::FiltEST_VTIFile, name::ASCIIString, vd::DataArray)
	if haskey(vti.voxeldata, "none")
		# replace empty DataArray
		vti.voxeldata = {name => vd}
	else
		# add DataArray to vector
		push!(vti.voxeldata, name, vd)
	end
end


# get number of voxels of each material
function get_material_count(vti::FiltEST_VTIFile)

	# access material data
	# FIXME
	#try
		materialdata = vti.voxeldata["Material"].data
	#catch e
	#	error(e)
	#end

	# count voxels
	materials = Dict()

	for voxel in materialdata
		name = MtName[voxel]
		if haskey(materials, name)
			materials[name] = materials[name] + 1
		else
			push!(materials, name, 1)
		end
	end

	# {"Fluid"=>12345, "Solid"=>4321}
	return materials	
end


# write FiltEST-XML-header
function write_FiltEST_header(vti::FiltEST_VTIFile, file::IOStream)

	materials = get_material_count(vti)

	write(file, "<FiltEST>
	<Date>$(strftime(time()))</Date>
	<TraverseDirection>$(vti.traversedirection)</TraverseDirection>
	<Units>$(vti.units)</Units>
	<Dimension>$(vti.dimension)</Dimension>
	<MaterialCount>$(length(materials))</MaterialCount>\n")
	
	totalvoxelcount = 0

	for material in materials
		write(file, "\t<Material>
		<Source>..</Source>
		<Type>$(material[1])</Type>
		<ID>$(Mt[material[1]])</ID>
		<VoxelCount>$(material[2])</VoxelCount>
	</Material>\n")
		totalvoxelcount += material[2]
	end

	write(file, "\t<TotalVoxelCount>$(totalvoxelcount)</TotalVoxelCount>\n</FiltEST>\n")

end


# write (compressed) voxel data
function write_data(vti::FiltEST_VTIFile, file::IOStream, zip::Bool)

	materialdata = vti.voxeldata["Material"].data
	wholeextent = size(materialdata)

	write(file, """<ImageData WholeExtent="0 $(wholeextent[1]) 0 $(wholeextent[2]) 0 $(wholeextent[3])" Origin="$(join(vti.origin, " "))" Spacing="$(strip(string(vti.spacing," ")^3))">
	<Piece Extent="0 $(wholeextent[1]) 0 $(wholeextent[2]) 0 $(wholeextent[3])">
		<CellData>\n""")

	for dataarray in vti.voxeldata
		offset = 0
		write(file, """\t\t\t<DataArray type="$(dataarray[2].datatype)" Name="$(dataarray[1])" NumberOfComponents="$(dataarray[2].components)" format="appended" offset="$offset"/>\n""")
	end

	write(file, """\t\t</CellData>
		</Piece>
	</ImageData>
	<AppendedData encoding="raw">""")

	for dataarray in vti.voxeldata
		# TODO: http://docs.julialang.org/en/release-0.3/manual/strings/#id3
		# http://docs.julialang.org/en/release-0.3/stdlib/base/#strings
		# base64encodepipe
		write(file, dataarray[2].data)
	end

	write(file, "\n</AppendedData>")
end


# write FiltEST-VTIFile
function write_file(vti::FiltEST_VTIFile, filename, zip::Bool)

	# TODO: try/catch

	# open file
	file = open(filename, "w")

	# write VTI-XML-header
	write(file, """<?xml version="1.0"?>\n
		<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian" """
		* (zip ? """compressor="vtkZLibDataCompressor">\n""" : ">\n"))

	# write FiltEST-XML-header
	write_FiltEST_header(vti, file)

	# write data
	write_data(vti, file, zip::Bool)

	# finish VTI-XML
	write(file, "\n</VTKFile>")

	# close file
	close(file)

end

end