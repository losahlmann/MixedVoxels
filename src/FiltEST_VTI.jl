module FiltEST_VTI

import Zlib

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

# FIXME: UInt16 in *.vti als Dataarray type string nötig?
typealias Material Uint16

Mt = Dict{ASCIIString, Material}({
	"Solid" => 0,
	"Fluid" => 5,
	"Porous" => 10,
	"Inflow" => 60,
	"Outflow" => 80,
	"Wall" => 101
}) #FIXME: Symmetry, 0 no slip?

MtName = Dict{Material, ASCIIString}({
	0 => "Solid",
	5 => "Fluid",
	10 => "Porous",
	60 => "Inflow",
	80 => "Outflow",
	101 => "Wall"
})


# DataArray
# TODO: Template for DataType
type DataArray
	name :: ASCIIString
	datatype :: DataType
	components :: Int64
	data
	datacompressed :: Array{Any,1}
	compressedblocksizes :: Array{Uint32,1}
	numberofblocks
	lastblocksize
end

DataArray() = DataArray("", Any, 0, [], [], [], 0, 0)
DataArray(datatype::DataType, components::Int64, data) = DataArray("", datatype, components, data, [], [], 0, 0)


# FiltEST-VTIFile
type FiltEST_VTIFile
	traversedirection
	units
	dimension
	origin
	spacing
	voxeldata :: Array{DataArray, 1} #Dict{ASCIIString, DataArray}
end

FiltEST_VTIFile() = FiltEST_VTIFile("ZYX", "mm", 3, [0.0,0.0,0.0], 0.0, DataArray[])


# TODO: name::UTF8 möglich?
function add_data(vti::FiltEST_VTIFile, name::ASCIIString, vd::DataArray)
		# set name of DataArray
		vd.name = name

		# add DataArray to vector
		push!(vti.voxeldata, vd)
end

function get_data(vti::FiltEST_VTIFile, name::ASCIIString)
	# return data of first DataArray in vti named "name"
	return filter(x->x.name==name, vti.voxeldata)[1].data
end

# get number of voxels of each material
function get_material_count(vti::FiltEST_VTIFile)

	# access material data
	# TODO: use try/catch for file handling
	#try
		materialdata = get_data(vti, "Material")
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

	materialdata = get_data(vti, "Material")
	wholeextent = size(materialdata)

	write(file, """<ImageData WholeExtent="0 $(wholeextent[1]) 0 $(wholeextent[2]) 0 $(wholeextent[3])" Origin="$(join(vti.origin, " "))" Spacing="$(strip(string(vti.spacing," ")^3))">
	<Piece Extent="0 $(wholeextent[1]) 0 $(wholeextent[2]) 0 $(wholeextent[3])">
		<CellData>\n""")

	blocksize = 32768
	offset = 0

	for dataarray in vti.voxeldata

		# write tag
		write(file, """\t\t\t<DataArray type="$(dataarray.datatype)" Name="$(dataarray.name)" NumberOfComponents="$(dataarray.components)" format="appended" offset="$offset"/>\n""")

		# calculate data size in bytes
		datasize = length(dataarray.data)*sizeof(dataarray.datatype)

		if zip
			# prepare data to be written compressed

			# clear
			dataarray.datacompressed = []
			dataarray.compressedblocksizes = []

			# chop data into blocks
			dataarray.numberofblocks = floor(datasize/blocksize)

			# remainder
			dataarray.lastblocksize = datasize % blocksize

			if dataarray.lastblocksize != 0
				# count last block as well
				dataarray.numberofblocks += 1
			else
				# there is only one block
				dataarray.lastblocksize = blocksize
			end

			# compress each block
			for i in 1:dataarray.numberofblocks
				# convert values into bytes
				if i == dataarray.numberofblocks
					# last block
					block = reinterpret(Uint8, vec(dataarray.data)[(i-1)*blocksize+1:end])
				else	
					block = reinterpret(Uint8, vec(dataarray.data)[(i-1)*blocksize+1:i*blocksize])
				end

				# compress with standard compression level
				compressedblock = Zlib.compress(block, 6)

				# note size of compressed block in bytes
				push!(dataarray.compressedblocksizes, length(compressedblock))
				
				# count bytes for offset
				offset += length(compressedblock)

				# save compressed block
				push!(dataarray.datacompressed, compressedblock)
			end
		end

		
		# count bytes of head for offset
		offset += (3+length(dataarray.compressedblocksizes))*sizeof(Uint32)
	end

	write(file, """\t\t</CellData>
		</Piece>
	</ImageData>
	<AppendedData encoding="raw">_""")

	for dataarray in vti.voxeldata
		if zip
			# header in UInt32 (==unsigned int)
			## datasize = num_elements*elementsize = NVoxel*sizeof(UInt16)

			# write number of blocks
			write(file, uint32(dataarray.numberofblocks))

			# write block size before compression
			write(file, uint32(blocksize))

			# write last block size
			write(file, uint32(dataarray.lastblocksize))

			# write sizes of compressed blocks
			write(file, dataarray.compressedblocksizes)

			# write compressed data
			write(file, dataarray.datacompressed)
		else
			datasize = length(dataarray.data)*sizeof(dataarray.datatype)
			write(file, uint32(datasize))
			write(file, dataarray.data)
		end
	end

	write(file, "\n\t</AppendedData>")
end


# write FiltEST-VTIFile
function write_file(vti::FiltEST_VTIFile, filename, zip::Bool)

	# TODO: try/catch

	# open file
	file = open(filename, "w")

	# write VTI-XML-header
	# FIXME: trailing space when data not zipped
	write(file, """<?xml version="1.0"?>
		<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian\""""
		* (zip ? """ compressor="vtkZLibDataCompressor">\n""" : ">\n"))

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