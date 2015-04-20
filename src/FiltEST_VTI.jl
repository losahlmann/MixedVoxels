module FiltEST_VTI

import Zlib

export Mt, Material, DataArray, FiltEST_VTIFile, add_data, write_file, read_file


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

Mt = Dict{String, Material}({
	"Solid" => 0,
	"Fluid" => 5,
	"Porous" => 10,
	"Inflow" => 60,
	"Outflow" => 80,
	"Wall" => 101
}) #FIXME: Symmetry, 0 no slip?

MtName = Dict{Material, String}({
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
	datatype :: DataType
	components :: Int64
	data
	datacompressed :: Array{Any, 1}
	compressedblocksizes :: Array{Uint32, 1}
	numberofblocks :: Int64
	lastblocksize :: Int64
end

DataArray() = DataArray(Any, 0, [], [], [], 0, 0)
DataArray(datatype::DataType, components::Int64, data) = DataArray(datatype, components, data, [], [], 0, 0)


# FiltEST-VTIFile
type FiltEST_VTIFile
	traversedirection :: String
	units :: String
	dimension :: Int64
	origin :: Array{Float64, 1}
	spacing :: Float64
	voxeldata :: Dict{String, DataArray}
	dataorder :: Array{String, 1}
end

FiltEST_VTIFile() = FiltEST_VTIFile("ZYX", "mm", 3, [0.0,0.0,0.0], 0.0, (String => DataArray)[], [])


# TODO: name::UTF8 möglich?
# TODO: !add_data ?
function add_data(vti::FiltEST_VTIFile, name::String, vd::DataArray)

		# add names in order of addition
		push!(vti.dataorder, name)

		# add DataArray to dict
		push!(vti.voxeldata, name, vd)
end

# get number of voxels of each material
function get_material_count(vti::FiltEST_VTIFile)

	# access material data
	# TODO: use try/catch for file handling
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

function g(t::DataType)
	(t==Uint16)? "UInt16" : t
end


# write (compressed) voxel data
function write_data(vti::FiltEST_VTIFile, file::IOStream, zipVTI::Bool)

	# get extents of data
	primarydata = vti.voxeldata[vti.dataorder[1]]
	wholeextent = size(primarydata.data)

	# write dimensions of data
	write(file, """<ImageData WholeExtent="0 $(wholeextent[1]) 0 $(wholeextent[2]) 0 $(wholeextent[3])" Origin="$(join(vti.origin, " "))" Spacing="$(strip(string(vti.spacing," ")^3))">
	<Piece Extent="0 $(wholeextent[1]) 0 $(wholeextent[2]) 0 $(wholeextent[3])">
		<CellData>\n""")

	# in bytes
	const blocksize = 32768
	offset = 0

	# write tags for all DataArrays
	for name in vti.dataorder

		dataarray = vti.voxeldata[name]

		# write tag
		write(file, """\t\t\t<DataArray type="$(g(dataarray.datatype))" Name="$(name)" NumberOfComponents="$(dataarray.components)" format="appended" offset="$offset"/>\n""")

		# calculate data size in bytes
		datasize = length(dataarray.data)*sizeof(dataarray.datatype)

		# prepare data to be written compressed
		if zipVTI
			# clear
			dataarray.datacompressed = []
			dataarray.compressedblocksizes = []

			# TODO: Int64

			# chop data into blocks
			dataarray.numberofblocks = ifloor(datasize/blocksize)

			# remainder
			dataarray.lastblocksize = int64(datasize % blocksize)

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
					block = reinterpret(Uint8, vec(dataarray.data))[(i-1)*blocksize+1:end]
				else
					block = reinterpret(Uint8, vec(dataarray.data))[(i-1)*blocksize+1:i*blocksize]
				end

				# compress with standard compression level
				compressedblock = Zlib.compress(block, 6)

				# note size of compressed block in bytes
				push!(dataarray.compressedblocksizes, length(compressedblock))
				
				# count bytes of compressed data block for offset
				offset += length(compressedblock)

				# save compressed block
				push!(dataarray.datacompressed, compressedblock)
			end
		
			# count bytes of head for offset
			offset += (3+length(dataarray.compressedblocksizes))*sizeof(Uint32)
		else
			offset += datasize + 1
		end
	end

	# close tags
	write(file, """\t\t</CellData>
		</Piece>
	</ImageData>
	<AppendedData encoding="raw">_""")

	# write data of all DataArrays
	for name in vti.dataorder
		dataarray = vti.voxeldata[name]
		if zipVTI
			# write head in UInt32 (==unsigned int)
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
			# write uncompressed
			datasize = length(dataarray.data)*sizeof(dataarray.datatype)
			write(file, uint32(datasize))
			write(file, dataarray.data)
		end
	end

	# close tag
	write(file, "\n\t</AppendedData>")
end


# write FiltEST-VTIFile
function write_file(vti::FiltEST_VTIFile, filename, zipVTI::Bool)

	# TODO: try/catch

	# open file
	file = open(filename, "w")

	# write VTI-XML-header
	# FIXME: trailing space when data not zipped
	write(file, """<?xml version="1.0"?>
		<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian\""""
		* (zipVTI ? """ compressor="vtkZLibDataCompressor">\n""" : ">\n"))

	# write FiltEST-XML-header
	write_FiltEST_header(vti, file)

	# write data
	write_data(vti, file, zipVTI::Bool)

	# finish VTI-XML
	write(file, "\n</VTKFile>")

	# close file
	close(file)

end

# TODO: if no match
function readtagpair(line, tag)
	m = match(Regex("<$tag>(.+)<\/$tag>"), line)
	if m == nothing
		return false
	else
		return m.captures[1]
	end
end

function readopentag(line, tag)
	#print(line)
	m = match(Regex("<$(tag)(?: (.+))?>"), line)
	if m == nothing
		return false
	elseif m.captures[1] == nothing
		return true
	else
		return m.captures[1]
	end
end

function readclosetag(line, tag)
	return ismatch(Regex("<\/$tag>"), line)
end

function readentity(entities, entity)
	m = match(Regex("$(entity)=\"([^\"]+)\""), entities)
	if m == nothing
		return false
	else
		return m.captures[1]
	end
end

function read_FiltEST_VTI(file, vti)
	#header = readuntil(file, delim)
	#line = readline(file)
	#println(line)
	if readline(file) != """<?xml version="1.0"?>\n"""
		error("error reading vti file!")
	elseif readline(file) != """<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n"""
		error("error reading vti file! maybe data not zipped!")
	elseif !readopentag(readline(file), "FiltEST")
		error("No <FiltEST>-tag")
	elseif readtagpair(readline(file), "Date") == false
		# no date: ignore
	end
	vti.traversedirection = readtagpair(readline(file), "TraverseDirection")
	vti.units = readtagpair(readline(file), "Units")
	vti.dimension = int64(readtagpair(readline(file), "Dimension"))

	readuntil(file, "</FiltEST>\n")

	entities = readopentag(readline(file), "ImageData")
	vti.origin = float64(split(readentity(entities, "Origin")))
	vti.spacing = float64(split(readentity(entities, "Spacing")))[1]
	extents = [int(split(readentity(entities, "WholeExtent")))[i] for i=(2,4,6)]

	readopentag(readline(file), "Piece")
	readopentag(readline(file), "CellData")

	typefromstring = {"UInt16" => Uint16,
	"Float64" => Float64}

	line = readline(file)
	while !readclosetag(line, "CellData")
		entities = readopentag(line, "DataArray")
		datatype = typefromstring[readentity(entities, "type")]
		dataarray = DataArray(datatype,
			int64(readentity(entities, "NumberOfComponents")),
			(datatype)[])
		add_data(vti, readentity(entities, "Name"), dataarray)

		line = readline(file)
	end

	if !readclosetag(readline(file), "Piece")
		error("No <Piece>-tag")
	end
	if !readclosetag(readline(file), "ImageData")
		error("No <ImageData>-tag")
	end

	# skip <AppendedData encoding="raw">_
	readuntil(file, """<AppendedData encoding="raw">_""")

	for name in vti.dataorder
		dataarray = vti.voxeldata[name]

		# read number of blocks
		dataarray.numberofblocks = int(read(file, Uint32))

		# read block size before compression
		blocksize = int(read(file, Uint32))

		# read last block size
		dataarray.lastblocksize = int(read(file, Uint32))

		# read sizes of compressed blocks
		dataarray.compressedblocksizes = int(read(file, Uint32, 1))

		# read compressed data
		for s in dataarray.compressedblocksizes
			push!(dataarray.datacompressed, readbytes(file, s))
		end

		# set type of data-array
		#dataarray.data = (dataarray.datatype)[]

		for compressedblock in dataarray.datacompressed
			append!(dataarray.data, reinterpret(dataarray.datatype, Zlib.decompress(compressedblock)))
		end

		dataarray.data = reshape(dataarray.data, extents...)


	end
end


# read FiltEST-VTIFile
function read_file(filename)

	# 
	vti = FiltEST_VTIFile()

	# open file
	file = open(filename, "r")
	#println(readall(file))
	# read VTI-XML-header
	read_FiltEST_VTI(file, vti)

	# close file
	close(file)

	return vti

end

end