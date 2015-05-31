module FiltEST_VTI

import Zlib

export Mt, Material, DataArray, FiltEST_VTIFile
export add_data, write_file, read_file
export readtagpair, readopentag, readclosetag, readentity


typefromstring = Dict{String, DataType}({
	"Int8" => Int8,
	"Int16" => Int16,
	"UInt16" => Uint16,
	"Float32" => Float32,
	"Float64" => Float64
})

# VTI writing of data types can differ from Julia
stringfromtype = Dict{DataType, String}({
	Int8 => "Int8",
	Int16 => "Int16",
	Uint16 => "UInt16",
	Float32 => "Float32",
	Float64 => "Float64"
})


typealias Material Uint16

Mt = Dict{String, Material}({
	"Solid" => 0, # no-slip BC
	"Fluid" => 5,
	"Porous" => 10,
	"Inflow" => 60,
	"Outflow" => 80,
	"Symmetry" => 101
})

MtName = Dict{Material, String}({
	0 => "Solid",
	5 => "Fluid",
	10 => "Porous",
	60 => "Inflow",
	80 => "Outflow",
	101 => "Symmetry"
})


# DataArray
# TODO: Template for DataType
type DataArray
	datatype :: DataType
	components :: Int
	data 
	datacompressed :: Array{Any, 1}
	compressedblocksizes :: Array{Uint32, 1}
	numberofblocks :: Int
	lastblocksize :: Int
end

DataArray() = DataArray(Any, 0, [], [], [], 0, 0)
DataArray(datatype::DataType, components::Int, data) = DataArray(datatype, components, data, [], [], 0, 0)


# FiltEST-VTIFile
type FiltEST_VTIFile
	traversedirection :: String
	units :: String
	dimension :: Int
	origin :: Array{Float64, 1}
	spacing :: Float64
	voxeldata :: Dict{String, DataArray}
	dataorder :: Array{String, 1}
end

FiltEST_VTIFile() = FiltEST_VTIFile("ZYX", "mm", 3, [0.0,0.0,0.0], 0.0, (String => DataArray)[], [])

# add DataArray to FiltEST_VTI-File
function add_data(vti::FiltEST_VTIFile, name::String, vd::DataArray)

		# remember names in order of addition
		push!(vti.dataorder, name)

		# add DataArray to dict
		push!(vti.voxeldata, name, vd)
end

# get number of voxels of each material
function get_material_count(vti::FiltEST_VTIFile)

	# access material data
	#try
		materialdata = vti.voxeldata["Material"].data
	#catch
	#	warn("No Material data!")
	#	return Dict()
	#end

	# create list of materials and their number of voxels
	materials = Dict()

	# count voxels
	for voxel in materialdata
		name = MtName[voxel]
		if haskey(materials, name)
			# material is already in list
			materials[name] +=  1
		else
			# add new material to list
			push!(materials, name, 1)
		end
	end

	# returns {"Fluid"=>12345, "Solid"=>4321}
	return materials	
end


# write FiltEST-XML-header
function write_FiltEST_header(vti::FiltEST_VTIFile, file::IOStream)

	# get list of materials and their voxel count
	materials = get_material_count(vti)

	write(file, "<FiltEST>
	<Date>$(strftime(time()))</Date>
	<TraverseDirection>$(vti.traversedirection)</TraverseDirection>
	<Units>$(vti.units)</Units>
	<Dimension>$(vti.dimension)</Dimension>
	<MaterialCount>$(length(materials))</MaterialCount>\n")
	
	# count total number of voxels
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
function write_data(vti::FiltEST_VTIFile, file::IOStream, zipVTI::Bool)

	# first data to be written (usually "Material")
	primarydata = vti.voxeldata[vti.dataorder[1]]

	# get extents of data
	wholeextent = size(primarydata.data)

	# write dimensions of data
	write(file, """<ImageData WholeExtent="0 $(wholeextent[1]) 0 $(wholeextent[2]) 0 $(wholeextent[3])" Origin="$(join(vti.origin, " "))" Spacing="$(strip(string(vti.spacing," ")^3))">
	<Piece Extent="0 $(wholeextent[1]) 0 $(wholeextent[2]) 0 $(wholeextent[3])">
		<CellData>\n""")

	# block size for compression
	const blocksize = 32768 # bytes

	# offset in data
	offset = 0 # bytes

	# write tags for all DataArrays
	for name in vti.dataorder

		# get DataArray
		dataarray = vti.voxeldata[name]

		# write tag
		write(file, """\t\t\t<DataArray type="$(stringfromtype[dataarray.datatype])" Name="$(name)" NumberOfComponents="$(dataarray.components)" format="appended" offset="$offset"/>\n""")
		# calculate data size in bytes
		datasize = length(dataarray.data) * dataarray.components * sizeof(dataarray.datatype)

		# prepare data to be written compressed
		if zipVTI
			# clear
			dataarray.datacompressed = []
			dataarray.compressedblocksizes = []

			# chop data into blocks
			dataarray.numberofblocks = ifloor(datasize/blocksize)

			# remainder
			dataarray.lastblocksize = int(datasize % blocksize)

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
					block = reinterpret(Uint8, [dataarray.data...])[(i-1)*blocksize+1:end]
				else
					block = reinterpret(Uint8, [dataarray.data...])[(i-1)*blocksize+1:i*blocksize]
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
			# count uncompressed data bytes
			offset += datasize + 1
		end
	end

	# close tags
	write(file, """\t\t</CellData>
		</Piece>
	</ImageData>
	<AppendedData encoding="raw">_""")
	# the underscore at the beginning of the data is necessary

	# write data of all DataArrays
	for name in vti.dataorder

		# get DataArray
		dataarray = vti.voxeldata[name]

		if zipVTI
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
function write_file(vti::FiltEST_VTIFile, filename::String, zipVTI::Bool, header::Bool=true)

	# TODO: try/catch for file handling

	# open file
	file = open(filename, "w")

	# write VTI-XML-header
	# FIXME: trailing space when data not zipped
	write(file, """<?xml version="1.0"?>
		<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian\""""
		* (zipVTI ? """ compressor="vtkZLibDataCompressor">\n""" : ">\n"))

	# write FiltEST-XML-header
	if header
		write_FiltEST_header(vti, file)
	end

	# write data
	write_data(vti, file, zipVTI::Bool)

	# finish VTI-XML
	write(file, "\n</VTKFile>")

	# close file
	close(file)

end

# TODO: if no match
function readtagpair(line::String, tag::String)
	m = match(Regex("<$tag>(.*)<\/$tag>"), line)
	if m == nothing
		return false
	else
		return m.captures[1]
	end
end

function readopentag(line::String, tag::String)
	m = match(Regex("<$(tag)(?: (.+))?>"), line)
	if m == nothing
		return false
	elseif m.captures[1] == nothing
		return true
	else
		return m.captures[1]
	end
end

function readclosetag(line::String, tag::String)
	return ismatch(Regex("<\/$tag>"), line)
end

function readentity(entities::String, entity::String)
	m = match(Regex("$(entity)=\"([^\"]*)\""), entities)
	if m == nothing
		return false
	else
		return m.captures[1]
	end
end

function parse_FiltEST_VTI(file::IOStream, vti::FiltEST_VTIFile)

	# TODO: use try/error in parse methods instead of if-construction
	if readline(file) != """<?xml version="1.0"?>\n"""
		error("error reading vti file!")
	elseif readline(file) != """<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n"""
		error("error reading vti file! maybe data not zipped!")
	end
	line = readline(file)
	if readopentag(line, "FiltEST")
		# if no date: ignore
		readtagpair(readline(file), "Date")
		
		# read header
		vti.traversedirection = readtagpair(readline(file), "TraverseDirection")
		vti.units = readtagpair(readline(file), "Units")
		vti.dimension = int(readtagpair(readline(file), "Dimension"))

		# skip material section
		readuntil(file, "</FiltEST>\n")
		line = readline(file)
	else
		warn("No <FiltEST>-tag")
	end

	# read data parameters
	entities = readopentag(line, "ImageData")
	vti.origin = float64(split(readentity(entities, "Origin")))
	vti.spacing = float64(split(readentity(entities, "Spacing")))[1]
	extents = [int(split(readentity(entities, "WholeExtent")))[i] for i=(2,4,6)]

	# skip tags
	readopentag(readline(file), "Piece")
	readopentag(readline(file), "CellData")

	line = readline(file)
	# read all DataArrays
	while !readclosetag(line, "CellData")
		# read entities
		entities = readopentag(line, "DataArray")

		# get data type as DataType object
		datatype = typefromstring[readentity(entities, "type")]

		# create empty DataArray
		dataarray = DataArray(datatype,
			int(readentity(entities, "NumberOfComponents")),
			(datatype)[])

		# add DataArray to FiltEST-VTI
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

	# read compressed data
	for name in vti.dataorder
		dataarray = vti.voxeldata[name]

		# read number of blocks
		dataarray.numberofblocks = int(read(file, Uint32))

		# read block size before compression
		blocksize = int(read(file, Uint32))

		# read last block size
		dataarray.lastblocksize = int(read(file, Uint32))

		# read sizes of compressed blocks
		dataarray.compressedblocksizes = int(read(file, Uint32, dataarray.numberofblocks))

		# read compressed data
		for s in dataarray.compressedblocksizes
			push!(dataarray.datacompressed, readbytes(file, s))
		end

		# decompress data
		for compressedblock in dataarray.datacompressed
			append!(dataarray.data, reinterpret(dataarray.datatype, Zlib.decompress(compressedblock)))
		end

		# chop data in vectors
		if dataarray.components > 1
			b = (Array{dataarray.datatype, 1})[]
			for i = 1:length(dataarray.data)/dataarray.components
				push!(b, dataarray.data[(i-1)*dataarray.components+1:i*dataarray.components])
			end
			dataarray.data = b
		end

		# resize data vector to actual dimensions
		dataarray.data = reshape(dataarray.data, extents...)
	end

end


# read FiltEST-VTI-File from disk
function read_file(filename::String)

	# create empty FiltEST-VTI object to be filled
	vti = FiltEST_VTIFile()

	# open file to be read
	file = open(filename, "r")

	# parse FiltEST-VTI-File
	parse_FiltEST_VTI(file, vti)

	# close file
	close(file)

	# return filled FiltEST-VTI object
	return vti

end

end