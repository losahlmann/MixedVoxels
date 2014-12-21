MixedVoxels
===========

## FiltEST_VTI

* Type Material
	```baremodule Direction
			const NORTH = 0
			const EAST = 1
			const WEST = 2
			const SOUTH = 3
		end```

* Type FiltEST_VTIFile
	* `traversedirection = "ZYX"`
	* `units = "mm"`
	* `dimension = "3"`
	* `WholeExtent`
	* `origin = [0.0,0.0,0.0]`
	* `spacing = 0.0`
	* `voxeldata::DataArray (Vector)`
	* 
	* `set_extent`
	* `set_origin(o::Float(3))`
	* `set_spacing`
	* `add_data(VoxelData)`
	* `write_file(filename, zip::Bool)`
	* `write_FiltEST_header(file::IOStream)`
	* `write_data(file::IOStream, zip::Bool)`
	* `get_material_count`

htol(x)
Converts the endianness of a value from that used by the Host to Little-endian.

* Type DataArray
	* `data::Material` X x Y x Z 
	* `datatype`
	* `name`
	* `dimensions`
	* `format`
	* (offset)

### Tests
* Voxel-Domain 3 * 3 * 7
	* wall
	* filter at 2,2,3
	* inlet/outlet

