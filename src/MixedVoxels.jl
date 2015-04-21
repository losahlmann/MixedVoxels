module MixedVoxels

# TODO: clean up includings
include("utils.jl")

export intersection_points, volumefraction, polygon_area, sort_vertices
export Geometry, Plane, add!, discretise
export Vertices, Vertex

abstract GeometryElement

type Plane <: GeometryElement
	# outward pointing normal
	n_p :: Vector{Float64}
	d :: Float64
end

Plane() = Plane([0,0,0], 0)

# TODO: type arc

type Geometry
	geometry_elements :: Vector{GeometryElement}
	# extensions in housing coordinate sys, containing geometry
	extent :: Vector{Uint64}
	# discretization length
	h :: Float64
end

Geometry() = Geometry(GeometryElement[], [0,0,0], 0.0)
Geometry(extent::Array{Int64, 1}, h::Float64) = Geometry(GeometryElement[], extent, h)

function add!(geom::Geometry, element::GeometryElement)
	push!(geom.geometry_elements, element)
end

typealias Vertex Vector{Float64}
typealias Vertices Vector{Vertex}


# intersect plane with unit cube
function intersection_points(plane::Plane)
	
	intersectionpoints = Vertex[]

	# minimal distance of an intersection point from a cube corner
	const eps = 1e-7

	# no intersection possible if
	if abs(plane.d) > sqrt(3.0)
		return intersectionpoints
	end

	# consider each corner (x,y,z) in {0,1}^3
	for z in (0, 1), y in (0, 1), x in (0, 1)
		# corner
		corner = [x, y, z]

		# check if corner is in plane and plane inside the voxel
		if equalszero(dot(plane.n_p, corner) - plane.d) && abs(dot(plane.n_p, [0.5,0.5,0.5]) - plane.d) < 0.5-eps
			push!(intersectionpoints, corner)
		end

		# skip corner (1,1,1)
		if corner == [1, 1, 1] continue end

		# 3 possible edges for this corner
		for j = 1:3

			# edge direction
			edge = [0, 0, 0]

			# only consider edges in positive direction
			if corner[j] == 1
				continue
			else
				edge[j] = 1
			end

			# scalar product of plane normal and edge direction
			denom = dot(plane.n_p, edge)

			# check edge for intersection with plane
			if equalszero(denom)
				# edge coplanar with plane
				lambda = -1.0
			else
				# intersection point = c+lambda*e
				lambda = (plane.d - dot(plane.n_p, corner)) / denom
			end

			# test for intersection within cube
			if eps < lambda < 1.0-eps
				
				# intersection
				point = corner + lambda * edge

				# add intersection point
				push!(intersectionpoints, point)

				# there's six points of intersection at maximum
				if length(intersectionpoints) == 6
					return intersectionpoints
				end
			end

		end # continue with another edge

	end # continue with another corner

	return intersectionpoints
end

# sort vertices of a polygon in positive (counterclockwise) direction
# around the normal vector n
function sort_vertices(vertices::Vertices, n_p::Vector{Float64})

	# start vector
	v0 = vertices[1]

	# Bubble sort
	for i in 1:length(vertices)-1
		for j in length(vertices):-1:(i+1)

			v1 = vertices[j]
			v2 = vertices[j-1]

			# compare v1 and v2: v1 < v2 ?
			vv = cross(v0-v1, v0-v2)

			if dot(vv, n_p) > 0
				# swap v1 and v2
				vertices[j] = v2
				vertices[j-1] = v1
			end
		end
	end

	return vertices
end

# Gauss's shoelace formula for the area of a polygon
# points must be ordered in positive direction (counterclockwise) around normal vector
# use sort_vertices() for that
function polygon_area(vertices::Vertices)

	v = [0.0, 0.0, 0.0]

	for i in 1:length(vertices)
		v += cross(vertices[i], vertices[mod(i, length(vertices))+1])
	end

	return 0.5 * norm(v)
end


# calculate volume of polyhedron by calculating the flow
# of a certain vector field through the surface
function volumefraction(plane::Plane, vertices::Vertices)

	# if there's not intersection, it's not a mixed voxel
	if length(vertices) == 0
		return -1.0
	end

	volume = 0.0

	# cube normals
	const cubenormals = Vector{Float64}[[1, 0, 0], [-1, 0, 0], [0, 1, 0],
		[0, -1, 0], [0, 0, 1], [0, 0, -1]]

	# cube vertices
	const corners::Vertices = {[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
		[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]}

	# check all 6 cube sides
	for i in 1:6

		# normal vector of cube side
		n_i = cubenormals[i]

		# distance of cube side from origin
		mod(i, 2) == 0 ? d_i = 0.0 : d_i = 1.0

		# choose vertices laying in this side
		v_i = filter(x -> equalszero(dot(x, n_i) - d_i), vertices)

		# choose cube corner points in this side and under plane
		c_i = filter(x -> equalszero(dot(x, n_i) - d_i) && (dot(x, plane.n_p) - plane.d < 0), corners)

		# continue if it's a polygon
		if length(v_i) + length(c_i) > 2
			#c_i = map(float64,c_i)
			A = polygon_area(sort_vertices([v_i, c_i], n_i))
			volume += A * dot(n_i, c_i[1])#flow(n_i, A, n_p, d)
		end
	end

	# last surface: plane cut by the cube
	
	# calculate surface area
	A = polygon_area(sort_vertices(vertices, plane.n_p))
	
	# calculate flow through this polygon and add to integral
	volume += A * dot(plane.n_p, vertices[1]) #flow(n_p, A, n_p, d)

	return volume/3.0
end


function translate(plane::Plane, voxel::Vector{Uint64}, h::Float64)

	#  origin of voxel
	x0 = (voxel .- 1) * h

	# translate plane into unit cube system
	d = (dot(x0, plane.n_p) - plane.d) / h

	# "-" because plane is moved against normal direction
	return Plane(plane.n_p, -d)

end



function inside(plane::Plane, vertex::Vertex)
	
	if dot(vertex, plane.n_p) > plane.d
		return false
	else
		return true
	end
end

function discretise(geom::Geometry, NOffOrigin::Vector{Int64}, mixedvoxel::Function, data1, data2)
	Nx, Ny, Nz = geom.extent
	# for each interior voxel
	# TODO if iterior, compare to extension of FiltEST_VTI
	for z in 1:Nz, y in 1:Ny, x in 1:Nx
		cut=false
		voxel = [x, y, z]
		pos = voxel-NOffOrigin
		
		for element in geom.geometry_elements

			element_loc = translate(element, voxel, geom.h)
			intersectionpoints = intersection_points(element_loc)

			# treat mixed voxels
			if length(intersectionpoints) > 0
				# needs local plane
				xi = 1-volumefraction(element_loc, intersectionpoints)

				m, p = mixedvoxel(xi)
				data1.data[pos...] = m
				data2.data[pos...] = p	

				# only intersection with one geometry element!
				cut=true
				break
			end	
		end

		if cut
			continue
		end

		# IDEE: mache angeschnittene Voxel porös und gehe dann von porösem bis zu porösem Voxel und fülle zwischendrin mit porös

		# check if voxel is entirely contained in filter medium: full porous voxel
		
		# center of voxel
		center = (voxel .- 0.5) * geom.h
		#=if dot(center, n_p) - d1 < 0 && dot(center, n_p) - d2 > 0
			# voxel is entirely contained in filter medium: full porous voxel
			material.data[pos...] = Mt["Porous"]
			permeability.data[pos...] = K_0
		else
			# else voxel is fluid
			material.data[pos...] = Mt["Fluid"]
			permeability.data[pos...] = 1.0
		end=#
		for element in geom.geometry_elements
			if inside(element, center)
				m, p = mixedvoxel(0.0)
				data1.data[pos...] = m
				data2.data[pos...] = p
				continue
			else
				# else voxel is fluid
				m, p = mixedvoxel(1.0)
				data1.data[pos...] = m
				data2.data[pos...] = p
				break
			end
		end
	end
end

end
