module MixedVoxels

# TODO: clean up includings
include("utils.jl")

export intersection_points, volumefraction, polygon_area, sort_vertices
export Geometry, Plane, add!, discretise

abstract GeometryElement

type Plane <: GeometryElement
	n_p :: Array{Float64, 1} # outer normal
	d :: Float64
end

Plane() = Plane([0,0,0], 0)

#type arc

type Geometry
	geometry_elements :: Array{GeometryElement, 1}
	# extensions in housing coordinate sys, containing geometry
	# test here for intersection
	extent :: Array{Uint64, 1}
	h :: Float64
end

Geometry() = Geometry(GeometryElement[], [0,0,0], 0.0)
Geometry(extent::Array{Int64, 1}, h::Float64) = Geometry(GeometryElement[], extent, h)

function add!(geom::Geometry, element::GeometryElement)
	push!(geom.geometry_elements, element)
end

# intersect a plane with the einheits voxel
function intersection_points(plane::Plane)
	
	intersectionpoints = (Array{Float64,1})[]

	# minimal distance of an intersection point from a cube corner
	const eps = 1e-7

	# no intersection possible if
	if abs(plane.d) > sqrt(3.0)
		return intersectionpoints
	end

	# consider each corner (x,y,z) in {0,1}^3
	for z in (0, 1)
		for y in (0, 1)
			for x in (0, 1)
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

					# only consider edges in positive directions
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
			end
		end
	end # continue with another corner

	return intersectionpoints
end

# sort vertices of a polygon in positive (counterclockwise) direction
# around the normal vector n_p
#isless(::Array{Float64,1}, ::Array{Float64,1})
# TODO: multiple dispatch
function sort_vertices(vertices, n_p)

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
# use sort() for that
function polygon_area(vertices)

	v = [0.0, 0.0, 0.0]

	for i in 1:length(vertices)
		v += cross(vertices[i], vertices[mod(i, length(vertices))+1])
	end

	return 0.5 * norm(v)
end

#=function g_x(y, z, n_p, d)
	#
	x = (d - n_p[2] * y - n_p[3] * z) / n_p[1]

	if x > 1.0
		return 1.0
	elseif x < 0.0
		return 0.0
	else
		return x
	end
end


function g_z(x, y, n_p, d)
	#
	z = (d - n_p[1] * x - n_p[2] * y) / n_p[3]
	if z > 1.0
		return 1.0
	elseif z < 0.0
		return z
	else
		return z
	end
end

function Tx(n_p, d, y)
	integral = 0.0
	x0 = g_x(y, 0, n_p, d)
	x1 = g_x(y, 1, n_p, d)
	if x0 < x1
		c = g_z(0, y, n_p, d)
		m = -n_p[1]/n_p[3]
		integral = m/3 * (x1^3-x0^3) + 0.5*c*(x1^2-x0^2)
	end
	integral += -0.5*x1^2+0.5

	return -n_p[1]/n_p[2]*integral
end

# calculate flow of vector field v=(x+1/n_1, -n_1/n_2*x, 0) through plane <x,n_p>=d
# it is div(v)=1 and dot(v,n)=1
function flow(n_i, A, n_p, d)
	# consider different cases of plane's normal vector n_p
	# choose appropriate vector field

	# case n_p.y == 0 and n_p.z == 0
	if equalszero(n_p[2]) && equalszero(n_p[3])
		if (n_i == n_p)
			return 1+d
		else
			return 0.0
		end
		
	# case n_p.y == 0
	# v=(x+1/n_1, 0, -n_1/n_3*x)
	elseif equalszero(n_p[2])
		if		n_i == [1, 0, 0] return (1.0 + 1/n_p[1]) * A
		elseif n_i == [-1,0, 0] return -1/n_p[1] * A
		elseif n_i == [0, 1, 0] return 0.0
		elseif n_i == [0,-1, 0] return 0.0
		elseif n_i == [0, 0, 1]
			return -n_p[1]/n_p[3]*(0.5-0.5*g_x(0,1,n_p,d)^2)
		elseif n_i == [0, 0,-1]
			return n_p[1]/n_p[3]*(0.5-0.5*g_x(0,0,n_p,d)^2)
		elseif n_i == n_p return A
		else return 0.0
		end
	
	# case n_p.z == 0
	# v=(x+1/n_1, -n_1/n_2*x, 0)
	elseif equalszero(n_p[3])
		if		n_i == [1, 0, 0] return (1.0+1.0/n_p[1]) * A
		elseif n_i == [-1,0, 0] return -1.0/n_p[1] * A
		elseif n_i == [0, 1, 0]
			return -n_p[1]/n_p[2]*(0.5-0.5*g_x(1,0,n_p,d)^2)
		elseif n_i == [0,-1, 0]
			return n_p[1]/n_p[2]*(0.5-0.5*g_x(0,0,n_p,d)^2)
		elseif n_i == [0, 0, 1] return 0.0
		elseif n_i == [0, 0,-1] return 0.0
		elseif n_i == n_p return A
		else return 0.0
		end
	
	# general case, n_p.y != 0 and n_p.z != 0
	# FIXME: noch altes Vektorfeld v=(x, -n_1/n_2*x, 1/n_3), updaten auf obiges
	else
		if		n_i == [1, 0, 0] return A
		elseif n_i == [-1,0, 0] return 0
		elseif n_i == [0, 0, 1] return A/n_p[3]
		elseif n_i == [0, 0,-1] return -A/n_p[3]
		elseif n_i == [0, 1, 0] return Tx(n_p, d, 1) # integrate
		elseif n_i == [0,-1, 0] return -Tx(n_p, d, 0) # integrate
		elseif n_i == n_p return A
		else return 0.0
		end
	end
end=#


# calculate volume of polyhedron by calculating the flow
# of a certain vector field through the surface

function volumefraction(plane::Plane, vertices)

	# if there's not intersection, it's not a mixed voxel
	if length(vertices) == 0
		return -1.0
	end

	volume = 0.0

	# cube normals
	const cubenormals = {[1, 0, 0], [-1, 0, 0], [0, 1, 0],
		[0, -1, 0], [0, 0, 1], [0, 0, -1]}

	# cube vertices
	const corners = {[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
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


function intersect(plane::Plane, voxel::Array{Uint64, 1}, h::Float64)

	#  origin of voxel
	x0 = (voxel .- 1) * h

	# translate plane into unit cube system
	d = (dot(x0, plane.n_p) - plane.d) / h

	# "-" because plane is moved against normal direction
	return intersection_points(Plane(plane.n_p, -d))

end



function inside(plane::Plane, vertex::Array{Float64, 1})
	
	if dot(vertex, plane.n_p) > plane.d
		return false
	else
		return true
	end
end

function discretise(geom::Geometry, NOffOrigin::Array{Int64, 1}, mixedvoxel::Function, data1, data2)
	Nx, Ny, Nz = geom.extent
	# for each interior voxel
	# TODO if iterior, compare to extension of FiltEST_VTI
	for z in 1:Nz, y in 1:Ny, x in 1:Nx
		cut=false
		voxel = [x, y, z]
		pos = voxel-NOffOrigin
		
		for element in geom.geometry_elements

			intersectionpoints = intersect(element, voxel, geom.h)

			# treat mixed voxels
			if length(intersectionpoints) > 0
				xi = 1-volumefraction(element, intersectionpoints)

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
		if dot(center, n_p) - d1 < 0 && dot(center, n_p) - d2 > 0
			# voxel is entirely contained in filter medium: full porous voxel
			material.data[pos...] = Mt["Porous"]
			permeability.data[pos...] = K_0
		else
			# else voxel is fluid
			material.data[pos...] = Mt["Fluid"]
			permeability.data[pos...] = 1.0
		end
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
